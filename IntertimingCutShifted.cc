#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TPaveText.h>
#include <TText.h>
#include <TF1.h>
#include <TLine.h>
#include <cmath>
#include <numeric>
#include <vector>
#include <string>
#include <iostream>
#include "GOptionParser.hh"
#include "CEventRec.hh"
#include "GRecoHit.hh"
#include "GGeometryObject.hh"
#include "progressbar.hpp"
#include <TStyle.h>
#include <TLegend.h>

#include <format>

using std::vector;
using std::string;
using std::format;
using std::cout;
using std::endl;

vector<int> parseCommaSeparatedInts(const string& input) {
    vector<int> result;
    std::stringstream ss(input);
    string token;
    while (std::getline(ss, token, ',')) {
        result.push_back(std::stoi(token));
    }
    return result;
}

vector<float> parseCommaSeparatedFloats(const string& input) {
	vector<float> result;
	std::stringstream ss(input);
	string token;
	while (std::getline(ss, token, ',')) {
		result.push_back(std::stod(token));
	}
	return result;
}

int main(int argc, char* argv[]) {

    gStyle->SetOptStat(0);
    // initialize paddle nums, paddle ids, calculated offsets, and 'base' panel id
    GOptionParser* parser = GOptionParser::GetInstance();
    parser->AddProgramDescription("Computes the interpaddle time differences for adjacent TOF paddles");
    parser->AddCommandLineOption<string>("rec_path", "path to instrument data files", "./*", "i");
    parser->AddCommandLineOption<string>("out_file", "name of output root file", "out.root", "o");
    parser->AddCommandLineOption<string>("paddle_nums", "comma seperated paddle numbers to be considered by this script", "1,2,3", "n");
    parser->AddCommandLineOption<string>("paddle_ids", "comma seperated paddle volume id endings to be considered by this script", "4,5,6", "b");
    parser->AddCommandLineOption<int>("vol_id_base", "base of volume ids to be considered by this script", 000000000, "s");
    parser->AddCommandLineOption<string>("pid", "panel id, used in title of output .root file", "p0", "p");
    parser->AddCommandLineOption<string>("offsets", "comma seperated offsets calculated to shift t_diff distributions", "0, 0, 0", "f");
    parser->ParseCommandLine(argc, argv);
    parser->Parse();

    string paddle_nums_str = parser->GetOption<string>("paddle_nums");
    string paddle_ids_str = parser->GetOption<string>("paddle_ids");
    string offsets_str = parser->GetOption<string>("offsets");

    vector<int> paddle_nums = parseCommaSeparatedInts(paddle_nums_str);
    vector<int> paddle_ids_suffix = parseCommaSeparatedInts(paddle_ids_str);
    vector<float> offsets = parseCommaSeparatedFloats(offsets_str);  
	
    int vol_id_base = parser->GetOption<int>("vol_id_base");
    string panel_id = parser->GetOption<string>("pid");

    vector<int> paddle_ids(paddle_ids_suffix.size());
    for (size_t i = 0; i < paddle_ids_suffix.size(); i++) {
    	paddle_ids[i] = vol_id_base + paddle_ids_suffix[i];
    }

    //create time_diff vector
    //handle histogram naming based on paddle numbers, and placement on canvas
    vector<TH1D*> time_diffs;
    string hist_name, hist_title;
    for (uint i = 0; i < paddle_ids.size() - 1; i++) {
        hist_name = "tdiff_" + paddle_nums[i] + "_" + paddle_nums[i + 1];
        hist_title = "T_{" + paddle_nums[i] + "}" + " - " + "T_{" + paddle_nums[i + 1] + "}";
        time_diffs.push_back(new TH1D(hist_name.c_str(), hist_title.c_str(), 150, -5, 5));
    }
    // initialize vectors that hold space for each paddle time and raw time based on the length of the paddle_nums vector
    vector<double> paddle_times(paddle_nums.size());
    vector<double> paddle_times_raw(paddle_nums.size());
    // remove?
    // initialize inter_board_diff between adjacent RBs
    vector<double> inter_board_diffs;

    string data_path = parser->GetOption<string>("rec_path");
    std::string out_path = parser->GetOption<std::string>("pid") + "_" + parser->GetOption<std::string>("out_file");

    // initialize new TChain to store information, open CEventRec, store events in Instrument_Events
    TChain* Instrument_Events = new TChain("TreeRec");

    Instrument_Events->SetAutoDelete(true);
    CEventRec* Event = new CEventRec;
    Instrument_Events->SetBranchAddress("Rec", &Event);
    Instrument_Events->Add(data_path.c_str());

    // initialize volume id and n_hits, initialize progress bar
    int vol_id = 0, n_relevant_hits = 0;
    progressbar progress(Instrument_Events->GetEntries() / 1000);

    Instrument_Events->SetBranchAddress("Rec", &Event);
    for (uint i = 0; i < Instrument_Events->GetEntries(); i++) {
        Instrument_Events->GetEntry(i);
	if(i%1000==0){
        	progress.update();
	}
    if (Event->GetNTracks() != 1) continue;

    vector<int> tof_track_indices = Event->GetHitTrackIndex();
    bool skip_event = false;
    int first_idx = -2;

    for (size_t j=0; j < tof_track_indices.size(); j++) {
        int idx = tof_track_indices[j];
        if (idx == -1) {
            skip_event = true;
            break;
        }
        if (j==0) {
            first_idx = idx;
        } else if (idx != first_idx) {
            skip_event = true;
            break;
        }
    }

    if (skip_event) continue;
    
    bool is_outer_tof = false;
    bool is_inner_tof = false;
    n_relevant_hits = 0;

    fill(paddle_times.begin(), paddle_times.end(), -1);
    fill(paddle_times_raw.begin(), paddle_times_raw.end(), -1);

    for (GRecoHit hit : Event->GetHitSeries()) {
        vol_id = hit.GetVolumeId();
        if (GGeometryObject::IsTofVolume(vol_id)) {
            if (GGeometryObject::IsUmbrellaVolume(vol_id)) {
                is_outer_tof = true;
            }
            if (GGeometryObject::IsCubeVolume(vol_id)) {
                is_inner_tof = true;
            }
            for (uint j = 0; j < paddle_ids.size(); j++) {
                if (vol_id == paddle_ids[j]) {
                    paddle_times_raw[j] = hit.GetTime();
                    n_relevant_hits++;
                }
            }
        }
    }

    if (n_relevant_hits < 2) continue;

    if (!(is_outer_tof && is_inner_tof)) {
        continue;
    }
    
    for (uint j = 0; j < paddle_times.size(); j++) {
        if (paddle_times_raw[j] > 0) {
            paddle_times[j] = paddle_times_raw[j] - offsets[j];
        }
    }

    for (uint j = 0; j < paddle_ids.size() - 1; j++) {
        if ((paddle_times[j] > 0) && (paddle_times[j + 1] > 0)) {
            double t_diff = paddle_times[j] - paddle_times[j + 1];
            time_diffs[j]->Fill(t_diff);
        }
    }
}

cout << endl;   

    // === Output histograms ===
    //std::string root_filename = out_path + ".root";
    TFile out_file(out_path.c_str(), "recreate");
    out_file.cd();

    TCanvas* canvas;

    for (uint i = 0; i < time_diffs.size(); i++) {
        string canvas_name = "p" + paddle_nums[i] + "" + paddle_nums[i+1];
        canvas = new TCanvas(canvas_name.c_str(), canvas_name.c_str(), 200, 10, 900, 900);
        canvas->SetLeftMargin(0.11);
        canvas->SetTopMargin(0.08);
        canvas->SetRightMargin(0.04);
        //canvas->SetLogy();

        time_diffs[i]->Fit("gaus", "Q");
        TF1* fitted_func = time_diffs[i]->TH1::GetFunction("gaus");

        if (fitted_func) {
            fitted_func->SetLineColor(kRed);
            fitted_func->SetLineWidth(2);

            double par1 = fitted_func->GetParameter(1);
            double par2 = fitted_func->GetParameter(2);
	        int n_entries = time_diffs[i]->GetEntries();

            time_diffs[i]->SetTitle("");
            time_diffs[i]->SetXTitle("Time Difference [ns]");
            time_diffs[i]->SetYTitle("Number of Events");
            time_diffs[i]->SetLineColor(kBlack);
            time_diffs[i]->Draw("HIST");
            fitted_func->Draw("SAME");

	        TLegend* legend = new TLegend(0.72, 0.78, 0.9, 0.9);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            legend->SetTextFont(42);
            legend->SetTextSize(0.02);
            
            string evt_entry = "Events = " + n_entries;
            string mean_entry = format("#mu = {:.3f} ns", par1);
            string stdv_entry = format("#sigma = {:.3f} ns", par2);

	        legend->AddEntry(nullptr, evt_entry.c_str(), "");
            legend->AddEntry(time_diffs[i], "Data", "l");
            legend->AddEntry(fitted_func, "Gaussian Fit", "l");
            legend->AddEntry(nullptr, mean_entry.c_str(), "");
            legend->AddEntry(nullptr, stdv_entry.c_str(), "");
            legend->Draw("SAME");
        }

        else {
            time_diffs[i]->Draw("HIST");
            std::cerr <<"Warning: Fit failed for histogram " << time_diffs[i]->GetName() << std::endl;
        }
       
        //canvas->SaveAs(pdf_name.c_str());
        canvas->Write();
        time_diffs[i]->Write();
    }
    out_file.Write();
    out_file.Close();
}
