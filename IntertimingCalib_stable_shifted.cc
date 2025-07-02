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
#include <boost/format.hpp>
#include "GOptionParser.hh"
#include "CEventRec.hh"
#include "GRecoHit.hh"
#include "GGeometryObject.hh"
#include "progressbar.hpp"
#include <TStyle.h>
#include <TLegend.h>

using std::vector;
using std::string;
using std::cout;
using std::endl;
using boost::format;

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
    format hist_name_fmt = format("tdiff_%1%_%2%");
    format hist_title_fmt = format("T_{%1%} - T_{%2%}");
    for (uint i = 0; i < paddle_ids.size() - 1; i++) {
        hist_name_fmt % paddle_nums[i] % paddle_nums[i + 1];
        hist_title_fmt % paddle_nums[i] % paddle_nums[i + 1];
        time_diffs.push_back(new TH1D(hist_name_fmt.str().c_str(), hist_title_fmt.str().c_str(), 150, -5, 5));
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

        n_relevant_hits = 0;
        fill(paddle_times.begin(), paddle_times.end(), -1);
        fill(paddle_times_raw.begin(), paddle_times_raw.end(), -1);

        for (GRecoHit hit : Event->GetHitSeries()) {
            vol_id = hit.GetVolumeId();
            if (GGeometryObject::IsTofVolume(vol_id)) {
                for (uint j = 0; j < paddle_ids.size(); j++) {
                    if (vol_id == paddle_ids[j]) {
                        paddle_times_raw[j] = hit.GetTime();
                        n_relevant_hits++;
                    }
                }
            }
        }

        if (n_relevant_hits < 2) continue;

     
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
    TFile out_file(out_path.c_str(), "recreate");
    out_file.cd();

    format canvas_name_fmt = format("p%1%%2%canvas");
    format pdf_name_fmt = format("paddle_%1%_%2%_tdiff.pdf");
    TCanvas* canvas;

    for (uint i = 0; i < time_diffs.size(); i++) {
        canvas_name_fmt % paddle_nums[i] % paddle_nums[i + 1];
        pdf_name_fmt % paddle_nums[i] % paddle_nums[i + 1];
        canvas = new TCanvas(canvas_name_fmt.str().c_str(), canvas_name_fmt.str().c_str(), 200, 10, 900, 900);
        canvas->SetLeftMargin(0.11);
        canvas->SetTopMargin(0.08);
        canvas->SetRightMargin(0.04);
        //canvas->SetLogy();

        if (time_diffs[i]->GetEntries() > 0) {
            time_diffs[i]->Fit("gaus", "Q");
            TF1* fitted_func = time_diffs[i]->TH1::GetFunction("gaus");
            fitted_func->SetLineColor(kRed);
            fitted_func->SetLineWidth(2);
            fitted_func->SetName("Gaussian Fit");

            double par1 = fitted_func->GetParameter(1);
            double par2 = fitted_func->GetParameter(2);

            time_diffs[i]->SetTitle("");
            time_diffs[i]->SetXTitle("Time Difference [ns]");
            time_diffs[i]->SetYTitle("Number of Events");
            time_diffs[i]->SetLineColor(kBlack);
            time_diffs[i]->Draw("HIST");
            time_diffs[i]->SetName("Data");
            fitted_func->Draw("SAME");
        

            TLegend* legend = new TLegend(0.65, 0.75, 0.9, 0.9);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            legend->SetTextFont(42);
            legend->SetTextSize(0.03);

            legend->AddEntry(time_diffs[i], "Data", "l");
            legend->AddEntry(fitted_func, "Gaussian Fit", "l");
            legend->AddEntry((TObject*)0, (boost::format("#mu = %.3f ns") % par1).str().c_str(), "");
            legend->AddEntry((TObject*)0, (boost::format("#sigma = %.3f ns") % par2).str().c_str(), "");
            legend->Draw("SAME");
        } else {
            TLegend* legend = new TLegend(0.65, 0.75, 0.9, 0.9);
            legend->SetBorderSize(0);
            legend->SetFillStyle(0);
            legend->SetTextFont(42);
            legend->SetTextSize(0.03);
            legend->AddEntry(time_diffs[i], "Data (empty)", "l");
            legend->Draw("SAME")
        }

        if (time_diffs[i]->GetEntries() == 0) {
            std::cout << "WARNING: Histogram " << time_diffs[i]->GetName() << " is empty!" << std::endl;
        }


        canvas->SaveAs(pdf_name_fmt.str().c_str());
        canvas->Write(canvas->GetName());
        time_diffs[i]->Write(time_diffs[i]->GetName());
    }

    out_file.Close();
}
