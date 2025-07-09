
+6
-2
Lines changed: 6 additions & 2 deletions
Original file line number	Diff line number	Diff line change
@@ -1,185 +1,189 @@
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

    TH1D* tail_edep_hist = new TH1D("tail_edeps", "Energy Depositions of tail events;Edep;Counts", 150, -5, 5);
    TH1D* tail_edep_hist = new TH1D("tail_edeps", "Energy Depositions of tail events;Edep;Counts", 150, 0, 25);

    // initialize vectors that hold space for each paddle time and raw time based on the length of the paddle_nums vector
    vector<double> paddle_times(paddle_nums.size());
    vector<double> paddle_edeps(paddle_nums.size());
    vector<double> paddle_times_raw(paddle_nums.size());

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
    double edep = 0;
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
        fill(paddle_edeps.begin(), paddle_edeps.end(), -1);
        fill(paddle_times_raw.begin(), paddle_times_raw.end(), -1);


        for (GRecoHit hit : Event->GetHitSeries()) {
            vol_id = hit.GetVolumeId();
            edep = hit.GetEdep();
            edep = hit.GetTotalEnergyDeposition();
            if (GGeometryObject::IsTofVolume(vol_id)) {
                for (uint j = 0; j < paddle_ids.size(); j++) {
                    if (vol_id == paddle_ids[j]) {
                        paddle_times_raw[j] = hit.GetTime();
                        paddle_edeps[j] = edep;
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
                if (t_diff >= 0.8) {
                    if (paddle_edeps[j] > 0) tail_edep_hist->Fill(paddle_edeps[j]);
                    if (paddle_edeps[j + 1] > 0) tail_edep_hist->Fill(paddle_edeps[j + 1]);
                }
            }
        }
    }

   // === Output histogram ===
    TFile out_file(out_path.c_str(), "RECREATE");
    out_file.cd();

    TCanvas* canvas = new TCanvas("tail_edep_canvas", "tail_edep_canvas", 200, 10, 900, 900);
    canvas->SetLeftMargin(0.11);
    canvas->SetTopMargin(0.08);
    canvas->SetRightMargin(0.04);

    tail_edep_hist->SetLineColor(kBlack);
    tail_edep_hist->SetTitle("");
    tail_edep_hist->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    tail_edep_hist->GetYaxis()->SetTitle("Number of Events");
    tail_edep_hist->Draw("HIST");

    int n_entries = tail_edep_hist->GetEntries();

    TLegend* legend = new TLegend(0.65, 0.80, 0.9, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.02);
    legend->AddEntry((TObject*)0, (boost::format("Entries = %d") % n_entries).str().c_str(), "");
    legend->AddEntry(tail_edep_hist, "Data", "l");
    legend->Draw("SAME");

    std::string canvas_name = panel_id + "_tail_edep_canvas";
    std::string pdf_name = panel_id + "_tail_edep_hist.pdf";
    std::string hist_name = panel_id + "_tail_edeps";

    canvas->SaveAs(pdf_name.c_str());
    canvas->Write(canvas_name.c_str());
    tail_edep_hist->Write(hist_name.c_str());

    out_file.Close();
}
