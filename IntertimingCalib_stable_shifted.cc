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

using std::vector;
using std::string;
using std::cout;
using std::endl;
using boost::format;

int main(int argc, char* argv[]) {
    // initialize paddle nums, paddle ids, calculated offsets, and 'base' panel id
    vector<int> paddle_nums = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12};
    vector<int> paddle_ids = {0, 100, 200, 300, 400, 500, 600, 700, 800, 900, 1000, 1100};
    vector<double> offsets = {0.0, 0.393, -0.066, 0.276, 0.001, 0.300, -0.417, -0.290, -0.173,-0.400, -0.391, -0.427 }; 

    for (uint i = 0; i < paddle_ids.size(); i++) {
        paddle_ids[i] += 110000000;
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

    // initialize GOptionParser for command line interactions with script
    GOptionParser* parser = GOptionParser::GetInstance();
    parser->AddProgramDescription("Computes the interpaddle time differences for adjacent TOF paddles");
    parser->AddCommandLineOption<string>("rec_path", "path to instrument data files", "./*", "i");
    parser->AddCommandLineOption<string>("out_file", "name of output root file", "out.root", "o");
    parser->ParseCommandLine(argc, argv);
    parser->Parse();

    string data_path = parser->GetOption<string>("rec_path");
    string out_path = parser->GetOption<string>("out_file");

    // initialize new TChain to store information, open CEventRec, store events in Instrument_Events
    TChain* Instrument_Events = new TChain("TreeRec");
    CEventRec* Event = new CEventRec;
    Instrument_Events->SetBranchAddress("Rec", &Event);
    Instrument_Events->Add(data_path.c_str());

    // initialize volume id and n_hits, initialize progress bar
    int vol_id = 0, n_relevant_hits = 0;
    progressbar progress(Instrument_Events->GetEntries() / 1000);

    Instrument_Events->SetBranchAddress("Rec", &Event);
    for (uint i = 0; i < Instrument_Events->GetEntries(); i++) {
        Instrument_Events->GetEntry(i);
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

        time_diffs[i]->Fit("gaus");
        TF1* fitted_func = time_diffs[i]->TH1::GetFunction("gaus");
        fitted_func->SetLineColor(kRed);
        fitted_func->SetLineWidth(2);

        double par1 = fitted_func->GetParameter(1);
        double par2 = fitted_func->GetParameter(2);

        std::string text0 = (boost::format("%-16s %1.3f") % "#mu" % par1).str();
        std::string text1 = (boost::format("%-20s %1.3f") % "#sigma" % par2).str();

        TPaveText* t = new TPaveText(0.78, 0.64, 0.98, 0.78, "blNDC");
        TText* fitText = t->AddText("Fit");
        fitText->SetTextSize(0.038);
        TLine* line = t->AddLine(0.0, 0.72, 1.0, 0.72);
        line->SetLineWidth(1);
        t->AddText(text0.c_str());
        t->AddText(text1.c_str());
        t->SetBorderSize(1);
        t->SetTextFont(42);
        t->SetFillColor(0);
        t->SetTextSize(0.025);
        t->SetMargin(0.0009);

        time_diffs[i]->SetXTitle("Time Difference [ns]");
        time_diffs[i]->SetYTitle("Number of Events");
        time_diffs[i]->SetLineColor(kBlack);
        time_diffs[i]->Draw("HIST");
        fitted_func->Draw("SAME");
        t->Draw("SAME");

        canvas->SaveAs(pdf_name_fmt.str().c_str());
        canvas->Write(canvas->GetName());
        time_diffs[i]->Write(time_diffs[i]->GetName());
    }

    out_file.Close();
}