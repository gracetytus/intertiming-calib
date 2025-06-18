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
    vector<int> paddle_nums = {1, 2, 3, 4, 5, 6};
    vector<int> paddle_ids = {0, 100, 200, 300, 400, 500};
    vector<double> offsets = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0}; // your paddle-specific offsets here

    for (uint i = 0; i < paddle_ids.size(); i++) {
        paddle_ids[i] += 110000000;
    }

    vector<TH1D*> time_diffs;
    format hist_name_fmt = format("tdiff_%1%_%2%");
    format hist_title_fmt = format("T_{%1%} - T_{%2%}");
    for (uint i = 0; i < paddle_ids.size() - 1; i++) {
        hist_name_fmt % paddle_nums[i] % paddle_nums[i + 1];
        hist_title_fmt % paddle_nums[i] % paddle_nums[i + 1];
        time_diffs.push_back(new TH1D(hist_name_fmt.str().c_str(), hist_title_fmt.str().c_str(), 150, -5, 5));
    }

    vector<double> paddle_times(paddle_nums.size());
    vector<double> paddle_times_raw(paddle_nums.size());
    vector<double> inter_board_diffs;

    GOptionParser* parser = GOptionParser::GetInstance();
    parser->AddProgramDescription("Computes the interpaddle time differences for adjacent TOF paddles");
    parser->AddCommandLineOption<string>("rec_path", "path to instrument data files", "./*", "i");
    parser->AddCommandLineOption<string>("out_file", "name of output root file", "out.root", "o");
    parser->ParseCommandLine(argc, argv);
    parser->Parse();

    string data_path = parser->GetOption<string>("rec_path");
    string out_path = parser->GetOption<string>("out_file");

    TChain* Instrument_Events = new TChain("TreeRec");
    CEventRec* Event = new CEventRec;
    Instrument_Events->SetBranchAddress("Rec", &Event);
    Instrument_Events->Add(data_path.c_str());

    int vol_id = 0, n_relevant_hits = 0;
    progressbar progress(Instrument_Events->GetEntries() / 1000);

    // === Pass 1: Compute rb_offset from RAW hit times ===
    for (uint i = 0; i < Instrument_Events->GetEntries(); i++) {
        Instrument_Events->GetEntry(i);
        if (i % 1000 == 0) progress.update();

        if (Event->GetNTracks() != 1) continue;

        n_relevant_hits = 0;
        fill(paddle_times_raw.begin(), paddle_times_raw.end(), -1);

        for (GRecoHit hit : Event->GetHitSeries()) {
            vol_id = hit.GetVolumeId();
            if (GGeometryObject::IsTofVolume(vol_id)) {
                for (uint j = 0; j < paddle_ids.size(); j++) {
                    if (vol_id == paddle_ids[j]) {
                        paddle_times_raw[j] = hit.GetTime(); // <- raw times only here
                        n_relevant_hits++;
                    }
                }
            }
        }

        if (n_relevant_hits < 2) continue;

        // Only inter-board pairs
        std::vector<std::pair<int, int>> inter_board_pairs = {{0, 1}, {2, 3}, {4, 5}};
        for (auto [idx_rb8, idx_rb16] : inter_board_pairs) {
            if (paddle_times_raw[idx_rb8] > 0 && paddle_times_raw[idx_rb16] > 0) {
                inter_board_diffs.push_back(paddle_times_raw[idx_rb8] - paddle_times_raw[idx_rb16]);
            }
        }
    }

    double rb_offset = 0.0;
    if (!inter_board_diffs.empty()) {
        double sum = std::accumulate(inter_board_diffs.begin(), inter_board_diffs.end(), 0.0);
        double mean = sum / inter_board_diffs.size();

        double sq_sum = 0.0;
        for (double val : inter_board_diffs) sq_sum += val * val;
        double stddev = std::sqrt((sq_sum / inter_board_diffs.size()) - (mean * mean));

        std::vector<double> filtered;
        for (double val : inter_board_diffs) {
            if (std::abs(val - mean) <= 2 * stddev) filtered.push_back(val);
        }

        if (!filtered.empty()) {
            rb_offset = std::accumulate(filtered.begin(), filtered.end(), 0.0) / filtered.size();
        }

        std::cout << "Inter-readout board offset (after 2σ filtering): " << rb_offset << " ns" << std::endl;
    }

    // === Pass 2: Fill histograms using rb_offset + paddle-specific offsets ===
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

        // Apply rb_offset then apply paddle-specific offsets
        for (uint j = 0; j < paddle_times.size(); j++) {
            if (paddle_times_raw[j] > 0) {
                paddle_times[j] = paddle_times_raw[j] - offsets[j];
            }
        }

        for (uint j = 0; j < paddle_ids.size() - 1; j++) {
            if ((paddle_times[j] > 0) && (paddle_times[j + 1] > 0)) {
                double t_diff = paddle_times[j] - paddle_times[j + 1];
                
                // Apply rb_offset only if pair crosses readout boards
                if ((j == 0) || (j == 2) || (j == 4)) {
                    t_diff -= rb_offset;
                }
                
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
        canvas->SetLogy();

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
        time_diffs[i]->Draw("HIST");
        fitted_func->Draw("SAME");
        t->Draw("SAME");

        canvas->SaveAs(pdf_name_fmt.str().c_str());
        canvas->Write(canvas->GetName());
        time_diffs[i]->Write(time_diffs[i]->GetName());
    }

    out_file.Close();
}