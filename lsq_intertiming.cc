#include <TFile.h>
#include "CEventMc.hh"
#include "GOptionParser.hh"
#include "GFileIO.hh"
#include "CEventRec.hh"
#include "GGeometryObject.hh"
#include "TChain.h"
#include "TH1.h"
#include "TH2.h"
#include "TH3.h"
#include "TCanvas.h"
#include <algorithm>
#include <iostream>
#include <vector>
#include "TGraph.h"
#include "TColor.h"
#include "TStyle.h"
#include "progressbar.hpp"
#include <format>
#include <boost/format.hpp>
#include <TH1.h>
#include <TF1.h>
#include <TPaveText.h>
#include "TFile.h"
#include <TLine.h>
#include <TMatrixD.h>
#include <TVectorD.h>
#include <TDecompLU.h>
#include <TDecompSVD.h>

using std::vector, std::string, std::cout, std::endl;
using boost::format;

int main(int argc, char* argv[]) {
    vector<int> paddle_nums = {79, 80, 81, 82, 83, 84};
    vector<int> paddle_ids = {500, 400, 300, 200, 100, 0};
    for (uint i = 0; i < paddle_ids.size(); i++) {
        paddle_ids[i] = paddle_ids[i] + 100200000;
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
    vector<double> offsets(paddle_ids.size(), 0.0);
    vector<double> means(paddle_ids.size() - 1);

    progressbar progress(Instrument_Events->GetEntries() / 1000);

    for (uint i = 0; i < Instrument_Events->GetEntries(); i++) {
        Instrument_Events->GetEntry(i);
        if (i % 1000 == 0) progress.update();

        if (Event->GetNTracks() != 1) continue;

        n_relevant_hits = 0;
        fill(paddle_times.begin(), paddle_times.end(), -1);

        for (GRecoHit hit : Event->GetHitSeries()) {
            vol_id = hit.GetVolumeId();
            if (GGeometryObject::IsTofVolume(vol_id)) {
                for (uint j = 0; j < paddle_ids.size(); j++) {
                    if (vol_id == paddle_ids[j]) {
                        paddle_times[j] = hit.GetTime();
                        n_relevant_hits++;
                    }
                }
            }
        }

        if (n_relevant_hits < 2) continue;

        for (uint j = 0; j < paddle_ids.size() - 1; j++) {
            if ((paddle_times[j] > 0) && (paddle_times[j + 1] > 0)) {
                time_diffs[j]->Fill(paddle_times[j] - paddle_times[j + 1]);
            }
        }
    }
    cout << endl;

    // Fit histograms to extract means
    for (uint i = 0; i < time_diffs.size(); i++) {
        time_diffs[i]->Fit("gaus", "Q");
        TF1* fit = time_diffs[i]->GetFunction("gaus");
        means[i] = fit->GetParameter(1);
    }

    // Solve least squares: A * c = mu
    const int N = paddle_ids.size();
    TMatrixD A(N - 1, N);
    TVectorD mu(N - 1);

    for (int i = 0; i < N - 1; i++) {
        A[i][i] = 1;
        A[i][i + 1] = -1;
        mu[i] = means[i];
    }

    TMatrixD At(TMatrixD::kTransposed, A);
    TMatrixD AtA = At * A;
    TVectorD AtMu = At * mu;
    TDecompSVD svd(AtA);
    TVectorD c(N);

    Bool_t ok;
    TVectorD solution = svd.Solve(AtMu, ok);

    if (!ok) {
	std::cerr << "Failed to solve least squares system" << std::endl;
	return 1;
    }
    c = solution;

    // Shift offsets so that paddle 79 (index 0) is reference
    double ref = c[0];
    for (int i = 0; i < N; ++i) {
        offsets[i] = c[i] - ref;
    }

    std::cout << "\nFinal paddle timing offsets (relative to paddle 79):" << std::endl;
    for (int i = 0; i < N; ++i) {
        std::cout << "Paddle " << paddle_nums[i] << ": " << offsets[i] << " ns" << std::endl;
    }

    // Clear and refill histograms with corrected times
    for (auto* hist : time_diffs) hist->Reset();

    progressbar progress2(Instrument_Events->GetEntries() / 1000);
    for (uint i = 0; i < Instrument_Events->GetEntries(); i++) {
        Instrument_Events->GetEntry(i);
        if (i % 1000 == 0) progress2.update();

        if (Event->GetNTracks() != 1) continue;
        n_relevant_hits = 0;
        fill(paddle_times.begin(), paddle_times.end(), -1);

        for (GRecoHit hit : Event->GetHitSeries()) {
            vol_id = hit.GetVolumeId();
            if (GGeometryObject::IsTofVolume(vol_id)) {
                for (uint j = 0; j < paddle_ids.size(); j++) {
                    if (vol_id == paddle_ids[j]) {
                        paddle_times[j] = hit.GetTime() - offsets[j];
                        n_relevant_hits++;
                    }
                }
            }
        }

        if (n_relevant_hits < 2) continue;

        for (uint j = 0; j < paddle_ids.size() - 1; j++) {
            if ((paddle_times[j] > 0) && (paddle_times[j + 1] > 0)) {
                time_diffs[j]->Fill(paddle_times[j] - paddle_times[j + 1]);
            }
        }
    }
    cout << endl;

    // Save updated histograms
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
        TF1* fitted_func = time_diffs[i]->GetFunction("gaus");
        fitted_func->SetLineColor(kRed);
        fitted_func->SetLineWidth(2);

        double par1 = fitted_func->GetParameter(1);
        double par2 = fitted_func->GetParameter(2);

        std::string text0 = (boost::format("%-16s %1.3f") % "#mu" % par1).str();
        std::string text1 = (boost::format("%-20s %1.3f") % "#sigma" % par2).str();

        TPaveText* t = new TPaveText(0.78, 0.64, 0.98, 0.78, "blNDC");
        t->AddText("Fit")->SetTextSize(0.038);
        t->AddLine(0.0, 0.72, 1.0, 0.72)->SetLineWidth(1);
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

