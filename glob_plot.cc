#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TF1.h>
#include <TStyle.h>
#include <TDirectory.h>
#include <TKey.h>
#include <iostream>
#include <vector>
#include <string>

using std::vector;
using std::string;
using std::cout;
using std::endl;

// Helper function to check if object is a time diff histogram
bool isTimeDiffHist(const TObject* obj) {
    return obj->InheritsFrom("TH1") && std::string(obj->GetName()).find("tdiff_") == 0;
}

void combine_tdiffs() {
    gStyle->SetOptStat(0);
    
    // === Step 1: List of input files ===
    vector<string> input_files = {
        "p0_out.root",
        "p1_out.root",
        "p2a_out.root", 
        "p2b_out.root",
        "p3_out.root", 
        "p4_out.root", 
        "p5a_out.root",
        "p5b_out.root", 
        "p6_out.root", 
        "p7_out.root", 
        "p8_out.root". 
        "p9_out.root", 
        "p10_out.root", 
        "p11_out.root", 
        "p12_out.root", 
        "p13_out.root", 
        "p14_out.root", 
        "p15_out.root", 
        "p16_out.root", 
        "p17_out.root", 
        "p18_out.root", 
        "p19_out.root", 
        "p20_out.root", 
        "p21_out.root"
    };

    TH1D* combined_hist = nullptr;

    // === Step 2: Loop through files and histograms ===
    for (const auto& fname : input_files) {
        TFile* f = TFile::Open(fname.c_str(), "READ");
        if (!f || f->IsZombie()) {
            cout << "Warning: Cannot open file " << fname << endl;
            continue;
        }

        TIter next(f->GetListOfKeys());
        TKey* key;
        while ((key = (TKey*)next())) {
            TObject* obj = key->ReadObj();
            if (!isTimeDiffHist(obj)) continue;

            TH1D* h = (TH1D*)obj;
            if (!combined_hist) {
                combined_hist = (TH1D*)h->Clone("combined_tdiff");
                combined_hist->SetDirectory(0);  // Detach from file
            } else {
                combined_hist->Add(h);
            }
        }

        f->Close();
    }

    // === Step 3: Fit and plot ===
    if (!combined_hist) {
        cout << "No histograms found." << endl;
        return;
    }

    TCanvas* canvas = new TCanvas("combined_canvas", "Combined Time Difference", 800, 800);
    canvas->SetLeftMargin(0.12);
    canvas->SetTopMargin(0.08);
    canvas->SetRightMargin(0.04);

    combined_hist->Fit("gaus", "Q");
    TF1* fit_func = combined_hist->GetFunction("gaus");

    combined_hist->SetLineColor(kBlack);
    combined_hist->SetTitle("Combined Time Difference");
    combined_hist->GetXaxis()->SetTitle("Time Difference [ns]");
    combined_hist->GetYaxis()->SetTitle("Number of Events");
    combined_hist->Draw("HIST");

    if (fit_func) {
        fit_func->SetLineColor(kRed);
        fit_func->SetLineWidth(2);
        fit_func->Draw("SAME");

        double mean = fit_func->GetParameter(1);
        double sigma = fit_func->GetParameter(2);

        TLegend* legend = new TLegend(0.65, 0.75, 0.9, 0.9);
        legend->SetBorderSize(0);
        legend->SetFillStyle(0);
        legend->SetTextFont(42);
        legend->SetTextSize(0.03);
        legend->AddEntry(combined_hist, "Combined Data", "l");
        legend->AddEntry(fit_func, "Gaussian Fit", "l");
        legend->AddEntry((TObject*)0, Form("#mu = %.3f ns", mean), "");
        legend->AddEntry((TObject*)0, Form("#sigma = %.3f ns", sigma), "");
        legend->Draw("SAME");
    }

    // === Step 4: Save ===
    canvas->SaveAs("combined_tdiff.pdf");

    TFile out_file("combined_tdiff.root", "RECREATE");
    combined_hist->Write();
    canvas->Write();
    out_file.Close();

    cout << "Combined histogram and canvas saved to combined_tdiff.root and combined_tdiff.pdf" << endl;
}
