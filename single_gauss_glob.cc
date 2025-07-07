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
#include <boost/format.hpp>

using std::vector;
using std::string;
using std::cout;
using std::endl;
using boost::format;

// Helper function to check if object is a time diff histogram
bool isTimeDiffHist(const TObject* obj) {
    return obj->InheritsFrom("TH1") && std::string(obj->GetName()).find("tdiff_") == 0;
}

int main(int argc, char* argv[]) {
    gStyle->SetOptStat(0);

    vector<string> input_files = {
        "p1_out.root",
        "p7_out.root", 
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
        return 1;
    }
	
    int n_entries = combined_hist->GetEntries();
    TCanvas* canvas = new TCanvas("combined_canvas", "Combined Time Difference", 800, 800);
    canvas->SetLeftMargin(0.12);
    canvas->SetTopMargin(0.08);
    canvas->SetRightMargin(0.04);
    

    combined_hist->SetLineColor(kBlack);
    combined_hist->SetTitle("");
    combined_hist->GetXaxis()->SetTitle("Time Difference [ns]");
    combined_hist->GetYaxis()->SetTitle("Number of Events");
    combined_hist->Draw("HIST");

    combined_hist->Fit("gaus", "Q");
    TF1* fit_func = combined_hist->TH1::GetFunction("gaus");

    fit_func->SetLineColor(kRed);
    fit_func->SetLineWidth(2);
    fit_func->Draw("SAME");

    TLegend* legend = new TLegend(0.65, 0.75, 0.9, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->AddEntry((TObject*)0, (boost::format("Events = %d") % n_entries).str().c_str(), "");
    legend->AddEntry(combined_hist, "Combined Data", "l");
    legend->AddEntry(fit_func, "Gaussian Fit", "l");
    legend->AddEntry((TObject*)0, (boost::format("#mu = %.3f ns") % fit_func->GetParameter(1)).str().c_str(), "");
    legend->AddEntry((TObject*)0, (boost::format("#sigma = %.3f ns") % fit_func->GetParameter(2)).str().c_str(), "");
    legend->Draw("SAME");

    // === Step 4: Save ===
    TFile* out_file = new TFile("combined_tdiff_single_gauss.root", "RECREATE");
    out_file->cd();

    canvas->SaveAs("combined_tdiff_single_gauss.pdf");
    canvas->Write("combinedcanvas_single_gauss");
    combined_hist->Write(combined_hist->GetName());
    
    out_file->Close();
    cout << "Combined histogram and canvas saved to combined_tdiff_single_gauss.root and combined_tdiff_single_gauss.pdf" << endl;
    return 0;
}
