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
#include <cmath>

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
        "p2a_out.root",
        "p2b_out.root",
        "p7_out.root",
        "p8_out.root",
        "p9_out.root",
        "p10_out.root",
        "p11_out.root",
        "p12_out.root",
        "p13_out.root"
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
                combined_hist = (TH1D*)h->Clone("horiz_tdiff");
                combined_hist->SetDirectory(0);  // Detach from file
            } else {
                h->SetDirectory(0);
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
    TCanvas* canvas = new TCanvas("horiz_canvas", "Combined Time Difference", 800, 800);
    canvas->SetLeftMargin(0.12);
    canvas->SetTopMargin(0.08);
    canvas->SetRightMargin(0.04);

    // Initial Gaussian fit
    double fit_range_min = -5.0;
    double fit_range_max = 5.0;

    TF1* f = new TF1("f", "gaus", fit_range_min, fit_range_max);
    combined_hist->Fit(f, "RQ"); // quiet fit

    combined_hist->SetLineColor(kBlack);
    combined_hist->SetTitle("");
    combined_hist->GetXaxis()->SetTitle("Time Difference [ns]");
    combined_hist->GetYaxis()->SetTitle("Number of Events");
    combined_hist->Draw("HIST");

    f->SetLineColor(kRed);
    f->SetLineWidth(2);
    f->Draw("SAME");

    const double PI = 3.14159265358979323846;

    double amp = f->GetParameter(0);
    double mean = f->GetParameter(1);
    double sigma = f->GetParameter(2);
    //double area = amp * sigma * std::sqrt(2 * PI);

    TLegend* legend = new TLegend(0.62, 0.64, 0.9, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);

    legend->AddEntry((TObject*)0, (format("Events = %d") % n_entries).str().c_str(), "");
    legend->AddEntry(combined_hist, "Combined Data", "l");
    legend->AddEntry(f, "Gaussian Fit", "l");
    legend->AddEntry((TObject*)0, (format("#mu = %.3f ns") % mean).str().c_str(), "");
    legend->AddEntry((TObject*)0, (format("#sigma = %.3f ns") % sigma).str().c_str(), "");
    //legend->AddEntry((TObject*)0, (format("Area = %.0f") % area).str().c_str(), "");
    
    legend->Draw("SAME");

    canvas->SaveAs("horiz_combined_tdiff_single_gauss.pdf");

    TFile* out_file = new TFile("horiz_combined_tdiff_single_gauss.root", "RECREATE");
    out_file->cd();
    canvas->Write("horiz_canvas");
    combined_hist->Write("horiz_tdiff");
    out_file->Close();

    return 0;
}

