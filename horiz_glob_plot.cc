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
	"p7_out.root"
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
    TCanvas* canvas = new TCanvas("combined_canvas", "Combined Time Difference", 800, 800);
    canvas->SetLeftMargin(0.12);
    canvas->SetTopMargin(0.08);
    canvas->SetRightMargin(0.04);

    combined_hist->Fit("gaus", "Q");
    TF1* fit_func = combined_hist->GetFunction("gaus");
    

    combined_hist->SetLineColor(kBlack);
    combined_hist->SetTitle("");
    combined_hist->GetXaxis()->SetTitle("Time Difference [ns]");
    combined_hist->GetYaxis()->SetTitle("Number of Events");
    combined_hist->Draw("HIST");
    
    double fit_range_min = -5.0;
    double fit_range_max = 5.0;
    double hist_mean = combined_hist->GetMean();
    double hist_sigma = combined_hist->GetRMS();
    double hist_peak = combined_hist->GetBinContent(combined_hist->GetMaximumBin());

    double par[6] = {hist_peak, hist_mean, hist_sigma, 0.01*hist_peak, hist_mean, 2*hist_sigma};

    TF1* f = new TF1("f", "gaus(0)+gaus(3)", fit_range_min, fit_range_max);
    f->SetParameters(par);
    combined_hist->Fit(f, "RQ"); // final fit (quiet)

    f->SetLineColor(kRed);
    f->SetLineWidth(2);
    f->Draw("SAME");

    const double PI = 3.14159265358979323846;

    double amp1 = f->GetParameter(0);
    double sigma1 = f->GetParameter(2);
    double area1 = amp1 * sigma1 * std::sqrt(2 * PI);

    double amp2 = f->GetParameter(3);
    double sigma2 = f->GetParameter(5);
    double area2 = amp2 * sigma2 * std::sqrt(2 * PI);

    TLegend* legend = new TLegend(0.62, 0.64, 0.9, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);

    legend->AddEntry((TObject*)0, (boost::format("Events = %d") % n_entries).str().c_str(), "");
    legend->AddEntry(combined_hist, "Combined Data", "l");
    legend->AddEntry(fit_func, "Gaussian Fit", "l");
    legend->AddEntry((TObject*)0, (boost::format("#mu_{1} = %.3f ns") % f->GetParameter(1)).str().c_str(), "");
    legend->AddEntry((TObject*)0, (boost::format("#sigma_{1} = %.3f ns") % f->GetParameter(2)).str().c_str(), "");
    legend->AddEntry((TObject*)0, (boost::format("#mu_{2} = %.3f ns") % f->GetParameter(4)).str().c_str(), "");
    legend->AddEntry((TObject*)0, (boost::format("#sigma_{2} = %.3f ns") % f->GetParameter(5)).str().c_str(), "");
    legend->AddEntry((TObject*)0, (boost::format("Area_{1} = %.0f") % area1).str().c_str(), "");
    legend->AddEntry((TObject*)0, (boost::format("Area_{2} = %.0f") % area2).str().c_str(), "");
    
    legend->Draw("SAME");

    TFile* out_file = new TFile("combined_tdiff_double_gauss.root", "RECREATE");
    out_file->cd();
    canvas->Write("combined_canvas");
    combined_hist->Write("combined_hist");
    f->Write("double_gaussian_fit");
    out_file->Close();

    canvas->SaveAs("combined_tdiff_double_gauss.pdf");

    return 0;
}
