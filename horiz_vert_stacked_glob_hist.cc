#include <TFile.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <cmath>
#include <sstream>
#include <string>

void compare_tdiffs() {
    const double PI = 3.14159265358979323846;
    double fit_range_min = -5.0;
    double fit_range_max = 5.0;

    // Open ROOT files
    TFile* f_horiz = TFile::Open("horiz_combined_tdiff_double_gauss.root");
    TFile* f_vert  = TFile::Open("vert_combined_tdiff_double_gauss.root");

    if (!f_horiz || !f_vert) {
        std::cerr << "Error opening ROOT files." << std::endl;
        return;
    }

    // Get histograms
    TH1D* h_horiz = (TH1D*)f_horiz->Get("horiz_tdiff");
    TH1D* h_vert  = (TH1D*)f_vert->Get("vert_tdiff");

    if (!h_horiz || !h_vert) {
        std::cerr << "Could not retrieve histograms." << std::endl;
        return;
    }

    h_horiz->SetLineColor(kBlack);
    h_horiz->SetLineWidth(1);
    h_horiz->SetStats(0);

    h_vert->SetLineColor(kBlack);
    h_vert->SetLineWidth(1);
    h_vert->SetStats(0);
	
    TFile* fout = new TFile("tdiff_comparison.root", "RECREATE");
    TCanvas* c = new TCanvas("c", "Tdiff Comparison", 1000, 700);
    c->cd();



    // Create canvas
    //TCanvas* c = new TCanvas("c", "Tdiff Comparison", 1000, 700);
    //c->cd();
    

    double max_y = std::max(h_horiz->GetMaximum(), h_vert->GetMaximum());
    h_horiz->SetMaximum(1.2 * max_y);
    // Draw histograms
    h_horiz->Draw("HIST");
    h_vert->Draw("HIST SAME");

    // Define and fit double Gaussian for both histograms
    auto fit_double_gauss = [&](TH1D* hist, Color_t color) -> TF1* {
        double mean = hist->GetMean();
        double sigma = hist->GetRMS();
        double peak = hist->GetBinContent(hist->GetMaximumBin());

        double par[6] = {peak, mean, sigma, 0.1*peak, mean, 2*sigma};
        TF1* f = new TF1(Form("fit_%s", hist->GetName()), "gaus(0)+gaus(3)", fit_range_min, fit_range_max);
        f->SetParameters(par);
        f->SetLineColor(color);
        f->SetLineWidth(2);
        hist->Fit(f, "RQ"); // Quiet fit

        return f;
    };

    TF1* fit_horiz = fit_double_gauss(h_horiz, kPink+5);
    TF1* fit_vert  = fit_double_gauss(h_vert,  kAzure+5);

    fit_horiz->Draw("SAME");
    fit_vert->Draw("SAME");

    // Create legend
    TLegend* leg = new TLegend(0.62, 0.58, 0.88, 0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->SetTextFont(42);
    leg->SetTextSize(0.03);
    
    auto add_fit_info = [&](const char* label, TF1* f, Color_t color) {
        double area1 = f->GetParameter(0) * f->GetParameter(2) * std::sqrt(2 * PI);
        double area2 = f->GetParameter(3) * f->GetParameter(5) * std::sqrt(2 * PI);
        
	leg->AddEntry(f, label, "l");
        leg->AddEntry((TObject*)nullptr, Form("#mu_{1} = %.3f ns", f->GetParameter(1)), "");
        leg->AddEntry((TObject*)nullptr, Form("#sigma_{1} = %.3f ns", f->GetParameter(2)), "");
        leg->AddEntry((TObject*)nullptr, Form("#mu_{2} = %.3f ns", f->GetParameter(4)), "");
        leg->AddEntry((TObject*)nullptr, Form("#sigma_{2} = %.3f ns", f->GetParameter(5)), "");
        //leg->AddEntry((TObject*)nullptr, Form("Area_{1} = %.0f", area1), "");
        //leg->AddEntry((TObject*)nullptr, Form("Area_{2} = %.0f", area2), "");
    };

    add_fit_info("Horizontal", fit_horiz, kBlue+1);
    add_fit_info("Vertical", fit_vert, kRed+1);

    leg->Draw();

    // Save canvas
    c->SaveAs("tdiff_comparison.pdf");

    //TFile* fout = new TFile("tdiff_comparison.root", "RECREATE");
    fout->cd(); 
    c->Write();
    h_horiz->Write("horiz_tdiff");
    h_vert->Write("vert_tdiff");
    
    fit_horiz->Write();
    fit_vert->Write();
    c->Write("c");
    fout->Close();
}

int main() {
    compare_tdiffs();
    return 0;
}

