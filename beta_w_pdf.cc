#include <TFile.h>
#include <TH1.h>
#include <TCanvas.h>
#include <TF1.h>
#include <TMath.h>
#include <TApplication.h>
#include <iostream>

int main(int argc, char** argv) {
    TApplication app("app", &argc, argv);

    // Parameters
    double mu = 3.40;
    double tau = TMath::Sqrt2() * 0.273;

    // Open ROOT file
    TFile *f = new TFile("beta_all.root", "READ");
    if (!f || f->IsZombie()) {
        std::cerr << "Could not open beta_all.root" << std::endl;
        return 1;
    }

    // Get histogram (replace "hbeta" with correct hist name)
    TH1 *h = dynamic_cast<TH1*>(f->Get("beta_all"));
    if (!h) {
        std::cerr << "Histogram not found in file" << std::endl;
        return 1;
    }

    // Mirror x-axis
    TH1 *h_mirror = (TH1*)h->Clone("h_mirror");
    h_mirror->Reset();

    for (int i = 1; i <= h->GetNbinsX(); ++i) {
        double x = h->GetBinCenter(i);
        double y = h->GetBinContent(i);

        int bin_mirror = h_mirror->FindBin(-x);
        h_mirror->SetBinContent(bin_mirror, y);
    }

    // Normalize histogram to unit area
    h_mirror->Scale(1.0 / h_mirror->Integral("width"));

    // Define normalized PDF (already normalized analytically)
    TF1 *pdf = new TF1("pdf", [mu, tau](double *x, double *) {
        double beta = x[0];
        if (beta == 0) return 0.0; // avoid division by zero
        double val = mu / (TMath::Sqrt(2*TMath::Pi()) * tau * TMath::Power(beta,4));
        val *= TMath::Exp(-TMath::Power(mu/beta - mu, 2) / (2*tau*tau));
        return val;
    }, h_mirror->GetXaxis()->GetXmin(), h_mirror->GetXaxis()->GetXmax(), 0);

    // Draw
    TCanvas *c = new TCanvas("c", "Beta Histogram vs Normalized PDF", 800, 600);
    h_mirror->Draw("HIST");
    pdf->SetLineColor(kRed);
    pdf->SetLineWidth(2);
    pdf->Draw("SAME");

    c->SaveAs("beta_with_pdf.png");
    
    TFile *fout = new TFile("beta_with_pdf.root", "RECREATE");
    h_mirror->Write();   // save the mirrored/normalized histogram
    pdf->Write();        // save the TF1 PDF function
    c->Write();          // optional: save the canvas
    fout->Close();
    delete fout;

    app.Run();
    return 0;
}

