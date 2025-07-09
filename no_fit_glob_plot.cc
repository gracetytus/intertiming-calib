#include <TFile.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
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

int main() {
    gStyle->SetOptStat(0);

    vector<string> input_files = {
        "p1_edep_out.root",  "p2a_edep_out.root", "p2b_edep_out.root", "p3_edep_out.root",
        "p4_edep_out.root",  "p5a_edep_out.root", "p5b_edep_out.root", "p6_edep_out.root",
        "p7_edep_out.root",  "p8_edep_out.root",  "p9_edep_out.root",  "p10_edep_out.root",
        "p11_edep_out.root", "p12_edep_out.root", "p13_edep_out.root", "p14_edep_out.root",
        "p15_edep_out.root", "p16_edep_out.root", "p17_edep_out.root", "p18_edep_out.root",
        "p19_edep_out.root", "p20_edep_out.root", "p21_edep_out.root"
    };

    TH1D* combined_hist = nullptr;

    for (const auto& fname : input_files) {
        TFile* file = TFile::Open(fname.c_str(), "READ");
        if (!file || file->IsZombie()) {
            cout << "Warning: Cannot open file " << fname << endl;
            continue;
        }

        // Loop through all keys and grab the first TH1D
        TIter next(file->GetListOfKeys());
        TKey* key;
        bool added = false;

        while ((key = (TKey*)next())) {
            TObject* obj = key->ReadObj();
            TH1D* h = dynamic_cast<TH1D*>(obj);
            if (!h) continue;

            h->SetDirectory(0); // Detach from file
            if (!combined_hist) {
                combined_hist = (TH1D*)h->Clone("combined_edeps");
                combined_hist->SetDirectory(0);
            } else {
                combined_hist->Add(h);
            }

            added = true;
            break; // Stop after adding the first TH1D
        }

        if (!added) {
            cout << "Warning: No TH1D found in file " << fname << endl;
        }

        file->Close();
    }

    if (!combined_hist) {
        cout << "No histograms were combined. Exiting." << endl;
        return 1;
    }

    // === Plot ===
    int n_entries = combined_hist->GetEntries();
    TCanvas* canvas = new TCanvas("combined_canvas", "Combined Energy Depositions", 800, 800);
    canvas->SetLeftMargin(0.12);
    canvas->SetTopMargin(0.08);
    canvas->SetRightMargin(0.04);

    combined_hist->SetLineColor(kBlack);
    combined_hist->SetTitle("");
    combined_hist->GetXaxis()->SetTitle("Energy Deposition [MeV]");
    combined_hist->GetYaxis()->SetTitle("Number of Events");
    combined_hist->Draw("HIST");

    TLegend* legend = new TLegend(0.62, 0.64, 0.9, 0.9);
    legend->SetBorderSize(0);
    legend->SetFillStyle(0);
    legend->SetTextFont(42);
    legend->SetTextSize(0.03);
    legend->AddEntry((TObject*)nullptr, (format("Events = %d") % n_entries).str().c_str(), "");
    legend->AddEntry(combined_hist, "Combined Edep Data", "l");
    legend->Draw("SAME");

    // === Save output ===
    TFile* out_file = new TFile("combined_edeps_output.root", "RECREATE");
    canvas->Write("combined_canvas");
    combined_hist->Write("combined_edeps");
    out_file->Close();

    canvas->SaveAs("combined_edeps_output.pdf");

    return 0;
}

