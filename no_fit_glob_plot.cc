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
    return obj->InheritsFrom("TH1") && std::string(obj->GetName()).find("p") == 0;
}

int main(int argc, char* argv[]) {
    gStyle->SetOptStat(0);

    vector<string> input_files = {
        "p1_edep_out.root",
        "p2a_edep_out.root",
	"p2b_edep_out.root",
	"p3_edep_out.root",
	"p4_edep_out.root",
	"p5a_edep_out.root",
	"p5b_edep_out.root",
	"p6_edep_out.root",
	"p7_edep_out.root",
	"p8_edep_out.root", 
	"p9_edep_out.root", 
	"p10_edep_out.root",
	"p11_edep_out.root",
	"p12_edep_out.root", 
	"p13_edep_out.root", 
	"p14_edep_out.root", 
	"p15_edep_out.root", 
	"p16_edep_out.root", 
	"p17_edep_out.root", 
	"p18_edep_out.root", 
	"p19_edep_out.root", 
	"p20_edep_out.root", 
	"p21_edep_out.root"
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
                combined_hist = (TH1D*)h->Clone("combined_edeps");
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

    legend->AddEntry((TObject*)0, (boost::format("Events = %d") % n_entries).str().c_str(), "");
    legend->AddEntry(combined_hist, "Combined Edep Data", "l");
    
    legend->Draw("SAME");

    // === Step 4: Save to .root and .pdf ===
    TFile* out_file = new TFile("combined_edeps_output.root", "RECREATE");
    out_file->cd();

    canvas->Write("combined_canvas");               // Save the canvas
    combined_hist->Write("combined_edeps");         // Save the histogram

    out_file->Close();                              // Close the .root file

    canvas->SaveAs("combined_edeps_output.pdf");    // Save the canvas as a .pdf

    return 0;
}
