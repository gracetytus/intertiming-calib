#include <TFile.h>
#include <TCanvas.h>
#include <TH1.h>
#include <TF1.h>
#include <TLegend.h>
#include <TString.h>
#include <iostream>
#include <TTree.h>
#include <TSystemDirectory.h>
#include <TList.h>
#include <TSystemFile.h>
#include <TKey.h>
#include <TClass.h>

TH1* find_first_histogram(TFile& f) {
    f.cd();
    TIter nextkey(f.GetListOfKeys());
    TKey* key;
    while ((key = (TKey*)nextkey())) {
        TObject* obj = key->ReadObj();
        if (obj->InheritsFrom("TH1")) {
            return dynamic_cast<TH1*>(obj);
        }
    }
    return nullptr;
}

void fit_histogram(TH1* h, const TString& outrootname) {
    if (!h) return;

    double mean = h->GetMean();
    double sigma = h->GetRMS();
    double peak = h->GetBinContent(h->GetMaximumBin());
    double par[6] = {peak, mean, sigma, 0.1*peak, mean, 2*sigma};

    TF1* f = new TF1("f", "gaus(0)+gaus(3)", mean-3*sigma, mean+3*sigma);
    f->SetParameters(par);

    h->Fit(f, "RQ"); // Quiet fit

    TCanvas* c = new TCanvas("c", "Double Gaussian Fit", 800, 600);
    h->Draw();
    f->SetLineColor(kRed);
    f->Draw("SAME");

    TLegend* leg = new TLegend(0.6,0.65,0.88,0.88);
    leg->SetBorderSize(0);
    leg->SetFillStyle(0);
    leg->AddEntry(h, "Data", "l");
    leg->AddEntry(f, "Double Gaussian Fit", "l");
    leg->AddEntry((TObject*)0, TString::Format("#mu1=%.3f", f->GetParameter(1)), "");
    leg->AddEntry((TObject*)0, TString::Format("#sigma1=%.3f", f->GetParameter(2)), "");
    leg->AddEntry((TObject*)0, TString::Format("#mu2=%.3f", f->GetParameter(4)), "");
    leg->AddEntry((TObject*)0, TString::Format("#sigma2=%.3f", f->GetParameter(5)), "");
    leg->Draw("SAME");

    TFile fout(outrootname, "RECREATE");
    c->Write("canvas");
    h->Write(h->GetName());
    f->Write("double_gaus_fit");
    fout.Close();

    TString pdfname(outrootname);
    pdfname.ReplaceAll(".root", ".pdf");
    c->SaveAs(pdfname);

    delete c;
    delete f;
}

int main(int argc, char** argv) {
    if (argc < 2) {
        std::cerr << "Usage: " << argv[0] << " <directory>" << std::endl;
        return 1;
    }

    TString dirname(argv[1]);
    TSystemDirectory dir(dirname, dirname);
    TList* files = dir.GetListOfFiles();
    if (!files) return 0;

    TIter next(files);
    TSystemFile* file;
    while ((file = (TSystemFile*)next())) {
        TString fname = file->GetName();
        if (!file->IsDirectory() && fname.EndsWith(".root")) {
            TString fullpath = dirname + "/" + fname;
            TFile f(fullpath, "READ");
            if (f.IsZombie()) continue;

            TH1* h = find_first_histogram(f);
            if (!h) {
                std::cerr << "No histogram found in " << fullpath << std::endl;
                continue;
            }

            TString outrootname = fname;
            outrootname.ReplaceAll(".root", "_fit.root");

            fit_histogram(h, outrootname);

            f.Close();
        }
    }

    return 0;
}
