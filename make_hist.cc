#include <TH1.h>
#include <TF1.h>
#include <TFitResultPtr.h>
#include <TFitResult.h>

#include <boost/format.hpp>
#include <TPaveText.h>
#include "TFile.h"
#include <TLine.h>

TFile in_file ("/home/gtytus/intertiming-calib/outs/cbe_top.root")

TH1D* Ht_diff = (TH1D*)in_file->Get("t_diff_1_2");
in_file.close()

Ht_diff->SetLineWidth(3);
Ht_diff->SetLineColor(kBlue);
Ht_diff->SetXTitle("nanoseconds");
Ht_diff->SetYTitle("n");

//canvas for t_diff histogram
TCanvas* t_diff_canvas = new TCanvas("Ht_diff", "Ht_diff", 200, 10, 900, 900);
t_diff_canvas->SetLeftMargin(0.11);
t_diff_canvas->SetRightMargin(0.04);
t_diff_canvas->SetTopMargin(0.08);
t_diff_canvas->SetLogy();
Ht_diff->Draw("HIST");

//fitting
Ht_diff->Fit("gaus");
TF1 *fitted_func = Ht_diff->TH1::GetFunction("gaus");
fitted_func->SetLineColor(kRed); 
fitted_func->SetLineWidth(2); 
fitted_func->Draw("SAME");

double chi2 = fitted_func->GetChisquare();
double par0 = fitted_func->GetParameter(0);
double err0 = fitted_func->GetParError(0);
double par1 = fitted_func->GetParameter(1);
double err1 = fitted_func->GetParError(1); 
double par2 = fitted_func->GetParameter(2);
double err2 = fitted_func->GetParError(2);

std::string text0 = (boost::format("%-16s %1.3f") % "#mu" % par1).str();
std::string text1 = (boost::format("%-20s %1.3f") % "#sigma" % par2).str();

TPaveText* t = new TPaveText(0.78, 0.64, 0.98, 0.78, "blNDC");

TText* fitText = t->AddText("Fit");
fitText->SetTextSize(0.038);  

TLine* line = t->AddLine(0.0, 0.72, 1.0, 0.72);
line->SetLineWidth(1); 
t->AddText(text0.c_str());
t->AddText(text1.c_str());

t->SetBorderSize(1);
t->SetTextFont(42);
t->SetFillColor(0);
t->SetTextSize(0.025);
t->SetMargin(0.0009);
t->Draw("SAME");

//saving
std::string filename = (t_diff_1_2)
t_diff_canvas->SaveAs(filename.c_str());

TFile out_file("output.root", "recreate");
out_file.cd();
Ht_diff->Write("Ht_diff");
out_file.Close();