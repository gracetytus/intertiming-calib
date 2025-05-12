#include <TH1.h>
#include "CEventRec.hh"
#include <TChain.h>
#include "TCanvas.h"
#include <iostream>
#include <string>
#include <boost/format.hpp>

using std::vector, std::cout, std::endl;


int main(int argc, char * argv[])
{

  //prepare root files to be read
  char FilenameRoot[400];
  sprintf(FilenameRoot,"%s/%s*.root",argv[1], argv[2]); 

  //prepare reconstronstructed event
  CEventRec* Event = new CEventRec;
  TChain * TreeRec = new TChain("TreeRec");
  TreeRec->SetBranchAddress("Rec", &Event);
  TreeRec->Add(FilenameRoot);
  TreeRec->GetEntry(0);

  //constructing histogram
  TH1D* Ht_diff = new TH1D("Ht_diff", "t_66_t_667", 150,-5,5);
  Ht_diff->SetLineWidth(3);
  Ht_diff->SetLineColor(kBlue);
  Ht_diff->SetXTitle("nanoseconds");
  Ht_diff->SetYTitle("n");

  //declaring variables
  int primary_vid   = 100000500;
  int secondary_vid = 100000600;
  int n_rejected    = 0;
  int n_accepted    = 0;
  int n_histo       = 0;


  for (uint i = 0; i < TreeRec->GetEntries(); i++) {
    TreeRec->GetEntry(i); 
      
    bool accepted1 = false;
    bool accepted2 = false;
    double t_diff = 0;
    for (const auto vid: Event->GetVolumeId()) {
        if (vid == primary_vid) {
            accepted1 = true;
        }
        if (vid == secondary_vid) {
            accepted2 = true;
        }
    }
    if (!accepted1 || !accepted2) {
        n_rejected++;
        continue;
    }
    n_accepted++;

    double primary_time   = -1;
    double secondary_time = -1;
    if (Event->GetNTracks() == 1) {
      for(unsigned int k = 0; k < Event->GetTotalEnergyDeposition().size(); k++) {
        unsigned int VolumeId = Event->GetVolumeId().at(k);
        //cout << " -- volumeid " << VolumeId << endl;
        if (VolumeId >= 200000000) {
          continue;
        }
        //cout << " -- time " << Event->GetTime().at(k) << endl;
        if(VolumeId == primary_vid){
            primary_time = Event->GetTime().at(k);
        
          if(VolumeId == secondary_vid){
              secondary_time = Event->GetTime().at(k);
          }
          if (primary_time == -1 || secondary_time == -1) {
            continue;
          }  
          t_diff = primary_time - secondary_time; 
          //cout << "  ==> primary time " << primary_time << " secondary time " << secondary_time << " t diff " << t_diff << endl;
          Ht_diff->Fill(t_diff);
          n_histo++;
          break;
        }
      }
    }
  }
  TCanvas* t_diff_canvas = new TCanvas("Ht_diff", "Ht_diff", 200, 10, 900, 900);
  t_diff_canvas->SetLeftMargin(0.11);
  t_diff_canvas->SetRightMargin(0.04);
  t_diff_canvas->SetTopMargin(0.04);
  t_diff_canvas->SetLogy();
  Ht_diff->Draw("HIST");
  std::string filename = (boost::format("%1%_%2%_t_diff.pdf") % primary_vid % secondary_vid).str();
  t_diff_canvas->SaveAs(filename.c_str());

}
