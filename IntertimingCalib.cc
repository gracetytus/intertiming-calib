#include <TH1.h>
#include "CEventRec.hh"
#include <TChain.h>
#include "TCanvas.h"
#include <iostream>
#include <string>
#include <boost/format.hpp>
#include "CAnalysisManager.hh"


using std::vector, std::cout, std::endl;


int main(int argc, char * argv[])
{

  //prepare reconstronstructed event
  CEventRec* Event = new CEventRec;
  TChain * TreeRec = new TChain("TreeRec");
  TreeRec->SetBranchAddress("Rec", &Event);
  TreeRec->Add(argv[1]);
  TreeRec->GetEntry(0);

  //constructing t_diff histogram
  TH1D* Ht_diff = new TH1D("Ht_diff", "t_65_t_66", 150,-5,5);
  Ht_diff->SetLineWidth(3);
  Ht_diff->SetLineColor(kBlue);
  Ht_diff->SetXTitle("nanoseconds");
  Ht_diff->SetYTitle("n");

  //constructing NTracks histogram
  TH1D* HNTracks = new TH1D("HNTracks", "n tracks", 150, -5, 5);
  HNTracks->SetLineWidth(3);
  HNTracks->SetLineColor(kRed);
  HNTracks->SetXTitle("n tracks");
  HNTracks->SetYTitle("n");

  //declaring variables for cuts
  int primary_vid   = 100000400;
  int secondary_vid = 100000500;
  int n_rejected    = 0;
  int n_accepted    = 0;
  int n_histo       = 0;
  double t_diff     = 0;


  vector<int> volume_ids;

  //cuts
  for (uint i = 0; i < TreeRec->GetEntries(); i++) {
    TreeRec->GetEntry(i); 
    t_diff = 0;
    volume_ids = Event->GetVolumeId();
    //if (Event->GetNTracks() != 1) {
      //continue;
    //}
    if(find(volume_ids.begin(), volume_ids.end(), primary_vid) != volume_ids.end()){
      if(find(volume_ids.begin(), volume_ids.end(), secondary_vid) != volume_ids.end()){
        n_accepted++;
      }
      else{
        n_rejected++;
        continue;
      }
    }
    else{
      n_rejected++;
      continue;
    }

    /*
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
*/
    //populating t_diff histogram
    double primary_time   = -1;
    double secondary_time = -1;
    if (Event->GetNTracks() == 1) {
      HNTracks->Fill(Event->GetNTracks());
      for(unsigned int k = 0; k < Event->GetTotalEnergyDeposition().size(); k++) {
        unsigned int VolumeId = Event->GetVolumeId().at(k);
        //debugging
        //cout << " -- volumeid " << VolumeId << endl;
        if (VolumeId >= 200000000) {
          continue;
        }
        //debugging
        //cout << " -- time " << Event->GetTime().at(k) << endl;
        if(VolumeId == primary_vid){
            primary_time = Event->GetTime().at(k);
        }
        if(VolumeId == secondary_vid){
            secondary_time = Event->GetTime().at(k);
        }
        if (primary_time == -1 || secondary_time == -1) {
          continue;
        }  
        t_diff = primary_time - secondary_time; 
        //debugging
        //cout << "  ==> primary time " << primary_time << " secondary time " << secondary_time << " t diff " << t_diff << endl;
        Ht_diff->Fill(t_diff);
        n_histo++;
        break;
        
      }
    }
  }

  //canvas for t_diff histogram
  TCanvas* t_diff_canvas = new TCanvas("Ht_diff", "Ht_diff", 200, 10, 900, 900);
  t_diff_canvas->SetLeftMargin(0.11);
  t_diff_canvas->SetRightMargin(0.04);
  t_diff_canvas->SetTopMargin(0.08);
  t_diff_canvas->SetLogy();
  Ht_diff->Draw("HIST");
  std::string filename = (boost::format("%1%_%2%_t_diff.pdf") % primary_vid % secondary_vid).str();
  t_diff_canvas->SaveAs(filename.c_str());

  //canvas for NTracks histogram
  TCanvas* CNTracks = new TCanvas("HNTracks", "HNTracks", 200, 10, 900, 900);
  CNTracks->SetLeftMargin(0.11);
  CNTracks->SetRightMargin(0.04);
  CNTracks->SetTopMargin(0.08);
  CNTracks->SetLogy();
  HNTracks->Draw("HIST");
  std::string filename1 = (boost::format("%1%_%2%_ntracks.pdf") % primary_vid % secondary_vid).str();
  CNTracks->SaveAs(filename1.c_str());

  //checks
  cout << n_rejected << " events did not pass our cut!" << endl;
  cout << n_accepted << " passed our cut!" << endl;
  cout << n_histo << " events were histogrammed!" << endl;
  return EXIT_SUCCESS;

}
