#include <TH1D.h>
#include <TH2D.h>
#include <TF1.h>
#include <TStyle.h>
#include <TFile.h>
#include <TTree.h>
#include <TROOT.h>
#include <TLegend.h>
#include <TColor.h>
#include <TProfile.h>
#include <TMath.h>
#include <TChain.h>
#include <TGraph.h>
#include <TGraphAsymmErrors.h>
#include <TLatex.h>
#include <TMinuit.h>
#include <vector>
#include <map>
#include <string>
#include <stdlib.h>

//FIXME: does this work on mac?
#include <sys/stat.h>

#include "CEventMc.hh"
#include "CAnalysisManager.hh"
#include "GAnalysisIdentification.hh"
#include "GBasicTrigger.hh"
#include "GSimulationParameter.hh"
#include "GPreselection.hh"
#include "CraneConstants.hh"
#include "CraneLogging.hh"
#include "GPlottingTools.hh"
#include "CNet.hh"
#include "CBackpropagation.hh"

#include "GGeometry.hh"

#ifdef USE_BOOST_PROGRAM_OPTIONS
#include "GOptionParser.hh"
#include "GFileIO.hh"
#endif

using namespace std;
using namespace Crane::Analysis;
namespace ca = Crane::Analysis;

namespace cl = Crane::Common;
//using Crane::Calibration;

//time stamp
//vol id name tracker
//event rate

int main(int argc, char * argv[]) {
// silence the root output a bit
gErrorIgnoreLevel = 2000; // warning
cl::set_loglevel(cl::LOGLEVEL::info);
//cl::set_loglevel(cl::SEVERITY::debug);

//--------------------  
    
bool Print = false; 

//step size to sample over the files, helpful for testing to not run over the full files, which is potentially slow
int MainLoopScaleFactor = 1;

//------------------------------------

//cuts

//min tracker enery deposition for charged particle track definition
double TrackerCut = 0.4;

double TrackerMpvCut = 0.45;

double TofCutLow = 0;
double TofCutHigh = 500;

//single track definition
//number of off hit tracker hits must be <= 
double OffHitCtrCut = 1;
//number of on hit tracker hits must be >= 
double OnHitCtrCut = 1;

//TOF cut
//single track
int OuterCtrLowCut = 1;
int CubeCtrLowCut = 1;

//antinucleus trigger
int OuterCtrHighCut = 3;
int CubeCtrHighCut = 3;
int TotalCtrHighCut = 8;

//this is for running code in Hawaii to shift displayed times to UTC
unsigned int TimeOffset = 10*3600;
    
//------------------------------------

unsigned int NoThreadCalcColumnDensity = 1;
double SampleLengthCalcColumnDensity = 1;//mm
    
//------------------------------------
    
char FilenameRoot[400];
sprintf(FilenameRoot,"%s/%s*.root",argv[1], argv[2]);   

//prepare reconstronstructed event
CEventRec* Event = new CEventRec;
TChain * TreeRec = new TChain("TreeRec");
TreeRec->SetBranchAddress("Rec", &Event);
TreeRec->Add(FilenameRoot);
TreeRec->GetEntry(0);
Event->SetEventTime(double(Event->GetEventTime())/(1000./64.)+1631030675);//placeholder for fc conversion to unix time

//------------------------------------

//These object are managing the analysis functionality
CAnalysisManager AnalysisManagerRec;
AnalysisManagerRec.SetEvent(Event);

//--------------------------------------

string DirName;
if(argc == 3) DirName = argv[2];
if(argc == 4) DirName = argv[3];

//create directory for output
char SaveDir[600];
sprintf(SaveDir, "mkdir %s", DirName.c_str());

string Directory;
Directory.append(DirName.c_str());
string Prefix;
Prefix.append(DirName.c_str());

for(int i = 4; i < argc; i += 2) 
    {
    sprintf(SaveDir, "%s_%s", SaveDir, argv[i]);
    
    Directory.append("_");
    Directory.append(argv[i]);
    
    Prefix.append("_");
    Prefix.append(argv[i]);
    }

int success = system(SaveDir);
if (success == 0){std::cout << "Directory " << SaveDir <<" created!" << std::endl;};

//--------------------------------------

//This object is managing the plotting
ca::GPlottingTools Plotting;
char text[400];

TH1D * Ht_diff = Plotting.DefineTH1D("Ht_diff", 150, -5, 5, "nanoseconds", "t_66 - t_67", 0.5, 1e5);
//--------------------------------------

cout<<endl<<"Event loop"<<endl;
//cout<<"0%------50%-----100%"<<endl;
int PercentageStep = 5;

int n_rejected = 0;
int n_accepted = 0;
int n_histo    = 0;
for(unsigned int i = 0; i < TreeRec->GetEntries(); i+=MainLoopScaleFactor) {
//for(unsigned int i = 0; i < 100; i+=MainLoopScaleFactor) {
  TreeRec->GetEntry(i); 
  Event->SetEventTime(double(Event->GetEventTime())/(1000./64.)+1631030675);//placeholder for fc conversion to unix time
  
  //all events
  int HitsTof = 0;
  int HitsUmbrella = 0;
  int HitsCortina = 0;
  int HitsCube = 0;   
  bool accepted1 = false;
  bool accepted2 = false;
  double t_diff = 0;
  for (const auto vid: Event->GetVolumeId()) {
      if (vid == 100000500) {
          accepted1 = true;
      }
      if (vid == 100000600) {
          accepted2 = true;
      }
  }
  if (!accepted1 || !accepted2) {
      n_rejected++;
      continue;
  }
  n_accepted++;
  //if (n_accepted == 100) {
  //  break;
  //}
  //cout << "======== Event ========" << endl;
  //cout << " -- n edeps " << Event->GetTotalEnergyDeposition().size() << endl;
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
      if(VolumeId == 100000600){
          primary_time = Event->GetTime().at(k);
      }
      else if(VolumeId == 100000700){
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
  //cout<<"|"<<endl<<endl;
}
//plotting
//--------------------------------------

gROOT->Reset();
TStyle * plain = new TStyle("plain","plain");
plain->SetCanvasBorderMode(0);
plain->SetPadBorderMode(0);
plain->SetPadColor(0);
plain->SetCanvasColor(0);
plain->SetTitleColor(1);
plain->SetStatColor(0);
plain->SetTitleFillColor(0);

gROOT->SetStyle("plain");

int NRGBs = 4;
int NCont = 100;
double stops[] = { 0.00, 0.33, 0.67, 1.00 };

double red[] = {    1,  0,  1,  0 };
double green[] = {  1,  1,  0,  0 };
double blue[] = {   1,  1,  0,  0 };

TColor::CreateGradientColorTable(NRGBs, stops, red,
green, blue, NCont);
gStyle->SetNumberContours(NCont);

//--------------------------------------

TCanvas Ct_diff = TCanvas("Ct_diff", "Ct_diff", 200, 10, 900, 900);

Ct_diff.SetLeftMargin(0.11);
Ct_diff.SetRightMargin(0.04);
Ct_diff.SetTopMargin(0.04);

Ht_diff->SetLineColor(2);
Ht_diff->Draw("hist");

gPad->SetGridx(1);
gPad->SetGridy(1);
gPad->SetLogy(1);

sprintf(text, "%s/%s_t_diff.root", Directory.c_str(), Prefix.c_str());
Ct_diff.SaveAs(text);
sprintf(text, "%s/%s_t_diff.pdf", Directory.c_str(), Prefix.c_str());
Ct_diff.SaveAs(text);
sprintf(text, "%s/%s_t_diff.png", Directory.c_str(), Prefix.c_str());
Ct_diff.SaveAs(text);   

cout << n_rejected << " events did not pass our cut!" << endl;
cout << n_accepted << " events did pass our cut!" << endl;
cout << n_histo << " events were histogrammed!" << endl;
return EXIT_SUCCESS;
}
