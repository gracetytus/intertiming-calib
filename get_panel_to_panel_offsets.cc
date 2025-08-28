#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <cmath>
#include <numeric>
#include <vector>
#include <string>
#include <iostream>
#include <boost/format.hpp>
#include "GOptionParser.hh"
#include "CEventRec.hh"
#include "GRecoHit.hh"
#include "GGeometryObject.hh"
#include "progressbar.hpp"
#include <fstream>

using std::vector;
using std::string;
using std::cout;
using std::endl;
using boost::format;

//volume ids by paddle
std::vector<int> panel_1_vids  = {110000000, 110000100, 110000200, 110000300, 110000400, 110000500, 110000600, 110000700, 110000800, 110000900, 110001000, 110001100};
std::vector<int> panel_2a_vids = {111001000, 111000900, 111000800, 111000700};
std::vector<int> panel_2b_vids = {110000500, 110000400, 110000300, 110000200, 110000100, 110000000};
std::vector<int> panel_3_vids  = {112000700, 112000600, 112000500, 112000400, 112000300, 112000200, 112000100, 112000000};
std::vector<int> panel_4_vids  = {114000700, 114000600, 114000500, 114000400, 114000300, 114000200, 114000100, 114000000};
std::vector<int> panel_5a_vids = {113000700, 113000600, 110000500};
std::vector<int> panel_5b_vids = {113000200, 113000100, 113000000};
std::vector<int> panel_6_vids  = {11500000, 115000100, 115000200, 115000300, 115000400, 115000500, 115000600, 115000700};
std::vector<int> panel_57_vids = {116000000};
std::vector<int> panel_58_vids = {116200000};
std::vector<int> panel_59_vids = {116300000};
std::vector<int> panel_60_vids = {116100000};
std::vector<int> panel_7_vids  = {100000000, 100000100, 100000200, 100000300, 100000400, 100000500, 100000600, 100000700, 100000800, 100000900, 100001000, 100001100};
std::vector<int> panel_8_vids  = {100300500, 100300400, 100300300, 100300200, 100300100, 100300000};
std::vector<int> panel_9_vids  = {100200500, 100200400, 100200300, 100200200, 100200100, 100200000};
std::vector<int> panel_10_vids = {100400000, 100400100, 100400200, 100400300, 100400400, 100400500};
std::vector<int> panel_11_vids = {100600500, 100600400, 100600300, 100600200, 100600100, 100600000};
std::vector<int> panel_12_vids = {100100400, 100100300, 100100200, 100100100, 100100000};
std::vector<int> panel_13_vids = {100500500, 100500400, 100500300, 100500200, 100500100, 100500000};
std::vector<int> panel_14_vids = {102000900, 102000800, 102000700, 102000600, 102000500, 102000400, 102000300, 102000200, 102000100, 102000000};
std::vector<int> panel_15_vids = {104000000, 104000100, 104000200, 104000300, 104000400, 104000500, 104000600, 104000700, 104000800, 104000900};
std::vector<int> panel_16_vids = {103000900, 103000800, 103000700, 103000600, 103000500, 103000400, 103000300, 103000200, 103000100, 103000000};
std::vector<int> panel_17_vids = {105000900, 105000800, 105000700, 105000600, 105000500, 105000400, 105000300, 105000200, 105000100, 105000000};
std::vector<int> panel_18_vids = {106000200, 106000100, 106000000};
std::vector<int> panel_19_vids = {106200000, 106200100, 106200200};
std::vector<int> panel_20_vids = {106300000, 106300100, 106300200};
std::vector<int> panel_21_vids = {106100200, 106100100, 106100000};

//offsets by paddle
std::vector<double> panel_1_offsets  = {0.000,0.405,-0.064,0.290,0.003,0.601,-0.099,0.030,0.150,-0.074,-0.058,-0.087};
std::vector<double> panel_2a_offsets = {0.000,0.267,-0.228,0.250};
std::vector<double> panel_2b_offsets = {0.000,0.875,0.752,0.681,1.018,0.898};
std::vector<double> panel_3_offsets  = {0.000,-0.397,0.307,0.334,0.503,0.523,0.867,0.698};
std::vector<double> panel_4_offsets  = {0.000,0.386,0.406,0.561,0.229,0.937,0.605,0.910};
std::vector<double> panel_5a_offsets = {0.000,0.234,0.953};
std::vector<double> panel_5b_offsets = {0.000,0.465,0.487};
std::vector<double> panel_6_offsets  = {0.000,0.373,0.207,0.657,0.537,0.834,1.035,1.133};
std::vector<double> panel_7_offsets  = {0.000,0.082,0.158,0.177,0.283,0.360,-0.549,0.091,-0.584,0.037,-0.715,0.128};
std::vector<double> panel_8_offsets  = {0.000,0.006,-0.304,-0.084,-0.415,-0.004};
std::vector<double> panel_9_offsets  = {0.000,-0.228,-0.065,-0.148,0.251,0.166};
std::vector<double> panel_10_offsets = {0.000,0.142,-0.841,-0.319,-0.862,-0.107};
std::vector<double> panel_11_offsets = {0.000,-0.701,0.876,1.726,1.224,1.459};
std::vector<double> panel_12_offsets = {0.000,0.377,0.180,0.533,0.284};
std::vector<double> panel_13_offsets = {0.000,0.484,0.460,0.709,0.781,0.824};
std::vector<double> panel_14_offsets = {0.000,0.082,0.266,0.548,0.090,0.208,0.332,0.366,0.963,0.525};
std::vector<double> panel_15_offsets = {0.000,0.353,0.427,0.703,0.947,1.195,2.390,1.480,0.947,0.643};
std::vector<double> panel_16_offsets = {0.000,-0.019,0.207,0.845,0.688,0.643,0.889,0.669,-0.014,0.887};
std::vector<double> panel_17_offsets = {0.000,-0.442,0.211,-0.157,0.658,0.194,0.921,0.322,1.184,0.403};
std::vector<double> panel_18_offsets = {0.000,-0.026,-0.112};
std::vector<double> panel_19_offsets = {0.000,0.008,0.056};
std::vector<double> panel_20_offsets = {0.000,0.623,-0.173};
std::vector<double> panel_21_offsets = {0.000,-0.439,0.027};
std::vector<double> panel_57_offsets = {0.000};
std::vector<double> panel_58_offsets = {0.000};
std::vector<double> panel_59_offsets = {0.000};
std::vector<double> panel_60_offsets = {0.000};

struct HitInfo {
    double adj_time;
    TVector3 pos;
};

int main(int argc, char* argv[]) {
    GOptionParser* parser = GOptionParser::GetInstance();
    parser->AddProgramDescription("Computes the panel to panel timing offsets for TOF panels");
    parser->AddCommandLineOption<string>("rec_path", "path to instrument data files", "./*", "i");
    parser->AddCommandLineOption<string>("out_file", "name of output root file", "out.root", "o");
    parser->ParseCommandLine(argc, argv);
    parser->Parse();

    string data_path = parser->GetOption<string>("rec_path");
    string out_path = parser->GetOption<std::string>("out_file");

    struct PanelInfo {
        std::string panel;
        size_t index; // index within that panel's offsets vector
    };

    std::map<int, PanelInfo> volid_lookup;
    auto add_panel_mapping = [&](const std::string& name, const std::vector<int>& vids) {
    for (size_t idx = 0; idx < vids.size(); ++idx) {
        volid_lookup[vids[idx]] = {name, idx};
        }
    };

    add_panel_mapping("panel_1", panel_1_vids);
    add_panel_mapping("panel_2a", panel_2a_vids);
    add_panel_mapping("panel_2b", panel_2b_vids);
    add_panel_mapping("panel_3", panel_3_vids);
    add_panel_mapping("panel_4", panel_4_vids);
    add_panel_mapping("panel_5a", panel_5a_vids);
    add_panel_mapping("panel_5b", panel_5b_vids);
    add_panel_mapping("panel_6", panel_6_vids);
    add_panel_mapping("panel_7", panel_7_vids);
    add_panel_mapping("panel_8", panel_8_vids);
    add_panel_mapping("panel_9", panel_9_vids);
    add_panel_mapping("panel_10", panel_10_vids);
    add_panel_mapping("panel_11", panel_11_vids);
    add_panel_mapping("panel_12", panel_12_vids);
    add_panel_mapping("panel_13", panel_13_vids);
    add_panel_mapping("panel_14", panel_14_vids);
    add_panel_mapping("panel_15", panel_15_vids);
    add_panel_mapping("panel_16", panel_16_vids);
    add_panel_mapping("panel_17", panel_17_vids);
    add_panel_mapping("panel_18", panel_18_vids);
    add_panel_mapping("panel_19", panel_19_vids);
    add_panel_mapping("panel_20", panel_20_vids);
    add_panel_mapping("panel_21", panel_21_vids);
    add_panel_mapping("panel_57", panel_57_vids);
    add_panel_mapping("panel_58", panel_58_vids);
    add_panel_mapping("panel_59", panel_59_vids);
    add_panel_mapping("panel_60", panel_60_vids);

    std::map<std::string, std::vector<int>> panel_vids = {
	{"panel_2a", panel_2a_vids},
        {"panel_2b", panel_2b_vids},
        {"panel_3",  panel_3_vids},
	{"panel_4", panel_4_vids},
	{"panel_5a", panel_5a_vids},
	{"panel_5b", panel_5b_vids},
	{"panel_6", panel_6_vids},
	{"panel_7", panel_7_vids},
	{"panel_8", panel_8_vids},
	{"panel_9", panel_9_vids},
	{"panel_10", panel_10_vids},
	{"panel_11", panel_11_vids},
	{"panel_12", panel_12_vids},
	{"panel_13", panel_13_vids},
	{"panel_14", panel_14_vids},
	{"panel_15", panel_15_vids}, 
	{"panel_16", panel_16_vids},
	{"panel_17", panel_17_vids},
	{"panel_18", panel_18_vids},
	{"panel_19", panel_19_vids},
	{"panel_20", panel_20_vids},
	{"panel_21", panel_21_vids}, 
	{"panel_57", panel_57_vids},
	{"panel_58", panel_58_vids}, 
	{"panel_59", panel_59_vids},
	{"panel_60", panel_60_vids}
    };

    std::map<std::string, std::vector<double>*> panel_offsets = {
        {"panel_1", &panel_1_offsets},
        {"panel_2a", &panel_2a_offsets},
        {"panel_2b", &panel_2b_offsets},
        {"panel_3",  &panel_3_offsets},
        {"panel_4", &panel_4_offsets},
        {"panel_5a", &panel_5a_offsets},
        {"panel_5b", &panel_5b_offsets},
        {"panel_6", &panel_6_offsets},
        {"panel_7", &panel_7_offsets},
        {"panel_8", &panel_8_offsets},
        {"panel_9", &panel_9_offsets},
        {"panel_10", &panel_10_offsets},
        {"panel_11", &panel_11_offsets},
        {"panel_12", &panel_12_offsets},
        {"panel_13", &panel_13_offsets},
        {"panel_14", &panel_14_offsets},
        {"panel_15", &panel_15_offsets},
        {"panel_16", &panel_16_offsets},
        {"panel_17", &panel_17_offsets},
        {"panel_18", &panel_18_offsets},
        {"panel_19", &panel_19_offsets},
        {"panel_20", &panel_20_offsets},
        {"panel_21", &panel_21_offsets},
        {"panel_57", &panel_57_offsets},
        {"panel_58", &panel_58_offsets},
        {"panel_59", &panel_59_offsets},
        {"panel_60", &panel_60_offsets}
    };
    std::map<std::string, TH1D*> hists_offsets;

    for (auto &kv : volid_lookup) {
        const std::string &panel_name = kv.second.panel;
        if (panel_name == "panel_1") continue;

        // Only create the histogram if it doesn't already exist
        if (hists_offsets.find(panel_name) == hists_offsets.end()) {
            hists_offsets[panel_name] = new TH1D(
                Form("offset_%s", panel_name.c_str()), 
                Form("Panel offset for %s vs panel_1; Δs - Δt (mm/ns); Counts", panel_name.c_str()), 
                200, -500, 500
        );
    }
}
    TChain* Instrument_Events = new TChain("TreeRec");
    Instrument_Events->SetAutoDelete(true);
    CEventRec* Event = new CEventRec;
    Instrument_Events->SetBranchAddress("Rec", &Event);
    Instrument_Events->Add(data_path.c_str()); 

    // initialize volume id and n_hits, initialize progress bar
    int vol_id = 0, n_relevant_hits = 0;
    progressbar progress(Instrument_Events->GetEntries() / 1000);

    Instrument_Events->SetBranchAddress("Rec", &Event);

    //begin loop
    for (size_t i = 0; i < Instrument_Events->GetEntries(); i++) {
        Instrument_Events->GetEntry(i);
        if(i%1000==0){
                progress.update();
        }
        // single track requirement 
            if (Event->GetNTracks() != 1) continue;

        // requiring all TOF hits to be on the same track
        vector<int> tof_track_indices = Event->GetHitTrackIndex();
        bool skip_event = false;
        int first_idx = -2;
        for (size_t j = 0; j < tof_track_indices.size(); j++) {
                int idx = tof_track_indices[j];
                if (idx == -1) {
                        skip_event = true;
                    break;
                }
                if (j == 0) {
                    first_idx = idx;
                } else if (idx != first_idx) {
                    skip_event = true;
                        break;
                }
        }
        if (skip_event) continue;

        // requiring at least one hit on outer tof and one hit on inner tof
        // and also getting the times and positions since it requires opening the event to see the volume id anyway
        //
        std::map<std::string, HitInfo> hit_infos; // re-initialized in every event
        bool is_outer_tof = false;
        bool is_inner_tof = false;
        int n_relevant_hits = 0;

        for (const auto &hit : Event->GetHitSeries()) {
            int vol_id = hit.GetVolumeId();
            if (!GGeometryObject::IsTofVolume(vol_id)) continue;
            
            if (GGeometryObject::IsUmbrellaVolume(vol_id)) is_outer_tof = true; //umb + cortina
            if (GGeometryObject::IsCubeVolume(vol_id)) is_inner_tof = true; //cube
        
            auto it = volid_lookup.find(vol_id);
            if (it != volid_lookup.end()) {
                const std::string &panel_name = it->second.panel;
                size_t paddle_offset_index = it->second.index;

                double raw_time = hit.GetTime();
                double paddle_offset = panel_offsets[panel_name]->at(paddle_offset_index);
                
                double adj_time = raw_time - paddle_offset;
                TVector3 pos = hit.GetPosition();
                
                hit_infos[panel_name] = {adj_time, pos};
                    n_relevant_hits++;
        
            }	    
        }

        if (n_relevant_hits < 2) continue; // checking if there are at least 2 relevant hits to be consdiered for the analysis
        if (!(is_outer_tof && is_inner_tof)) continue; // checking if track has one hit on inner tof + one hit on outer tof for proper reconstruction

        if (hit_infos.find("panel_1") == hit_infos.end()) continue; // check if one of the hits is on panel 1

        double t_panel1 = hit_infos["panel_1"].adj_time;
        TVector3 pos_panel1 = hit_infos["panel_1"].pos;

        const double c_mm_per_ns = 299.705; // calculated @ McMurdo with a refractive index of 1.000305 

        for (const auto &kv : hit_infos) {
            if (kv.first == "panel_1") continue; //don't compare panel 1 to itself

            double t_other = kv.second.adj_time;
            TVector3 pos_other = kv.second.pos;

            double delta_t = std::abs(t_other - t_panel1); 
            if (delta_t == 0) continue; //avoid seg-fault from somehow dividing by 0

            TVector3 diff = pos_other - pos_panel1;
            double distance = diff.Mag();

            double inter_panel_offset = (distance/c_mm_per_ns) - delta_t;

            auto it = hists_offsets.find(kv.first);
            if (it != hists_offsets.end()) {
                it->second->Fill(inter_panel_offset);
            }
        }
    }

    std::map<std::string, double> mean_panel_offsets;

    std::ofstream outfile("panel_offsets.csv");
    if (!outfile.is_open()) {
        std::cerr << "Error: could not open panel_offsets.csv for writing" << std::endl;
    } else {
        outfile << "Panel,Offset\n"; // CSV header

        for (auto &kv : hists_offsets) {
            const std::string &panel = kv.first;
            TH1D *hist = kv.second;
            double mean_offset = hist->GetMean();
            mean_panel_offsets[panel] = mean_offset;

            // print to terminal
            std::cout << "Panel offset for " << panel
                    << " relative to panel_1 = "
                    << mean_offset << std::endl;

            // write to CSV
            outfile << panel << "," << mean_offset << "\n";
        }

        outfile.close();
    }
}
        
        


