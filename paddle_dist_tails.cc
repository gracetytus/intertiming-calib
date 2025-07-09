#include <TFile.h>
#include <TTree.h>
#include <TChain.h>
#include <TH1D.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TStyle.h>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
#include <numeric>
#include <boost/format.hpp>

#include "GOptionParser.hh"
#include "CEventRec.hh"
#include "GRecoHit.hh"
#include "GGeometryObject.hh"
#include "progressbar.hpp"

using std::vector;
using std::string;
using std::cout;
using std::endl;
using boost::format;

vector<int> parseCommaSeparatedInts(const string& input) {
    vector<int> result;
    std::stringstream ss(input);
    string token;
    while (std::getline(ss, token, ',')) {
        result.push_back(std::stoi(token));
    }
    return result;
}

vector<float> parseCommaSeparatedFloats(const string& input) {
    vector<float> result;
    std::stringstream ss(input);
    string token;
    while (std::getline(ss, token, ',')) {
        result.push_back(std::stof(token));
    }
    return result;
}

int main(int argc, char* argv[]) {
    gStyle->SetOptStat(0);

    GOptionParser* parser = GOptionParser::GetInstance();
    parser->AddProgramDescription("Counts leading paddle IDs when t_diff >= 0.8");
    parser->AddCommandLineOption<string>("rec_path", "Path to input .root files", "./*", "i");
    parser->AddCommandLineOption<string>("paddle_ids", "Comma-separated paddle volume ID suffixes", "4,5,6", "b");
    parser->AddCommandLineOption<int>("vol_id_base", "Base of paddle volume ID", 0, "s");
    parser->AddCommandLineOption<string>("offsets", "Comma-separated time offsets for paddles", "0,0,0", "f");
    parser->AddCommandLineOption<string>("out_file", "Output file name", "paddle_id_counts.root", "o");
    parser->AddCommandLineOption<string>("pid", "Panel ID prefix", "p0", "p");
    parser->ParseCommandLine(argc, argv);
    parser->Parse();

    string data_path = parser->GetOption<string>("rec_path");
    string panel_id = parser->GetOption<string>("pid");
    string out_file_name = panel_id + "_" + parser->GetOption<string>("out_file");

    vector<int> paddle_ids_suffix = parseCommaSeparatedInts(parser->GetOption<string>("paddle_ids"));
    vector<float> offsets = parseCommaSeparatedFloats(parser->GetOption<string>("offsets"));
    int vol_id_base = parser->GetOption<int>("vol_id_base");

    vector<int> paddle_ids(paddle_ids_suffix.size());
    for (size_t i = 0; i < paddle_ids_suffix.size(); i++) {
        paddle_ids[i] = vol_id_base + paddle_ids_suffix[i];
    }

    vector<double> paddle_times(paddle_ids.size(), -1);
    vector<double> paddle_times_raw(paddle_ids.size(), -1);

    TChain* chain = new TChain("TreeRec");
    CEventRec* event = new CEventRec;
    chain->SetBranchAddress("Rec", &event);
    chain->Add(data_path.c_str());

    // Create histogram: bin for each paddle ID
    int min_id = *std::min_element(paddle_ids.begin(), paddle_ids.end());
    int max_id = *std::max_element(paddle_ids.begin(), paddle_ids.end());
    TH1D* paddle_id_hist = new TH1D("leading_paddle_ids", "Leading Paddle IDs for t_{diff} >= 0.8;Paddle ID;Count", 
                                    max_id - min_id + 1, min_id - 0.5, max_id + 0.5);

    progressbar bar(chain->GetEntries() / 1000);

    for (Long64_t i = 0; i < chain->GetEntries(); ++i) {
        chain->GetEntry(i);
        if (i % 1000 == 0) bar.update();

        if (event->GetNTracks() != 1) continue;

        std::fill(paddle_times.begin(), paddle_times.end(), -1);
        std::fill(paddle_times_raw.begin(), paddle_times_raw.end(), -1);

        for (const GRecoHit& hit : event->GetHitSeries()) {
            int vol_id = hit.GetVolumeId();
            if (GGeometryObject::IsTofVolume(vol_id)) {
                for (size_t j = 0; j < paddle_ids.size(); ++j) {
                    if (vol_id == paddle_ids[j]) {
                        paddle_times_raw[j] = hit.GetTime();
                    }
                }
            }
        }

        for (size_t j = 0; j < paddle_ids.size(); ++j) {
            if (paddle_times_raw[j] > 0) {
                paddle_times[j] = paddle_times_raw[j] - offsets[j];
            }
        }

        for (size_t j = 0; j < paddle_ids.size() - 1; ++j) {
            if (paddle_times[j] > 0 && paddle_times[j + 1] > 0) {
                double t_diff = paddle_times[j] - paddle_times[j + 1];
                if (t_diff >= 0.8) {
                    paddle_id_hist->Fill(paddle_ids[j]);
                }
            }
        }
    }

    // Output
    TFile out_file(out_file_name.c_str(), "RECREATE");
    out_file.cd();

    TCanvas* canvas = new TCanvas("paddle_id_canvas", "paddle_id_canvas", 800, 600);
    paddle_id_hist->SetLineColor(kBlue + 2);
    paddle_id_hist->Draw("HIST");
    canvas->SaveAs((panel_id + "_leading_paddle_ids.pdf").c_str());
    canvas->Write(panel_id + "_paddle_id_canvas");
    paddle_id_hist->Write(panel_id + "_leading_paddle_ids");
    out_file.Close();

    return 0;
}

