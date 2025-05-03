#include "RawtoRoot.h"
#include "TGaxis.h"
#include "TGraph.h"
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TMath.h>
#include <TROOT.h>
#include <TSpectrum.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <algorithm>
#include <iostream>
#include <stdlib.h>
#include <string>
#include <unordered_map>
#include <vector>
using namespace std;
int main(int argc, char *argv[]) {
    double start = clock();
    raw2Root tw;
    if (argc < 3) {
        tw.PEDlist(argv[1]);
    } else {
        tw.PEDlist(argv[1], argv[2]);
    }
    double end = clock();
    cout << "end of pedestal : Time : " << (end - start) / CLOCKS_PER_SEC << endl;
    return 0;
}
int raw2Root::PEDlist(const string _list, const string ishittag) {
    ReadList(_list);
    TFile *fin, *fout;
    TTree *tin, *tout;

    fout = TFile::Open(TString(_list) + "_pedestal.root", "recreate");
    unordered_map<int, TH1D *> mped_high;
    unordered_map<int, TH1D *> mped_low;
    for (int layer = 0; layer < 40; layer++) {
        for (int chip = 0; chip < 9; chip++) {
            for (int channel = 0; channel < 36; channel++) {
                int cellid = layer * 1e5 + chip * 1e4 + channel;
                TString name = "highgain_" + TString(to_string(cellid).c_str());
                mped_high[cellid] = new TH1D(name, name, 500, 0, 1000);
                name = "lowgain_" + TString(to_string(cellid).c_str());
                mped_low[cellid] = new TH1D(name, name, 500, 0, 1000);
            }
        }
    }
    tout = new TTree("pedestal", "pedestal");
    double highpeak = 0, highrms = 0, lowpeak = 0, lowrms = 0;
    int cellid = 0;
    tout->Branch("cellid", &cellid);
    tout->Branch("highgain_peak", &highpeak);
    tout->Branch("highgain_rms", &highrms);
    tout->Branch("lowgain_peak", &lowpeak);
    tout->Branch("lowgain_rms", &lowrms);
    for_each(list.begin(), list.end(), [&](string tmp) {
        cout << "Reading: " << tmp << endl;
        fin = TFile::Open(TString(tmp), "read");
        tin = (TTree *) fin->Get("Raw_Hit");
        ReadTreeBranch(tin);
        for (int n = 0; n < tin->GetEntries(); n++) {
            tin->GetEntry(n);
            for (int i = 0; i < cellID->size(); i++) {
                cellid = cellID->at(i);
                if (cellid / 100 % 100 != 0) continue;
                if (ishittag == "1" && hitTag->at(i) == 1) continue;
                mped_high[cellid]->Fill(HG_Charge->at(i));
                mped_low[cellid]->Fill(LG_Charge->at(i));
            }
        }
        fin->Close();
    });
    fout->cd();
    // TString dir="";
    // fout->mkdir(dir);
    fout->mkdir("highgain");
    fout->mkdir("lowgain");
    TF1 *ffit = new TF1("ffit", "gaus", 0, 1000);
    for (int layer = 0; layer < 40; layer++) {
        TString slayer = "layer_" + TString(to_string(layer).c_str());
        cout << "fitting " << slayer << " ..." << endl;
        fout->mkdir("highgain/" + slayer);
        fout->mkdir("lowgain/" + slayer);
        for (int chip = 0; chip < 9; chip++) {
            // TString schip="chip"+TString(to_string(chip).c_str());
            // fout->mkdir(dir+"/highgain/"+slayer+"/"+schip);
            // fout->mkdir(dir+"/lowgain/"+slayer+"/"+schip);
            for (int channel = 0; channel < 36; channel++) {
                cellid = layer * 1e5 + chip * 1e4 + channel;
                fout->cd("highgain/" + slayer);
                double mean = mped_high[cellid]->GetMean();
                double sigma = mped_high[cellid]->GetRMS();
                mped_high[cellid]->Fit(ffit, "qL", "", mean - 2 * sigma, mean + 2 * sigma);
                highpeak = ffit->GetParameter(1);
                highrms = ffit->GetParameter(2);
                mped_high[cellid]->Write();

                fout->cd("lowgain/" + slayer);
                mean = mped_low[cellid]->GetMean();
                sigma = mped_low[cellid]->GetRMS();
                mped_low[cellid]->Fit(ffit, "qL", "", mean - 2 * sigma, mean + 2 * sigma);
                lowpeak = ffit->GetParameter(1);
                lowrms = ffit->GetParameter(2);
                mped_low[cellid]->Write();

                tout->Fill();
            }
        }
    }
    fout->cd();
    tout->Write();
    fout->Close();
    return 1;
}