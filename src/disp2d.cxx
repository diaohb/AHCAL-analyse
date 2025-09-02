#include "Global.h"
#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH2D.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TTree.h"
#include "TView.h"
#include <iostream>
#include <string>
#include <vector>
using namespace std;
const int cchannel[6][6] = {{2, 1, 0, 11, 12, 15},
                            {3, 4, 5, 10, 13, 16},
                            {6, 7, 8, 9, 14, 17},
                            {29, 28, 27, 26, 21, 18},
                            {32, 31, 30, 25, 22, 19},
                            {35, 34, 33, 24, 23, 20}};
double MIP_E = 0.461;
int main(int argc, char *argv[]) {
    string filename = argv[1];
    string nentries = argv[2];
    string nlayer = argv[3];
    TApplication theApp("app", &argc, argv);
    TFile *fin = TFile::Open(TString(filename), "READ");
    if (!fin) {
        cout << "can not open " << filename << endl;
        return -1;
    }
    int n = stoi(nentries);
    int layer = stoi(nlayer);
    cout << "layer:  " << layer << endl;
    TTree *tin = (TTree *) fin->Get("EventTree");
    vector<double> *Hit_Energy = 0;
    vector<int> *cellid = 0;
    vector<double> *Hit_X = 0;
    vector<double> *Hit_Y = 0;
    vector<double> *Hit_Z = 0;
    tin->SetBranchAddress("CellID", &cellid);
    tin->SetBranchAddress("Hit_Energy", &Hit_Energy);
    tin->SetBranchAddress("Hit_X", &Hit_X);
    tin->SetBranchAddress("Hit_Y", &Hit_Y);
    tin->SetBranchAddress("Hit_Z", &Hit_Z);

    TH2D *h_display = new TH2D("display", "display", 18, -360, 360, 18, -360, 360);
    // h_display->GetXaxis()->SetRangeUser(-200,200);
    // h_display->GetYaxis()->SetRangeUser(-200,200);
    TCanvas *c = new TCanvas("c", "c", 48, 130, 1000, 723);
    gStyle->SetOptStat(0);
    bool flag = 1;
    string tmp = "";
    while (flag) {
        tin->GetEntry(n);
        h_display->SetTitle("entries_" + TString(nentries) + "  layer_" + TString(nlayer));
        h_display->Reset();
        for (int i = 0; i < cellid->size(); i++) {
            if (int(cellid->at(i) / 100000) == layer && Hit_Energy->at(i) > 0.1 * MIP_E)
                h_display->Fill(Hit_X->at(i), Hit_Y->at(i), Hit_Energy->at(i));
        }
        h_display->Draw("colztext");
        gStyle->SetPaintTextFormat("2.1f");
        c->Update();
        c->Modified();
        // while(tmp!=nentries){
        cout << "next display entry no('q' for quit):" << endl;
        cin >> tmp;
        if (tmp == "q") {
            flag = 0;
            break;
        }
        cin >> nlayer;
        // cout << tmp << "  " << nlayer << endl;
        try {
            n = stoi(tmp);
            layer = stoi(nlayer);
        } catch (std::invalid_argument) {
            cout << "wrong input!!!" << endl;
        }
        // }
        nentries = tmp;
    }
    return 1;
}
