#include "Global.h"
#include "TApplication.h"
#include "TAxis.h"
#include "TCanvas.h"
#include "TFile.h"
#include "TH3D.h"
#include "TStyle.h"
#include "TTree.h"
#include <iostream>
#include <string>
#include <vector>
using namespace std;
double MIP_E = 0.461;
int main(int argc, char *argv[]) {
    string filename = argv[1];
    string nentries = argv[2];
    TApplication theApp("app", &argc, argv);
    TFile *fin = TFile::Open(TString(filename), "READ");
    if (!fin) {
        cout << "can not open " << filename << endl;
        return -1;
    }
    int n = stoi(nentries);
    TTree *tin = (TTree *) fin->Get("EventTree");
    vector<double> *Hit_Energy = 0;
    vector<double> *HG_Charge = 0;
    vector<int> *HitTag = 0;
    vector<int> *CellID = 0;
    int file_flag = 0;
    if (tin) {
        file_flag = 0;
        tin->SetBranchAddress("CellID", &CellID);
        tin->SetBranchAddress("Hit_Energy", &Hit_Energy);
    } else {
        file_flag = 1;
        tin = (TTree *) fin->Get("Raw_Hit");
        tin->SetBranchAddress("CellID", &CellID);
        tin->SetBranchAddress("HG_Charge", &HG_Charge);
        tin->SetBranchAddress("HitTag", &HitTag);
    }

    TH3D *h_display = new TH3D("display", "display", 400, 0, 1200, 18, -360, 360, 18, -360, 360);
    // h_display->GetXaxis()->SetRangeUser(0, 800);
    // h_display->GetYaxis()->SetRangeUser(-120,120);
    // h_display->GetZaxis()->SetRangeUser(-120,120);
    TCanvas *c = new TCanvas("c", "c", 48, 130, 1000, 723);
    gStyle->SetOptStat(0);
    // c->SetHighLightColor(2);
    // c->Range(-1.020015,-0.8450986,1.020015,0.8450986);
    // TView *view1 = TView::CreateView(1);
    // view1->SetRange(0,0,0,5,5,4);
    // c->SetFillColor(0);
    // c->SetBorderMode(0);
    // c->SetBorderSize(2);
    c->SetTheta(40);
    c->SetPhi(16);
    // c->SetFrameBorderMode(0);
    // gStyle->SetCanvasPreferGL(1);
    // h_display->GetXaxis()->SetNdivisions(40);
    // h_display->GetYaxis()->SetNdivisions(18);
    // h_display->GetZaxis()->SetNdivisions(18);
    bool flag = 1;
    string tmp = "";
    while (flag) {
        tin->GetEntry(n);
        h_display->SetTitle("entries_" + TString(nentries));
        h_display->Reset();
        for (int i = 0; i < CellID->size(); i++) {
            // if (Hit_Energy->at(i) > 0.5 * MIP_E) {
            double e = 0, x = 0, y = 0, z = 0;
            int cellid = CellID->at(i);
            x = Pos_X_1(cellid);
            y = Pos_Y_1(cellid);
            z = cellid / 1E5 * 30;
            if (file_flag && HitTag->at(i) == 1) {
                e = HG_Charge->at(i);
            } else if (Hit_Energy->at(i) > 0.5 * MIP_E) {
                e = Hit_Energy->at(i);
                e = log(e + 1) + 5;
            }
            h_display->Fill(z, x, y, e);
            // }
        }
        // h_display->SetMinimum(20);
        h_display->Draw("box2");
        c->Update();
        c->Modified();
        do {
            cout << "next display entry no('q' for quit):" << endl;
            cin >> tmp;
            if (tmp == "q") {
                flag = 0;
                break;
            } else {
                try {
                    n = stoi(tmp);
                } catch (std::invalid_argument) {
                    cout << "wrong input!!!" << endl;
                }
            }
        } while (tmp == nentries);
        nentries = tmp;
        // cout << nentries << endl;
    }
    return 1;
}
