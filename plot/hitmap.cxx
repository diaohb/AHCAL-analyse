#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <vector>
#include <sstream>
#include <unordered_map>
#include "TH1.h"
#include "TH2.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TF1.h"
using namespace std;
int hitmap()
{
    double energy = 40;
    stringstream ss;
    ss<<energy;
    TString senergy = ss.str();
    cout<<senergy<<endl;
    TFile *fin=TFile::Open("/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0116_25000126_1010_trigger40/Calib_nobirks_non/Calib_e-_40GeV.root","READ");
    TTree* tin=(TTree*)fin->Get("EventTree");
    vector<int>* cellid;
    vector<double>* hit_x;
    vector<double>* hit_y;
    vector<double>* hit_z;
    vector<double>* hit_e;
    tin->SetBranchAddress("CellID",&cellid);
    tin->SetBranchAddress("Hit_X",&hit_x);
    tin->SetBranchAddress("Hit_Y",&hit_y);
    tin->SetBranchAddress("Hit_Z",&hit_z);
    tin->SetBranchAddress("Hit_Energy",&hit_e);
    TCanvas* c=new TCanvas("c","c",3000,2000);
    unordered_map<int,TH2I*> hitmap;
    for(int i=0;i<40;i++){
        TString sname=senergy+"GeV_layer"+TString(to_string(i).c_str());
        hitmap[i]=new TH2I(sname,sname,18,-360,360,18,-360,360);
        hitmap[i]->GetXaxis()->SetTitle("x");
        hitmap[i]->GetXaxis()->CenterTitle();
        hitmap[i]->GetYaxis()->SetTitle("y");
        hitmap[i]->GetYaxis()->CenterTitle();
    }
    for(int i=0;i<tin->GetEntries();i++){
        tin->GetEntry(i);
        for(int j=0;j<cellid->size();j++){
            int layer=hit_z->at(j)/30;
            hitmap[layer]->Fill(hit_x->at(j),hit_y->at(j));
        }
    }
    for(int i=0;i<40;i++){
        c->Clear();
        gStyle->SetOptStat(0);
        TString sname=senergy+"GeV_layer"+TString(to_string(i).c_str());
        hitmap[i]->SetMinimum(0);
        hitmap[i]->SetContour(50);
        hitmap[i]->GetYaxis()->SetTitleOffset(1);
        hitmap[i]->Draw("colz");
        c->SaveAs("hitmap/MC/"+senergy+"GeV/"+sname+".png");
    }
    return 1;
}