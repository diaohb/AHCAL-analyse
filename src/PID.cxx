#include "HBase.h"
#include <TGraph.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <map>
#include <unordered_map>
#include <algorithm>
#include "Select.h"
#include "PID_analyse.h"
using namespace std;
int main(int argc, char *argv[])
{
    HBase *hbase = new HBase();
    hbase->ReadList(argv[1]);
    TString s = argv[1];
    TFile *fout = TFile::Open(s + "_PID.root", "RECREATE");
    TTree *tout = new TTree("PID", "PID");
    double beam_energy = 0;
    double FD = 0;
    double E_Hit = 0;
    double RMS = 0;
    double R[10] = {0};
    double alpha[10] = {0};
    for (int i = 0; i < 9; i++)
    {
        alpha[i] = log(i + 2);
    }
    alpha[9] = log(20);
    tout->Branch("beam_energy", &beam_energy);
    tout->Branch("FD", &FD);
    tout->Branch("E_Hit", &E_Hit);
    tout->Branch("RMS", &RMS);
    tout->Branch("R", &R, "R[10]/D");
    tout->Branch("alpha", &alpha, "alpha[10]/D");
    TGraph *g = new TGraph(20);
    for_each(hbase->list.begin(), hbase->list.end(), [&](string tmp)
    {
        string beamenergy = tmp;
		beamenergy = beamenergy.substr(0,beamenergy.find("eV")-1);
		beamenergy = beamenergy.substr(beamenergy.rfind('/')+1);
        beamenergy = beamenergy.substr(beamenergy.rfind('_')+1);
        beam_energy=stod(beamenergy);
        if(beam_energy==500)beam_energy=0.5;
        hbase->ReadTree(tmp,"EventTree");
        for (int ientry = 0; ientry < hbase->tin->GetEntries(); ientry++){
            // hbase->Clear();
            hbase->tin->GetEntry(ientry);
            // if(ientry%10000==0)cout<<ientry<<endl;
            PID_analyse *pid = new PID_analyse(hbase->Hit_X, hbase->Hit_Y, hbase->Hit_Z, hbase->Hit_Energy, beam_energy);
            FD = pid->GetFD();
            E_Hit = pid->GetE_Hit();
            RMS = pid->GetRMS();
            pid->GetR_alpha(R);
            tout->Fill();
            delete pid;
        }
        hbase->fin->Close(); 
    });
    fout->cd();
    tout->Write();
    fout->Close();
    return 0;
}