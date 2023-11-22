#include "HBase.h"
#include <TGraph.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <vector>
int main(int argc,char* argv[]){
    HBase *hbase=new HBase();
    hbase->ReadList(argv[1]);
    TF1* f=new TF1("f","gaus",0,2000000);
    vector<int> adc;
    vector<double> mean;
    for_each(hbase->list.begin(),hbase->list.end(),[&](string tmp){
        string beamenergy = tmp;
		beamenergy = beamenergy.substr(beamenergy.find_last_of('_')+1);
		beamenergy = beamenergy.substr(0,beamenergy.find("GeV"));
        int a=stoi(beamenergy);
        if(a==500)a=0.5;
		adc.push_back(a);
		hbase->ReadTree(tmp,"EventTree");
        TH1D *h=new TH1D("h","h",2000000,0,2000000);
        for(int ientry=0;ientry<hbase->tin->GetEntries();ientry++){
            hbase->tin->GetEntry(ientry);
            h->Fill(hbase->_Digi_Energy_HCAL);
        }
        double mmean=h->GetMean();
        double rms=h->GetRMS();
        double fs=mmean-2*rms;
        if(mmean<4*rms)fs=mmean/2;
        h->Fit(f,"","",fs,mmean+2*rms);
        mmean=f->GetParameter(1);
        rms=f->GetParameter(2);
        fs=mmean-2*rms;
        if(mmean<4*rms)fs=mmean/2;
        h->Fit(f,"","",fs,mmean+2*rms);
        mean.push_back(f->GetParameter(1));
        hbase->fin->Close();
    });
    TString s=argv[1];
    TFile *fnew=TFile::Open(s+".root","RECREATE");
    TTree *tout=new TTree("tout","tout");
    int a;
    double m;
    tout->Branch("beam_energy",&a);
    tout->Branch("mean",&m);
    for(int i=0;i<adc.size();i++){
        a=adc.at(i);
        m=mean.at(i);
        tout->Fill();
    }
    tout->Write();
}