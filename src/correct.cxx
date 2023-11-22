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
using namespace std;
int n_pixel=7284;
TString MIP="/cefs/higgs/shiyk/Beam_2023/DataBase/Calib/MIP/all_V2_Fit.root";
int main(int argc,char* argv[]){
    TFile* f;
    TTree* tree_in;
    int CellID=0;
    double MPV=0;
    int layer=0,chip=0,channel=0;
    double MIP[40][9][36]={0};
    fin = TFile::Open(MIP,"READ");
    if (!fin){
        cout<<"cant open "<<MIP<<endl;
        return 0;
    }
    tree_in = (TTree*)fin->Get("MIP_Calibration");
    if(!tree_in){
        cout<<"cant get MIP tree"<<endl;
        return 0;
    }
    tree_in->SetBranchAddress("CellID",&CellID);
    tree_in->SetBranchAddress("MPV",&MPV);
    //tree_in->SetBranchAddress("ID",&CellID);
    //tree_in->SetBranchAddress("mpv",&MPV);
    for (int i = 0; i < tree_in->GetEntries(); ++i){
        tree_in->GetEntry(i);
        decode_cellid(CellID,layer,chip,channel);
        if(MPV<100)MPV=510;
        MIP[layer][chip][channel]=MPV;
    }
    tree_in->Delete();
    fin->Close();

    HBase *hbase=new HBase();
    hbase->ReadList(argv[1]);
    TF1 *fcor=new TF1("fcor","-[0]*log(1-x/[0])",0,5000);
    fcor->SetParameter(0,n_pixel);
    for_each(hbase->list.begin(),hbase->list.end(),[&](string tmp){
		hbase->ReadTree(tmp,"EventTree","digi");
        string s=tmp;
        s=s.substr(s.find_last_of('/')+1);
        TString fname="correct_"+TString(s);
        TFile *fout=TFile::Open(fname,"RECREATE");
        TTree *tout=new TTree("EventTree","EventTree");
        int   _Event_No=0;
        double _Digi_Energy_ECAL=0;
        double _Digi_Energy_HCAL=0;
	    vector< int > *_detectorID=0;
	    vector< int > *_cellID=0;
	    vector< double > *_Digi_Hit_Energy=0;
	    vector< double > *_Hit_Time=0;
	    vector< double > *_Hit_X=0;
	    vector< double > *_Hit_Y=0;
	    vector< double > *_Hit_Z=0;
        // tout->Branch("EventNum",&_Event_No);
	    // tout->Branch("DetectorID",&_detectorID);
        tout->Branch("CellID",&_cellID);
        tout->Branch("Digi_Hit_Energy",&_Digi_Hit_Energy);
	    // tout->Branch("Digi_Energy_ECAL",&_Digi_Energy_ECAL);
        tout->Branch("Digi_Energy_HCAL",&_Digi_Energy_HCAL);
        tout->Branch("Hit_Time",&_Hit_Time);
        tout->Branch("Hit_X",&_Hit_X);
        tout->Branch("Hit_Y",&_Hit_Y);
        tout->Branch("Hit_Z",&_Hit_Z);
        for(int ientry=0;ientry<hbase->tin->GetEntries();ientry++){
            hbase->tin->GetEntry(ientry);
            _Event_No=0;
            _Digi_Energy_ECAL=0;
            // _detectorID->clear();
            _cellID->clear();
            _Digi_Hit_Energy->clear();
            _Hit_Time->clear();
            _Hit_X->clear();
            _Hit_Y->clear();
            _Hit_Z->clear();
            _Digi_Energy_HCAL=0;          
            for(int icell=0;icell<hbase->_cellID->size();icell++){
                _cellID->push_back(hbase->_cellID->at(icell));
                _Hit_Time->push_back(hbase->_Hit_Time->at(icell));
                double x=hbase->_Hit_X->at(icell);
                double y=hbase->_Hit_Y->at(icell);
                _Hit_X->push_back(x);
                _Hit_Y->push_back(y);
                _Hit_Z->push_back(hbase->_Hit_Z->at(icell));
                int i=(x+360)/40;
                int j=(360-y)/40;
                double e=hbase->_Digi_Hit_Energy->at(icell)*f[i][j];
                _Digi_Hit_Energy->push_back(e);
                _Digi_Energy_HCAL+=e;
                // std::cout<<e<<" "<<_Digi_Energy_HCAL<<" "<<f[i][j]<<std::endl;
            }
            if(ientry%10000==0)std::cout<<ientry<<" / "<<hbase->tin->GetEntries()<<std::endl;
            tout->Fill();         
        }
        hbase->fin->Close();
        fout->cd();
        tout->Write();
    });
}