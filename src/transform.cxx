#include "HBase_data.h"
#include "Global.h"
#include <vector>
#include <TFile.h>
#include <TTree.h>
#include <iostream>
#include <TString.h>
int main(int argc,char* argv[]){
    HBase_data *hbase_data=new HBase_data();
    hbase_data->ReadList(argv[1]);
    string fname=argv[1];
    fname=fname.substr(fname.find_last_of('/')+1);
    fname=fname.substr(0,fname.find("_"));
    TFile *fout=TFile::Open("run/transform/e-/e-_"+TString(fname)+".root","RECREATE");

    
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
    tout->Branch("EventNum",&_Event_No);
	tout->Branch("DetectorID",&_detectorID);
	tout->Branch("CellID",&_cellID);
	tout->Branch("Digi_Hit_Energy",&_Digi_Hit_Energy);
	tout->Branch("Digi_Energy_ECAL",&_Digi_Energy_ECAL);
	tout->Branch("Digi_Energy_HCAL",&_Digi_Energy_HCAL);
	tout->Branch("Hit_Time",&_Hit_Time);
	tout->Branch("Hit_X",&_Hit_X);
	tout->Branch("Hit_Y",&_Hit_Y);
	tout->Branch("Hit_Z",&_Hit_Z);
    for_each(hbase_data->list.begin(),hbase_data->list.end(),[&](string tmp){
        // string beamenergy = tmp;
		// beamenergy = beamenergy.substr(beamenergy.find_last_of('_')+1);
		// beamenergy = beamenergy.substr(0,beamenergy.find("GeV"));
        // int a=stoi(beamenergy);
       
		hbase_data->ReadTree(tmp,"Calib_Hit");
        for(int ientry=0;ientry<hbase_data->tin->GetEntries();ientry++){
            hbase_data->tin->GetEntry(ientry);
            _Event_No=ientry;
            _detectorID->clear();
            _Digi_Hit_Energy->clear();
            _cellID->clear();
            _Digi_Energy_ECAL=0;
            _Digi_Energy_HCAL=0;
            _Hit_Time->clear();
            _Hit_X->clear();
            _Hit_Y->clear();
            _Hit_Z->clear();
            double en=0;
            for(int icell=0;icell<(hbase_data->_cellID->size());icell++){
                int cellid=hbase_data->_cellID->at(icell);
                _cellID->push_back(cellid);
                // if(hbase_data->_hitTag->at(icell)==0)continue;
                en=hbase_data->_hitenergy->at(icell);
                _Digi_Energy_HCAL+=en;
                _Digi_Hit_Energy->push_back(en);
                // _Hit_Time->push_back(hbase_data->_Hit_Time->at(icell));
                double x=Pos_X_1(cellid),y=Pos_Y_1(cellid);
                // std::cout<<"-------------------"<<x<<" "<<y<<std::endl;
                _Hit_X->push_back(x);
                _Hit_Y->push_back(y);
                _Hit_Z->push_back(cellid/100000 *30);
            }
            tout->Fill();
            if(ientry%1000==0)std::cout<<ientry<<std::endl;
        }
        hbase_data->fin->Close();
    });
    // tout->SetDirectory(fout);
    fout->cd();
    // tout->Write();
    tout->AutoSave();
    fout->Write();
}