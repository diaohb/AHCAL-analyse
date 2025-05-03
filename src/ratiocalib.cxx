#include "TString.h"
#include <unordered_map>
#include "TH1.h"
#include "TH2.h"
#include "TFile.h"
#include "TTree.h"
#include <iostream>
#include <fstream>
#include "TF1.h"
using namespace std;
int ratiocalib(){
	double MIP_E=0.461;
	TString pedestal_dir="/cefs/higgs/diaohb/CEPC2023/PS/calibration/pedestal.root";
	std::unordered_map<int,TH2S*> hratio;
	for(int i=0;i<40;i++){
		for(int j=0;j<36;j++){
			for(int k=0;k<9;k++){
				int cellid=i*1e5+k*1e4+j;
				TString scellid=to_string(cellid).c_str();
				// hratio[cellid]=new TH2D("hdac_"+scellid,"hdac_"+scellid,1200,-50,3550,150,-50,250);
				hratio[cellid]=new TH2S("hdac_"+scellid,"hdac_"+scellid,900,100,1000,1000,2000,4000);
			}
		}
	}
	TFile* fin;
	TTree* tin;
	fin=TFile::Open(pedestal_dir,"READ");
	tin=(TTree*)fin->Get("pedestal");
	int cellId=0;
	double highp=0;
	double lowp=0;
	std::unordered_map<int,double> hpedestal; 
	std::unordered_map<int,double> lpedestal;
	tin->SetBranchAddress("cellid",&cellId);
	tin->SetBranchAddress("highgain_peak",&highp);
	tin->SetBranchAddress("lowgain_peak",&lowp);
	for(int i=0;i<tin->GetEntries();i++){
		tin->GetEntry(i);
		hpedestal[cellId]=highp;
		lpedestal[cellId]=lowp;
	}

	int   _Run_No;
	int   _cycleID;
	int   _triggerID;
	unsigned int   _Event_Time;
	vector< int > *_cellID;
	vector< int > *_bcid;
	vector< int > *_hitTag;
	vector< int > *_gainTag;
	vector< int > *_cherenkov;
	vector< double > *_HG_Charge;
	vector< double > *_LG_Charge;
	vector< double > *_Hit_Time;
	ifstream in("e-list");
	//fill histo
	while(!in.eof()){
		TString dat_temp;
		in >> dat_temp;
		std::cout << dat_temp << std::endl;
		if (in.fail())
		{
			break;
		}
		fin=TFile::Open(dat_temp,"READ");
		tin=(TTree*)fin->Get("Raw_Hit");
		_cellID=0;_hitTag=0;_HG_Charge=0;_LG_Charge=0;
		// _cellID=0;_bcid=0;_hitTag=0;_gainTag=0;_cherenkov=0;_HG_Charge=0;_LG_Charge=0;_Hit_Time=0;
		// tin->SetBranchAddress("Run_Num",&_Run_No);
		// tin->SetBranchAddress("Event_Time",&_Event_Time);
		// tin->SetBranchAddress("CycleID",&_cycleID);
		// tin->SetBranchAddress("TriggerID",&_triggerID);
		tin->SetBranchAddress("CellID",&_cellID);
		// tin->SetBranchAddress("BCID",&_bcid);
		tin->SetBranchAddress("HitTag",&_hitTag);
		// tin->SetBranchAddress("GainTag",&_gainTag);
		tin->SetBranchAddress("HG_Charge",&_HG_Charge);
		tin->SetBranchAddress("LG_Charge",&_LG_Charge);
		// tin->SetBranchAddress("Hit_Time",&_Hit_Time);
		// tin->SetBranchAddress("Cherenkov",&_cherenkov);
		
		for(int n=0;n<tin->GetEntries();n++){
			tin->GetEntry(n);
			for(int i=0;i<_cellID->size();i++){
				// for(int j=0;j<9;j++){
					int cid=_cellID->at(i);
					if(_hitTag->at(i)==1&&((cid/100)%100==0)){
						hratio[_cellID->at(i)]->Fill(_LG_Charge->at(i)-lpedestal[_cellID->at(i)],_HG_Charge->at(i)-hpedestal[_cellID->at(i)]);
					}
				// }
			}
			
	if(n%10000==0)std::cout<<n<<"  "<<tin->GetEntries()<<std::endl;
		}
	}
	std::cout<<"filling done"<<std::endl;

	TFile* fout=new TFile("e-calib_HG.root","RECREATE");
	TTree* tout=new TTree("dac","dac");
	float ratio=0;
	float plat=0;
	int Cellid=0;
	int entries=0;
	float intercept=0;
	tout->Branch("cellid",&Cellid);
	tout->Branch("slope",&ratio);
	tout->Branch("plat",&plat);
	tout->Branch("entries",&entries);
	tout->Branch("intercept",&intercept);

	TF1* f1=new TF1("f1","[0]*x+[1]",2);
	// TCanvas *c=new TCanvas("c","c",1500,1200);
	// c->Divede(3,3);
	TString histo="histogram";
	fout->mkdir(histo);
	for(int i=0;i<40;i++){
		TString slayer=histo+TString("/layer_"+to_string(i));
		fout->mkdir(slayer);
		fout->cd(slayer);
		for(int k=0;k<9;k++){
			TString schip=slayer+TString("/chip_"+to_string(k));
			fout->mkdir(schip);
			fout->cd(schip);
			for(int j=0;j<36;j++){
				Cellid=i*1e5+k*1e4+j;
				TH1D* tmp=hratio[Cellid]->ProjectionX("tmp");
				entries=tmp->Integral(1,tmp->GetNbinsX());
				if(entries<10){
					tout->Fill();
					continue;
				}
				plat=tmp->GetBinCenter(tmp->GetNbinsX());
				double n=0;
				// for(int left=tmp->GetNbinsX();left>0;left--){
				// 	n+=tmp->GetBinContent(left);
				// 	if(n>entries/200.){
				// 		plat=tmp->GetBinCenter(left);
				// 		break;
				// 	}
				// }
				// hratio[Cellid]->GetXaxis()->SetRangeUser(500,plat-500);
				// hratio[Cellid]->GetYaxis()->SetRangeUser(-20,90);
				// hratio[Cellid]->GetXaxis()->SetRangeUser(250,1000);
				hratio[Cellid]->Fit(f1,"qL","",250,1000);
				// hratio[Cellid]->GetXaxis()->UnZoom();
				// hratio[Cellid]->GetYaxis()->UnZoom();
				hratio[Cellid]->Write();
				// ratio=1/f1->GetParameter(0);
				// intercept=-ratio * f1->GetParameter(1);
				ratio=f1->GetParameter(0);
				intercept=f1->GetParameter(1);
				tout->Fill();
			}
		}
	}
	fout->cd();
	tout->Write();
	fout->Close();
	return 1;
}
