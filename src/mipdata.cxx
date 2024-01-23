#include "RawtoRoot.h"
#include "Global.h"
#include "langaus.h"
#include <fstream>
#include <iostream>
#include <stdio.h>
#include <algorithm>
#include <stdlib.h>
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include <TROOT.h>
#include <TSystem.h>
#include <TMath.h>
#include <string>
#include <TF1.h>
#include <TH1.h>
#include <TStyle.h>
#include "TMath.h"
#include "TRandom.h"
#include <TSpectrum.h>
#include "TLegend.h"
#include "TLine.h"
#include "TVirtualFitter.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TGaxis.h"
#include <unordered_map>
using namespace std;
int main(int argc,char *argv[]){
    double start = clock();
    raw2Root tw;
    tw.MIPlist(argv[1],argv[2]);
    double end = clock();
    cout<<"end of mip : Time : "<<(end-start)/CLOCKS_PER_SEC<<endl;
    return 0;
}
int raw2Root::MIPlist(const string _list,string pedestal){
    ReadList(_list);
    TFile *fin,*fout;
    TTree *tin,*tout;
    fin = TFile::Open(pedestal.c_str(),"READ");
    if (!fin){
        cout<<"cant open "<<pedestal<<endl;
        return 0;
    }
    tin = (TTree*)fin->Get("pedestal");
    if(!tin){
        cout<<"cant get tree pedestal"<<endl;
        return 0;
    }
    int CellID=0;
    double pedestal_high=0,pedestal_low=0;
    tin->SetBranchAddress("cellid",&CellID);
    tin->SetBranchAddress("highgain_peak",&pedestal_high);
    unordered_map<int, double> ped_high;
    for (int i = 0; i < tin->GetEntries(); ++i){
        tin->GetEntry(i);
        ped_high[CellID]=pedestal_high;
    }

    fout=TFile::Open(TString(_list)+"_mip.root","recreate");
    unordered_map<int, TH1D*> mmip;
    unordered_map<int, TH1D*> mmip_chip;
    for(int layer=0;layer<40;layer++){
        TString slayer="layer"+TString(to_string(layer).c_str());
        for(int chip=0;chip<9;chip++){
            TString schip="chip"+TString(to_string(chip).c_str());
            TString name_chip ="MIP Spectrum "+slayer+" "+schip;
            mmip_chip[layer*10+chip]=new TH1D(name_chip,name_chip,500,0,1000);
            for(int channel=0;channel<36;channel++){
                TString schannel="channel"+TString(to_string(channel).c_str());
                TString name ="MIP Spectrum "+slayer+" "+schip+" "+schannel;
                int cellid=layer*1e5+chip*1e4+channel;
                mmip[cellid]=new TH1D(name,name,500,0,1000);
            }
        }
    }
    tout=new TTree("mip","mip");
    double MPV=0,width=0,gaus_sigma=0;
    int cellid=0,entries=0;
    tout->Branch("MPV",&MPV);
    tout->Branch("width",&width);
    tout->Branch("gaus_sigma",&gaus_sigma);
    tout->Branch("cellid",&cellid);
    tout->Branch("entries",&entries);
    for_each(list.begin(),list.end(),[&](string tmp){
        cout<<"Reading: "<<tmp<<endl;
        fin=TFile::Open(TString(tmp),"read");
        tin=(TTree*)fin->Get("Raw_Hit");
        ReadTreeBranch(tin);
        for(int n=0;n<tin->GetEntries();n++){
            tin->GetEntry(n);
            if(MIP(cellID,hitTag)){
                for(int i=0;i<cellID->size();i++){
                    if(cellID->at(i)/100%100==0){
                        int cid=cellID->at(i);
                        mmip_chip[cid/10000]->Fill(HG_Charge->at(i)-ped_high[cid]);
                        mmip[cid]->Fill(HG_Charge->at(i)-ped_high[cid]);
                    }
                }
            }
        }    
        fin->Close();
    });
    fout->cd();
    TString dir="histogram";
    fout->mkdir(dir);
    double fr[2];
	double sv[4], pllo[4], plhi[4], fps[4], fpe[4];
	double chisqr;
	int ndf;
    for(int layer=0;layer<40;layer++){
        TString slayer="layer"+TString(to_string(layer).c_str());
        cout<<"fitting "<<slayer<<" ..."<<endl;
        fout->mkdir(dir+"/"+slayer);
        for(int chip=0;chip<9;chip++){
            TString schip="chip"+TString(to_string(chip).c_str());
            fout->mkdir(dir+"/"+slayer+"/"+schip);
            fout->cd(dir+"/"+slayer+"/"+schip);
            fr[0] = 150;
            fr[1] = 900;
            pllo[0] = 0;
            pllo[1] = 100;
            pllo[2] = 1;
            pllo[3] = 0;
            plhi[0] = 200;
            plhi[1] = 700;
            plhi[2] = mmip_chip[layer*10+chip]->GetEntries()*10;
            plhi[3] = 200;
            sv[0] = 20;
            sv[1] = 344;
            sv[2] = mmip_chip[layer*10+chip]->GetEntries();
            sv[3] = 10;
            langaufit(mmip_chip[layer*10+chip], fr, sv, pllo, plhi, fps, fpe, &chisqr, &ndf);
            sv[0] = fps[0];
			sv[1] = fps[1];
			sv[3] = fps[3];
            mmip_chip[layer*10+chip]->Write();
            for(int channel=0;channel<36;channel++){
                cellid=layer*1e5+chip*1e4+channel;
                entries=mmip[cellid]->GetEntries();
                sv[2] = mmip[cellid]->GetEntries();
                if(entries<100){
                    tout->Fill();
                    mmip[cellid]->Write();
                    continue;
                }
                langaufit(mmip[cellid], fr, sv, pllo, plhi, fps, fpe, &chisqr, &ndf);
                mmip[cellid]->Write();
                MPV=fps[1];
                width=fps[0];
                gaus_sigma=fps[3];
                tout->Fill();
            }
        }
    }
    fout->cd();
    tout->Write();
    fout->Close();
    return 1;
}
// int raw2Root::ReadList(const string _list){
//     ifstream data(_list);
// 	while(!data.eof())
// 	{
// 			string temp;
// 			data>>temp;
// 			if(temp=="")continue;
//             if (data.fail())break;
// 			list.push_back(temp);
// 	}
//     return 1;
// }
int raw2Root::MIP(vector<int> *_cellid,vector<int> *hitTag){
    int layerlen[40]={0};
    for(int i=0;i<_cellid->size();i++){
        if(hitTag->at(i)==1){
            int layer=_cellid->at(i)/100000;
            layerlen[layer]++;
        }
    }
    for(int i=0;i<40;i++){
        if(layerlen[i]>2){
            return 0;
        }
    }
    return 1;
}