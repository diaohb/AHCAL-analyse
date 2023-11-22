#define info(x) std::cout<<"==>"<<__LINE__<<" "<<#x<<"  | "<<(x)<<" |"<<std::endl;
#include "HBase.h"
#include <TGraph.h>
#include <TH1D.h>
#include <TF1.h>
#include <TCanvas.h>
#include <iostream>
#include <TFile.h>
#include <TTree.h>
#include <vector>
#include "Global.h"
#include "Select.h"
#include <algorithm>
int main(int argc,char* argv[]){
    HBase *hbase=new HBase();
    hbase->ReadList(argv[1]);
    TString s=argv[1];
    string mode;
    if(argc==3){
        mode=argv[2];
    }
    else{
        mode="digi";
    }

    TFile *fout=TFile::Open(s+"_maxEnCell.root","RECREATE");
    TTree *tout=new TTree("tout","tout");
    double a;
    int event_no=0;
    vector<int> cellid(10,0);
    vector<int> energy_no(10,0);
    vector<double> energy(10,0);
    vector<vector<double>> neighbor(10,vector<double>(8));
    
    tout->Branch("beam_energy",&a);
    tout->Branch("event_no",&event_no);
    tout->Branch("cellid",&cellid);
    tout->Branch("energy",&energy);
    tout->Branch("neighbor",&neighbor);
    tout->Branch("energy_no",&energy_no);

    double max[10]={0};
    for_each(hbase->list.begin(),hbase->list.end(),[&](string tmp){
        string beamenergy = tmp;
        beamenergy = beamenergy.substr(0,beamenergy.find("eV")-1);
		beamenergy = beamenergy.substr(beamenergy.rfind('/')+1);
        beamenergy = beamenergy.substr(beamenergy.rfind('_')+1);
        a=stoi(beamenergy);
        if(a==500)a=0.5;
		hbase->ReadTree(tmp,"EventTree",mode);
        int flag=0;
        for(int ientry=0;ientry<hbase->tin->GetEntries();ientry++){
            hbase->tin->GetEntry(ientry);
            // info(hbase->_cellID->size());
            // info(hbase->_Digi_Hit_Energy->size());
            // info(hbase->_Hit_X->size());
            // info(hbase->_Hit_Y->size());
            // info(hbase->_Hit_Z->size());
            double totalenergy=0;
            vector<double> vtmp={1};
            vector<double>* layer_energy=&vtmp;
            int hitlayer=0;
            int hitno=0;
            Select *select =new Select(hbase->_Hit_X,hbase->_Hit_Y,hbase->_Hit_Z,hbase->_Digi_Hit_Energy,a);
            if(select->Result(totalenergy,layer_energy,hitlayer,hitno)!=6){
                flag++;
                continue;
            }
            event_no=ientry;
            fill(energy.begin(),energy.end(),0);
            fill(neighbor.begin(),neighbor.end(),vector<double>(8,0));
            for(int icell=0;icell<hbase->_Digi_Hit_Energy->size();icell++){
                auto minp=min_element(energy.begin(),energy.end());
                double tmp=hbase->_Digi_Hit_Energy->at(icell);
                if(*minp<tmp){
                    *minp=tmp;
                    double x=hbase->_Hit_X->at(icell);
                    double y=hbase->_Hit_Y->at(icell);
                    double z=hbase->_Hit_Z->at(icell);
                    int layer=z/30;
                    int chip=0,channel=0;
                    inverse(x,y,chip,channel);
                    int cid=layer*1e5+chip*1e4+channel;
                    cellid[minp-energy.begin()]=cid;
                }
            }
            // cout<<"-------"<<endl;
            for(int i=0;i<energy_no.size();i++){
                energy_no.at(i)=i+1;
            }
            double tmp=0;
            double no=0;
            for(int i=0;i<energy.size();i++){
                for(int j=i;j<energy.size();j++){
                    if(energy.at(i)<energy.at(j)){
                        tmp=energy.at(i);
                        energy.at(i)=energy.at(j);
                        energy.at(j)=tmp;
                        no=cellid[i];
                        cellid[i]=cellid[j];
                        cellid[j]=no;
                    }
                }
            }
            for(int icell=0;icell<hbase->_Hit_X->size();icell++){
                double x=hbase->_Hit_X->at(icell);
                double y=hbase->_Hit_Y->at(icell);
                double z=hbase->_Hit_Z->at(icell);
                int layer=z/30;
                int chip=0,channel=0;
                inverse(x,y,chip,channel);
                int cid=layer*1e5+chip*1e4+channel;
                for(int i=0;i<10;i++){
                    int icid=cellid[i];
                    if(layer!=icid/100000)continue;
                    double x0=Pos_X_1(icid);int deltaX=round((x-x0)/40);
                    double y0=Pos_Y_1(icid);int deltaY=round((y-y0)/40);
                    int pos=deltaX-3*deltaY+4;
                    if(pos==4)continue;
                    if(pos>=5)pos--;
                    // cout<<pos<<endl;
                    if(pos>=0&&pos<=7)neighbor[i][pos]=hbase->_Digi_Hit_Energy->at(icell);
                    // if(deltaX==1&&deltaY==1)neighbor[i][0]=hbase->_Digi_Hit_Energy->at(icell);
                    // if(deltaX==1&&deltaY==0)neighbor[i][1]=hbase->_Digi_Hit_Energy->at(icell);
                    // if(deltaX==1&&deltaY==-1)neighbor[i][2]=hbase->_Digi_Hit_Energy->at(icell);
                    // if(deltaX==0&&deltaY==-1)neighbor[i][3]=hbase->_Digi_Hit_Energy->at(icell);
                    // if(deltaX==-1&&deltaY==-1)neighbor[i][4]=hbase->_Digi_Hit_Energy->at(icell);
                    // if(deltaX==-1&&deltaY==0)neighbor[i][5]=hbase->_Digi_Hit_Energy->at(icell);
                    // if(deltaX==-1&&deltaY==1)neighbor[i][6]=hbase->_Digi_Hit_Energy->at(icell);
                    // if(deltaX==0&&deltaY==1)neighbor[i][7]=hbase->_Digi_Hit_Energy->at(icell);
                    // if((deltaX==1||deltaX==-1||deltaX==0)&&(deltaY==1||deltaY==-1||deltaY==0)&&(deltaX!=0||deltaY!=0)){
                    // // if(deltaX==0&&deltaY==0){
                    //     // cout<<deltaX<<"  "<<deltaY<<"  "<<x<<"  "<<x0<<"  "<<y<<"  "<<y0<<"  "<<cid<<"  "<<icid<<endl;
                    //     neighbor[i]+=hbase->_Digi_Hit_Energy->at(icell);
                    // }
                }
            }
            // for(int i=0;i<10;i++){
            //     cout<<energy_no[i]<<"  "<<energy[i]<<"  "<<cellid[i]<<"  "<<neighbor[i]<<endl;
            // }
            // cout<<"---------------------------------------------------------"<<endl; 
            delete select;
            tout->Fill();
        }
        cout<<"\e[44m"<<flag<<" / "<<hbase->tin->GetEntries()<<"\e[0m"<<endl;
        hbase->fin->Close();
    });
    
    tout->AutoSave();
    fout->Write();
}