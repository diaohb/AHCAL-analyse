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
using namespace std;
// analyse list digi/mc select gaus/cb
// default digi yes cb
int main(int argc, char *argv[])
{
    HBase *hbase = new HBase();
    hbase->ReadList(argv[1]);
    TString s = argv[1];
    string mode = "mc";
    bool isselect = true;
    string func = "cb";
    cout << "\e[44mmode: " << mode << " , selection: " << isselect << " , function: " << func << "\e[0m" << endl;
    TFile *fout = TFile::Open(s + "_analyse.root", "RECREATE");
    TTree *tanalyse = new TTree("tanalyse", "tanalyse");
    double beam_energy;
    double total_energy = 0;
    int event_no = 0;
    int hitlayer = 0;
    int hitno = 0;
    int big_hitno=0;
    double zeta = 0;
    double rms = 0;
    int shower_start = 0;
    int shower_end = 0;
    double f_15 = 0; // 15th, index=14
    double f_20 = 0;
    int select_flag = 0;
    double FD=0;
    double E_Hit=0;
    double R[10] = {0};

    vector<double> *layer_energy = 0;
    vector<double> *layer_hitno = 0;
    vector<int> *layer = 0;
    vector<double> *layer_hit_x = 0;
    vector<double> *layer_hit_y = 0;
    vector<int> *cellid = 0;
    vector<double> *hit_energy = 0;
    vector<double> *layer_max_energy = 0;
    vector<double> *layer_max_x =0 ;
    vector<double> *layer_max_y =0 ;

    tanalyse->Branch("beam_energy", &beam_energy);
    tanalyse->Branch("total_energy", &total_energy);
    tanalyse->Branch("layer", &layer);
    tanalyse->Branch("event_no", &event_no);
    tanalyse->Branch("hitlayer", &hitlayer);
    tanalyse->Branch("hitno", &hitno);
    tanalyse->Branch("big_hitno", &big_hitno);
    tanalyse->Branch("zeta", &zeta);
    tanalyse->Branch("rms", &rms);
    tanalyse->Branch("shower_start", &shower_start);
    tanalyse->Branch("shower_end", &shower_end);
    tanalyse->Branch("f_15", &f_15);
    tanalyse->Branch("f_20", &f_20);
    tanalyse->Branch("layer_energy", &layer_energy);
    tanalyse->Branch("layer_hitno", &layer_hitno);
    tanalyse->Branch("layer_hit_x", &layer_hit_x);
    tanalyse->Branch("layer_hit_y", &layer_hit_y);
    tanalyse->Branch("layer_max_energy", &layer_max_energy);
    tanalyse->Branch("layer_max_x", &layer_max_x);
    tanalyse->Branch("layer_max_y", &layer_max_y);
    tanalyse->Branch("select_flag", &select_flag);
    tanalyse->Branch("cellid", &cellid);
    tanalyse->Branch("hit_energy", &hit_energy);
    tanalyse->Branch("FD", &FD);
    tanalyse->Branch("E_Hit", &E_Hit);
    tanalyse->Branch("R", &R, "R[10]/D");

    for (int l = 0; l < 40; l++)
    {
        layer->push_back(l);
    }
    for_each(hbase->list.begin(), hbase->list.end(), [&](string tmp){
        string beamenergy = tmp;
		// beamenergy = beamenergy.substr(0,beamenergy.find_last_of('/'));
		beamenergy = beamenergy.substr(0,beamenergy.find("eV")-1);
		beamenergy = beamenergy.substr(beamenergy.rfind('/')+1);
        beamenergy = beamenergy.substr(beamenergy.rfind('_')+1);
        beam_energy=stod(beamenergy);
        if(beam_energy==500)beam_energy=0.5;
        hbase->Clear();
		hbase->ReadTree(tmp,"EventTree");

        int flag[10]={0};
        for(int ientry=0;ientry<hbase->tin->GetEntries();ientry++){
            hbase->tin->GetEntry(ientry);

            event_no=ientry;
            layer_energy->clear();
            layer_hitno->clear();
            layer_hit_x->clear();
            layer_hit_y->clear();
            layer_max_energy->clear();
            layer_max_x->clear();
            layer_max_y->clear();
            cellid->clear();
            hit_energy->clear();
            f_15=0;
            f_20=0;
            // for (int i_hit = 0;i_hit<hbase->Hit_Energy->size();i_hit++){
            //     hbase->Hit_Energy->at(i_hit) += hbase->SiPM_Energy->at(i_hit)/10 * 1e6 / 3.6 /12 *0.461;
            // }

            Select select(hbase->Hit_X, hbase->Hit_Y, hbase->Hit_Z, hbase->Hit_Energy, beam_energy);
            int result = select.Result(total_energy, layer_energy, layer_hitno, hitlayer, hitno, big_hitno);

            zeta = select.GetZeta();
            rms = select.GetRMS();
            FD = select.GetFD();
            E_Hit = select.GetE_Hit();
            select.GetR_alpha(R);
            shower_start = select.GetShower_Start();
            shower_end = select.GetShower_End();
            select.GetCenter(layer_hit_x, layer_hit_y);
            select.GetHitEnergy(cellid, hit_energy);
            select.GetMax(layer_max_energy, layer_max_x, layer_max_y);
            flag[result]++;
            select_flag = result;
            if(result!=7){
                f_15=layer_energy->at(14)/total_energy;
                f_20=layer_energy->at(19)/total_energy;
            }
            tanalyse->Fill();
        }
        cout<<"----------------------------------------------\e[45m"<<beam_energy<<"GeV\e[0m----------------------------------------------------"<<endl;
        if(isselect==true)cout<<"\e[44mcorner: "<<flag[7]<<"  <0.5MIP: "<<flag[0]<<"  hitno: "<<flag[1]<<"  hitlayer: "<<flag[2]<<"  s start: "<<flag[3]<<"  lastratio: "<<flag[4]<<"  n straight: "<<flag[5]<<"  preshower: "<<flag[8]<<"  hadron: "<<flag[9]<<"  ramain: "<<flag[6]<<"  total eneries: "<<hbase->tin->GetEntries()<<"\e[0m"<<endl;
        cout<<"--------------------------------------------------------------------------------------------------------"<<endl;
        hbase->fin->Close(); });

    fout->cd();
    tanalyse->Write("", TObject::kOverwrite);
    fout->Close();
    return 1;
}