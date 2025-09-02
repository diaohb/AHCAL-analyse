#include "Global.h"
#include "RawtoRoot.h"
#include <TF1.h>
#include <TFile.h>
#include <TH1.h>
#include <TMath.h>
#include <TROOT.h>
#include <TSpectrum.h>
#include <TStyle.h>
#include <TSystem.h>
#include <TTree.h>
#include <iostream>
#include <stdio.h>
#include <stdlib.h>
#include <string>
#include <vector>
int layer, chip, channel;
using namespace std;
char char_tmp[200];
Int_t main(int argc, char *argv[]) {
    double start = clock();
    raw2Root tw;
    tw.EnergyCalib(argv[1], argv[2], argv[3], argv[4], argv[5]);
    double end = clock();
    cout << "end of Reconstruction : Time : " << (end - start) / CLOCKS_PER_SEC << endl;
    return 0;
}
int raw2Root::EnergyCalib(string str_dat, string str_ped, string str_dac, string str_MIP, string output_file) {
    string str_out = output_file;
    TFile *fin, *fout;
    TTree *tree_in, *tree_out;
    double ped_time[Layer_No][chip_No][channel_No];
    double ped_charge[Layer_No][chip_No][channel_No];
    double gain_ratio[Layer_No][chip_No][channel_No];
    double gain_plat[Layer_No][chip_No][channel_No];
    double gain_intercept[Layer_No][chip_No][channel_No];
    double MIP[Layer_No][chip_No][channel_No];
    double hitE;
    double Edep = 0;
    double MIP_E = 0.461;//MeV
    double SwitchPoint = 500;
    const double ref_ped_time = 390;
    const double ref_ped_charge = 384;
    const double ref_MIP = 344.3;
    const double ref_gain_plat = 2927;
    const double ref_gain_ratio = 26;
    // const int lowgain_plat = 2000;
    int CellID = 0;
    float slope = 0;
    float intercept = 0;
    double MPV = 0;
    int Tag = 0;
    float plat = 0;
    double pedestal_time = 0;
    double pedestal_charge = 0;
    TH1D *h_Edep;
    TH2D *h2_HitMap;
    sprintf(char_tmp, "Energy deposition");
    h_Edep = new TH1D(char_tmp, char_tmp, 1000, 0, 3);
    sprintf(char_tmp, "Hit Map");
    h2_HitMap = new TH2D(char_tmp, char_tmp, 18, -HBU_X * 3 / 2., HBU_X * 3 / 2., 18, -HBU_Y / 2., HBU_Y / 2.);
    //read ped
    fin = TFile::Open(str_ped.c_str(), "READ");
    if (!fin) {
        cout << "cant open " << str_ped << endl;
        return 0;
    }
    tree_in = (TTree *) fin->Get("pedestal");
    if (!tree_in) {
        cout << "cant get tree pedestal" << endl;
        return 0;
    }
    tree_in->SetBranchAddress("cellid", &CellID);
    //tree_in->SetBranchAddress("time_peak",&pedestal_time);
    //tree_in->SetBranchAddress("charge_peak",&pedestal_charge);
    tree_in->SetBranchAddress("highgain_peak", &pedestal_time);
    tree_in->SetBranchAddress("lowgain_peak", &pedestal_charge);
    for (int i_layer = 0; i_layer < Layer_No; ++i_layer) {
        for (int i_chip = 0; i_chip < chip_No; ++i_chip) {
            for (int i_chan = 0; i_chan < channel_No; ++i_chan) {
                ped_time[i_layer][i_chip][i_chan] = -1;
                ped_charge[i_layer][i_chip][i_chan] = -1;
                MIP[i_layer][i_chip][i_chan] = -1;
                gain_ratio[i_layer][i_chip][i_chan] = -1;
            }
        }
    }
    for (int i = 0; i < tree_in->GetEntries(); ++i) {
        tree_in->GetEntry(i);
        decode_cellid(CellID, layer, chip, channel);
        ped_time[layer][chip][channel] = pedestal_time;
        ped_charge[layer][chip][channel] = pedestal_charge;
    }
    tree_in->Delete();
    fin->Close();
    //read dac
    fin = TFile::Open(str_dac.c_str(), "READ");
    if (!fin) {
        cout << "cant open " << str_dac << endl;
        return 0;
    }
    //tree_in = (TTree*)fin->Get("calib");
    tree_in = (TTree *) fin->Get("dac");
    if (!tree_in) {
        cout << "cant get tree dac" << endl;
        return 0;
    }
    tree_in->SetBranchAddress("cellid", &CellID);
    tree_in->SetBranchAddress("slope", &slope);
    tree_in->SetBranchAddress("plat", &plat);
    tree_in->SetBranchAddress("intercept", &intercept);
    for (int i = 0; i < tree_in->GetEntries(); ++i) {
        tree_in->GetEntry(i);
        decode_cellid(CellID, layer, chip, channel);
        //cout<<layer<<" "<<chip<<" "<<channel<<" "<<slope<<" "<<endl;
        gain_ratio[layer][chip][channel] = slope;
        gain_plat[layer][chip][channel] = plat;
        gain_intercept[layer][chip][channel] = intercept;
        if (gain_ratio[layer][chip][channel] < 10 || gain_ratio[layer][chip][channel] > 50)
            cout << CellID << " abnormal gain ratio " << layer << " " << chip << " " << channel << " " << slope << endl;
    }
    tree_in->Delete();
    fin->Close();
    //read MIP
    fin = TFile::Open(str_MIP.c_str(), "READ");
    if (!fin) {
        cout << "cant open " << str_MIP << endl;
        return 0;
    }
    tree_in = (TTree *) fin->Get("MIP_Calibration");
    if (!tree_in) {
        cout << "cant get MIP tree" << endl;
        return 0;
    }
    tree_in->SetBranchAddress("CellID", &CellID);
    tree_in->SetBranchAddress("MPV", &MPV);
    tree_in->SetBranchAddress("Tag", &Tag);
    for (int i = 0; i < tree_in->GetEntries(); ++i) {
        tree_in->GetEntry(i);
        decode_cellid(CellID, layer, chip, channel);
        MIP[layer][chip][channel] = MPV;
        if (Tag != 1)
            MIP[layer][chip][channel] = ref_MIP;
    }
    tree_in->Delete();
    fin->Close();
    //read dat
    fin = TFile::Open(str_dat.c_str(), "READ");
    if (!fin) {
        cout << "cant open " << str_dat << endl;
        return 0;
    }
    cout << "Read TTree " << endl;
    tree_in = (TTree *) fin->Get("Raw_Hit");
    cout << "Read TTree Over" << endl;
    ReadTreeBranch(tree_in);
    for (int i_layer = 0; i_layer < Layer_No; ++i_layer) {
        for (int i_chip = 0; i_chip < chip_No; ++i_chip) {
            for (int i_chan = 0; i_chan < channel_No; ++i_chan) {
                // ped_time[i_layer][i_chip][i_chan] += gain_intercept[i_layer][i_chip][i_chan];
                // ped_charge[i_layer][i_chip][i_chan] -= gain_intercept[i_layer][i_chip][i_chan]/gain_ratio[i_layer][i_chip][i_chan];
                if (ped_time[i_layer][i_chip][i_chan] < 0) ped_time[i_layer][i_chip][i_chan] = ref_ped_time;
                if (ped_charge[i_layer][i_chip][i_chan] < 0) ped_charge[i_layer][i_chip][i_chan] = ref_ped_charge;
                if (MIP[i_layer][i_chip][i_chan] < 100 || MIP[i_layer][i_chip][i_chan] > 600) MIP[i_layer][i_chip][i_chan] = ref_MIP;
                if (gain_ratio[i_layer][i_chip][i_chan] < 20 || gain_ratio[i_layer][i_chip][i_chan] > 40) gain_ratio[i_layer][i_chip][i_chan] = ref_gain_ratio;
                if (gain_plat[i_layer][i_chip][i_chan] < 0)
                    gain_plat[i_layer][i_chip][i_chan] = ref_gain_plat;
                //cout<<i_layer<<" "<<i_chip<<" "<<i_chan<<endl;
                //cout<<ped_time[i_layer][i_chip][i_chan]<<" "<<ped_time[i_layer][i_chip][i_chan]<<" "<<MIP[i_layer][i_chip][i_chan]<<" "<<gain_ratio[i_layer][i_chip][i_chan]<<endl;
            }
        }
    }
    fout = TFile::Open(str_out.c_str(), "RECREATE");
    if (!fout) {
        cout << "cant open " << str_out << endl;
        return 0;
    }
    tree_out = new TTree("EventTree", "Hits afer energy calibration");
    SetTreeBranch(tree_out);
    for (int i = 0; i < tree_in->GetEntries(); ++i) {
        Edep = 0;
        BranchClear();
        if ((i % 1000) == 0) cout << i << " out of " << tree_in->GetEntries() << endl;
        tree_in->GetEntry(i);
        _Event_No = _triggerID;
        _Detector_ID = 1;
        for (int i_hit = 0; i_hit < cellID->size(); ++i_hit) {
            if ((hitTag->at(i_hit)) == 0) continue;
            _cellID.push_back(cellID->at(i_hit));
            decode_cellid(cellID->at(i_hit), layer, chip, channel);
            if ((HG_Charge->at(i_hit)) - ped_time[layer][chip][channel] < gain_plat[layer][chip][channel] - SwitchPoint) {
                hitE = (HG_Charge->at(i_hit) - ped_time[layer][chip][channel]) * MIP_E / MIP[layer][chip][channel];
            }
            // else
            // {
            //     if(LG_Charge->at(i_hit)<800)
            //         hitE=( LG_Charge->at(i_hit) - ped_charge[layer][chip][channel] )*gain_ratio[layer][chip][channel]*MIP_E/MIP[layer][chip][channel];
            //     else
            //     {
            //         hitE = -100;
            //         // hitE = (800 - ped_charge[layer][chip][channel]) * gain_ratio[layer][chip][channel] * MIP_E / MIP[layer][chip][channel];
            //     }
            // }
            else
                hitE = ((LG_Charge->at(i_hit) - ped_charge[layer][chip][channel]) * gain_ratio[layer][chip][channel] + gain_intercept[layer][chip][channel]) * MIP_E / MIP[layer][chip][channel];
            double x = Pos_X(channel, chip);
            double y = Pos_Y(channel, chip);
            _Hit_E.push_back(hitE);
            _Hit_X.push_back(x);
            _Hit_Y.push_back(y);
            _Hit_Z.push_back(layer * 30);
            // _Hit_Time.push_back(Hit_Time->at(i_hit));
            Edep += hitE;
            h2_HitMap->Fill(x, y);
            if (hitE > 500 * MIP_E) {
                cout << hitE / MIP_E << " high energy alert " << layer << " " << chip << " " << channel << endl;
                cout << HG_Charge->at(i_hit) << " " << MIP[layer][chip][channel] << " " << gain_ratio[layer][chip][channel] << endl;
            }
        }
        h_Edep->Fill(Edep);
        for (int i_c = 0; i_c < cherenkov->size(); ++i_c) {
            if ((cherenkov->size()) != 2) cout << "abnormal cherenkov " << i << " " << cherenkov->size() << endl;
            _cherenkov.push_back(cherenkov->at(i_c));
        }
        _Digi_Energy = Edep;
        tree_out->Fill();
    }
    fout->cd();
    tree_out->Write();
    h2_HitMap->Write();
    h_Edep->Write();
    fout->Close();
    return 1;
}
