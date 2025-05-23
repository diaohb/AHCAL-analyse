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
#include "TRandom3.h"
#include <TSpectrum.h>
#include "TLegend.h"
#include "TLine.h"
#include "TVirtualFitter.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TGaxis.h"
#include <thread>
#include <mutex>
#include <unordered_map>
using namespace std;
void fit(TH1D *his, Double_t *fitrange, Double_t *startvalues, Double_t *parlimitslo, Double_t *parlimitshi, Double_t *fitparams, Double_t *fiterrors, Double_t *ChiSqr, Int_t *NDF);
int main(int argc, char *argv[])
{
    double start = clock();
    raw2Root tw;
    tw.MIPlist(argv[1], argv[2]);
    double end = clock();
    cout << "end of mip : Time : " << (end - start) / CLOCKS_PER_SEC << endl;
    return 0;
}
int raw2Root::MIPlist(const string _list, string pedestal) {
    ReadList(_list);
    TFile *fin, *fout1, *fout2;
    TTree *tin, *tout1, *tout2;
    gRandom = new TRandom3(0);
    fin = TFile::Open(pedestal.c_str(), "READ");
    if (!fin) {
        cout << "cant open " << pedestal << endl;
        return 0;
    }
    tin = (TTree *) fin->Get("pedestal");
    if (!tin) {
        cout << "cant get tree pedestal" << endl;
        return 0;
    }
    int CellID = 0;
    double pedestal_high = 0, pedestal_low = 0;
    tin->SetBranchAddress("cellid", &CellID);
    tin->SetBranchAddress("highgain_peak", &pedestal_high);
    unordered_map<int, double> ped_high;
    for (int i = 0; i < tin->GetEntries(); ++i)
    {
        tin->GetEntry(i);
        ped_high[CellID] = pedestal_high;
    }

    fout1 = TFile::Open(TString(_list) + "_mip1.root", "recreate");
    fout2 = TFile::Open(TString(_list) + "_mip2.root", "recreate");
    unordered_map<int, TH1D *> mmip1;
    unordered_map<int, TH1D *> mmip_chip1;
    unordered_map<int, TH1D *> mmip2;
    unordered_map<int, TH1D *> mmip_chip2;
    for (int layer = 0; layer < 40; layer++) {
        TString slayer = "layer" + TString(to_string(layer).c_str());
        for (int chip = 0; chip < 9; chip++) {
            TString schip = "chip" + TString(to_string(chip).c_str());
            TString name_chip = "MIP Spectrum " + slayer + " " + schip;
            mmip_chip1[layer * 10 + chip] = new TH1D(name_chip, name_chip, 500, 0, 2000);
            mmip_chip2[layer * 10 + chip] = new TH1D(name_chip, name_chip, 500, 0, 2000);
            for (int channel = 0; channel < 36; channel++) {
                TString schannel = "channel" + TString(to_string(channel).c_str());
                TString name = "MIP Spectrum " + slayer + " " + schip + " " + schannel;
                int cellid = layer * 1e5 + chip * 1e4 + channel;
                mmip1[cellid] = new TH1D(name, name, 500, 0, 2000);
                mmip2[cellid] = new TH1D(name, name, 500, 0, 2000);
            }
        }
    }
    double MPV = 0, width = 0, gaus_sigma = 0, max_x = 0, FWHM = 0, chi2_ndf = 0;
    int cellid = 0, entries = 0;
    int tag = 0;
    for_each(list.begin(), list.end(), [&](string tmp) {
        cout<<"Reading: "<<tmp<<endl;
        fin=TFile::Open(TString(tmp),"read");
        tin=(TTree*)fin->Get("Raw_Hit");
        ReadTreeBranch(tin);
        for (int n = 0; n < tin->GetEntries(); n++)
        {
            tin->GetEntry(n);
            if (MIP(cellID, hitTag, HG_Charge)) {
                int l = cellID->size();
                for (int i = 0; i < l; i++)
                {
                    int cid = cellID->at(i);
                    if (gRandom->Uniform() < 0.5) {
                        mmip_chip1[cid / 10000]->Fill(HG_Charge->at(i) - ped_high[cid]);
                        mmip1[cid]->Fill(HG_Charge->at(i) - ped_high[cid]);
                    } else {
                        mmip_chip2[cid / 10000]->Fill(HG_Charge->at(i) - ped_high[cid]);
                        mmip2[cid]->Fill(HG_Charge->at(i) - ped_high[cid]);
                    }
                }
            }
        }
        fin->Close(); });
    TString dir = "histogram";
    double fr[2];
    double sv[4], pllo[4], plhi[4], fps[4], fpe[4];
    double chisqr;
    int ndf;
    int cut = 0;
    // mutex g_mutex;
    auto ffit = [&](int layer, unordered_map<int, TH1D *> mmip_chip, unordered_map<int, TH1D *> mmip, TTree *tout, TFile *fout) {
        TString slayer = "layer" + TString(to_string(layer).c_str());
        cout << "fitting " << slayer << " ..." << endl;
        fout->mkdir(dir + "/" + slayer);
        for (int chip = 0; chip < 9; chip++) {
            TString schip = "chip" + TString(to_string(chip).c_str());
            fout->mkdir(dir + "/" + slayer + "/" + schip);
            fout->cd(dir + "/" + slayer + "/" + schip);
            fr[0] = 150;
            fr[1] = 900;
            sv[0] = 40;
            sv[1] = 344;
            sv[2] = mmip_chip[layer * 10 + chip]->Integral(50, 500);
            sv[3] = 80;
            pllo[0] = 10;
            pllo[1] = 100;
            pllo[2] = sv[2] / 50;
            pllo[3] = 10;
            plhi[0] = 100;
            plhi[1] = 700;
            plhi[2] = sv[2] * 10;
            plhi[3] = 200;

            if (sv[2] >= 1000) {
                for (cut = 75; cut > 0; cut--) {
                    if (mmip_chip[layer * 10 + chip]->GetBinContent(cut) < sv[2] / 400) {
                        break;
                    }
                }
                // cut += 5;
                fr[0] = mmip_chip[layer * 10 + chip]->GetBinCenter(cut);
                for (cut = 100; cut < 250; cut++) {
                    if ((mmip_chip[layer * 10 + chip]->GetBinContent(cut) + mmip_chip[layer * 10 + chip]->GetBinContent(cut + 1)) < sv[2] / 500) {
                        break;
                    }
                }
                cut -= 10;
                fr[1] = mmip_chip[layer * 10 + chip]->GetBinCenter(cut);
                // g_mutex.lock();
                langaufit(mmip_chip[layer * 10 + chip], fr, sv, pllo, plhi, fps, fpe, &chisqr, &ndf);
                // g_mutex.unlock();
            }
            sv[0] = fps[0];
            sv[1] = fps[1];
            sv[3] = fps[3];
            mmip_chip[layer * 10 + chip]->Write();
            for (int channel = 0; channel < 36; channel++) {
                cellid = layer * 1e5 + chip * 1e4 + channel;
                entries = mmip[cellid]->Integral(50, 500);
                sv[2] = entries;
                pllo[2] = sv[2] / 50;
                plhi[2] = sv[2] * 10;
                if (sv[2] < 500) {
                    tag = -2;
                    tout->Fill();
                    mmip[cellid]->Write();
                    continue;
                }
                // for (cut = 75; cut > 25; cut--)
                // {
                //     if (mmip1[cellid]->GetBinContent(cut) < sv[2] / 400)
                //     {
                //         break;
                //     }
                // }
                // cut += 5;
                // fr[0] = mmip1[cellid]->GetBinCenter(cut);
                for (cut = 100; cut < 250; cut++) {
                    if ((mmip[cellid]->GetBinContent(cut) + mmip[cellid]->GetBinContent(cut + 1)) < sv[2] / 500) {
                        break;
                    }
                }
                cut -= 10;
                fr[1] = mmip[cellid]->GetBinCenter(cut);
                // g_mutex.lock();
                langaufit(mmip[cellid], fr, sv, pllo, plhi, fps, fpe, &chisqr, &ndf);
                // g_mutex.unlock();
                langaupro(fps, max_x, FWHM);
                double x = mmip[cellid]->FindBin(max_x);
                mmip[cellid]->Write();
                MPV = fps[1];
                width = fps[0];
                gaus_sigma = fps[3];
                tag = 1;
                chi2_ndf = chisqr / double(ndf);
                if (chi2_ndf > 30 || x > fr[1] - 40 || gaus_sigma < 30 || gaus_sigma > 160 || width < 15 || width > 80) {
                    tag = -3;
                }
                tout->Fill();
            }
        }
        cout << "fitting " << slayer << " done" << endl;
    };
    auto write = [&](TFile *fout, TTree *tout, unordered_map<int, TH1D *> mmip_chip, unordered_map<int, TH1D *> mmip) {
        tout = new TTree("MIP_Calibration", "MIP_Calibration");
        tout->Branch("MPV", &MPV);
        tout->Branch("width", &width);
        tout->Branch("gaus_sigma", &gaus_sigma);
        tout->Branch("CellID", &cellid);
        tout->Branch("Entries", &entries);
        tout->Branch("max_x", &max_x);
        tout->Branch("FWHM", &FWHM);
        tout->Branch("chi2_ndf", &chi2_ndf);
        tout->Branch("Tag", &tag);
        fout->mkdir(dir);
        for (int i = 0; i < 40; i++) {
            ffit(i, mmip_chip, mmip, tout, fout);
        }
        fout->cd();
        tout->Write();
        fout->Close();
    };
    write(fout1, tout1, mmip_chip1, mmip1);
    write(fout2, tout2, mmip_chip2, mmip2);
    return 1;
}
int raw2Root::MIP(vector<int> *_cellid, vector<int> *hitTag, vector<double> *highgain) {
    for (int i = 0; i < _cellid->size(); i++) {
        int cid = _cellid->at(i);
        if (hitTag->at(i) == 0 || (hitTag->at(i) && cid / 100 % 100 != 0)) {
            _cellid->erase(_cellid->begin() + i);
            highgain->erase(highgain->begin() + i);
            hitTag->erase(hitTag->begin() + i);
            i--;
        }
    }
    int l = _cellid->size();
    // cout<<l<<endl;
    if (l < 30)
        return 0;
    double hit_x = 0, hit_y = 0;
    for (int i = 0; i < l; i++) {
        int cid = _cellid->at(i);
        hit_x += Pos_X_1(cid);
        hit_y += Pos_Y_1(cid);
    }
    hit_x /= l;
    hit_y /= l;
    for (int i = 0; i < _cellid->size(); i++)
    {
        int cid = _cellid->at(i);
        if (abs(Pos_X_1(cid) - hit_x) > 60 || abs(Pos_Y_1(cid) - hit_y) > 60)
        {
            _cellid->erase(_cellid->begin() + i);
            highgain->erase(highgain->begin() + i);
            i--;
        }
    }
    int layerlen[40] = {0};
    for (int i = 0; i < _cellid->size(); i++) {
        int layer = _cellid->at(i) / 100000;
        layerlen[layer]++;
    }
    // for (int i = 0; i < 40; i++)
    // {
    //     if (layerlen[i] > 2)
    //     {
    //         return 0;
    //     }
    // }
    int hitlayer = 0;
    for (int i = 0; i < 40; i++) {
        if (layerlen[i])
            hitlayer++;
    }
    if (hitlayer < 30)
        return 0;
    int flag03 = 0, flag358 = 0;
    int cid0 = 0, cid38 = 0;
    // int h0=0,h38=0;
    for (int i = 0; i < _cellid->size(); i++) {
        int cid = _cellid->at(i);
        if (flag03 == 0 && int(cid / 1e5) <= 2) {
            flag03++;
            cid0 = cid % 100000;
            // break;
            // h0 = highgain->at(i);
        }
        if (flag358 == 0 && int(cid / 1e5) >= 35 && int(cid / 1e5) < 38) {
            flag358++;
            cid38 = cid % 100000;
            // break;
            // h38 = highgain->at(i);
        }
    }
    if (flag03 * flag358 == 0 || abs(Pos_X_1(cid0) - Pos_X_1(cid38)) > 41 || abs(Pos_Y_1(cid0) - Pos_Y_1(cid38)) > 41)
        return -1;
    return 1;
}
