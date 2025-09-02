#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TPad.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include <iostream>
#include <sstream>
#include <vector>
using namespace std;

int layer() {
    // double energy[6] = {0.5, 1, 2, 3, 4, 5};
    // double energy_array[9] = {10, 20, 30, 40, 50, 60, 70, 80, 100};
    // int rms_limit[6] = {100, 200, 300, 400, 400, 500};
    double energy[11] = {1, 2, 3, 4, 5, 6, 7, 8, 10, 12, 15};
    double limit_2345[11] = {7, 15, 20, 25, 35, 40, 45, 50, 60, 70, 90};
    int rms_limit[11] = {300, 600, 1000, 1200, 1500, 1700, 2000, 2200, 2500, 2700, 3000};
    const int n_total = 11;
    TString pdfname = "pi-figure/layer_profile40_PS_energy_3008_cherenkov2345_noise2.pdf";
    TCanvas *c = new TCanvas("c", "c", 4096, 4096);
    TString fname[2] = {"/home/diaohb/CEPC/AHCAL-simulation/run/0422_30080126_trigger40/list_PS_noise2/calib_analyse.root", "/home/diaohb/CEPC/ihep/PS/analyse/intercept/pi-list_analyse.root"};
    TString _fname[2] = {"/home/diaohb/CEPC/AHCAL-simulation/run/0422_30080126_trigger40/list_PS_noise2/calib_analyse.root", "/home/diaohb/CEPC/ihep/PS/analyse/intercept/pi-list_analyse.root"};
    // TString fname[2] = {"/home/diaohb/CEPC/AHCAL-simulation/run/0405_30080126_trigger40/list_PS9/calib_analyse.root", "/home/diaohb/CEPC/ihep/PS/analyse/intercept/e-list_analyse.root"};
    // TString _fname[2] = {"/home/diaohb/CEPC/AHCAL-simulation/run/0405_30080126_trigger40/list_PS9/calib_analyse.root", "/home/diaohb/CEPC/ihep/PS/analyse/intercept/e-list_analyse.root"};
    int color[2] = {4, 2};
    TString name[2] = {"sim", "data"};
    double parameters[6][2] = {{0.2, -0.3}, {0.35, -0.3}, {0.45, -0.25}, {0.5, -0.2}, {0.6, -0.2}, {0.6, -0.15}};
    TF1 *fun = new TF1("fun", "[0]+[1]*x", 0, 10);
    TFile *fin, *_fin;
    TTree *tin, *_tin;
    TH1D *hlayer[n_total][40][2];
    int n1[n_total][2];
    int n2[n_total][2];
    fill(&n1[0][0], &n1[0][0] + n_total * 2, 0);
    fill(&n2[0][0], &n2[0][0] + n_total * 2, 0);
    double temp = 0;
    for (int mode = 0; mode < 2; mode++) {
        for (int i = 0; i < 40; i++) {
            TString s = "layer_" + TString(to_string(i).c_str());
            for (int n_e = 0; n_e < n_total; n_e++) {
                TString sn = TString::Format("MC %dGeV  ", int(energy[n_e])) + s + "  " + name[mode];
                hlayer[n_e][i][mode] = new TH1D(sn, sn, 30030, 0, 1000);
            }
        }
        fin = TFile::Open(fname[mode], "READ");
        tin = (TTree *) fin->Get("tanalyse");
        double beam_energy = 0;
        vector<int> *layer = 0;
        vector<double> *layer_energy = 0;
        vector<double> *layer_hitno = 0;
        int select_flag = 0;
        vector<double> *layer_hit_x = 0;
        vector<double> *layer_hit_y = 0;
        double rms = 0;
        double fd = 0;
        double e_hit = 0;
        int cherenkov = 0;
        tin->SetBranchAddress("beam_energy", &beam_energy);
        tin->SetBranchAddress("layer", &layer);
        tin->SetBranchAddress("layer_energy", &layer_energy);
        tin->SetBranchAddress("layer_hitno", &layer_hitno);
        tin->SetBranchAddress("select_flag", &select_flag);
        tin->SetBranchAddress("layer_hit_x", &layer_hit_x);
        tin->SetBranchAddress("layer_hit_y", &layer_hit_y);
        tin->SetBranchAddress("rms", &rms);
        tin->SetBranchAddress("FD", &fd);
        tin->SetBranchAddress("E_Hit", &e_hit);
        tin->SetBranchAddress("cherenkov", &cherenkov);
        _fin = TFile::Open(_fname[mode], "READ");
        _tin = (TTree *) _fin->Get("tanalyse");
        double _beam_energy = 0;
        vector<int> *_layer = 0;
        vector<double> *_layer_energy = 0;
        vector<double> *_layer_hitno = 0;
        int _select_flag = 0;
        vector<double> *_layer_hit_x = 0;
        vector<double> *_layer_hit_y = 0;
        double _rms = 0;
        double _fd = 0;
        double _e_hit = 0;
        int _cherenkov = 0;
        _tin->SetBranchAddress("beam_energy", &_beam_energy);
        _tin->SetBranchAddress("layer", &_layer);
        _tin->SetBranchAddress("layer_energy", &_layer_energy);
        _tin->SetBranchAddress("layer_hitno", &_layer_hitno);
        _tin->SetBranchAddress("select_flag", &_select_flag);
        _tin->SetBranchAddress("layer_hit_x", &_layer_hit_x);
        _tin->SetBranchAddress("layer_hit_y", &_layer_hit_y);
        _tin->SetBranchAddress("rms", &_rms);
        _tin->SetBranchAddress("FD", &_fd);
        _tin->SetBranchAddress("E_Hit", &_e_hit);
        _tin->SetBranchAddress("cherenkov", &_cherenkov);

        for (int nentries = 0; nentries < _tin->GetEntries(); nentries++) {
            _tin->GetEntry(nentries);
            tin->GetEntry(nentries);
            for (int n_e = 0; n_e < n_total; n_e++) {
                fun->SetParameters(parameters[n_e]);
                if (_beam_energy == energy[n_e] && _select_flag == 6 && _rms < rms_limit[n_e]) {
                    if (_cherenkov == 1) continue;
                    n1[n_e][mode]++;
                    // int a = 0;
                    // for (int i = 0; i < _layer->size(); i++)
                    // {
                    //     if (_layer->at(i) >= 2 && _layer->at(i) < 3 && _layer_hit_x->at(i) != -500 && _layer_hit_y->at(i) != -500 && (abs(_layer_hit_x->at(i) + 20) > 15 || abs(_layer_hit_y->at(i) - 20) > 15))
                    //     {
                    //         a = 1;
                    //         break;
                    //     }
                    // }
                    // if (a == 1)
                    //     break;
                    // int flag = 1;
                    // for (int layer = 3; layer < 6; layer++)
                    // {
                    //     if (layer_energy->at(layer) < 2)
                    //         flag = 0;
                    // }
                    // if (flag == 0)
                    //     continue;

                    double energy2345 = 0;
                    // int flag = 1;
                    for (int layer = 2; layer <= 5; layer++) {
                        energy2345 += layer_energy->at(layer);
                    }
                    if (energy2345 > limit_2345[n_e]) continue;
                    n2[n_e][mode]++;
                    for (int i = 0; i < layer_energy->size(); i++) {
                        if (layer_energy->at(i) > 0.)
                            hlayer[n_e][layer->at(i)][mode]->Fill(layer_energy->at(i));
                    }
                    break;
                }
            }
        }
    }
    cout << "filling done" << endl;

    gStyle->SetOptTitle(1);
    gStyle->SetOptStat(0);
    for (int n_e = 0; n_e < n_total; n_e++) {
        stringstream ss;
        ss << energy[n_e];
        TString senergy = ss.str();
        cout << senergy << endl;
        c->Clear();
        c->Divide(5, 4);
        cout << n2[n_e][0] << "   " << n2[n_e][1] << endl;
        for (int i = 0; i < 20; i++) {
            c->cd(i + 1);
            TPad *pad1 = new TPad("pad1", "pad1", 0, 0, 1, 1);
            pad1->SetBottomMargin(0.3);
            pad1->Draw();
            pad1->SetFillStyle(0);
            pad1->SetLogy(1);
            // pad1->SetLogx(1);
            pad1->cd();
            double range = hlayer[n_e][i][0]->GetMean() + 5 * hlayer[n_e][i][0]->GetRMS();
            for (int mode = 0; mode < 2; mode++) {
                hlayer[n_e][i][mode]->Scale(10000 / double(n2[n_e][mode]));
                hlayer[n_e][i][mode]->Rebin((range + 2) / 4);
                hlayer[n_e][i][mode]->GetXaxis()->SetRangeUser(0., range);
                temp = hlayer[n_e][i][mode]->GetXaxis()->GetLabelSize();
                hlayer[n_e][i][mode]->GetXaxis()->SetLabelSize(0);
                hlayer[n_e][i][mode]->SetMarkerStyle(6);
                hlayer[n_e][i][mode]->SetMarkerColor(color[mode]);
                hlayer[n_e][i][mode]->SetLineColor(color[mode]);
            }
            hlayer[n_e][i][0]->SetTitle("");
            hlayer[n_e][i][1]->DrawCopy("hist p");
            hlayer[n_e][i][0]->DrawCopy("same hist");

            TPad *pad2 = new TPad("pad2", "pad2", 0, 0, 1, 1);
            pad2->SetTopMargin(0.7);
            pad2->Draw();
            pad2->SetFillStyle(0);
            // pad2->SetLogx(1);
            pad2->cd();
            hlayer[n_e][i][0]->Add(hlayer[n_e][i][1], -1);
            hlayer[n_e][i][0]->Divide(hlayer[n_e][i][1]);
            hlayer[n_e][i][0]->SetLineColor(4);
            hlayer[n_e][i][0]->GetYaxis()->SetRangeUser(-1, 1);
            hlayer[n_e][i][0]->GetXaxis()->SetTitle("energy/MeV");
            hlayer[n_e][i][0]->GetXaxis()->SetLabelSize(temp);
            hlayer[n_e][i][0]->DrawCopy("hist");
        }

        TString s_suffix = "";
        if (n_e == 0)
            s_suffix = "(";
        if (n_e == n_total - 1)
            s_suffix = ")";
        c->Print(pdfname + s_suffix, "Title:" + senergy + "GeV pi-");
    }
    return 1;
}
