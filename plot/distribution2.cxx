#include "/home/diaohb/plotstyle/plotstyle.c"
#include "TCanvas.h"
#include "TF1.h"
#include "TFile.h"
#include "TH1.h"
#include "TPaveStats.h"
#include "TString.h"
#include "TStyle.h"
#include "TTree.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <vector>
using namespace std;
int distribution2() {
    SetStyle();
    SetPrelimStyle();
    // int nbins[9] = {200, 200, 200, 200, 200, 200, 200, 200, 200};
    // double energy[9]={10,20,30,40,50,60,70,80,100};
    // double ranges[9] = {0, 0, 0, 0, 0, 0,0,0,0};
    // double rangee[9] = {300, 550, 800, 1000, 1300, 1500,2000,2500,3000};
    // double ranges[6] = {0, 40, 60, 80, 100, 100};
    // double rangee[6] = {100, 140, 160, 180, 200, 220};
    int nbins[6] = {200, 200, 200, 200, 200, 200};
    double energy[6] = {0.5, 1, 2, 3, 4, 5};
    double ranges[6] = {0, 0, 10, 30, 40, 50};
    double rangee[6] = {25, 40, 70, 80, 120, 130};
    int rms_limit[6] = {100, 200, 300, 400, 400, 500};
    TString pdf_name = "figure/trigger40_PS_energy_3008_0.pdf";
    string txtname = "txt/energy_PS_3008_0.txt";
    TString filename[3] = {"/home/diaohb/CEPC/ihep/PS/analyse/intercept/e-list_analyse.root",
                           //    "/home/diaohb/CEPC/ihep/PS/analyse/split/list2_analyse.root",
                           "/home/diaohb/CEPC/AHCAL-simulation/run/0405_30080126_trigger40/list_PS/calib_analyse.root",
                           "/home/diaohb/CEPC/AHCAL-simulation/run/0405_30080126_trigger40/list_PS5/truth_analyse.root"};
    // TString filename_a[3] = {"/cefs/higgs/diaohb/SIM/v2code/select_run/data/PS_newmip/e-list_analyse.root",
    //                         "/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0312_30050126_trigger40/list_PS_lrms_spe_mpv/calib_analyse.root",
    //                         "/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0312_30050126_trigger40/list_PS_lrms/truth_analyse.root"};
    TString histoname[3] = {"data", "data 1", "data 2"};
    double parameters[6][2] = {{0.2, -0.3}, {0.35, -0.3}, {0.45, -0.25}, {0.5, -0.2}, {0.6, -0.2}, {0.6, -0.15}};
    TF1 *fun = new TF1("fun", "[0]+[1]*x", 0, 10);
    int color[3] = {2, 1, 3};
    const int n_total = 6;
    TH1D *hen[3][n_total];
    int n_energy[n_total] = {0};
    int _n_energy[n_total] = {0};
    for (int i = 0; i < 3; i++) {
        for (int n = 0; n < n_total; n++) {
            TString sn = TString::Format("%d GeV  ", int(energy[n]));
            hen[i][n] = new TH1D(sn + "hen", sn + "hen", nbins[n], ranges[n], rangee[n]);
        }
    }
    TFile *fin, *_fin;
    TTree *tin, *_tin;
    TCanvas *c = new TCanvas();
    for (int i = 0; i < 1; i++) {
        cout << histoname[i] << endl;
        fin = TFile::Open(filename[i], "READ");
        tin = (TTree *) fin->Get("tanalyse");
        double beam_energy = 0;
        int select_flag = 0;
        vector<int> *layer = 0;
        vector<double> *layer_hitno = 0;
        vector<double> *layer_energy = 0;
        vector<double> *layer_hit_x = 0;
        vector<double> *layer_hit_y = 0;
        double total_energy = 0;
        double rms = 0;
        int hitno = 0;
        double fd = 0;
        double e_hit = 0;
        tin->SetBranchAddress("beam_energy", &beam_energy);
        tin->SetBranchAddress("select_flag", &select_flag);
        tin->SetBranchAddress("layer", &layer);
        tin->SetBranchAddress("layer_energy", &layer_energy);
        tin->SetBranchAddress("layer_hitno", &layer_hitno);
        tin->SetBranchAddress("layer_hit_x", &layer_hit_x);
        tin->SetBranchAddress("layer_hit_y", &layer_hit_y);
        tin->SetBranchAddress("total_energy", &total_energy);
        tin->SetBranchAddress("hitno", &hitno);
        tin->SetBranchAddress("rms", &rms);
        tin->SetBranchAddress("FD", &fd);
        tin->SetBranchAddress("E_Hit", &e_hit);
        double _beam_energy = 0;
        int _select_flag = 0;
        vector<int> *_layer = 0;
        vector<double> *_layer_hitno = 0;
        vector<double> *_layer_energy = 0;
        vector<double> *_layer_hit_x = 0;
        vector<double> *_layer_hit_y = 0;
        double _total_energy = 0;
        int _hitno = 0;
        double _rms = 0;
        double _fd = 0;
        double _e_hit = 0;
        _fin = TFile::Open(filename[i], "READ");
        _tin = (TTree *) _fin->Get("tanalyse");
        _tin->SetBranchAddress("beam_energy", &_beam_energy);
        _tin->SetBranchAddress("select_flag", &_select_flag);
        _tin->SetBranchAddress("layer", &_layer);
        _tin->SetBranchAddress("layer_energy", &_layer_energy);
        _tin->SetBranchAddress("layer_hitno", &_layer_hitno);
        _tin->SetBranchAddress("layer_hit_x", &_layer_hit_x);
        _tin->SetBranchAddress("layer_hit_y", &_layer_hit_y);
        _tin->SetBranchAddress("total_energy", &_total_energy);
        _tin->SetBranchAddress("hitno", &_hitno);
        _tin->SetBranchAddress("rms", &_rms);
        _tin->SetBranchAddress("FD", &_fd);
        _tin->SetBranchAddress("E_Hit", &_e_hit);
        for (int n_entries = 0; n_entries < tin->GetEntries(); n_entries++) {
            _tin->GetEntry(n_entries);
            tin->GetEntry(n_entries);
            for (int n_e = 0; n_e < n_total; n_e++) {
                fun->SetParameters(parameters[n_e]);
                if (_beam_energy == energy[n_e] && _select_flag == 6 && _rms < rms_limit[n_e]) {
                    n_energy[n_e]++;
                    if (layer_energy->size() == 0)
                        break;
                    double tote = 0;
                    for (int layer = 0; layer < 20; layer++) {
                        tote += layer_energy->at(layer);
                    }
                    hen[0][n_e]->Fill(tote);
                    break;
                }
            }
        }
        for (int n_entries = 0; n_entries < tin->GetEntries(); n_entries++) {
            _tin->GetEntry(n_entries);
            tin->GetEntry(n_entries);
            for (int n_e = 0; n_e < n_total; n_e++) {
                fun->SetParameters(parameters[n_e]);
                if (_beam_energy == energy[n_e] && _select_flag == 6 && _rms < rms_limit[n_e]) {
                    _n_energy[n_e]++;
                    if (layer_energy->size() == 0)
                        break;
                    double tote = 0;
                    for (int layer = 0; layer < 20; layer++) {
                        tote += layer_energy->at(layer);
                    }
                    if (_n_energy[n_e   ] < n_energy[n_e] / 2) {
                        hen[1][n_e]->Fill(tote);
                    } else {
                        hen[2][n_e]->Fill(tote);
                    }
                    break;
                }
            }
        }
    }
    ofstream result(txtname, ios::out);
    gStyle->SetOptTitle(1);
    gStyle->SetFrameLineWidth(3);
    gStyle->SetPadTopMargin(0.15);
    for (int n_e = 0; n_e < n_total; n_e++) {
        int n_entries = hen[0][n_e]->GetEntries();
        cout << n_entries << endl;
        TF1 *fun = new TF1("fun", "gaus", ranges[n_e], rangee[n_e]);
        TLegend *tl = new TLegend(0.62, 0.64, 0.8, 0.84);
        TPad *pad1 = new TPad("p1", "p1", 0, 0, 1, 1);
        pad1->Draw();
        pad1->cd();
        TH1D *hframe = new TH1D("hframe", "hframe", nbins[n_e], ranges[n_e], rangee[n_e]);
        hframe->SetTitle(TString::Format("%d GeV e-", int(energy[n_e])));
        NameAxis(hframe, "Energy [MeV]", "Event Number");
        FormatData(hframe);
        double ymax = hen[0][n_e]->GetBinContent(hen[0][n_e]->GetMaximumBin());
        hframe->GetYaxis()->SetRangeUser(0, 1.5 * ymax);
        hframe->GetYaxis()->SetMaxDigits(3);
        hframe->DrawCopy();
        for (int i = 0; i < 3; i++) {
            double mean = hen[i][n_e]->GetMean();
            double rms = hen[i][n_e]->GetRMS();
            // cout<<mean<<"  "<<rms<<endl;
            fun->SetRange(mean - 1.5 * rms, mean + 1.5 * rms);
            hen[i][n_e]->SetTitle(histoname[i]);
            FormatData(hen[i][n_e], color[i], 20);
            hen[i][n_e]->Scale(n_entries / double(hen[i][n_e]->GetEntries()));
            cout << histoname[i] + " entries: " << hen[i][n_e]->GetEntries() << endl;
            auto f = hen[i][n_e]->Fit(fun, "0rS", "");
            if (histoname[i] == "data") {
                hen[i][n_e]->Draw("sames hist p");
                tl->AddEntry(hen[i][n_e], "", "p");
            } else {
                hen[i][n_e]->Draw("sames hist");
                tl->AddEntry(hen[i][n_e], "", "l");
            }
            // fun->DrawCopy("same");
            mean = f->Parameter(1);
            rms = f->Parameter(2);
            result << i << "  " << int(energy[n_e]) << "  " << mean << "  " << rms << "\n";
        }
        tl->Draw();
        tl->SetFillColor(0);
        tl->SetTextSize(0.04);
        TString s_suffix = "";
        if (n_e == 0)
            s_suffix = "(";
        if (n_e == n_total - 1)
            s_suffix = ")";
        c->SaveAs(pdf_name + s_suffix);
    }
    result.close();
    return 1;
}
