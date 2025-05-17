#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include "TH1.h"
#include "TH2.h"
#include "TPaveStats.h"
#include "TProfile.h"
#include "TStyle.h"
#include "TF1.h"
#include "/cefs/higgs/diaohb/plotstyle/plotstyle.c"
using namespace std;
int maxfit()
{
    SetStyle();
    SetPrelimStyle();
    int n_layer = 6;
    int nbins[9] = {300, 300, 300, 300, 300, 300,300,300,300};
    double energy[9]={10,20,30,40,50,60,70,80,100};
    double ranges[9] = {0};
    double rangee[9] = {100, 100, 100, 150, 200, 200,250,250,250};
    // double ranges[9] = {-100,-100,-100,-100,-100,-100,-100,-100,-100};
    // double rangee[9] = { 100, 100, 100, 100, 100, 100, 100, 100, 100};
    TString filename[3] = {"/cefs/higgs/diaohb/SIM/v2code/select_run/data/origin/e-list_analyse.root",
                           "/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0115_25000126_1010_trigger100/list/calib_analyse.root",
                           "/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0115_25000126_1010_trigger100/list/truth_analyse.root"};
    TString histoname[3] = {"data","MC digi", "MC truth"};
    int color[3] = {2,1,3};
    int n_total = 9;
    TH1D *h[3][n_total];
    int entries[3][n_total];
    for (int i = 0; i < 3; i++)
    {
        for (int n = 0; n < n_total; n++){
            TString sn = TString::Format("%d GeV  ", int(energy[n])) + histoname[i];
            h[i][n] = new TH1D(sn, sn, nbins[n],ranges[n],rangee[n]);
            entries[i][n] = 0;
        }
    }
    TFile *fin;
    TTree *tin;
    TCanvas *c = new TCanvas();
    for (int i_file = 0; i_file < 3; i_file++){
        cout << histoname[i_file] << endl;
        fin = TFile::Open(filename[i_file], "READ");
        tin = (TTree *)fin->Get("tzprofile");
        double beam_energy = 0;
        int select_flag = 0;
        vector<int> *layer = 0;
        vector<double> *layer_hitno = 0;
        vector<double> *layer_energy = 0;
        vector<double> *layer_hit_x = 0;
        vector<double> *layer_hit_y = 0;
        vector<double> *layer_max_energy = 0;
        vector<double> *layer_max_x = 0;
        vector<double> *layer_max_y = 0;
        double total_energy = 0;
        int hitno = 0;
        tin->SetBranchAddress("beam_energy", &beam_energy);
        tin->SetBranchAddress("select_flag", &select_flag);
        tin->SetBranchAddress("layer", &layer);
        tin->SetBranchAddress("layer_energy", &layer_energy);
        tin->SetBranchAddress("layer_hitno", &layer_hitno);
        tin->SetBranchAddress("layer_hit_x", &layer_hit_x);
        tin->SetBranchAddress("layer_hit_y", &layer_hit_y);
        tin->SetBranchAddress("layer_max_energy", &layer_max_energy);
        tin->SetBranchAddress("layer_max_x", &layer_max_x);
        tin->SetBranchAddress("layer_max_y", &layer_max_y);
        tin->SetBranchAddress("total_energy", &total_energy);
        tin->SetBranchAddress("hitno", &hitno);
        for (int n_entries = 0; n_entries < tin->GetEntries(); n_entries++)
        {
            tin->GetEntry(n_entries);
            for (int n_e = 0; n_e < n_total;n_e++){
                if (beam_energy == energy[n_e] && select_flag == 6)
                {
                    int a = 0;
                    for (int i = 0; i < layer->size(); i++)
                    {
                        if (layer->at(i) >= 4 && layer->at(i) < 5 && layer_hit_x->at(i) != -500 && layer_hit_y->at(i) != -500 && (abs(layer_hit_x->at(i) - 20) > 10 || abs(layer_hit_y->at(i) - 20) > 10))
                        {
                            a = 1;
                            break;
                        }
                    }
                    if (a == 1)
                        continue;
                    h[i_file][n_e]->Fill(layer_energy->at(n_layer)-layer_max_energy->at(n_layer));
                    entries[i_file][n_e]++;
                    break;
                }
            }
        }
    }
    cout << "fill done" << endl;
    gStyle->SetOptTitle(1);
    gStyle->SetFrameLineWidth(3);
    gStyle->SetPadTopMargin(0.15);
    TF1* fun=new TF1("fun","gaus",0,200);
    fun->SetNpx(2000);
    double mean = 0, sigma = 0;
    for (int n_e = 0; n_e < n_total; n_e++)
    {
        cout << n_e << endl;
        TLegend *tl = new TLegend(0.52, 0.64, 0.8, 0.84);
        TPad *pad1 = new TPad("p1", "p1", 0, 0, 1, 1);
        pad1->Draw();
        pad1->cd();
        TH1D *hframe = new TH1D("hframe", "hframe", nbins[n_e],ranges[n_e],rangee[n_e]);
        hframe->SetTitle(TString::Format("%d GeV e-", int(energy[n_e])));
        NameAxis(hframe, "layer", "Max Energy [%]");
        FormatData(hframe);
        double y_max=h[0][n_e]->GetBinContent(h[0][n_e]->GetMaximumBin());
        hframe->GetYaxis()->SetRangeUser(0, y_max*1.5);
        hframe->GetYaxis()->SetMaxDigits(3);
        hframe->DrawCopy();
        for (int i = 0; i < 3; i++)
        {
            // h[i][n_e]->SetTitle(histoname[i]);
            FormatData(h[i][n_e], color[i], 20);
            cout << entries[i][n_e] << endl;
            h[i][n_e]->Scale(entries[0][n_e] / double(entries[i][n_e]));
            mean=h[i][n_e]->GetMean(1);
            sigma=h[i][n_e]->GetRMS(1);
            fun->SetParameter(1,mean);
            fun->SetParameter(2,sigma);
            h[i][n_e]->Fit(fun,"q0","",mean - 0.8 * sigma, mean + 0.8 * sigma);
            mean = fun->GetParameter(1);
            sigma = fun->GetParameter(2);
            h[i][n_e]->Fit(fun, "q0", "", mean -  0.8 * sigma, mean + 0.8 * sigma);
            mean = fun->GetParameter(1);
            sigma = fun->GetParameter(2);
            h[i][n_e]->Fit(fun, "q0", "", mean -  0.8 * sigma, mean + 0.8 * sigma);
            mean = fun->GetParameter(1);
            sigma = fun->GetParameter(2);
            fun->SetRange(mean -  0.8 * sigma, mean + 0.8 * sigma);
            h[i][n_e]->Fit(fun, "q0", "", mean -  0.8 * sigma, mean + 0.8 * sigma);
            mean = fun->GetParameter(1);
            sigma = fun->GetParameter(2);
            h[i][n_e]->SetTitle(histoname[i]+TString::Format(": %4.1lf #pm %4.1lf",mean,sigma));
            if (histoname[i] == "data")
            {
                h[i][n_e]->Draw("same p");
                tl->AddEntry(h[i][n_e], "", "p");
            }
            else
            {
                h[i][n_e]->Draw("same hist");
                tl->AddEntry(h[i][n_e], "", "p");
            }
            fun->SetLineColor(kBlue);
            fun->DrawCopy("same");
        }
        tl->Draw();
        tl->SetFillColor(0);
        tl->SetTextSize(0.04);
        TString s_suffix = "";
        if (n_e == 0)
            s_suffix = "(";
        if (n_e == n_total - 1)
            s_suffix = ")";
        c->SaveAs("layer6-max_trigger100_2500_pos.pdf" + s_suffix);
    }
    return 1;
}
