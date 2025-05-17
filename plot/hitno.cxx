#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include "TCanvas.h"
#include <vector>
#include <iostream>
#include <sstream>
#include <fstream>
#include "TH1.h"
#include "TPaveStats.h"
#include "TStyle.h"
#include "TF1.h"
#include "/cefs/higgs/diaohb/plotstyle/plotstyle.c"
using namespace std;
int hitno()
{
    SetStyle();
    SetPrelimStyle();
    // int nbins[9] = {200, 200, 200, 200, 200, 200, 200, 200, 200};
    // double energy[9]={10,20,30,40,50,60,70,80,100};
    // double ranges[9] = {0, 0, 0, 0, 0, 0,0,0,0};
    // double rangee[9] = {300, 550, 800, 1000, 1300, 1500,2000,2500,3000};
    // double ranges[6] = {0, 40, 60, 80, 100, 100};
    // double rangee[6] = {100, 140, 160, 180, 200, 220};
    double energy[6] = {0.5, 1, 2, 3, 4, 5};
    double ranges[6] = {0, 0, 0, 0, 0, 0};
    double rangee[6] = {20, 30, 40, 60, 60, 80};
    int nbins[6] = {20, 30, 40, 60, 60, 80};
    int rms_limit[6]={100,200,300,400,400,500};
    TString pdf_name = "trigger40_PS_energy_2500_hitno_2mip_.pdf";
    string txtname = "hitno_PS_2mip.txt";
    TString filename[3] = {"/cefs/higgs/diaohb/SIM/v2code/select_run/data/PS_newmip/e-list_analyse.root",
                           "/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0119_25000126_trigger40/list_PS_newmip/calib_analyse.root",
                           "/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0119_25000126_trigger40/list_PS/truth_analyse.root"};
    // TString filename_a[3] = {"/cefs/higgs/diaohb/SIM/v2code/select_run/data/PS/e-list_analyse.root", "/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0119_25000126_trigger40/list_PS_newmip/calib_analyse.root", "/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0119_25000126_trigger40/list_PS/truth_analyse.root"};
    TString histoname[3] = {"data","MC digi", "MC truth"};
    int color[3] = {2,1,3};
    int n_total = 6;
    TH1D *hen[3][n_total];
    for (int i = 0; i < 3; i++)
    {
        for (int n = 0; n < n_total; n++)
        {
            TString sn = TString::Format("%d GeV  ", int(energy[n]));
            hen[i][n] = new TH1D(sn + "hen", sn + "hen", nbins[n], ranges[n], rangee[n]);
        }
    }
    TFile *fin,*_fin;
    TTree *tin,*_tin;
    TCanvas *c = new TCanvas();
    for (int i = 0;   i < 3; i++){
        cout << histoname[i] << endl;
        fin = TFile::Open(filename[i], "READ");
        tin = (TTree *)fin->Get("tzprofile");
        double beam_energy = 0;
        int select_flag = 0;
        vector<int> *layer = 0;
        vector<double> *layer_hitno = 0;
        vector<double> *layer_energy = 0;
        vector<double> *layer_hit_x = 0;
        vector<double> *layer_hit_y = 0;
        double total_energy = 0;
        double rms=0;
        int hitno = 0;
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
            
            for (int n_entries = 0; n_entries < tin->GetEntries(); n_entries++)
            {
                // tin->GetEntry(n_entries);
                tin->GetEntry(n_entries);
                for (int n_e = 0; n_e < n_total; n_e++)
                {
                    if (beam_energy == energy[n_e] && select_flag == 6 && rms < rms_limit[n_e])
                    {
                        int a = 0;
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
                        if (layer_energy->size() == 0)
                            break;

                        double tote = 0;
                        for (int layer = 0; layer < 20; layer++)
                        {
                            // if (layer == 6)
                            //     continue;
                            tote += layer_hitno->at(layer);
                        }
                        hen[i][n_e]->Fill(tote);
                        break;
                    }
                }
        }
    }
    ofstream result(txtname, ios::out);
    gStyle->SetOptTitle(1);
    gStyle->SetFrameLineWidth(3);
    gStyle->SetPadTopMargin(0.15);
    for (int n_e = 0; n_e < n_total; n_e++)
    {
        int n_entries = hen[0][n_e]->GetEntries();
        cout << n_entries << endl;
        TF1 *fun = new TF1("fun", "gaus", ranges[n_e], rangee[n_e]);
        TLegend *tl = new TLegend(0.62, 0.64, 0.8, 0.84);
        TPad *pad1 = new TPad("p1", "p1", 0, 0, 1, 1);
        pad1->Draw();
        pad1->cd();
        TH1D *hframe = new TH1D("hframe", "hframe", nbins[n_e], ranges[n_e], rangee[n_e]);
        hframe->SetTitle(TString::Format("%d GeV e-", int(energy[n_e])));
        NameAxis(hframe, "Hitno", "Event Number");
        FormatData(hframe);
        double ymax = hen[0][n_e]->GetBinContent(hen[0][n_e]->GetMaximumBin());
        hframe->GetYaxis()->SetRangeUser(0, 1.5 * ymax);
        hframe->GetYaxis()->SetMaxDigits(3);
        hframe->DrawCopy();
        for (int i = 0; i < 3; i++)
        {
            double mean = hen[i][n_e]->GetMean();
            double rms = hen[i][n_e]->GetRMS();
            // cout<<mean<<"  "<<rms<<endl;
            fun->SetRange(mean - 1.5 * rms, mean + 1.5 * rms);
            hen[i][n_e]->SetTitle(histoname[i]);
            FormatData(hen[i][n_e], color[i], 20);
            hen[i][n_e]->Scale(n_entries / double(hen[i][n_e]->GetEntries()));
            cout<<histoname[i]+" entries: "<<hen[i][n_e]->GetEntries()<<endl;
            hen[i][n_e]->Fit(fun, "0r", "");
            if (histoname[i] == "data")
            {
                hen[i][n_e]->Draw("sames hist p");
                tl->AddEntry(hen[i][n_e], "", "p");
            }
            else
            {
                hen[i][n_e]->Draw("sames hist");
                tl->AddEntry(hen[i][n_e], "", "l");
            }
            mean = fun->GetParameter(1);
            rms = fun->GetParameter(2);
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
