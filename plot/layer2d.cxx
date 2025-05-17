#include "TString.h"
#include "TFile.h"
#include "TTree.h"
#include <vector>
#include <sstream>
#include "TH1.h"
#include "TH2.h"
#include "TCanvas.h"
#include "TStyle.h"
#include "TPad.h"
#include "TF1.h"
#include <iostream>
using namespace std;

int layer2d()
{
    // double energy_array[9] = {10, 20, 30, 40, 50, 60, 70, 80, 100};
    double energy_array[6] = {0.5, 1, 2, 3, 4, 5};
    double range[6] = {15, 20, 25, 30, 35, 40};
    int rms_limit[6]={100,200,300,400,400,500};
    int n_total = 6;
    TString pdfname="layer34.pdf";
    TCanvas *c = new TCanvas("c", "c", 4096, 3072);
    TFile *fsim = TFile::Open("/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0119_25000126_trigger40/list_PS_newmip/calib_analyse.root", "READ");
    TFile *_fsim = TFile::Open("/cefs/higgs/diaohb/SIM/cern-testbeam-simulation-for-scecal-and-ahcal/run/0119_25000126_trigger40/list_PS_newmip/calib_analyse.root", "READ");
    // TFile *fsim = TFile::Open("/cefs/higgs/diaohb/SIM/v2code/select_run/data/PS_newmip/e-list_analyse.root", "READ");
    // TFile *_fsim = TFile::Open("/cefs/higgs/diaohb/SIM/v2code/select_run/data/PS_newmip/e-list_analyse.root", "READ");
    double parameters[6][2]={{0.2,-0.3},{0.35,-0.3},{0.45,-0.25},{0.5,-0.2},{0.6,-0.2},{0.6,-0.15}};
    TF1 *fun = new TF1("fun", "[0]+[1]*x", 0, 10);
    TTree *tsim = (TTree *)fsim->Get("tanalyse");
    double sim_beam_energy = 0;
    vector<int> *sim_layer = 0;
    vector<double> *sim_layer_energy = 0;
    vector<double> *sim_layer_hitno = 0;
    int sim_flag = 0;
    vector<double> *simlayer_hit_x = 0;
    vector<double> *simlayer_hit_y = 0;
    double sim_rms=0;
    double fd=0;
    double e_hit=0;
    tsim->SetBranchAddress("beam_energy", &sim_beam_energy);
    tsim->SetBranchAddress("layer", &sim_layer);
    tsim->SetBranchAddress("layer_energy", &sim_layer_energy);
    tsim->SetBranchAddress("layer_hitno", &sim_layer_hitno);
    tsim->SetBranchAddress("select_flag", &sim_flag);
    tsim->SetBranchAddress("layer_hit_x", &simlayer_hit_x);
    tsim->SetBranchAddress("layer_hit_y", &simlayer_hit_y);
    tsim->SetBranchAddress("rms", &sim_rms);
    tsim->SetBranchAddress("FD", &fd);
    tsim->SetBranchAddress("E_Hit", &e_hit);

    TTree *_tsim = (TTree *)_fsim->Get("tanalyse");
    double _sim_beam_energy = 0;
    vector<int> *_sim_layer = 0;
    vector<double> *_sim_layer_energy = 0;
    vector<double> *_sim_layer_hitno = 0;
    int _sim_flag = 0;
    vector<double> *_simlayer_hit_x = 0;
    vector<double> *_simlayer_hit_y = 0;
    double _sim_rms=0;
    double _fd=0;
    double _e_hit=0;
    _tsim->SetBranchAddress("beam_energy", &_sim_beam_energy);
    _tsim->SetBranchAddress("layer", &_sim_layer);
    _tsim->SetBranchAddress("layer_energy", &_sim_layer_energy);
    _tsim->SetBranchAddress("layer_hitno", &_sim_layer_hitno);
    _tsim->SetBranchAddress("select_flag", &_sim_flag);
    _tsim->SetBranchAddress("layer_hit_x", &_simlayer_hit_x);
    _tsim->SetBranchAddress("layer_hit_y", &_simlayer_hit_y);
    _tsim->SetBranchAddress("rms", &_sim_rms);
    _tsim->SetBranchAddress("FD", &_fd);
    _tsim->SetBranchAddress("E_Hit", &_e_hit);

    TH2D *hlayer45[n_total];
    for (int n_e = 0; n_e < n_total; n_e++)
    {
        double energy = energy_array[n_e];
        stringstream ss;
        ss << energy;
        TString name =ss.str()+" GeV layer34";
        hlayer45[n_e] = new TH2D(name, name, 200, 0, range[n_e], 200, 0, range[n_e]);
    }
    // cout<<"test"<<endl;
    for (int nentries = 0; nentries < _tsim->GetEntries(); nentries++)
    {
        _tsim->GetEntry(nentries);
        tsim->GetEntry(nentries);
        for(int n_e=0;n_e<n_total;n_e++)
        {
            fun->SetParameters(parameters[n_e]);
            if (_sim_beam_energy == energy_array[n_e] && _sim_flag == 6 && _sim_rms < rms_limit[n_e] && fun->Eval(_e_hit) < _fd)
            {
                int a = 0;
                // for (int n_e = 0; n_e < _sim_layer->size(); n_e++)
                // {
                //     if (_sim_layer->at(n_e) >= 2 && _sim_layer->at(n_e) < 3 && _simlayer_hit_x->at(n_e) != -500 && _simlayer_hit_y->at(n_e) != -500 && (abs(_simlayer_hit_x->at(n_e) + 20) > 15 || abs(_simlayer_hit_y->at(n_e) - 20) > 15))
                //     {
                //         a = 1;
                //         break;
                //     }
                // }
                // if (a == 1)
                //     break;
                if (sim_layer_energy->at(3) != 0 && sim_layer_energy->at(4)!=0)
                    hlayer45[n_e]->Fill(sim_layer_energy->at(3), sim_layer_energy->at(4));
                // cout<<"test"<<endl;
                break;
            }
        }
    }
    cout << "filling done" << endl;

    for (int n_e = 0; n_e < n_total; n_e++)
    {
        double energy=energy_array[n_e];
        stringstream ss;
        ss<<energy;
        TString senergy = ss.str();
        cout << senergy << endl;
        cout << hlayer45[n_e]->GetEntries() << endl;
        c->Clear();
        cout << hlayer45[n_e]->GetEntries() << endl;
        gPad->SetLogz();
        hlayer45[n_e]->Draw("colz");
        cout << senergy << endl;

        TString s_suffix = "";
        if (n_e == 0)
            s_suffix = "(";
        if (n_e == n_total - 1)
            s_suffix = ")";
        c->Print(pdfname+s_suffix, "Title:" + senergy + "GeV e-");
    }
    return 1;
}
