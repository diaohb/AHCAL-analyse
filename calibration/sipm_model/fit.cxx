#include "TFile.h"
#include "TTree.h"
#include "TGraph.h"
#include "TF1.h"
#include "TPad.h"
#include "TAxis.h"
#include </cefs/higgs/diaohb/plotstyle/plotstyle.c>
int fit(){
    SetStyle();
    TFile* fin=TFile::Open("sipm_model_0.0267xt.root","READ");
    TTree* tin=(TTree*)fin->Get("sipm_fire");
    int n_photon=0;
    double mean=0,sigma=0;
    tin->SetBranchAddress("n_photon",&n_photon);
    tin->SetBranchAddress("mean",&mean);
    tin->SetBranchAddress("sigma",&sigma);
    TGraph* gr=new TGraph();
    for(int i=0;i<tin->GetEntries();i++){
        tin->GetEntry(i);
        // if(mean>0)
            // gr->AddPoint(i,mean);
		gr->AddPoint(i,sigma);
    }
    // TF1 *fcor=new TF1("fcor","[0]*(1-e^(-x/[1]))",0,4900);
    TF1 *fcor=new TF1("fcor","[0]*(x+[1])^0.5",0,4900);
	fcor->SetParameters(7000,7);
    gr->Fit(fcor);
    TString s;
    s.Form("Fit Function: %.2f (1-e^{-#frac{x}{%.2f}})",fcor->GetParameter(0),fcor->GetParameter(1));
    fcor->SetTitle(s);
    gr->SetTitle("Model Sampling");
    FormatData(gr,kBlack,20);
    gr->SetMarkerSize(0.1);
    NameAxis(gr,"Incident Photons","Sigma of N_fired");
	// gr->Fit(fcor);
    gr->Draw("");
    gStyle->SetOptFit(0);
    gStyle->SetOptStat(0);
	fcor->Draw("same");
    gPad->BuildLegend(0.45,0.3,0.8,0.5);
    return 1;
}
