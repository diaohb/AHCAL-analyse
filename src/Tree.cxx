#include "RawtoRoot.h"
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
#include "TRandom.h"
#include <TSpectrum.h>
#include "TLegend.h"
#include "TLine.h"
#include "TVirtualFitter.h"
#include "TGraphErrors.h"
#include "TGraph.h"
#include "TGaxis.h"
// #include "Event.h"
using namespace std;
void raw2Root::Reset()
{
    _Run_No = 0;
    _Event_Time = 0;
    _layer_id = 0;
    _triggerID = 0;
    _Event_No = 0;
    _Detector_ID = 1;
    _cycleID = 0;
    _buffer_v.clear();
    //_chip_v=0;
    _cellID.clear();
    _cherenkov.clear();
    _Hit_E.clear();
    _Hit_X.clear();
    _Hit_Y.clear();
    _Hit_Z.clear();
    _Digi_Energy = 0;
    _Hit_Time.clear();
    _time.clear();
    _temp.clear();
    cellID = 0;
    data_cellid = 0;
    bcid = 0;
    hitTag = 0;
    gainTag = 0;
    cherenkov = 0;
    HG_Charge = 0;
    LG_Charge = 0;
    Hit_Time = 0;
    Hit_E = 0;
    Hit_X = 0;
    Hit_Y = 0;
    Hit_Z = 0;
    Digi_Energy = 0;
}
void raw2Root::ReadTreeBranch(TTree *tree)
{
    Reset();
    tree->SetBranchAddress("Run_Num", &_Run_No);
    tree->SetBranchAddress("Event_Time", &_Event_Time);
    tree->SetBranchAddress("CycleID", &_cycleID);
    tree->SetBranchAddress("TriggerID", &_triggerID);
    tree->SetBranchAddress("CellID", &cellID);
    // tree ->SetBranchAddress("layerIDs",&_layer_id);
    tree->SetBranchAddress("BCID", &bcid);
    tree->SetBranchAddress("HitTag", &hitTag);
    tree->SetBranchAddress("GainTag", &gainTag);
    tree->SetBranchAddress("HG_Charge", &HG_Charge);
    tree->SetBranchAddress("LG_Charge", &LG_Charge);
    tree->SetBranchAddress("Hit_Time", &Hit_Time);
    tree->SetBranchAddress("Cherenkov", &cherenkov);
}
void raw2Root::SetTreeBranch(TTree *tree)
{
    tree->Branch("Run_Num", &_Run_No);
    tree->Branch("Event_Time", &_Event_Time);
    tree->Branch("Event_Num", &_Event_No);
    tree->Branch("DetectorID", &_Detector_ID);
    tree->Branch("CellID", &_cellID);
    tree->Branch("Energy_HCAL", &_Digi_Energy);
    tree->Branch("Hit_Energy", &_Hit_E);
    tree->Branch("Hit_X", &_Hit_X);
    tree->Branch("Hit_Y", &_Hit_Y);
    tree->Branch("Hit_Z", &_Hit_Z);
    tree->Branch("Hit_Time", &_Hit_Time);
    tree->Branch("Cherenkov", &_cherenkov);
}
void raw2Root::ReadCalibTree(TTree *tree)
{
    Reset();
    // tree ->SetBranchAddress("Run_Num",&_Run_No);
    // tree ->SetBranchAddress("Event_Time",&_Event_Time);
    // tree->SetBranchAddress("Event_Num", &_Event_No);
    // tree ->SetBranchAddress("DetectorID",&_Detector_ID);
    tree->SetBranchAddress("CellID", &cellID);
    tree->SetBranchAddress("Hit_Energy", &Hit_E);
    tree->SetBranchAddress("Hit_X", &Hit_X);
    tree->SetBranchAddress("Hit_Y", &Hit_Y);
    tree->SetBranchAddress("Hit_Z", &Hit_Z);
    // tree ->SetBranchAddress("Digi_Energy_HCAL",&Digi_Energy);
    tree->SetBranchAddress("Hit_Time", &Hit_Time);
    // tree ->SetBranchAddress("Cherenkov",&cherenkov);
}
// void raw2Root::ReadMCTree(TTree *tree){
//    Reset();
//    tree->SetBranchAddress("EventNum",&_Event_No);
//    tree->SetBranchAddress("CellID",&cellID);
//    tree->SetBranchAddress("Hit_X",&Hit_X);
//    tree->SetBranchAddress("Hit_Energy",&Hit_E);
//    tree->SetBranchAddress("Hit_Time",&Hit_Time);
//    //tree->SetBranchAddress("EvtID",&_Event_No);
//    //tree->SetBranchAddress("vecHcalCellID",&cellID);
//    //tree->SetBranchAddress("vecHcalVisibleEdepCell",&Hit_E);
//    //tree->SetBranchAddress("vecHcalHitTimeCell",&Hit_Time);
// }

void raw2Root::ReadMCTree(TTree *tree)
{
    Reset();
    // tree ->SetBranchAddress("Run_Num",&_Run_No);
    // tree ->SetBranchAddress("Event_Time",&_Event_Time);
    tree->SetBranchAddress("EventNum", &_Event_No);
    // tree ->SetBranchAddress("DetectorID",&_Detector_ID);
    tree->SetBranchAddress("CellID", &cellID);
    tree->SetBranchAddress("Hit_Energy", &Hit_E);
    tree->SetBranchAddress("Hit_X", &Hit_X);
    tree->SetBranchAddress("Hit_Y", &Hit_Y);
    tree->SetBranchAddress("Hit_Z", &Hit_Z);
    // tree ->SetBranchAddress("Digi_Energy_HCAL",&Digi_Energy);
    tree->SetBranchAddress("Hit_Time", &Hit_Time);
    // tree ->SetBranchAddress("Cherenkov",&cherenkov);
}
void raw2Root::SetDataTree(TTree *tree)
{
    Reset();
    tree->Branch("Run_Num", &_Run_No);
    tree->Branch("Event_Time", &_Event_Time);
    tree->Branch("CycleID", &_cycleID);
    tree->Branch("TriggerID", &_triggerID);
    tree->Branch("CellID", &data_cellid);
    // tree ->SetBranchAddress("layerIDs",&_layer_id);
    tree->Branch("BCID", &bcid);
    tree->Branch("HitTag", &hitTag);
    tree->Branch("GainTag", &gainTag);
    tree->Branch("HG_Charge", &HG_Charge);
    tree->Branch("LG_Charge", &LG_Charge);
    tree->Branch("Hit_Time", &Hit_Time);
    tree->Branch("Cherenkov", &cherenkov);
}
void raw2Root::CalibBranchClear()
{
    cellID->clear();
    Hit_E->clear();
    Hit_X->clear();
    Hit_Y->clear();
    Hit_Z->clear();
    Hit_Time->clear();
    cherenkov->clear();
}
void raw2Root::BranchClear()
{
    _cellID.clear();
    _Hit_E.clear();
    _Hit_X.clear();
    _Hit_Y.clear();
    _Hit_Z.clear();
    _Hit_Time.clear();
    _cherenkov.clear();
    if(data_cellid){
        HG_Charge->clear();
        LG_Charge->clear();
        Hit_Time->clear();
        data_cellid->clear();
        hitTag->clear();
    }
}
int raw2Root::ReadList(const string _list){
    ifstream data(_list);
	while(!data.eof())
	{
			string temp;
			data>>temp;
			if(temp=="")continue;
            if (data.fail())break;
			list.push_back(temp);
	}
    return 1;
}
