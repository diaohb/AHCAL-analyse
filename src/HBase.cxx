#include "HBase.h"

using namespace std;

HBase::HBase() : fin(0), fout(0), tin(0), tout(0) {
    list.clear();
    cout << "HBase class instance initialized." << endl;
}

HBase::~HBase() {
    cout << "Base destructor called" << endl;
    tin->Delete();
    fin->Close();
    //fin->Close();
}

void HBase::Init(const TString &_outname) {
    CreateFile(_outname);
}

void HBase::CreateFile(const TString &_outname) {
    fout = new TFile(TString(_outname), "RECREATE");
}


void HBase::ReadList(const string &_list) {
    ifstream data(_list);
    while (!data.eof()) {
        string temp;
        data >> temp;
        if (temp == "") continue;
        if (data.fail()) break;
        list.push_back(temp);
    }
}

void HBase::ReadTree(const TString &fname, const TString &tname) {
    cout << "Reading file " << fname << endl;
    fin = TFile::Open(TString(fname), "READ");
    tin = (TTree *) fin->Get(TString(tname));

    // tin->SetBranchAddress("EventNum",&_Event_No);
    // tin->SetBranchAddress("DetectorID",&_detectorID);
    // tin->SetBranchAddress("CellID",&_cellID);
    tin->SetBranchAddress("Hit_Energy", &Hit_Energy);
    // tin->SetBranchAddress("Hit_Time",&_Hit_Time);
    tin->SetBranchAddress("Hit_X", &Hit_X);
    tin->SetBranchAddress("Hit_Y", &Hit_Y);
    tin->SetBranchAddress("Hit_Z", &Hit_Z);
    tin->SetBranchAddress("Cherenkov", &Cherenkov);
    // tin->SetBranchAddress("SiPM_Energy", &SiPM_Energy);
}
void HBase::Clear() {
    Event_No = 0;
    Energy_HCAL = 0;
    if (cellID) cellID->clear();
    if (detectorID) detectorID->clear();
    if (Hit_Energy) Hit_Energy->clear();
    if (Hit_Time) Hit_Time->clear();
    if (Hit_X) Hit_X->clear();
    if (Hit_Y) Hit_Y->clear();
    if (Hit_Z) Hit_Z->clear();
    if (Cherenkov) Cherenkov->clear();
    if (SiPM_Energy) SiPM_Energy->clear();
}
