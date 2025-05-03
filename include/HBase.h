#ifndef HBASE_HH
#define HBASE_HH

#include "TFile.h"
#include "TTree.h"
#include <TH2D.h>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

using namespace std;

class HBase {
public:
    //Constructor, destructor and instance of base class
    HBase();
    virtual ~HBase();

    // protected:
    //Protected member functions
    virtual void ReadTree(const TString &fname, const TString &tname);//Read TTree from ROOT files
    virtual void ReadList(const string &_list);                       // Read the file list and save to the protected vector
    virtual void CreateFile(const TString &_outname);                 // Create output file
    virtual void Init(const TString &_outname);                       // Initialize derived members
                                                                      // protected:
    // Protected member variables
    vector<string> list;
    TFile *fin; //TFile pointer to open files
    TTree *tin; //Read Tree from fin
    TFile *fout;//Create output files
    TTree *tout;//Create output trees
    int Event_No = 0;
    double Energy_HCAL = 0;
    vector<int> *detectorID = 0;
    vector<int> *cellID = 0;
    vector<double> *Hit_Energy = 0;
    vector<double> *Hit_Time = 0;
    vector<double> *Hit_X = 0;
    vector<double> *Hit_Y = 0;
    vector<double> *Hit_Z = 0;
    vector<double> *SiPM_Energy = 0;
    vector<int> *Cherenkov = 0;
    void Clear();
};

#endif
