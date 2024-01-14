#ifndef HBASE_HH
#define HBASE_HH

#include <TH2D.h>
#include <vector>
#include "TFile.h"
#include "TTree.h"
#include <fstream>
#include <iostream>
#include <string>

using namespace std;

class HBase{
		public:
				//Constructor, destructor and instance of base class
				HBase();
				virtual ~HBase();

		// protected:
				//Protected member functions
				virtual void ReadTree(const TString &fname,const TString &tname,const string mode); //Read TTree from ROOT files
				virtual void ReadList(const string &_list); // Read the file list and save to the protected vector
				virtual void CreateFile(const TString &_outname); // Create output file
				virtual void Init(const TString &_outname);// Initialize derived members
		// protected:
				// Protected member variables
				vector<string>	list;
				TFile *fin;	//TFile pointer to open files
				TTree *tin;	//Read Tree from fin
				TFile *fout;	//Create output files
				TTree *tout;	//Create output trees
				int   _Event_No;
				double _Digi_Energy_ECAL;
				double _Digi_Energy_HCAL;
				vector< int > *_detectorID;
				vector< int > *_cellID;
				vector< double > *_Digi_Hit_Energy;
				vector< double > *_Hit_Time;
				vector< double > *_Hit_X;
				vector< double > *_Hit_Y;
				vector< double > *_Hit_Z;
				void Clear();

};

#endif
