#include "HBase_data.h"

using namespace std;

HBase_data::HBase_data() : fin(0),fout(0),tin(0),tout(0)
{
		list.clear();
		cout<<"HBase_data class instance initialized."<<endl;
}

HBase_data::~HBase_data()
{
		cout<<"Base_data destructor called"<<endl;
		fout->Close();
		//fin->Close();
}

void HBase_data::Init(const TString &_outname)
{
		CreateFile(_outname);
}

void HBase_data::CreateFile(const TString &_outname)
{
		fout = new TFile(TString(_outname),"RECREATE");
}


void HBase_data::ReadList(const string &_list)
{
		ifstream data(_list);
		while(!data.eof())
		{
				string temp;
				data>>temp;
				if(temp=="")continue;
				list.push_back(temp);
		}
}

void HBase_data::ReadTree(const TString &fname,const TString &tname)
{
		cout<<"Reading tree "<<fname<<endl;
		fin = TFile::Open(TString(fname),"READ");
		tin = (TTree*)fin->Get(TString(tname));
		_cellID=0;_bcid=0;_hitTag=0;_gainTag=0;_cherenkov=0;_HG_Charge=0;_LG_Charge=0;_Hit_Time=0;_hitenergy=0;
		// tin->SetBranchAddress("Run_Num",&_Run_No);
		// tin->SetBranchAddress("Event_Time",&_Event_Time);
		// tin->SetBranchAddress("CycleID",&_cycleID);
		// tin->SetBranchAddress("TriggerID",&_triggerID);
		tin->SetBranchAddress("CellID",&_cellID);
		tin->SetBranchAddress("Hit_Energy",&_hitenergy);
		// tin->SetBranchAddress("BCID",&_bcid);
		// tin->SetBranchAddress("HitTag",&_hitTag);
		// tin->SetBranchAddress("GainTag",&_gainTag);
		// tin->SetBranchAddress("HG_Charge",&_HG_Charge);
		// tin->SetBranchAddress("LG_Charge",&_LG_Charge);
		// tin->SetBranchAddress("Hit_Time",&_Hit_Time);
		// tin->SetBranchAddress("Cherenkov",&_cherenkov);
}
