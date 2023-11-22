#include "HBase.h"

using namespace std;

HBase::HBase() : fin(0),fout(0),tin(0),tout(0)
{
		list.clear();
		cout<<"HBase class instance initialized."<<endl;
}

HBase::~HBase()
{
		cout<<"Base destructor called"<<endl;
		fout->Close();
		//fin->Close();
}

void HBase::Init(const TString &_outname)
{
		CreateFile(_outname);
}

void HBase::CreateFile(const TString &_outname)
{
		fout = new TFile(TString(_outname),"RECREATE");
}


void HBase::ReadList(const string &_list)
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

void HBase::ReadTree(const TString &fname,const TString &tname,const string mode="digi")
{
		cout<<"Reading file "<<fname<<endl;
		fin = TFile::Open(TString(fname),"READ");
		tin = (TTree*)fin->Get(TString(tname));
		_Event_No=0;_detectorID=0;_cellID=0;_Digi_Hit_Energy=0;_Digi_Energy_ECAL=0;_Digi_Energy_HCAL=0;_Hit_Time=0;_Hit_X=0;_Hit_Y=0;_Hit_Z=0;
		// tin->SetBranchAddress("EventNum",&_Event_No);
		// tin->SetBranchAddress("DetectorID",&_detectorID);
		// tin->SetBranchAddress("CellID",&_cellID);
		if(mode=="digi"){
			tin->SetBranchAddress("Digi_Hit_Energy",&_Digi_Hit_Energy);
		}
		else if(mode=="mc"){
			
			tin->SetBranchAddress("Hit_Energy",&_Digi_Hit_Energy);
		}
		// tin->SetBranchAddress("Digi_Energy_ECAL",&_Digi_Energy_ECAL);
		// tin->SetBranchAddress("Digi_Energy_HCAL",&_Digi_Energy_HCAL);
		// tin->SetBranchAddress("Hit_Time",&_Hit_Time);
		tin->SetBranchAddress("Hit_X",&_Hit_X);
		tin->SetBranchAddress("Hit_Y",&_Hit_Y);
		tin->SetBranchAddress("Hit_Z",&_Hit_Z);
}
