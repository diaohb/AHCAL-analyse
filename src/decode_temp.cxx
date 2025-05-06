#include <TFile.h>
#include <TTree.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <stdio.h>
#include <vector>
#include <TMath.h>
#include <string>
#include <bitset>
using namespace std;
extern char char_tmp[200];
vector<double> *temperature = 0;
vector<int> *layer = 0;
vector<int> *channel = 0;
int day;
int hour;
int minute;
int second;
int CatchEventBag(ifstream &f_in, vector<int> &buffer_v);
int DecodeAEvent(vector<int> &buffer_v);
void SetTreeBranch(TTree *tree);
void BranchClear();
int main(int argc, char *argv[])
{
    string input_list=argv[1];
    string output_dir;
    if(argc==3){
        output_dir=argv[2];
    }
    else{
        output_dir=".";
    }
	int layer_id;
	int cycleID;
	int triggerID;
    vector<int> _buffer_v;
    TFile *fout;
    TTree *tree;
    ifstream data(input_list);
	ifstream f_in;
	while(!data.eof())
	{
		string tmp_string;
        data >> tmp_string;
        if (tmp_string == "")
            continue;
        if (data.fail())break;
        f_in.open(tmp_string, ios::in);
        if (!f_in)
        {
            cout<<"cant open "<<tmp_string<<endl;
            return 0;
        }
        cout << "Reading " << tmp_string << " ..." << endl;
        tmp_string = tmp_string.substr(tmp_string.find_last_of('/') + 1);
        tmp_string = tmp_string.substr(0, tmp_string.find_last_of('.')); //remove .dat
        fout = TFile::Open(TString(output_dir + "/" + tmp_string + ".root"), "recreate");
        tree = new TTree("temperature", "data from binary file");
        SetTreeBranch(tree);
        while (!(f_in.eof()))
        {
            //if(Event_No%1000==0)cout<<"Event_No: "<<Event_No<<" Bag_No "<<Bag_No<<endl;
            _buffer_v.clear();
            BranchClear();
            // cout<<"catching a bag"<<endl;
            CatchEventBag(f_in,_buffer_v);
            // cout<<"catch a bag"<<endl;
            if(DecodeAEvent(_buffer_v))
                tree->Fill();
        }
        f_in.close();
        fout->cd();
        tree->Write();
        fout->Close();
    }
	return 1;
}
int CatchEventBag(ifstream &f_in, vector<int> &buffer_v)
{
    int buffer = 0;
    int tmp = 0;
    buffer_v.clear();
    while (f_in.read((char *)&buffer, 1))
    {
        // cout<<hex<<buffer<<" ";
        buffer_v.push_back(buffer);
        if(buffer_v.size()>8){
            buffer_v.erase(buffer_v.begin());
        }
        if (tmp == 0x00 && buffer == 0x00)
            break;
        tmp = buffer;
    }
    // cout<<hex<<tmp<<"  "<<buffer<<endl;
    while (f_in.read((char *)&buffer, 1))
    {
        buffer_v.push_back(buffer);
        if (tmp == 0x2f && buffer == 0x27)
            break;
        tmp = buffer;
    }
    int int_tmp = buffer_v.size();
    // int_tmp = buffer_v.size();
    // cout<<int_tmp<<endl;
    // cout<<tmp<<"  "<<buffer<<"  "<<buffer_v.at(int_tmp-2)<<"  "<<buffer_v.at(int_tmp-1)<<endl;
    return 1;
}
int DecodeAEvent(vector<int> &buffer_v)
{
    int size = buffer_v.size();
    if (size != 7684)
    {
        cout << "wrong bag size " << buffer_v.size() << endl;
        return 0;
    }
    day = buffer_v.at(0);
    hour = buffer_v.at(1);
    minute = buffer_v.at(2);
    second = buffer_v.at(3);
    for (int i = 0; i < 40 * 48; i++)
    {
        temperature->push_back((buffer_v.at(4 + 4 * i) * 256 + buffer_v.at(5 + 4 * i)) / 128.);
        channel->push_back(buffer_v.at(6 + 4 * i));
        layer->push_back(buffer_v.at(7 + 4 * i));
    }
    return 1;
}
void SetTreeBranch(TTree *tree){
    tree->Branch("temperature", &temperature);
    tree->Branch("layer", &layer);
    tree->Branch("channel", &channel);
    tree->Branch("day", &day);
    tree->Branch("hour", &hour);
    tree->Branch("minute", &minute);
    tree->Branch("second", &second);
}
void BranchClear(){
    temperature->clear();
    layer->clear();
    channel->clear();
    day = 0;
    hour = 0;
    minute = 0;
    second = 0;
}
