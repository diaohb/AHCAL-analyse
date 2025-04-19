#include <TFile.h>
#include <TMath.h>
#include <TTree.h>
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <vector>
using namespace std;
extern char char_tmp[200];
int int_tmp = 0;
int flag = 0;
int channel_No = 36;
double lowgain = 0;
double highgain = 0;
int channel = 0;
int CatchEventBag(ifstream& f_in, vector<int>& buffer_v) {
    int buffer = 0;
    int tmp = 0;
    buffer_v.clear();
    while (f_in.read((char*)&buffer, 1)) {
        // cout<<hex<<buffer<<" ";
        if (tmp == 0xfa && buffer == 0x5a)
            break;
        tmp = buffer;
    }
    // cout<<hex<<tmp<<"  "<<buffer<<endl;
    while (f_in.read((char*)&buffer, 1)) {
        buffer_v.push_back(buffer);
        if (tmp == 0xfe && buffer == 0xee)
            break;
        tmp = buffer;
    }
    int_tmp = buffer_v.size();
    // cout<<int_tmp<<endl;
    // cout<<hex<<tmp<<"  "<<buffer<<"  "<<buffer_v.at(int_tmp-2)<<"
    // "<<buffer_v.at(int_tmp-1)<<endl;
    return 1;
}
int DecodeAEvent(vector<int>& chip_v) {
    int size = chip_v.size();
    if (size != 150) {
        cout << "wrong chip size " << chip_v.size() << endl;
        return 0;
    }
    vector<int> buffer_v;
    for (int i = 0; i < size / 2; i++) {
        // int t=chip_v.at(2*i);
        // int tt=chip_v.at(2*i+1);
        // cout<<(bitset<sizeof(t)*8>(t))<<"  "<<(bitset<sizeof(t)*8>(tt))<<endl;
        buffer_v.push_back(chip_v.at(2 * i) * 256 + chip_v.at(2 * i + 1));
        // t=buffer_v.at(i);
        // cout<<(bitset<sizeof(t)*8>(t))<<endl;
    }
    // chip_v.erase(chip_v.begin()+size/2,chip_v.end());
    for (int i_ch = 0; i_ch < channel_No; ++i_ch) {  // 72+BCID+Chip
        int hit = buffer_v.at(i_ch) & 0x1000;
        if (hit < 1)
            continue;
        // int gain = chip_v.at(i_ch+channel_No)&0x2000;
        int tdc = buffer_v.at(i_ch) & 0x0fff;
        // int gainTag_tdc  = chip_v.at(i_ch)&0x2000;
        int adc = buffer_v.at(i_ch + channel_No) & 0x0fff;
        highgain = tdc;
        lowgain = adc;
        channel = channel_No - 1 - i_ch;
        // cout<<hex<<chip_v.at(i_ch)<<dec<<"  "<<channel<<endl;
        // cout<<hex<<chip_v.at(i_ch)<<"  "<<tdc<<endl;
    }
    return 1;
}
int main(int argc, char* argv[]) {
    string input_list = argv[1];
    string output_file = argv[2];
    int layer_id;
    int cycleID;
    int triggerID;
    vector<int> _buffer_v;
    TFile* fout = TFile::Open(output_file.c_str(), "RECREATE");
    TTree* tree = new TTree("e_calibrate", "data from binary file");
    int voltage;
    tree->Branch("volatage", &voltage);
    tree->Branch("lowgain", &lowgain);
    tree->Branch("highgain", &highgain);
    tree->Branch("channel", &channel);
    ifstream data(input_list);
    ifstream f_in;
    while (!data.eof()) {
        // cout<<"safdsaf"<<endl;
        string tmp_string;
        data >> tmp_string;
        if (tmp_string == "")
            continue;
        if (data.fail())
            break;
        f_in.open(tmp_string, ios::in);
        if (!f_in) {
            cout << "cant open " << tmp_string << endl;
            return 0;
        }
        cout << "Reading " << tmp_string << " ..." << endl;
        tmp_string = tmp_string.substr(tmp_string.find_last_of('/') + 1);
        // tmp_string=tmp_string.substr(0,tmp_string.find_last_of('.'));
        // string str_out=output_file+"/"+tmp_string+".root";
        tmp_string = tmp_string.substr(0, tmp_string.find("mV"));
        stringstream geek(tmp_string);
        geek >> voltage;
        while (!(f_in.eof())) {
            // if(Event_No%1000==0)cout<<"Event_No: "<<Event_No<<" Bag_No "<<Bag_No<<endl;
            _buffer_v.clear();
            // cout<<"catching a bag"<<endl;
            CatchEventBag(f_in, _buffer_v);
            // cout<<"catch a bag"<<endl;
            DecodeAEvent(_buffer_v);
            tree->Fill();
        }
        f_in.close();
    }
    tree->Write();
    fout->Close();
    return 1;
}
