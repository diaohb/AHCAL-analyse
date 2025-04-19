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
int lowgain = 0;
int highgain = 0;
int channel = 0;
int col_id = 0;
int chip_id = 0;
int bx_id = 0;
int CatchEventBag(ifstream &f_in, vector<int> &buffer_v) {
    int buffer = 0;
    int tmp = 0;
    buffer_v.clear();
    while (f_in.read((char *) &buffer, 1)) {
        if (tmp == 0xfa && buffer == 0x5a)
            break;
        tmp = buffer;
    }
    while (f_in.read((char *) &buffer, 1)) {
        buffer_v.push_back(buffer);
        if (tmp == 0xfe && buffer == 0xee)
            break;
        tmp = buffer;
    }
    // int_tmp = buffer_v.size();
    return 1;
}
int CatchChipBag(vector<int> &buffer_v, vector<int> &chip_v) {
    int size = buffer_v.size();
    if (size < 296) {
        // cout << "empty buffer" << endl;
        return 0;
    }
    if (size % 294 != 2) {
        // cout << "wrong chip size " << buffer_v.size() << endl;
        return 0;
    }
    chip_v.clear();
    for (int i = 0; i < 147; i++) {
        chip_v.push_back(buffer_v.at(2 * i) * 256 + buffer_v.at(2 * i + 1));
    }
    buffer_v.erase(buffer_v.begin(), buffer_v.begin() + 293);
    return 1;
}
int DecodeAEvent(vector<int> &chip_v) {
    int size = chip_v.size();
    size = chip_v.size();
    if (size < 75) {
        // cout << "wrong size " << size << endl;
        return 0;
    }
    int col_num = size / 73;
    highgain = chip_v.at(channel_No - channel - 1) & 0x0fff;
    lowgain = chip_v.at(channel_No * 2 - channel - 1) & 0x0fff;
    chip_id = chip_v.at(size - 1);
    bx_id = chip_v.at(size - col_num - 1);
    chip_v.erase(chip_v.begin(), chip_v.begin() + 72);
    return 1;
}
int main(int argc, char *argv[]) {
    string input_list = argv[1];
    string output_file = argv[2];
    vector<int> _buffer_v;
    vector<int> _chip_v;
    TFile *fout = TFile::Open(output_file.c_str(), "RECREATE");
    TTree *tree = new TTree("e_calibrate", "data from binary file");
    int voltage;
    tree->Branch("dac", &voltage);
    tree->Branch("lowgain", &lowgain);
    tree->Branch("highgain", &highgain);
    tree->Branch("channel", &channel);
    tree->Branch("chip_id", &chip_id);
    tree->Branch("bx_id", &bx_id);
    tree->Branch("col_id", &col_id);
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
        tmp_string = tmp_string.substr(tmp_string.find("chn") + 3);
        tmp_string = tmp_string.substr(0, tmp_string.find("DAC"));
        string str_chn = tmp_string.substr(0, tmp_string.find("_"));
        string str_dac = tmp_string.substr(tmp_string.find("_") + 1);
        stringstream geek(str_dac);
        geek >> voltage;
        stringstream geek_chn(str_chn);
        geek_chn >> channel;
        cout << "channel: " << str_chn << " voltage: " << voltage << endl;
        while (!(f_in.eof())) {
            // if(Event_No%1000==0)cout<<"Event_No: "<<Event_No<<" Bag_No "<<Bag_No<<endl;
            _buffer_v.clear();
            CatchEventBag(f_in, _buffer_v);
            // cout << "catch a bag" << endl;
            while (CatchChipBag(_buffer_v, _chip_v)) {
                col_id = 0;
                // cout << "decoding" << endl;
                while (DecodeAEvent(_chip_v)) {
                    tree->Fill();
                    col_id++;
                }
            }
        }
        cout << "decode " << tree->GetEntries() << " events" << endl;
        f_in.close();
    }
    tree->Write("", TObject::kOverwrite);
    fout->Close();
    return 1;
}
