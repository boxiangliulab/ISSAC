#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<map>
#include<set>
#include<sstream>
#include<algorithm>
#include "junctions_extractor.h"
using namespace std;

vector<string> split(const string& str,const string& sep) {
    vector<string> result;
    string r1 = str;
    for(int i=0;i<=8;i++){
        string item=r1.substr(0,r1.find("\t"));
        string r2=r1.substr(r1.find("\t"));
        r1=r2.substr(1);
        result.push_back(item);
    }
    result.push_back(r1);
    return result;
}

int juncstat(int argc,char* argv[]){
    string barcode_file, junc_file, out_file;
    int c;
    while ((c = getopt(argc, argv, "hb:j:o:")) != -1) {
    switch (c) {
        case 'h':
            cout << "Usage: " << argv[0] << " [options]\n\n"
                 << "Options:\n"
                 << "  -h          Show this help message and exit\n"
                 << "  -b <file>   Barcode file\n"
                 << "  -j <file>   Junction file\n"
                 << "  -o <file>   Output file\n";
            exit(0);
        case 'b':
            barcode_file = string(optarg);
            break;
        case 'j':
            junc_file = string(optarg);
            break;
        case 'o':
            out_file = string(optarg);
            break;
        case '?':
        default:
            throw runtime_error("Error parsing inputs!\n\n");
    }}
    if (barcode_file.empty() || junc_file.empty() || out_file.empty()) {
    cerr << "Error: options -b, -j, -o are all required.\n"
         << "Run with -h for usage.\n";
    return 1;}
    //read in barcode info
    ifstream fin;
    fin.open(barcode_file);
    if (!fin.is_open()) {
    cerr << "Failed to open the barcode file." << endl;
    return 1; // Or handle the error as appropriate
    }
    string line;
    vector<string> barcodelist;

    while(getline(fin,line)){
        barcodelist.push_back(line);
    }
    fin.close();
    //read in junc files
    ifstream fin1;
    fin1.open(junc_file);
    //fin1.open("test.junc");
    if (!fin1.is_open()) {
    cerr << "Failed to open the junc file." << endl;
    return 1; // Or handle the error as appropriate
    }

    //define a junc_cell ~ the numbers of UMI
    map<string,set<string>> dict;
    while(getline(fin1,line)){
        string sep="\t";
        vector<string> splitline = split(line,sep);
        auto it = find(barcodelist.begin(),barcodelist.end(),splitline[8]);
        if(it!=barcodelist.end()){
            string tmp_junc = splitline[0]+":"+splitline[4]+":"+splitline[5]+":"+splitline[6];
            string CB_UMI = splitline[8]+":"+splitline[9];
            dict[tmp_junc].insert(CB_UMI);
        }
    }
    fin1.close();

    ofstream fout(out_file);
    for (const auto& pair:dict){
        fout<<pair.first<<"\t"<<pair.second.size()<<endl;
    }
    fout.close();
    return 0;
}