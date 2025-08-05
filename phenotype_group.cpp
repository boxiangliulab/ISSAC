#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<map>
#include<set>
#include<sstream>
#include<algorithm>
#include<cmath>
#include<vector>
#include<map>
#include<numeric>
#include<utility>
#include "pheno.h"

using namespace std;

bool haveCommonElements(const std::set<string>& set1, const std::set<string>& set2) {
    // Create an empty set to store the intersection
    std::set<string> intersection;

    // Use std::set_intersection to find common elements
    std::set_intersection(set1.begin(), set1.end(),
                          set2.begin(), set2.end(),
                          std::inserter(intersection, intersection.begin()));

    // If the intersection set is not empty, the sets share common elements
    return !intersection.empty();
}

map<string,vector<string>> cluster(set<string> intron_set){
    map<string,map<int,set<string>>> cluster_site;
    map<string,vector<string>> cluster_intron;
    string chrom,strand,donor,acceptor;
    string tmp1,tmp2,tmp3;

    for(const string & elem : intron_set){
        chrom = elem.substr(0,elem.find(":"));
        tmp1 = elem.substr(elem.find(":")+1);
        strand = tmp1.substr(0,tmp1.find(":"));
        tmp2 = tmp1.substr(tmp1.find(":")+1);
        donor = tmp2.substr(0,tmp2.find(":"));
        acceptor = tmp2.substr(tmp2.find(":")+1);
        tmp3 = chrom + ":" + strand;
        if(cluster_site.find(tmp3)!=cluster_site.end()){
            bool if_find = false;
            for(const auto& [key,valueSet] : cluster_site[tmp3]){
                for(const auto&value : valueSet){
                    if((value==donor)||(value==acceptor)){
                        cluster_site[tmp3][key].insert(donor);
                        cluster_site[tmp3][key].insert(acceptor);
                        string clu_name = tmp3 + ":" + to_string(key);
                        cluster_intron[clu_name].push_back(elem);
                        if_find = true;
                        break;
                    }
                }
                if(if_find == true){
                    break;
                }
            }
            if(if_find == false){
                int seq = cluster_site[tmp3].size()+1;
                cluster_site[tmp3][seq].insert(donor);
                cluster_site[tmp3][seq].insert(acceptor);
                string clu_name = tmp3 + ":" + to_string(seq);
                vector<string> tmp_vec;
                tmp_vec.push_back(elem);
                cluster_intron[clu_name] = tmp_vec;
            }
        }
        else{
            map<int,set<string>> tmp_set;
            tmp_set[1].insert(donor);
            tmp_set[1].insert(acceptor);
            cluster_site[tmp3] = tmp_set;
        }

    }
    // merge all the clusters with common site
    for(const auto & outer_pair : cluster_site){
        const string & outer_key = outer_pair.first;
        const auto& inner_map = outer_pair.second;
        //
        cout<<outer_key<<endl;
        for(const auto & inner_pair1 : inner_map){
            int inner_key1=inner_pair1.first;
            cout<<"inner_key1"<<inner_key1<<endl;
            const auto & string_set1 = inner_pair1.second;
            if(string_set1.size()>0){
            for(const auto & inner_pair2 : inner_map){
                int inner_key2=inner_pair2.first;
                if(inner_key2!=inner_key1){
                const auto & string_set2 = inner_pair2.second;
                cout<<"loaded"<<endl;
                if(string_set2.size()>0){
                    cout<<inner_key1<<"\t"<<inner_key2<<"\t"<<string_set1.size()<<"\t"<<string_set2.size()<<endl;
                    if(haveCommonElements(string_set1,string_set2)){
                        //add elements of set2 to set1 and clear set2
                        cluster_site[outer_key][inner_key1].insert(string_set2.begin(),string_set2.end());
                        cluster_site[outer_key][inner_key2].clear();
                        string clu_name1=outer_key + ":" + to_string(inner_key1);
                        string clu_name2=outer_key + ":" + to_string(inner_key2);
                        cluster_intron[clu_name1].insert(cluster_intron[clu_name1].end(),cluster_intron[clu_name2].begin(),cluster_intron[clu_name2].end());
                        cluster_intron.erase(clu_name2);
                    }
                }
                cout<<"final"<<endl;
            }}
        }}
    }
    return cluster_intron;
}

void print(map<string,vector<string>> cluster_intron, map<string,vector<int>> intron_read_map,string out_file){
    ofstream fout(out_file);
    string intron_name;
    for(const auto & [value,valuevec]: cluster_intron){
        for(int i=0;i<cluster_intron[value].size();i++){
            intron_name = cluster_intron[value][i];
            fout<<value<<" "<<cluster_intron[value][i];
            for(int j=0;j<intron_read_map[intron_name].size();j++){
                fout<<" "<<intron_read_map[intron_name][j];
            }
            fout<<endl;
        }
    }
    fout.close();
}

void print_total(map<string,vector<string>> cluster_intron, map<string,vector<int>> intron_read_map, string refined_output){
    ofstream fout(refined_output);
    string intron_name, chrom, strand, donor, acceptor;
    string tmp_1,tmp_2;
    int total_read;
    for(const auto & [value,valuevec]: cluster_intron){
        fout<<value;
        for(int i=0;i<cluster_intron[value].size();i++){
            intron_name = cluster_intron[value][i];
            chrom = intron_name.substr(0,intron_name.find(":"));
            tmp_1 = intron_name.substr(intron_name.find(":")+1);
            strand = tmp_1.substr(0,tmp_1.find(":"));
            tmp_2 = tmp_1.substr(tmp_1.find(":")+1);
            donor = tmp_2.substr(0,tmp_2.find(":"));
            acceptor = tmp_2.substr(tmp_2.find(":")+1);
            total_read = accumulate(intron_read_map[intron_name].begin(),intron_read_map[intron_name].end(),0);
            fout<<" "<<donor<<":"<<acceptor<<":"<<to_string(total_read);
        }
        fout<<endl;
    }
}

pair<map<int,set<string>>,map<int,set<string>>> site_anno(set<int> site, map<int,set<int>> site_corres){
    int i=0;
    int end=site.size();
    pair<map<int,set<string>>,map<int,set<string>>> cluster_site_inclu_exclu;
    map<int,set<string>> cluster_include;
    map<int,set<string>> cluster_exclude;
    for(const auto&elem:site){
        set<string> included;
        set<string> excluded;
        for(const auto&elem1:site_corres[elem]){
            //check 5' splice site
            if(i==0){
                for(const auto&elem2:site_corres[elem1]){
                    bool competive=true;
                    for(const auto&elem3:site_corres[elem]){
                        if(elem2>elem3){
                            competive=false;
                        }
                    }
                    if(elem2==elem){
                        competive=false;
                    }
                    if(competive==true){
                        string intron=to_string(elem2)+":"+to_string(elem1);
                        excluded.insert(intron);
                    }
                }
            }
            else if(i==end-1){
                for(const auto&elem2:site_corres[elem1]){
                    bool competive=true;
                    for(const auto&elem3:site_corres[elem]){
                        if(elem2<elem3){
                            competive=false;
                        }
                    }
                    if(elem2==elem){
                        competive=false;
                    }
                    if(competive==true){
                        string intron=to_string(elem1)+":"+to_string(elem2);
                        excluded.insert(intron);
                    }
                }
            }
            if(elem1>elem){
                string intron=to_string(elem)+":"+to_string(elem1);
                included.insert(intron);
            }
            else{
                string intron=to_string(elem1)+":"+to_string(elem);
                included.insert(intron);
            }
        }
        if((i>0)&&(i<end-1)){
            for(const auto& [value,valueset]:site_corres){
                for(const auto&elem1:site_corres[value]){
                    if((value>elem)&&(elem>elem1)){
                        string intron = to_string(elem1)+":"+to_string(value);
                        excluded.insert(intron);
                    }
                    if((value<elem)&&(elem<elem1)){
                        string intron = to_string(value)+":"+to_string(elem1);
                        excluded.insert(intron);
                    }
                }
            }
        }
        cluster_include[elem] = included;
        cluster_exclude[elem] = excluded;
        i++;
    }
    cluster_site_inclu_exclu.first = cluster_include;
    cluster_site_inclu_exclu.second = cluster_exclude;
    return cluster_site_inclu_exclu;
}


void site_detection(map<string,vector<string>> cluster_intron, map<string,vector<int>> intron_read_map, vector<string> sample_name, string inclu_exclu_out, string site_readcount){
    string intron_name, chrom, strand, donor, acceptor;
    string tmp_1,tmp_2;
    int total_read;
    ofstream fout(inclu_exclu_out);
    ofstream fout1(site_readcount);
    // write in sample name
    fout1<<sample_name[0];
    for(int i=1;i<sample_name.size();i++){
        fout1<<" "<<sample_name[i];
    }
    fout1<<endl;
    for(const auto & [value,valuevec]: cluster_intron){
        set<int> site;
        map<int,set<int>> site_corres;
        if(cluster_intron[value].size()>1){
        for(int i=0;i<cluster_intron[value].size();i++){
            intron_name = cluster_intron[value][i];
            chrom = intron_name.substr(0,intron_name.find(":"));
            tmp_1 = intron_name.substr(intron_name.find(":")+1);
            strand = tmp_1.substr(0,tmp_1.find(":"));
            tmp_2 = tmp_1.substr(tmp_1.find(":")+1);
            donor = tmp_2.substr(0,tmp_2.find(":"));
            acceptor = tmp_2.substr(tmp_2.find(":")+1);
            site.insert(stoi(donor));
            site.insert(stoi(acceptor));
            site_corres[stoi(donor)].insert(stoi(acceptor));
            site_corres[stoi(acceptor)].insert(stoi(donor));
        }
        pair<map<int,set<string>>,map<int,set<string>>> clu_ex_in = site_anno(site,site_corres);
        for(const auto &elem:site){
            fout<<value<<" "<<elem<<endl;
            fout<<"included ";
            for(const auto &elem1:clu_ex_in.first[elem]){
                fout<<elem1<<" ";
            }
            fout<<endl;
            fout<<"excluded ";
            for(const auto &elem2:clu_ex_in.second[elem]){
                fout<<elem2<<" ";
            }
            fout<<endl;
        }

        for(const auto&elem:site){
            if(clu_ex_in.second[elem].size()>0){
                chrom=value.substr(0,value.find(":"));
                tmp_1=value.substr(value.find(":")+1);
                strand=tmp_1.substr(0,tmp_1.find(":"));
                fout1<<chrom<<":"<<strand<<":"<<to_string(elem);
                cout<<chrom<<":"<<strand<<":"<<to_string(elem)<<endl;
                for(int sam_num=0;sam_num<sample_name.size();sam_num++){
                    int include_num=0;
                    int total_num=0;
                    for(const auto &included_intron:clu_ex_in.first[elem]){
                        intron_name=chrom+":"+strand+":"+included_intron;
                        include_num+=intron_read_map[intron_name][sam_num];
                        total_num+=intron_read_map[intron_name][sam_num];
                    }
                    for(const auto &excluded_intron:clu_ex_in.second[elem]){
                        intron_name=chrom+":"+strand+":"+excluded_intron;
                        total_num+=intron_read_map[intron_name][sam_num];
                    }
                    fout1<<" "<<to_string(include_num)<<":"<<to_string(total_num);
                }
                fout1<<endl;
            }
        }
        }
    }
    fout.close();
    fout1.close();
}



int pheno_group(int argc,char* argv[]){
    if(argc!=6){
        cerr<<"Usage: "<<argv[0]<<endl;
        return 1;
    }
    string sample_file=argv[1];
    string junc_pos=argv[2];
    string out_file_prefix=argv[3];
    string threshold_read=argv[4]; //100
    string log_out=argv[5]; //100
    //output
    string out_file=out_file_prefix+".intron.out";
    string refined_output=out_file_prefix+".refined";
    string site_readcount=out_file_prefix+".site";
    string inclu_exclu_out=out_file_prefix+".inclu_exclu";
    //output
    vector<string> sample_name;
    ifstream fin;
    fin.open(sample_file);
    string line;
    while(getline(fin,line)){
        sample_name.push_back(line);
    }
    fin.close();
    map<string,vector<int>> intron_read_map;
    string junc_name;
    cout<<"Start construct junction cluster"<<endl;
    string intron,num_read;
    int int_num_read;
    int length = sample_name.size();
    set<string> intron_set;
    ofstream fout;
    fout.open(log_out);
    for(int i=0;i<length;i++){
        ifstream fin1;
        junc_name = junc_pos + "/" + sample_name[i] + ".stat";
        fin1.open(junc_name);
        fout<<i<<endl;
        while(getline(fin1,line)){
            intron = line.substr(0,line.find("\t"));
            intron_set.insert(intron);
            num_read = line.substr(line.find("\t")+1);
            int_num_read = stoi(num_read);
            if(intron_read_map.find(intron)!=intron_read_map.end()){
                //intron_read_map[intron][i] = int_num_read;
                intron_read_map[intron].push_back(int_num_read);
            }
            else{
                intron_read_map[intron] = vector<int>(1,0);
                intron_read_map[intron][0] = int_num_read;
            }
        }
        fin1.close();
    }
    //filter intron total read count>100;
    set<string> filter_intron;
    for(const auto& element:intron_set){
        int sum = accumulate(intron_read_map[element].begin(),intron_read_map[element].end(),0);
        if(sum>stoi(threshold_read)){
            filter_intron.insert(element);
        }
    }
    // cluster each intron
    fout<<"size of intron\t"<<filter_intron.size()<<endl;

    map<string,vector<string>> clu_intron = cluster(filter_intron);
    fout<<"Site detection"<<endl;
    //rebuild intron_read_map
    intron_read_map.clear();
    for(const auto& element:filter_intron){
        intron_read_map[element] = vector<int>(length,0);
    }
    //
    for(int i=0;i<length;i++){
        ifstream fin1;
        junc_name = junc_pos + "/" + sample_name[i] + ".stat";
        fin1.open(junc_name);
        fout<<i<<endl;
        while(getline(fin1,line)){
            intron = line.substr(0,line.find("\t"));
            intron_set.insert(intron);
            num_read = line.substr(line.find("\t")+1);
            int_num_read = stoi(num_read);
            if((int_num_read>10000)||(int_num_read<-10000))fout<<"Exception"<<sample_name[i]<<endl;
            if(filter_intron.find(intron)!=filter_intron.end()){
                intron_read_map[intron][i] = int_num_read;
                //intron_read_map[intron].push_back(int_num_read);
            }
        }
        fin1.close();
    }

    site_detection(clu_intron, intron_read_map,sample_name,inclu_exclu_out,site_readcount);
    //
    print(clu_intron, intron_read_map,out_file);
    print_total(clu_intron, intron_read_map,refined_output);
    return 0;
}