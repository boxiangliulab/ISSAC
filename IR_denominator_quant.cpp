#include<iostream>
#include<fstream>
#include<string>
#include<vector>
#include<map>
#include<set>
#include<sstream>
#include<algorithm>
#include <iomanip>
#include <sstream>
#include <stdexcept>
#include "common.h"
#include "htslib/sam.h"
#include "htslib/hts.h"
#include "htslib/faidx.h"
#include "htslib/kstring.h"
#include <unordered_map>
#include <unordered_set>
#include "IR_denominator_quant.h"
#include <getopt.h>
#include "gtf_parser.h"

using namespace std;

//Parse the options passed to this tool
int USRExtractor::parse_options(int argc, char *argv[]) {
    optind = 1; //Reset before parsing again.
    int c;
    stringstream help_ss;
    while((c = getopt(argc, argv, "h:o:s:b:t:a:")) != -1) {
        switch(c) {
            case 'h':
                usage(help_ss);
                throw common::cmdline_help_exception(help_ss.str());
            case 'o':
                output_file_ = string(optarg);
                break;
            case 'a':
                bam_ = string(optarg);
                break;
            case 's':
                if (string(optarg).compare("XS") == 0){
                    strandness_ = 0;
                } else if (string(optarg).compare("RF") == 0) {
                    strandness_ = 1;
                } else if (string(optarg).compare("FR") == 0) {
                    strandness_ = 2;
                } else {
                    throw runtime_error("Unrecognized strandness argument!\n\n");
                }
                break;
            case 'b':
                barcodes_file_ = string(optarg);
                break;
            case 't':
                sites_file_ = string(optarg);
                break;
            case '?':
            default:
                usage();
                throw runtime_error("Error parsing inputs!(1)\n\n");
        }
    }
    if(argc - optind >= 1) {
        bam_ = string(argv[optind++]);
    }

    if(optind < argc || bam_ == "NA") {
        usage();
        throw runtime_error("Error parsing inputs!(2)\n\n");
    }
    if(strandness_ == -1){
        usage();
        throw runtime_error("Please supply strandness mode with '-s' option!\n\n");
    }

    cerr << "Alignment: " << bam_ << endl;
    cerr << "Barcode: " << barcodes_file_ << endl;
    cerr << "Site: " << sites_file_ << endl;
    cerr << "Output file: " << output_file_ << endl;
    cerr << "Input bam file: " << bam_ << endl;
    cerr << endl;
    return 0;
}

int USRExtractor::usage(ostream& out) {
    out << "Usage:" 
        << "\t\t" << "ISSAC unspliced read count extract" << endl;
    out << "Options:" << endl;
    out << "\t\t" << "-o FILE\tThe file to write output to. [STDOUT]" << endl;
    out << "\t\t" << "-s INT\tStrandness mode \n"
        << "\t\t\t " << "XS, use XS tags provided by aligner; RF, first-strand; FR, second-strand. REQUIRED" << endl;
    out << "\t\t" << "-t FILE\tInterested splice site" << endl;
    out << "\t\t" << "-a FILE\tbam file" << endl;
    out << "\t\t" << "-b BARCODE\tThe file containing the barcodes of interest for single cell data." << endl;
        
    out << endl;
    return 0;
}

//Get the BAM filename
string USRExtractor::get_bam() {
    return bam_;
}


//Name the junction based on the number of junctions
// in the map.
string USRExtractor::get_new_USR_name() {
    int index = USR_.size() + 1;
    stringstream name_ss;
    name_ss << "USR" << setfill('0') << setw(8) << index;
    return name_ss.str();
}


void USRExtractor::set_USR_strand_flag(bam1_t *aln, USR& j1) {
    uint32_t flag = (aln->core).flag;
    int reversed      = (flag >> 4) & 1;
    int mate_reversed = (flag >> 5) & 1;
    int first_in_pair = (flag >> 6) & 1;
    int second_in_pair= (flag >> 7) & 1;
    int bool_strandness = strandness_ - 1;  // 0 for RF, 1 for FR
    int strand_flag;
    if (first_in_pair) {
    strand_flag = !bool_strandness ^ reversed;
    } else if (second_in_pair) {
    strand_flag = !bool_strandness ^ !mate_reversed;  // second read is opposite
    } else {
    // unpaired read
    strand_flag = !bool_strandness ^ reversed;
    }
    j1.strand = string(1, strand_flag ? '+' : '-');
    return;
}

void USRExtractor::set_site(USR& j1){
    //search for the site possible covered by j1
    string chrom_strand=j1.chrom+":"+j1.strand;
    int read_len=500000;
    int start_tmp=start_pos[chrom_strand];
    for(int i=start_tmp;i<sitelist[chrom_strand].size();i++){
        if(j1.start>=sitelist[chrom_strand][i]){
            start_pos[chrom_strand] = i;
        }
        if((j1.start<sitelist[chrom_strand][i])&&(sitelist[chrom_strand][i]<j1.start+read_len)){
            j1.site_pos.insert(sitelist[chrom_strand][i]);
        }
        if(j1.start+read_len <= sitelist[chrom_strand][i]){
            break;
        }
    }
}

int USRExtractor::parse_alignment_into_unspliced_read(bam_hdr_t *header, bam1_t *aln) {
    int n_cigar = aln->core.n_cigar;
    USR j1;
    uint8_t *cb = bam_aux_get(aln, "CB");  // Get the cell barcode
    uint8_t *ub = bam_aux_get(aln, "UB");  // Get the UMI
    if(!(cb&&ub))return 0;
    j1.CB = bam_aux2Z(cb); 
    j1.UB = bam_aux2Z(ub);
    if(!((aln->core.tid)&&(aln->core.pos)))return 0;
    int chr_id = aln->core.tid;
    if((chr_id>22)||(chr_id<1))return 0;
    int read_pos = aln->core.pos;
    string chr(header->target_name[chr_id]);
    if(chr.find("chr")==string::npos)return 0;
    j1.chrom = chr;
    j1.start = read_pos; 
    set_USR_strand_flag(aln, j1); 
    set_site(j1);
    
    if(j1.site_pos.size()==0)return 0;
    
    auto it = find(barcodelist.begin(),barcodelist.end(),j1.CB);
    
    if(it==barcodelist.end()){
        return 0;
    }
    
    if(j1.site_pos.size()==0){
        return 0;
    }
    
    uint32_t *cigar = bam_get_cigar(aln);
    if(!(cigar))return 0;
    set<int> possible_site;
    for (int i = 0; i < n_cigar; ++i) {
        char op =
               bam_cigar_opchr(cigar[i]);
        int len =
               bam_cigar_oplen(cigar[i]);
        switch(op) {
            case 'N':  //for competitive read count
                j1.start = j1.start + len;
                break;
            case '=':
                for(set<int>::iterator it = j1.site_pos.begin();it!=j1.site_pos.end();++it){
                    if((*it>j1.start+1)&&(*it<j1.start+len-1)){
                        possible_site.insert(*it);
                    }
                }
                j1.start = j1.start + len;
                break;
            case 'M':
                for(set<int>::iterator it = j1.site_pos.begin();it!=j1.site_pos.end();++it){
                    if((*it>j1.start+1)&&(*it<j1.start+len-1)){
                        possible_site.insert(*it);
                    }
                }
                j1.start = j1.start + len;
                break;
            //No mismatches allowed in anchor
            case 'D':
                j1.start = j1.start + len;
                break;
            case 'X':
                j1.start = j1.start +len;
                break;
            case 'I':
            case 'S':
            case 'H':
                break;
            default:
                cerr << "Unknown cigar " << op;
                break;
        }
    }
    j1.site_pos = possible_site;
    if(j1.site_pos.size()>0){
        add_USR(j1);
    }
    return 0;
}


int USRExtractor::identify_USR_from_BAM(){
    if(!bam_.empty()) {
        //read in barcode list && sitelist
        ifstream fin;
        fin.open(barcodes_file_);
        if (!fin.is_open()) {
            cerr << "Failed to open the barcode file." << endl;
            return 1; // Or handle the error as appropriate
            }
        string line;

        while(getline(fin,line)){
        barcodelist.push_back(line);
        }
        fin.close();
        ifstream fin1;
        fin1.open(sites_file_);
        if (!fin1.is_open()) {
            cerr << "Failed to open the site file." << endl;
            return 1; // Or handle the error as appropriate
            }
        string chrom, tmp_1, strand, tmp_2, pos;
        string chrom_strand;
        while(getline(fin1,line)){
            chrom = line.substr(0,line.find(":"));
            tmp_1 = line.substr(line.find(":")+1);
            strand = tmp_1.substr(0,tmp_1.find(":"));
            tmp_2 = tmp_1.substr(tmp_1.find(":")+1);
            pos = tmp_2.substr(0,tmp_2.find(":"));
            int pos_int = stoi(pos);
            chrom_strand = chrom + ":" + strand;
            sitelist[chrom_strand].push_back(pos_int);
            sort(sitelist[chrom_strand].begin(),sitelist[chrom_strand].end());
        }
        fin1.close();
        for(int i=0;i<sitelist["chr1:+"].size();i++){
            cout<<sitelist["chr1:+"][i]<<endl;
        }

        for(int i=1;i<=22;i++){
            chrom_strand="chr"+to_string(i)+":+";
            start_pos[chrom_strand]=0;
            chrom_strand="chr"+to_string(i)+":-";
            start_pos[chrom_strand]=0;
        }
        //open BAM for reading
        samFile *in = sam_open(bam_.c_str(), "r");
        //
        cout<<"test1"<<endl;
        if(in == NULL) {
            throw runtime_error("Unable to open BAM/SAM file.\n\n");
        }
        //Load the index
        hts_idx_t *idx = sam_index_load(in, bam_.c_str());
        if(idx == NULL) {
            throw runtime_error("Unable to open BAM/SAM index."
                                " Make sure alignments are indexed\n\n");
        }
        //Get the header
        bam_hdr_t *header = sam_hdr_read(in);
        //Initialize iterator
        hts_itr_t *iter = NULL;
        //Move the iterator to the region we are interested in
        iter  = sam_itr_querys(idx, header, region_.c_str());
        if(header == NULL || iter == NULL) {
            sam_close(in);
            throw runtime_error("Unable to iterate to region within BAM.\n\n");
        }
        cout<<"test2"<<endl;
        //Initiate the alignment record
        bam1_t *aln = bam_init1();
        cout<<"test2"<<endl;
        while(sam_itr_next(in, iter, aln) >= 0) {
            parse_alignment_into_unspliced_read(header, aln);
        }
        cout<<"test3"<<endl;
        hts_itr_destroy(iter);
        cout<<"test4"<<endl;
        hts_idx_destroy(idx);
        cout<<"test5"<<endl;
        bam_destroy1(aln);
        cout<<"test6"<<endl;
        //bam_hdr_destroy(header);
        sam_close(in);
    }
    return 0;
}

//Add a junction to the junctions map
//The read_count field is the number of reads supporting the junction.
int USRExtractor::add_USR(USR j1) {

    for (int element : j1.site_pos){
        string key = j1.chrom + string(":") + to_string(element) + ":" + j1.strand + ":" + j1.CB + ":" + j1.UB;  
        USR_[key] +=1;
        cout<<key<<"\t"<<USR_[key]<<endl;
    }
    
    return 0;
}


//Print all the junctions - this function needs work
void USRExtractor::print_all_USR(ostream& out) {
    ofstream fout;
    if(output_file_ != string("NA")) {
        fout.open(output_file_.c_str());
    }
    for(map<string,int> :: iterator it = USR_.begin();
        it != USR_.end(); it++) {
        fout<<it->first<<"\t"<<it->second<<endl;
    }
    if(fout.is_open())
        fout.close();
}




//Run 'junctions extract'
int IR_extract(int argc, char *argv[]) {
    USRExtractor extract;
    try {
        extract.parse_options(argc, argv);
        extract.identify_USR_from_BAM();
        extract.print_all_USR();
    } catch(const common::cmdline_help_exception& e) {
        cerr << e.what() << endl;
        return 0;
    } catch(const runtime_error& error) {
        cerr << error.what() << endl;
        return 1;
    }
    return 0;
}

int IR_usage() {
    cout << "Usage:\t\t" << "extract nonsplit reads crossing splice sites <command> [options]" << endl;
    cout << "Command:\t" << "extract\t\tIdentify UMI-based nonsplit reads from alignments." << endl;
    cout << endl;
    return 0;
}

//Parse out subcommands under junctions
int IR_main(int argc, char *argv[]) {
    if(argc > 1) {
        string subcmd(argv[1]);
        if(subcmd == "extract") {
            return IR_extract(argc - 1, argv + 1);
        }
    }
    return IR_usage();
}





