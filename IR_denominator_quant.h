#ifndef IR_DENOMINATOR_QUANT_H
#define IR_DENOMINATOR_QUANT_H

#include <algorithm>
#include <iomanip>
#include <iostream>
#include "bedFile.h"
#include "htslib/sam.h"
#include <unordered_map>
#include <unordered_set>

using namespace std;

//Format of an unspliced read
struct USR : BED {
    //Number of reads supporting the junction
    unsigned int read_count;
    string chrom;
    //Has the junction been added to the map
    bool added; 
    //Does the junction satisfy the min anchor requirement
    //Minimum anchor on either side can be supported by different
    //reads, only junctions anchored on both sides are reported.
    int start;
    int end;
    string CB;
    string UB;
    string strand;
    string name;
    set<int> site_pos;

    USR() {
        chrom = "NA";
        start = 0;
        end = 0;
        read_count = 0;
        strand = "NA";
        name = "NA";
        added = false;
        CB = "NA";
        UB = "NA";
    }
    USR(string chrom1, CHRPOS start1, CHRPOS end1,
             string strand1) {
        chrom = chrom1;
        start = start1;
        end = end1;
        read_count = 0;
        strand = strand1;
        added = false;
        CB = "NA";
        UB = "NA";
    }
    //Print junction
    void print(ostream& out) const {
        out << chrom <<
            "\t" << start << "\t" << end <<
            "\t" << name << "\t" << read_count << "\t" << strand <<
            "\t" << CB << "\t" << UB << "\t" << endl;
    }
};


//The class that deals with creating the unspliced read count
class USRExtractor {
    private:
        //Alignment file
        string bam_;
        //The key is "chr:site_pos:strand:CB:UB"
        //The value is an object of type USR(see above)
        map<string, int> USR_;
        //File to write output to - optional, write to STDOUT by default
        string output_file_;
        //File to barcodes - optional
        string barcodes_file_;
        //File to sites - optional
        string sites_file_;
        //strandness of data; 0 = unstranded, 1 = RF, 2 = FR
        int strandness_;
        //tag used in BAM to denote strand, default "XS"
        string region_;
        string strand_tag_;
        vector<string> barcodelist;
        map<string,vector<int>> sitelist;
        map<string,int> start_pos;
    public:
        //Default constructor
        USRExtractor() {
            strandness_ = -1;
            strand_tag_ = "XS";
            bam_ = "NA";
            output_file_ = "NA";
            region_ = ".";
        }
        USRExtractor(string bam1, int strandness1, string strand_tag1, string region1) : 
            bam_(bam1), strandness_(strandness1), strand_tag_(strand_tag1){
            output_file_ = "NA";
        }
        //Name the junction based on the number of junctions
        // in the map.
        string get_new_USR_name();  
        //Parse command-line options for this tool
        int parse_options(int argc, char *argv[]); 
        //Print default usage
        int usage(ostream& out = cerr);  
        //Identify exon-exon junctions
        int identify_USR_from_BAM();
        //Print all the junctions
        void print_all_USR(ostream& out = cout);  
        //Get the BAM filename
        string get_bam();
        //Parse the alignment into the junctions map
        int parse_alignment_into_unspliced_read(bam_hdr_t *header, bam1_t *aln);
        //Add a junction to the junctions map
        int add_USR(USR j1);  
        //Get the strand from bitwise flag
        void set_site(USR& j1);
        //set the interested site for each read count
        void set_USR_strand_flag(bam1_t *aln, USR& j1);
};


int IR_extract(int argc, char *argv[]);
int IR_main(int argc, char *argv[]);

int IR_usage();

#endif
