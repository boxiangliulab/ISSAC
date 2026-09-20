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
    int first_in_pair = (flag >> 6) & 1;
    int second_in_pair= (flag >> 7) & 1;
    int bool_strandness = strandness_ - 1;  // 0 for RF, 1 for FR
    int strand_flag;
    if (first_in_pair) {
    strand_flag = !bool_strandness ^ reversed;
    } else if (second_in_pair) {
    strand_flag = !bool_strandness ^ !reversed;  // second read is opposite
    } else {
    // unpaired read
    strand_flag = !bool_strandness ^ reversed;
    }
    j1.strand = string(1, strand_flag ? '+' : '-');
    return;
}

void USRExtractor::set_site(USR& j1, int alignment_end) {

    string chrom_strand = j1.chrom + ":" + j1.strand;

    int start_tmp = start_pos[chrom_strand];

    for (int i = start_tmp;
         i < sitelist[chrom_strand].size();
         ++i) {

        int site = sitelist[chrom_strand][i];

        if (j1.start >= site) {
            start_pos[chrom_strand] = i;
        }

        if ((j1.start < site) &&
            (site <= alignment_end)) {
            j1.site_pos.insert(site);
        }

        if (alignment_end < site) {
            break;
        }
    }
}

int USRExtractor::parse_alignment_into_unspliced_read(
    bam_hdr_t *header,
    bam1_t *aln) {
    int n_cigar = aln->core.n_cigar;
    USR j1;

    // Skip unmapped reads
    if (aln->core.flag & BAM_FUNMAP) return 0;

    // Get cell barcode and UMI
    uint8_t *cb = bam_aux_get(aln, "CB");
    uint8_t *ub = bam_aux_get(aln, "UB");
    if (!(cb && ub)) return 0;

    const char *cb_str = bam_aux2Z(cb);
    const char *ub_str = bam_aux2Z(ub);
    if (!(cb_str && ub_str)) return 0;

    j1.CB = cb_str;
    j1.UB = ub_str;

    // Filter barcode early
    if (barcodeset.find(j1.CB) == barcodeset.end()) {
        return 0;
    }

    int chr_id = aln->core.tid;
    if (chr_id < 0 || chr_id >= header->n_targets) return 0;

    int read_pos = aln->core.pos;
    string chr(header->target_name[chr_id]);

    // Require chromosome names beginning with "chr"
    if (chr.rfind("chr", 0) != 0) return 0;

    // Restrict to autosomes chr1-chr22
    string chr_suffix = chr.substr(3);

    try {
        int chr_num = stoi(chr_suffix);

        // Ensure the entire suffix is numeric
        if (chr_suffix != to_string(chr_num)) return 0;

        if (chr_num < 1 || chr_num > 22) return 0;
    } catch (const std::exception&) {
        return 0;
    }

    j1.chrom = chr;
    j1.start = read_pos;

    set_USR_strand_flag(aln, j1);

    // Actual reference end of this alignment
    int alignment_end = bam_endpos(aln);

    // Find candidate sites only within the alignment span
    set_site(j1, alignment_end);

    if (j1.site_pos.empty()) return 0;

    uint32_t *cigar = bam_get_cigar(aln);
    if (!cigar) return 0;

    set<int> possible_site;

    for (int i = 0; i < n_cigar; ++i) {

        char op = bam_cigar_opchr(cigar[i]);
        int len = bam_cigar_oplen(cigar[i]);

        switch (op) {

            case 'N':
                // Reference skip
                j1.start += len;
                break;

            case '=':
            case 'M':
                for (set<int>::iterator it = j1.site_pos.begin();
                     it != j1.site_pos.end(); ++it) {

                    if ((*it > j1.start + 1) &&
                        (*it < j1.start + len - 1)) {
                        possible_site.insert(*it);
                    }
                }

                j1.start += len;
                break;

            // No mismatches/deletions allowed in anchor
            case 'D':
            case 'X':
                j1.start += len;
                break;

            // Do not consume reference
            case 'I':
            case 'S':
            case 'H':
            case 'P':
                break;

            default:
                cerr << "Unknown CIGAR operation: " << op << endl;
                break;
        }
    }

    j1.site_pos = possible_site;

    if (!j1.site_pos.empty()) {
        add_USR(j1);
    }

    return 0;
}


int USRExtractor::identify_USR_from_BAM() {

    if (bam_.empty()) {
        return 0;
    }

    // ============================================================
    // 1. Read barcode list
    // ============================================================
    fstream fin(barcodes_file_);
    if (!fin.is_open()) {
    cerr << "Failed to open the barcode file: "
         << barcodes_file_ << endl;
    return 1;}
    string line;
    while (getline(fin, line)) {
    if (!line.empty()) {
        barcodelist.push_back(line);
    }}
    fin.close();

// Build barcode set for fast lookup during BAM processing
    barcodeset.insert(barcodelist.begin(), barcodelist.end());

    // ============================================================
    // 2. Read splice-site list
    // ============================================================
    ifstream fin1(sites_file_);

    if (!fin1.is_open()) {
        cerr << "Failed to open the site file: "
             << sites_file_ << endl;
        return 1;
    }

    string chrom;
    string strand;
    string pos;
    string chrom_strand;

    while (getline(fin1, line)) {

        if (line.empty()) {
            continue;
        }

        // Expected format:
        // chr:strand:position[:...]
        //
        // Example:
        // chr1:+:123456

        size_t p1 = line.find(':');

        if (p1 == string::npos) {
            cerr << "Invalid site entry: " << line << endl;
            continue;
        }

        size_t p2 = line.find(':', p1 + 1);

        if (p2 == string::npos) {
            cerr << "Invalid site entry: " << line << endl;
            continue;
        }

        size_t p3 = line.find(':', p2 + 1);

        chrom = line.substr(0, p1);
        strand = line.substr(p1 + 1, p2 - p1 - 1);

        if (p3 == string::npos) {
            pos = line.substr(p2 + 1);
        } else {
            pos = line.substr(p2 + 1, p3 - p2 - 1);
        }

        // Only accept valid strands
        if (strand != "+" && strand != "-") {
            cerr << "Invalid strand in site entry: "
                 << line << endl;
            continue;
        }

        int pos_int;

        try {
            pos_int = stoi(pos);
        }
        catch (const std::exception&) {
            cerr << "Invalid genomic position in site entry: "
                 << line << endl;
            continue;
        }

        chrom_strand = chrom + ":" + strand;

        sitelist[chrom_strand].push_back(pos_int);
    }

    fin1.close();


    // ============================================================
    // 3. Sort site lists once after all sites have been loaded
    // ============================================================
    for (auto &entry : sitelist) {
        sort(entry.second.begin(), entry.second.end());
    }


    // ============================================================
    // 4. Initialize site-search positions
    // ============================================================
    for (int i = 1; i <= 22; ++i) {

        chrom_strand = "chr" + to_string(i) + ":+";
        start_pos[chrom_strand] = 0;

        chrom_strand = "chr" + to_string(i) + ":-";
        start_pos[chrom_strand] = 0;
    }


    // ============================================================
    // 5. Open BAM/SAM file
    // ============================================================
    samFile *in = sam_open(bam_.c_str(), "r");

    if (in == NULL) {
        throw runtime_error(
            "Unable to open BAM/SAM file: " + bam_ + "\n"
        );
    }


    // ============================================================
    // 6. Load BAM index
    // ============================================================
    hts_idx_t *idx = sam_index_load(in, bam_.c_str());

    if (idx == NULL) {

        sam_close(in);

        throw runtime_error(
            "Unable to open BAM/SAM index. "
            "Make sure the alignment file is coordinate-sorted "
            "and indexed.\n"
        );
    }


    // ============================================================
    // 7. Read BAM header
    // ============================================================
    bam_hdr_t *header = sam_hdr_read(in);

    if (header == NULL) {

        hts_idx_destroy(idx);
        sam_close(in);

        throw runtime_error(
            "Unable to read BAM/SAM header.\n"
        );
    }


    // ============================================================
    // 8. Create iterator for requested region
    // ============================================================
    hts_itr_t *iter =
        sam_itr_querys(idx, header, region_.c_str());

    if (iter == NULL) {

        bam_hdr_destroy(header);
        hts_idx_destroy(idx);
        sam_close(in);

        throw runtime_error(
            "Unable to iterate over region: "
            + region_ + "\n"
        );
    }


    // ============================================================
    // 9. Allocate alignment record
    // ============================================================
    bam1_t *aln = bam_init1();

    if (aln == NULL) {

        hts_itr_destroy(iter);
        bam_hdr_destroy(header);
        hts_idx_destroy(idx);
        sam_close(in);

        throw runtime_error(
            "Unable to allocate BAM alignment record.\n"
        );
    }


    // ============================================================
    // 10. Process alignments
    // ============================================================
    while (sam_itr_next(in, iter, aln) >= 0) {

        parse_alignment_into_unspliced_read(
            header,
            aln
        );
    }


    // ============================================================
    // 11. Clean up
    // ============================================================
    bam_destroy1(aln);
    hts_itr_destroy(iter);
    bam_hdr_destroy(header);
    hts_idx_destroy(idx);
    sam_close(in);


    return 0;
}

//Add a junction to the junctions map
//The read_count field is the number of reads supporting the junction.
int USRExtractor::add_USR(const USR& j1) {

    for (int element : j1.site_pos) {

        string key =
            j1.chrom + ":" +
            to_string(element) + ":" +
            j1.strand + ":" +
            j1.CB + ":" +
            j1.UB;

        USR_[key] += 1;
    }

    return 0;
}


void USRExtractor::print_all_USR() {

    if (output_file_ == "NA") {
        return;
    }

    ofstream fout(output_file_);

    if (!fout.is_open()) {
        throw runtime_error(
            "Unable to open output file: " + output_file_
        );
    }

    for (const auto& entry : USR_) {
        fout << entry.first
             << "\t"
             << entry.second
             << "\n";
    }
}


int IR_extract(int argc, char *argv[]) {

    USRExtractor extract;

    try {

        extract.parse_options(argc, argv);

        int status = extract.identify_USR_from_BAM();

        if (status != 0) {
            return status;
        }

        extract.print_all_USR();

    }
    catch (const common::cmdline_help_exception& e) {
        cerr << e.what() << endl;
        return 0;
    }
    catch (const runtime_error& error) {
        cerr << error.what() << endl;
        return 1;
    }

    return 0;
}


int IR_main(int argc, char *argv[]) {
    if(argc > 1) {
        string subcmd(argv[1]);
        if(subcmd == "extract") {
            return IR_extract(argc - 1, argv + 1);
        }
    }
    return IR_usage();
}

int IR_usage() {
    cout << "Usage:\t\t"
         << "extract nonsplit reads crossing splice sites <command> [options]"
         << endl;
    cout << "Command:\t"
         << "extract\t\tIdentify UMI-based nonsplit reads from alignments."
         << endl;
    cout << endl;
    return 0;
}

