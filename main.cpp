#include <iostream>
#include <getopt.h>
#include <stdexcept>
#include "common.h"
#include "gtf_parser.h"
#include "junctions_extractor.h"
#include "eigen_decom.h"
#include "IR_denominator_quant.h"
#include "pheno.h"
using namespace std;

//Usage for junctions subcommands
int junctions_usage(ostream &out = cout) {
    out << "Usage:\t\t" << "junctools junctions <command> [options]" << endl;
    out << "Command:\t" << "extract\t\tIdentify UMI-based exon-exon junctions from alignments." << endl;
    out << endl;
    return 0;
}

//Run 'junctions extract'
int junctions_extract(int argc, char *argv[]) {
    JunctionsExtractor extract;
    try {
        extract.parse_options(argc, argv);
        extract.identify_junctions_from_BAM();
        extract.print_all_junctions();
    } catch(const common::cmdline_help_exception& e) {
        cerr << e.what() << endl;
        return 0;
    } catch(const runtime_error& error) {
        cerr << error.what() << endl;
        return 1;
    }
    return 0;
}

//Parse out subcommands under junctions
int junctions_main(int argc, char *argv[]) {
    if(argc > 1) {
        string subcmd(argv[1]);
        if(subcmd == "extract") {
            return junctions_extract(argc - 1, argv + 1);
        }
    }
    return junctions_usage();
}


int usage() {
    cerr << "Usage:" 
        << "\t\t" << "ISSAC <command> [options]" << endl;
    cerr << "Command:\t" << "Integrative single-cell splicing analysis and QTL caller" << endl;
    cerr << endl;
    return 0;
}

//Everything starts here
int main(int argc, char* argv[]) {
    if(argc > 1) {
        string subcmd(argv[1]);
        if(subcmd == "junctools") {
            return junctions_main(argc - 1, argv + 1);
        }
        else if(subcmd == "juncstat") {
            return juncstat(argc - 1, argv + 1);
        }
        else if(subcmd == "IR") {
            return IR_main(argc - 1, argv + 1);
        }
        else if(subcmd == "IR_combine") {
            return IR_combine(argc - 1, argv + 1);
        }
        else if(subcmd == "model") {
            return model(argc - 1, argv + 1);
        }
        else if(subcmd == "QTL"){
            return QTL(argc - 1, argv+1);
        }
        else if(subcmd == "trans_QTL"){
            return trans_QTL(argc - 1, argv+1);
        }
        else if(subcmd == "DS"){
            return DS(argc - 1, argv+1);
        }
        else if(subcmd == "pheno_group"){
            return pheno_group(argc - 1, argv+1);
        }
        else if(subcmd == "pheno_output"){
            return pheno_output(argc - 1, argv+1);
        }
        else if(subcmd == "rare_variant"){
            return rare_mapping(argc - 1, argv+1);
        }
    }
    return usage();
}
