import os
 
os.environ['R_HOME']='/home/users/nus/e0950183/ISSAC_env/bin/R'
import sys
import scReadSim.Utility as Utility
#import scReadSim.GenerateSyntheticCount as GenerateSyntheticCount
import scReadSim.scATAC_GenerateBAM as scATAC_GenerateBAM
#import pkg_resources

import pandas as pd
import numpy as np
from scIsoSim import scIsoSim

INPUT_cells_barcode_file = "/data/projects/11003054/e0950183/compare_site_intron_based/JP_RIK_H002_CD4+_T_naive_simulation/tmp_scisosim/JP_RIK_H002.CD4+_T_naive.barcodes.tsv"
filename = "10X_demo"
INPUT_bamfile = "/data/projects/11003054/e0950183/compare_site_intron_based/JP_RIK_H002_CD4+_T_naive_simulation/tmp_scisosim/JP_RIK_H002.CD4+_T_naive.bam"
samtools_directory="/home/users/nus/e0950183/ISSAC_env/bin" 
bedtools_directory="/home/users/nus/e0950183/ISSAC_env/bin"
seqtk_directory="/home/users/nus/e0950183/ISSAC_env/bin"
gffread_dir = "/home/users/nus/e0950183/ISSAC_env/bin"

outdirectory = "/data/projects/11003054/e0950183/compare_site_intron_based/JP_RIK_H002_CD4+_T_naive_simulation/tmp_scisosim"
INPUT_genome_size_file = "/data/projects/11003054/e0950183/compare_site_intron_based/common_file/hg38.chrom.sizes"
INPUT_genome_annotation = "/data/projects/11003054/e0950183/compare_site_intron_based/common_file/gencode.v32.primary_assembly.annotation.gff3"


gene_bedfile = "/data/projects/11003054/e0950183/compare_site_intron_based/JP_RIK_H002_CD4+_T_naive_simulation/tmp_scisosim/scReadSim.Gene.bed"

filename="CD4+_T_naive"

UMI_gene_count_mat_filename = "%s.gene.countmatrix" % filename

isoform_annotation_file = INPUT_genome_annotation


referenceGenome_name = "GRCh38.primary_assembly.genome"
referenceGenome_dir = "/data/projects/11003054/e0950183/compare_site_intron_based/common_file"  # may use users' own path
referenceGenome_file = "%s/%s.fa" % (referenceGenome_dir, referenceGenome_name)

UMI_count_mat_filename = "CD4+_T_naive.gene.countmatrix"
synthetic_cell_label_file = "/data/projects/11003054/e0950183/compare_site_intron_based/JP_RIK_H002_CD4+_T_naive_simulation/tmp_scisosim/" + UMI_count_mat_filename + ".scDesign2Simulated.CellTypeLabel.txt"
OUTPUT_cells_barcode_file = "synthetic_cell_barcode.txt"

synthetic_countmat_file="CD4+_T_naive.gene.countmatrix.scDesign2Simulated.txt"
UMI_count_mat_df = pd.read_csv("%s" % (outdirectory + "/" + synthetic_countmat_file), header=None, delimiter="\t")
UMI_count_mat = UMI_count_mat_df.iloc[:,1:].to_numpy() # remove feature names
random_cellbarcode_list = pd.read_csv(outdirectory + "/" + OUTPUT_cells_barcode_file, header=None, delimiter="\t").to_numpy()

effect="1_1"
out_isoform="/data/projects/11003054/e0950183/compare_site_intron_based/effect_size/fq/"+effect


for i in range(1,401):
    random_intergars=np.random.randint(1,145,size=100)
    sub_cellbarcode = random_cellbarcode_list[random_intergars]
    sub_UMI_count_mat = UMI_count_mat[:,random_intergars]
    read_bedfile_prename = "Isoform"+str(i)
    isoform_proportion_file = "/data/projects/11003054/e0950183/compare_site_intron_based/effect_size/iso_prop/"+effect+"/iso_prop_2_S"+str(i)+".txt"
    scIsoSim.generateIsoformBED(bed_file=gene_bedfile,
                                sub_count_mat=sub_UMI_count_mat,
                                sub_cell_barcode=sub_cellbarcode,
                                INPUT_bamfile=INPUT_bamfile,
                                isoform_proportion_file= isoform_proportion_file,
                                trancriptome_file=outdirectory + "/" + "transcripts.fa",
                                outdirectory=out_isoform,
                                read_bedfile_prename=read_bedfile_prename,
                                UMI_tag='UB:Z',
                                sequencing_protocol = "10x",
                                library_type="5primePairedEndReads",
                                read_len=150,
                                UMI_len=10)
    
    scATAC_GenerateBAM.scATAC_BED2FASTQ(bedtools_directory=bedtools_directory, 
                              seqtk_directory=seqtk_directory, 
                              referenceGenome_file=outdirectory + "/" + "transcripts.fa", 
                              outdirectory=out_isoform, 
                              BED_filename_combined=read_bedfile_prename, 
                              synthetic_fastq_prename = read_bedfile_prename)

###Convert bed file to fastq files
    scIsoSim.incorporate_CB_UMI_TSO(input_file = out_isoform + "/" + read_bedfile_prename + ".read1.bed2fa.sorted.fq")
