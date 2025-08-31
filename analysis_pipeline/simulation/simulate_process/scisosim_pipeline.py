import sys, os
os.environ['R_HOME']='/home/e0950183/miniconda3/envs/my_java_env/lib/R'

import scReadSim.Utility as Utility
import scReadSim.GenerateSyntheticCount as GenerateSyntheticCount
import scReadSim.scRNA_GenerateBAM as scRNA_GenerateBAM
import pkg_resources

import pandas as pd
import numpy as np
from scIsoSim import scIsoSim

INPUT_cells_barcode_file = "/home/e0950183/project/compare_splice_based_intron_based/simulated_data/intermediate_data/JP_RIK_H002.CD4+_T_naive.barcodes.tsv"
filename = "10X_demo"
INPUT_bamfile = "/data/projects/11003054/e0950183/compare_site_intron_based/JP_RIK_H002_CD4+_T_naive_simulation/tmp_scisosim/JP_RIK_H002.CD4+_T_naive.bam"
samtools_directory="/home/e0950183/miniconda3/envs/my_java_env/bin" 
bedtools_directory="/home/e0950183/miniconda3/envs/my_java_env/bin"
seqtk_directory="/home/e0950183/miniconda3/envs/my_java_env/bin"
gffread_dir = "/home/e0950183/miniconda3/envs/my_java_env/bin"

outdirectory = "/home/e0950183/project/compare_splice_based_intron_based/simulated_data/intermediate_data"
INPUT_genome_size_file = "/home/e0950183/project/compare_splice_based_intron_based/reference_data/scisosim_commonfile/hg38.chrom.sizes"
INPUT_genome_annotation = "/home/e0950183/project/compare_splice_based_intron_based/reference_data/scisosim_commonfile/gencode.v32.annotation.gff3"

###generate features
Utility.scRNA_CreateFeatureSets(INPUT_bamfile=INPUT_bamfile,
    samtools_directory=samtools_directory,
    bedtools_directory=bedtools_directory, outdirectory=outdirectory,
    genome_annotation=INPUT_genome_annotation, genome_size_file=INPUT_genome_size_file)




###generate real count matrices

gene_bedfile = outdirectory + "/" + "scReadSim.Gene.bed"

filename="CD4+_T_naive"

UMI_gene_count_mat_filename = "%s.gene.countmatrix" % filename

Utility.scRNA_bam2countmat_paral(cells_barcode_file=INPUT_cells_barcode_file,
                                 bed_file=gene_bedfile, INPUT_bamfile=INPUT_bamfile,
                                 outdirectory=outdirectory, count_mat_filename=UMI_gene_count_mat_filename,
                                 UMI_modeling=True, UMI_tag = "UB:Z", n_cores=8)


###Sythetic count matrix simulation
GenerateSyntheticCount.scRNA_GenerateSyntheticCount(count_mat_filename=UMI_gene_count_mat_filename,
                                                    directory=outdirectory,
                                                    outdirectory=outdirectory)

###Prepare ground-truth isoform annotations
isoform_annotation_file = INPUT_genome_annotation

referenceGenome_name = "GRCh38.primary_assembly.genome"
referenceGenome_dir = "/home/e0950183/project/compare_splice_based_intron_based/reference_data/scisosim_commonfile"  # may use users' own path
referenceGenome_file = "%s/%s.fa" % (referenceGenome_dir, referenceGenome_name)

###Prepare trascriptome
UMI_count_mat_filename = "CD4+_T_naive.gene.countmatrix"
synthetic_cell_label_file = outdirectory+ "/" + UMI_count_mat_filename + ".scDesign2Simulated.CellTypeLabel.txt"
OUTPUT_cells_barcode_file = "synthetic_cell_barcode.txt"
scIsoSim.prepare_scIsoSim(outdirectory=outdirectory,
                        synthetic_cell_label_file=synthetic_cell_label_file,
                        gffread_dir=gffread_dir,
                        referenceGenome_file=referenceGenome_file,
                        isoform_annotation_file=isoform_annotation_file,
                        OUTPUT_cells_barcode_file=OUTPUT_cells_barcode_file,
                        CB_len=16)


###Prepare ground-truth isoform proportions
isoform_proportion_file = "IsoformProportion.txt"
###Genes
gene_selected=pd.read_csv("/data/zhangyuntian/project/scSplice/data/output/selected_Genes.txt",header=None,delimiter="\t")
gene_ind_selected=gene_selected.iloc[:,1].to_numpy()
gene_ind_selected=gene_ind_selected -1
scIsoSim.generateIsoformProp(bed_file=gene_bedfile,
                    isoform_annotation_file=isoform_annotation_file,
                    outdirectory=outdirectory,
                    isoform_proportion_file=isoform_proportion_file,
                    gene_ind=gene_ind_selected)


###Output synthetic read
synthetic_countmat_file="CD4+_T_naive.gene.countmatrix.scDesign2Simulated.txt"
UMI_count_mat_df = pd.read_csv("%s" % (outdirectory + "/" + synthetic_countmat_file), header=None, delimiter="\t")
UMI_count_mat = UMI_count_mat_df.iloc[:,1:].to_numpy() # remove feature names
random_cellbarcode_list = pd.read_csv(outdirectory + "/" + OUTPUT_cells_barcode_file, header=None, delimiter="\t").to_numpy()

read_bedfile_prename = "IsoformSyntheticRead"


scIsoSim.generateIsoformBED(bed_file=gene_bedfile,
                                sub_count_mat=UMI_count_mat,
                                sub_cell_barcode=random_cellbarcode_list,
                                INPUT_bamfile=INPUT_bamfile,
                                isoform_proportion_file=outdirectory + "/" + isoform_proportion_file,
                                trancriptome_file=outdirectory + "/" + "transcripts.fa",
                                outdirectory=outdirectory,
                                read_bedfile_prename=read_bedfile_prename,
                                UMI_tag='UB:Z',
                                sequencing_protocol = "10x",
                                library_type="5prime",
                                read_len=126,
                                UMI_len=10)


###Convert bed file to fastq files

scRNA_GenerateBAM.scRNA_BED2FASTQ(bedtools_directory=bedtools_directory,
                                  seqtk_directory=seqtk_directory,
                                  referenceGenome_file=outdirectory + "/" + "transcripts.fa",
                                  outdirectory=outdirectory,
                                  BED_filename_combined=read_bedfile_prename,
                                  synthetic_fastq_prename = read_bedfile_prename)
