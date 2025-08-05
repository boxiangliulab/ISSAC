ISSAC=/home/users/nus/e0950183/code/ISSAC_combine_tools_v2/QTL_mapping/build/ISSAC

sample_num=400
cell=100
title=sample${sample_num}
isoform=2
effect="1_7"
ls /home/users/nus/e0950183/scratch/new_simulation_result/junc/iso${isoform}/${effect}_cell${cell}/*.stat | cut -d '/' -f 11 | cut -d '.' -f 1 > ${cell}_${title}_${effect}.txt

output_prefix=/home/users/nus/e0950183/scratch/new_simulation_result/phenotype/sample400_iso2/${title}_iso${isoform}_effect${effect}

log_out=${title}.log

readnum=0
read_threshold=`echo $readnum/20 | bc -l`

junc_pos=/home/users/nus/e0950183/scratch/new_simulation_result/junc/iso${isoform}/${effect}_cell${cell}
sample_file=${cell}_${title}_${effect}.txt
$ISSAC pheno_group $sample_file $junc_pos $output_prefix $read_threshold $log_out

out_file=${output_prefix}.filtered
prop_file=${output_prefix}.prop
sd=0.001
na_prop=1

splice_read_count=${output_prefix}.site
$ISSAC pheno_output $splice_read_count $out_file $prop_file $sd $na_prop



####cell 20 sample 50,100,200

ISSAC=/home/users/nus/e0950183/code/ISSAC_combine_tools_v2/QTL_mapping/build/ISSAC

sample_num=50
cell=20
title=sample${sample_num}
isoform=2
effect="1_5"

readnum=0
read_threshold=`echo $readnum/20 | bc -l`
output_prefix=/data/projects/11003054/e0950183/compare_site_intron_based/effect_size/phenotype/ISSAC/${title}_iso${isoform}_effect${effect}_cell${cell}

junc_pos=/home/users/nus/e0950183/scratch/new_simulation_result/metacell/iso${isoform}_meta2_cell${cell}/${effect}
sample_file=/home/users/nus/e0950183/scratch/new_simulation_result/sample${sample_num}_meta.list
log_out=${title}.log
$ISSAC pheno_group $sample_file $junc_pos $output_prefix $read_threshold $log_out

out_file=${output_prefix}.filtered
prop_file=${output_prefix}.prop
sd=0.001
na_prop=1

splice_read_count=${output_prefix}.site
$ISSAC pheno_output $splice_read_count $out_file $prop_file $sd $na_prop
