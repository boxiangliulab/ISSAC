effect="1_1"
sample=400


/data/projects/11003054/e0950183/software/QTLtools cis --vcf /data/projects/11003054/e0950183/compare_site_intron_based/effect_size/genotype/pseudobulk_geno.vcf.gz \
                      --bed /home/users/nus/e0950183/scratch/new_simulation_result/phenotype/sample400_iso2/sample400_iso2_effect${effect}.bed.gz \
                      --cov /data/projects/11003054/e0950183/compare_site_intron_based/effect_size/PC/effect1_5.PCs \
                      --nominal 1 \
                      --normal \
                      --out /home/users/nus/e0950183/scratch/new_simulation_result/result/QTLtools/${effect}_${sample}_nominal


e=("1_1" "1_3" "1_5" "1_7" "1_9")
sample=(50 100 200)

for effect in ${e[@]};do
    bgzip /home/users/nus/e0950183/scratch/new_simulation_result/phenotype/sample400_iso2/site_based_sample400_iso2_effect${effect}.bed
    tabix -p bed /home/users/nus/e0950183/scratch/new_simulation_result/phenotype/sample400_iso2/site_based_sample400_iso2_effect${effect}.bed.gz
    /data/projects/11003054/e0950183/software/QTLtools cis --vcf /data/projects/11003054/e0950183/compare_site_intron_based/effect_size/genotype/pseudobulk_geno.vcf.gz \
                      --bed /home/users/nus/e0950183/scratch/new_simulation_result/phenotype/sample400_iso2/site_based_sample400_iso2_effect${effect}.bed.gz \
                      --cov /data/projects/11003054/e0950183/compare_site_intron_based/effect_size/PC/effect1_5.PCs \
                      --nominal 1 \
                      --normal \
                      --out /home/users/nus/e0950183/scratch/new_simulation_result/result/QTLtools/site_based_${effect}_${sample}_nominal
done
