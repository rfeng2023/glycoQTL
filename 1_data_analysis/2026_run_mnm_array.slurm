
REGION_FILE=TADB_enhanced_cis.coding_nochrX.glycoproject.bed

# # Skip header line and get gene_id for this task
# REGION_NAME=$(awk -v line=$((SLURM_ARRAY_TASK_ID + 2)) 'NR==line {print $4}' $REGION_FILE)

# Calculate how many regions to process per job
TOTAL_REGIONS=$(tail -n +2 $REGION_FILE | wc -l)
REGIONS_PER_JOB=$(( (TOTAL_REGIONS + 199) / 200 ))   # 200 arrays
START_LINE=$((SLURM_ARRAY_TASK_ID * REGIONS_PER_JOB + 2))
END_LINE=$(( (SLURM_ARRAY_TASK_ID + 1) * REGIONS_PER_JOB + 1 ))

for LINE_NUM in $(seq $START_LINE $END_LINE); do
    REGION_NAME=$(awk -v line=$LINE_NUM 'NR==line {print $4}' $REGION_FILE)
    [ -z "$REGION_NAME" ] && break  # Stop if past end of file
    
    echo "[$(date)] Processing region: $REGION_NAME (line $LINE_NUM)"
    
    # Run your command
    sos run xqtl-protocol/code/mnm_analysis/mnm_methods/mnm_regression.ipynb susie_twas \
        --name ROSMAP_Bennett_Klein_glycoQTL_twas \
        --genoFile genotype_meta/ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.genotype_by_chrom_files.txt \
        --phenoFile \
            phenotype/protein_RegionList.txt \
            phenotype/unadjusted_RegionList.txt \
            phenotype/adjusted_RegionList.txt \
        --phenoIDFile \
            phenotype/protein_PhenotypeID.txt \
            phenotype/unadjusted_PhenotypeID.txt \
            phenotype/adjusted_PhenotypeID.txt \
        --covFile  \
            covariate/protein.sorted.covariates_protein.ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.pQTL.plink_qc.prune.pca.Marchenko_PC.gz  \
            covariate/unadj_glyco.sorted.covariates_unadjusted_glyco.ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.pQTL.plink_qc.prune.pca.Marchenko_PC.gz  \
            covariate/adj_glyco.sorted.covariates_adjusted_glyco.ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.pQTL.plink_qc.prune.pca.Marchenko_PC.gz  \
        --customized_association_windows resource/TADB_enhanced_cis.coding.bed \
        --phenotype-names DLPFC_Bennett_pQTL DLPFC_Klein_gpQTL_unadjusted DLPFC_Klein_gpQTL_adjusted \
        --keep_variants GWAS/AD/AD_GWAS_NO_QC_NO_imputed_concatenate_sumstat/AD_Bellenguez_2022.sumstats.tsv.gz \
        --max-cv-variants 8000 \
        --min-twas-maf 0.05 \
        --skip_analysis_pip_cutoff 0 -1 -1 \
        --ld_reference_meta_file resource/LD_EUR/ld_meta_file.tsv \
        --cwd ${home_dir}/susie_twas/ \
        -s force \
        --region-name $REGION_NAME
done