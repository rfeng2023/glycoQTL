
REGION_NAME=${1:?region name required}
CWD=${2:-colocboost}

# 1. xQTL colocalization: snuc-eQTL (6 cell types) + unadjusted glycoQTL, no GWAS
sos run colocboost.ipynb colocboost \
    --name eQTL_gpQTL \
    --genoFile ./genotype/analysis_ready/geno_by_chrom/ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.CHR.bed \
    --phenoFile \
        ./DLPFC_Mic/analysis_ready/phenotype_preprocessing/Mic.log2cpm.region_list.txt \
        ./DLPFC_Ast/analysis_ready/phenotype_preprocessing/Ast.log2cpm.region_list.txt \
        ./DLPFC_Oli/analysis_ready/phenotype_preprocessing/Oli.log2cpm.region_list.txt \
        ./DLPFC_OPC/analysis_ready/phenotype_preprocessing/OPC.log2cpm.region_list.txt \
        ./DLPFC_Exc/analysis_ready/phenotype_preprocessing/Exc.log2cpm.region_list.txt \
        ./DLPFC_Inh/analysis_ready/phenotype_preprocessing/Inh.log2cpm.region_list.txt \
        ./gpQTL/unadjusted/analysis_ready/phenotype_preprocessing/unadjusted_RegionList.txt \
    --covFile \
        ./DLPFC_Mic/analysis_ready/covariate_preprocessing/Mic.log2cpm.Mic.rosmap_cov.ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.snuc_pseudo_bulk.related.plink_qc.extracted.pca.projected.Marchenko_PC.gz \
        ./DLPFC_Ast/analysis_ready/covariate_preprocessing/Ast.log2cpm.Ast.rosmap_cov.ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.snuc_pseudo_bulk.related.plink_qc.extracted.pca.projected.Marchenko_PC.gz \
        ./DLPFC_Oli/analysis_ready/covariate_preprocessing/Oli.log2cpm.Oli.rosmap_cov.ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.snuc_pseudo_bulk.related.plink_qc.extracted.pca.projected.Marchenko_PC.gz \
        ./DLPFC_OPC/analysis_ready/covariate_preprocessing/OPC.log2cpm.OPC.rosmap_cov.ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.snuc_pseudo_bulk.related.plink_qc.extracted.pca.projected.Marchenko_PC.gz \
        ./DLPFC_Exc/analysis_ready/covariate_preprocessing/Exc.log2cpm.Exc.rosmap_cov.ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.snuc_pseudo_bulk.related.plink_qc.extracted.pca.projected.Marchenko_PC.gz \
        ./DLPFC_Inh/analysis_ready/covariate_preprocessing/Inh.log2cpm.Inh.rosmap_cov.ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.snuc_pseudo_bulk.related.plink_qc.extracted.pca.projected.Marchenko_PC.gz \
        ./gpQTL/unadjusted/analysis_ready/covariate_preprocessing/unadj_glyco.sorted.covariates_unadjusted_glyco.ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.pQTL.plink_qc.prune.pca.Marchenko_PC.gz \
    --phenoIDFile \
        ./unused_phenotype_id.txt \
        ./unused_phenotype_id.txt \
        ./unused_phenotype_id.txt \
        ./unused_phenotype_id.txt \
        ./unused_phenotype_id.txt \
        ./unused_phenotype_id.txt \
        ./gpQTL/unadjusted/analysis_ready/phenotype_preprocessing/unadjusted_PhenotypeID.txt \
    --xqtl_coloc \
    --no-separate_gwas \
    --customized_association_windows TADB_enhanced_cis.coding.bed \
    --phenotype-names ROSMAP_Mic ROSMAP_Ast ROSMAP_Oli ROSMAP_OPC ROSMAP_Exc ROSMAP_Inh DLPFC_Klein_gpQTL_unadjusted \
    --maf 0.05 \
    --skip_analysis_pip_cutoff 0 0 0 0 0 0 -1 \
    --region-name $REGION_NAME \
    --cwd $CWD \
    -s build

# 2. xQTL + unadjusted glycoQTL + AD GWAS (Bellenguez 2022)
sos run colocboost.ipynb colocboost \
    --name unadj_glyco_AD \
    --genoFile ./genotype/analysis_ready/geno_by_chrom/ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.CHR.bed \
    --phenoFile \
        ./pQTL/ROSMAP/DLPFC/analysis_ready/phenotype_preprocessing/protein_RegionList.msk.txt \
        ./gpQTL/unadjusted/analysis_ready/phenotype_preprocessing/unadjusted_RegionList.msk.txt \
    --covFile \
        ./pQTL/ROSMAP/DLPFC/analysis_ready/covariate_preprocessing/protein.sorted.covariates_protein.ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.pQTL.plink_qc.prune.pca.Marchenko_PC.gz \
        ./gpQTL/unadjusted/analysis_ready/covariate_preprocessing/unadj_glyco.sorted.covariates_unadjusted_glyco.ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.pQTL.plink_qc.prune.pca.Marchenko_PC.gz \
    --phenoIDFile \
        ./pQTL/ROSMAP/DLPFC/analysis_ready/phenotype_preprocessing/protein_PhenotypeID.txt \
        ./gpQTL/unadjusted/analysis_ready/phenotype_preprocessing/unadjusted_PhenotypeID.txt \
    --ld-meta-data ./ADSP_R4_EUR_LD/ld_meta_file.tsv \
    --gwas-meta-data ./ADGWAS/GWAS_sumstat_meta_AD_Bellenguez_tabix_msk.tsv \
    --xqtl_coloc \
    --separate_gwas \
    --customized_association_windows TADB_enhanced_cis.coding.bed \
    --phenotype-names DLPFC_Bennett_pQTL DLPFC_Klein_gpQTL_unadjusted \
    --maf 0.05 \
    --skip_analysis_pip_cutoff 0 -1 \
    --region-name $REGION_NAME \
    --cwd $CWD \
    -s build

# 3. xQTL + adjusted glycoQTL + AD GWAS (Bellenguez 2022)
sos run colocboost.ipynb colocboost \
    --name adj_glyco_AD \
    --genoFile ./genotype/analysis_ready/geno_by_chrom/ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.CHR.bed \
    --phenoFile \
        ./pQTL/ROSMAP/DLPFC/analysis_ready/phenotype_preprocessing/protein_RegionList.msk.txt \
        ./gpQTL/adjusted/analysis_ready/phenotype_preprocessing/adjusted_RegionList.msk.txt \
    --covFile \
        ./pQTL/ROSMAP/DLPFC/analysis_ready/covariate_preprocessing/protein.sorted.covariates_protein.ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.pQTL.plink_qc.prune.pca.Marchenko_PC.gz \
        ./gpQTL/adjusted/analysis_ready/covariate_preprocessing/adj_glyco.sorted.covariates_adjusted_glyco.ROSMAP_NIA_WGS.leftnorm.bcftools_qc.plink_qc.pQTL.plink_qc.prune.pca.Marchenko_PC.gz \
    --phenoIDFile \
        ./pQTL/ROSMAP/DLPFC/analysis_ready/phenotype_preprocessing/protein_PhenotypeID.txt \
        ./gpQTL/adjusted/analysis_ready/phenotype_preprocessing/adjusted_PhenotypeID.txt \
    --ld-meta-data ./ADSP_R4_EUR_LD/ld_meta_file.tsv \
    --gwas-meta-data ./ADGWAS/GWAS_sumstat_meta_AD_Bellenguez_tabix_msk.tsv \
    --xqtl_coloc \
    --separate_gwas \
    --customized_association_windows TADB_enhanced_cis.coding.bed \
    --phenotype-names DLPFC_Bennett_pQTL DLPFC_Klein_gpQTL_adjusted \
    --maf 0.05 \
    --skip_analysis_pip_cutoff 0 -1 \
    --region-name $REGION_NAME \
    --cwd $CWD \
    -s build
