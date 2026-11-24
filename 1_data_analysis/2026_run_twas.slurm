
export PATH="$PIXI_ENV/bin:$PATH"  # PIXI_ENV: pixi environment providing pecotmr and SoS
source ~/.bashrc

TMPDIR=tmp/slurm_${SLURM_JOB_ID}
mkdir -p $TMPDIR
export TMPDIR

trap "rm -rf $TMPDIR" EXIT

# Region list: header line, then one LD block name per line in column 1 (e.g. chr7_11259571_12595836)
REGION_FILE=missing_regions1.txt

LINES_PER_JOB=5
START_LINE=$((SLURM_ARRAY_TASK_ID * LINES_PER_JOB + 1))
END_LINE=$((START_LINE + LINES_PER_JOB - 1))

for LINE_NUM in $(seq $START_LINE $END_LINE); do
    # Skip header line, so add 1 to LINE_NUM
    REGION_NAME=$(awk -v line=$((LINE_NUM + 1)) 'NR==line {print $1}' $REGION_FILE)
    [ -z "$REGION_NAME" ] && break  # Stop if past end of file

    echo "[$(date)] Processing region: $REGION_NAME"
    sos run xqtl-protocol/code/pecotmr_integration/twas_ctwas.ipynb twas \
        --name ROSMAP_Bennett_Klein_glycoQTL \
        --gwas_meta_data RSS_QC_nonimputed_result/gwas_meta.tsv \
        --ld_meta_data resource/LD_EUR/ld_meta_file.tsv \
        --xqtl_meta_data ROSMAP_gpQTL.meta_twas.tsv \
        --xqtl_type_table data_type_table.txt \
        --rsq_pval_cutoff 0.05 --rsq_cutoff 0.01 \
        --cwd RSS_QC_nonimputed_result/ -s build \
        --regions resource/EUR_LD_blocks.bed \
        --region-name $REGION_NAME
done
