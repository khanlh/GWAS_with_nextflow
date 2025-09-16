nextflow.enable.dsl=2

params.outdir = "./results"
params.p_threshold = 5e-8 // P-value threshold to filter SNPs of interest (adjustable)

// ----------------------------
// Simulate genotype
// ----------------------------
process simulate_genotype {
    publishDir "${params.outdir}/genotype", mode: 'copy'

    output:
    path "geno.*", emit: geno

    script:
    """
    plink2 --dummy 50 100 0.1 --make-bed --out geno
    """
}

// ----------------------------
// Genotype QC
// ----------------------------
process qc_genotype {
    publishDir "${params.outdir}/qc", mode: 'copy'

    input:
    path geno_files

    output:
    path "geno_qc.*", emit: qc

    script:
    """
    plink2 --bfile geno --geno 0.1 --mind 0.1 --maf 0.05 --make-bed --out geno_qc
    """
}

// ----------------------------
// Association test
// ----------------------------
process assoc_test {
    publishDir "${params.outdir}/assoc", mode: 'copy'

    input:
    path qc_files

    output:
    path "assoc_results.*", emit: gwas

    script:
    """
    awk '{ if(NR>1) print \$1, \$2, rand() }' geno_qc.fam > pheno.txt

    plink2 --bfile geno_qc \\
           --pheno pheno.txt \\
           --glm allow-no-covars \\
           --out assoc_results
    """
}

// ----------------------------
// Plot Manhattan
// ----------------------------
process plot_manhattan {
    container 'gwas-demo:latest'
    publishDir "${params.outdir}/manhattan", mode: 'copy'
    input:
        path glm_files
    output:
        path "*.manhattan.png"
    script:
    """
    set -euo pipefail

    for f in *.glm.linear; do
        echo "Plotting Manhattan for \$f"

        cat <<EOF > plot_manhattan.R
library(qqman)
dat <- read.table(file="\$f", header=TRUE, sep="\\t", check.names=FALSE, comment.char="")
if ("#CHROM" %in% colnames(dat)) colnames(dat)[colnames(dat) == "#CHROM"] <- "CHR"

if ("P" %in% colnames(dat)) {
    dat\\\$PVAL <- as.numeric(dat\\\$P)
} else if ("T_STAT" %in% colnames(dat)) {
    dat\\\$T_STAT <- as.numeric(dat\\\$T_STAT)
    dat <- dat[!is.na(dat\\\$T_STAT), ]
    dat\\\$PVAL <- 2 * pnorm(-abs(dat\\\$T_STAT))
} else {
    stop("T_STAT or P column not found")
}

png(paste0(gsub(".*/", "", "\$f"), ".manhattan.png"), width=1000, height=600)
manhattan(dat, chr="CHR", bp="POS", snp="ID", p="PVAL", main="GWAS Manhattan Plot", suggestiveline = -log10(1e-5), genomewideline = -log10(${params.p_threshold}))
dev.off()
EOF

        Rscript plot_manhattan.R
    done
    """
}

// ----------------------------
// Filter SNPs of interest
// ----------------------------
process filter_snp_of_interest {
    publishDir "${params.outdir}/snp_of_interest", mode: 'copy'

    input:
    path glm_files

    output:
    path "snp_of_interest.txt", emit: snp_interest

    script:
    """
    # Filter SNPs based on p-value threshold
    awk -v threshold=${params.p_threshold} 'NR>1 && \$15<threshold {print \$0}' ${glm_files} > snp_of_interest.txt
    """
}

// ----------------------------
// Run VEP
// ----------------------------
process run_vep {
    container 'ensemblorg/ensembl-vep:latest'
    publishDir "${params.outdir}/vep", mode: 'copy'

    input:
    path snp_file

    output:
    path "vep_output.txt"

    script:
    """
    # Check input file
    if [ ! -s ${snp_file} ]; then
        echo "Error: ${snp_file} is empty or does not exist" >&2
        exit 1
    fi

    # Copy input file
    cp ${snp_file} temp.txt

    # Run VEP online with --database
    vep --input_file temp.txt \
        --output_file vep_output.txt \
        --species homo_sapiens \
        --assembly GRCh38 \
        --database \
        --force_overwrite \
        --tab
    """
}

// ----------------------------
// Workflow
// ----------------------------
workflow {
    geno_ch = simulate_genotype()
    qc_ch   = qc_genotype(geno_ch)
    gwas_ch = assoc_test(qc_ch)
    manhattan_ch = plot_manhattan(gwas_ch)
    snp_interest_ch = filter_snp_of_interest(gwas_ch)
    run_vep(snp_interest_ch)
}
