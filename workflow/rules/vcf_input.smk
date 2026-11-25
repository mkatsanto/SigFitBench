# ==============================================================================
# VCF INPUT (data_type == "samples")
# ==============================================================================
# Turns one VCF into a table of mutation counts. The tools do not read this file: to get it,
# ask Snakemake for it by name, e.g. `snakemake input/samples/<sample>_mutations.dat`.
# To fit real VCFs with the tools, use a real_datasets entry with vcf_dir (below).
#
# Config:
#   samples:
#     vcf_dir: "input/vcf"                     # {identifier}.vcf(.gz) files live here
#     reference_fasta: "input/ref/genome.fa"   # needs a matching genome.fa.fai

# SigVCF is run locally from its source folder (no installation needed): config 'sigvcf_source' (relative to the folder
# snakemake is started from), or by default the SigVCF folder next to SigFitBench
SIGVCF_SOURCE = os.path.abspath(config.get("sigvcf_source") or os.path.join(workflow.basedir, "..", "..", "SigVCF"))
if not os.path.exists(os.path.join(SIGVCF_SOURCE, "sigvcf", "cli.py")):
    logger.warning(f"SigVCF not found in {SIGVCF_SOURCE}: set 'sigvcf_source' in the config to build catalogs from VCFs")
SIGVCF = f"PYTHONPATH='{SIGVCF_SOURCE}' python -m sigvcf.cli"


rule generate_matrix_from_vcf:
    """Count the mutations of one VCF per mutation type (96 types), with SigVCF.
    SigVCF reads the VCF and the genome bit by bit, so it needs little memory even for
    whole genomes."""
    input:
        vcf=lambda wc: _resolve_sample_vcf(wc.identifier),
        fasta=lambda wc: config["samples"]["reference_fasta"],
        fai=lambda wc: config["samples"]["reference_fasta"] + ".fai",
    output:
        data=os.path.join(config.get("sample_input_dir", "input/samples"), "{identifier}_mutations.dat")
    params:
        output_dir=config.get("sample_input_dir", "input/samples"),
        matrix_order_arg=f"--matrix-order '{config['matrix_order']}'" if config.get("matrix_order") else "",
    conda:
        os.path.join(workflow.basedir, "envs", "sigvcf.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 1000 * attempt,
    log:
        stderr=os.path.join(config["log_dir"], "generate_matrix_from_vcf.{identifier}.stderr.log"),
        stdout=os.path.join(config["log_dir"], "generate_matrix_from_vcf.{identifier}.stdout.log")
    shell:
        f"({SIGVCF} "
        "--vcf {input.vcf:q} --sample-names {wildcards.identifier} "
        "--fasta {input.fasta:q} "
        "--output-dir {params.output_dir:q} "
        "--output-prefix {wildcards.identifier} "
        "--output-formats mutations "
        "{params.matrix_order_arg}) "
        "> {log.stdout:q} 2> {log.stderr:q}"


def _resolve_sample_vcf(identifier):
    """First existing {vcf_dir}/{identifier}.vcf[.gz], in that order."""
    vcf_dir = config.get("samples", {}).get("vcf_dir", "input/vcf")
    for ext in (".vcf", ".vcf.gz"):
        candidate = os.path.join(vcf_dir, f"{identifier}{ext}")
        if os.path.exists(candidate):
            return candidate
    # neither exists: Snakemake then reports the .vcf as missing
    return os.path.join(vcf_dir, f"{identifier}.vcf")


# ==============================================================================
# CATALOG FROM A FOLDER OF VCFs (data_type == "real", real_datasets entry with vcf_dir)
# ==============================================================================
# Config:
#   real_datasets:
#     MyCohort:
#       vcf_dir: "path/to/vcfs"                  # one {sample}.vcf(.gz) per sample
#       reference_fasta: "input/ref/genome.fa"   # optional, default samples.reference_fasta; needs a .fai index

import glob
import re

VCF_DATASETS = [name for name, ds in config.get("real_datasets", {}).items() if (ds or {}).get("vcf_dir") and not ds.get("catalog")]


def _dataset_vcfs(dataset):
    """Every .vcf / .vcf.gz in the dataset's vcf_dir, sorted (one sample each)."""
    vcf_dir = config["real_datasets"][dataset]["vcf_dir"]
    vcfs = sorted(glob.glob(os.path.join(vcf_dir, "*.vcf")) + glob.glob(os.path.join(vcf_dir, "*.vcf.gz")))
    if not vcfs:
        raise ValueError(f"no .vcf or .vcf.gz files in {vcf_dir} (real_datasets entry {dataset!r})")
    return vcfs


def _dataset_fasta(dataset):
    return config["real_datasets"][dataset].get("reference_fasta") or config.get("samples", {}).get("reference_fasta")


if VCF_DATASETS:
    rule catalog_from_vcfs:
        """Count the mutations per mutation type (96 types) of every VCF in a real dataset's
        vcf_dir, with SigVCF. Each VCF becomes one sample, named after the file."""
        input:
            vcfs=lambda wc: _dataset_vcfs(wc.dataset),
            fasta=lambda wc: _dataset_fasta(wc.dataset),
            fai=lambda wc: _dataset_fasta(wc.dataset) + ".fai",
        output:
            os.path.join(config["output_dir"], "internal", "catalogs_from_vcfs", "{dataset}_sbs96.dat")
        wildcard_constraints:
            dataset="|".join(re.escape(name) for name in VCF_DATASETS)
        params:
            output_dir=os.path.join(config["output_dir"], "internal", "catalogs_from_vcfs"),
            matrix_order_arg=f"--matrix-order '{config['matrix_order']}'" if config.get("matrix_order") else "",
        conda:
            os.path.join(workflow.basedir, "envs", "sigvcf.yaml")
        threads: 1
        resources:
            mem_mb=lambda wildcards, attempt: 2000 * attempt,
        log:
            stderr=os.path.join(config["log_dir"], "internal", "catalogs_from_vcfs", "{dataset}.stderr.log"),
            stdout=os.path.join(config["log_dir"], "internal", "catalogs_from_vcfs", "{dataset}.stdout.log")
        shell:
            f"({SIGVCF} "
            "--vcf {input.vcfs:q} "
            "--fasta {input.fasta:q} "
            "--output-dir {params.output_dir:q} "
            "--output-prefix {wildcards.dataset} "
            "--output-formats sbs96 "
            "{params.matrix_order_arg}) "
            "> {log.stdout:q} 2> {log.stderr:q}"
