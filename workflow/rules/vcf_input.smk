# ==============================================================================
# VCF INPUT (data_type == "samples")
# ==============================================================================
# Lets a raw .vcf/.vcf.gz stand in for a pre-computed counts_*.dat matrix.
# get_input_data() (common.smk) already resolves data_type=="samples" to
# "{sample_input_dir}/{identifier}_mutations.dat" *before* falling back to a
# tool-specific format or a configured input_catalog -- this rule just
# supplies that file whenever it doesn't already exist but a matching VCF
# does, so every existing tool rule picks it up unchanged. No wiring of
# `rule all` or the data_type=="samples" aggregation is included here --
# that path is otherwise unbuilt in this workflow; target this rule's
# output (or a downstream tool_contribution.dat under data_type=samples)
# directly with `snakemake <path>` until/unless that aggregation exists.
#
# Config (see config/config_test.yaml for the wired example):
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
    """Build an SBS96 mutation-count matrix directly from a VCF via SigVCF.

    Streams the VCF and (indexed) FASTA rather than loading either into
    memory, so this is safe to run on WGS-scale VCFs -- see ../../SigVCF.
    """
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
        mem_mb=lambda wildcards, attempt: 1000 * attempt,  # streaming: flat, low memory
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
    # Let Snakemake raise its normal "missing input" error against the .vcf
    # path rather than failing early here with a less specific message.
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
        """One SBS96 catalog (mutation types x samples) from every VCF in a real_datasets entry's vcf_dir, via SigVCF.
        Each VCF becomes one column, named after the file (without .vcf / .vcf.gz)."""
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
            mem_mb=lambda wildcards, attempt: 2000 * attempt,   # streaming: flat, low memory
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
