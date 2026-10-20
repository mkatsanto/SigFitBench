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
        "(sigvcf "
        "--vcf {input.vcf} --sample-names {wildcards.identifier} "
        "--fasta {input.fasta} "
        "--output-dir {params.output_dir} "
        "--output-prefix {wildcards.identifier} "
        "--output-formats mutations "
        "{params.matrix_order_arg}) "
        "> {log.stdout} 2> {log.stderr}"


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
