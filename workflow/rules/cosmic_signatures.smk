"""Snakemake rules for generating format-specific COSMIC signature files

This module converts base COSMIC signature files into various formats required by
different signature analysis tools, supporting both standard and relevant (subset) modes.
"""

import os


# ==============================================================================
# WILDCARD CONSTRAINTS
# ==============================================================================

wildcard_constraints:
    format=r"(tab_with_index|csv_no_index|csv_with_index|csv_quoted)",
    data_type=r"(synthetic|realsim)",
    num_muts=r"(\d+|none)"


# ==============================================================================
# MAIN COSMIC CONVERSION RULES
# ==============================================================================

rule cosmic_signatures_by_format:
    """
    Generate format-specific COSMIC files for FITTING with out-of-reference exclusions.
    These files are given to tools during fitting and may have some signatures excluded.
    Simulation uses the base_file directly which has ALL signatures.
    """
    input:
        base_cosmic=config.get("cosmic_signatures", {}).get("base_file", "input/COSMIC_v3.4_SBS_GRCh38.txt"),
        script=os.path.join(workflow.basedir, "..", "scripts", "convert_cosmic_format.py")
    output:
        os.path.join(config["output_dir"], "cosmic_signatures", "COSMIC_v3_SBS_GRCh38_{format}.txt")
    params:
        exclude_sigs=lambda wc: " -e '" + ",".join(config["out_of_reference_sigs"]) + "'" if config["out_of_reference_sigs"] else ""
    log:
        stderr=os.path.join(config["log_dir"], "cosmic_signatures_by_format_{format}.stderr.log"),
        stdout=os.path.join(config["log_dir"], "cosmic_signatures_by_format_{format}.stdout.log")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    shell:
        "(python '{input.script}' \
        -i {input.base_cosmic} \
        -o {output} \
        -f {wildcards.format}{params.exclude_sigs}) > {log.stdout} 2> {log.stderr}"


rule cosmic_signatures_by_format_relevant:
    """
    Generate RELEVANT (subset) versions for specific data_type and num_muts.
    Allows automatic detection from each run's active signatures.
    """
    input:
        base_cosmic=config["cosmic_signatures"]["base_file"],
        script=os.path.join(workflow.basedir, "..", "scripts", "convert_cosmic_format.py"),
        active_sigs=lambda wc: (
            os.path.join(config['output_dir'], "relevant_signatures",
                        f"{wc.data_type}_{wc.num_muts}_active_signatures.txt")
            if config.get('relevant_signatures') == "auto" else []
        )
    output:
        os.path.join(config["output_dir"], "cosmic_signatures",
                    "COSMIC_v3_SBS_GRCh38_{format}_RELEVANT_{data_type}_{num_muts}.txt")
    params:
        relevant_sigs=lambda wc: get_relevant_signatures_for_cosmic(wc)
    log:
        stderr=os.path.join(config["log_dir"],
                    "cosmic_signatures_by_format_relevant_{format}_{data_type}_{num_muts}.stderr.log"),
        stdout=os.path.join(config["log_dir"],
                    "cosmic_signatures_by_format_relevant_{format}_{data_type}_{num_muts}.stdout.log")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    shell:
        "(python '{input.script}'  -i {input.base_cosmic} -o {output} -f {wildcards.format} -r '{params.relevant_sigs}') > {log.stdout} 2> {log.stderr}"
