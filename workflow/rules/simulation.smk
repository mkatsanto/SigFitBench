
"""Simulation rules for generating synthetic mutational catalogs"""

import os


def _custom_mode_args(wildcards):
    """Return CLI args string for the active custom mode.

    Combinations (combinations_file key):
        → --signature_contributions_file path/to/file.dat
          Each sample randomly draws one row from the table (same logic as tissue data).
    """
    ds = config["custom_datasets"][wildcards.cancer_type]

    if "combinations_file" not in ds:
        raise ValueError(
            f"custom_datasets[{wildcards.cancer_type}] must have 'combinations_file'"
        )
    return f"--signature_contributions_file {ds['combinations_file']}"

def _custom_out_of_reference_weights(wildcards):
    ds = config["custom_datasets"][wildcards.cancer_type]
    return ds.get("out_of_reference_weights", [])  # no unknowns by default for custom datasets

def _custom_out_sigs(wildcards):
    ds = config["custom_datasets"][wildcards.cancer_type]
    return ds.get("out_of_reference_sigs", [])

# Data generation rules
rule ms_create_synthetic_catalog:
    input:
        # Signature contributions file for specific cancer type
        signature_contributions_file=lambda wildcards: os.path.join(
            config["signature_contributions_dir"],
            f"signature_contributions_{config['wgs_or_wes']}_{wildcards.cancer_type}.dat"
        ),
        script=os.path.join(workflow.basedir, "..", "scripts", "MS_create_synthetic.py")
    output:
        weights=os.path.join(config["output_dir"], "generated_data", "synthetic", "weights_{cancer_type}-w_{rep}-{num_muts}.dat"),
        # Dynamic outputs based on selected tools - at least one file must exist
        data_files=expand(
            os.path.join(config["output_dir"], "generated_data", "synthetic", "counts_{{cancer_type}}-w_{{rep}}-{{num_muts}}_{format}.dat"),
            format=get_required_output_formats()
        )
    params:
        num_samples=get_param("ms_create_synthetic", "num_samples"),
        input_sigs=get_param("ms_create_synthetic", "input_sigs"),
        epsilon=get_param("ms_create_synthetic", "epsilon"),
        random_seed=get_param("ms_create_synthetic", "random_seed", "rng_seed"),
        output_formats=get_required_output_formats(),
        wgs_or_wes=get_param("ms_create_synthetic", "wgs_or_wes"),
        out_of_reference_weights=get_param("ms_create_synthetic", "out_of_reference_weights"),
        out_sigs=config["out_of_reference_sigs"]
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    log:
        stderr=os.path.join(config["log_dir"], "ms_create_synthetic_catalog.{cancer_type}.{rep}.{num_muts}.stderr.log"),
        stdout=os.path.join(config["log_dir"], "ms_create_synthetic_catalog.{cancer_type}.{rep}.{num_muts}.stdout.log")
    shell:
        "(python '{input.script}'  \
        --num_samples {params.num_samples} \
        --num_muts {wildcards.num_muts} \
        --signature_contributions_file {input.signature_contributions_file} \
        --output_dir $(dirname {output.weights}) \
        --cancer_type {wildcards.cancer_type} \
        --realization {wildcards.rep} \
        --input_sigs {params.input_sigs} \
        --epsilon {params.epsilon} \
        --random_seed {params.random_seed} \
        --wgs_or_wes {params.wgs_or_wes} \
        --output_formats {params.output_formats} \
        --out_of_reference_weights {params.out_of_reference_weights} \
        --out_sigs {params.out_sigs}) > {log.stdout} 2> {log.stderr}"


# Real samples: subset a real mutation catalog to one cancer type, no downsampling
rule create_realsim_catalog:
    """Select real samples for one cancer type from a PCAWG-style catalog, at their
    native mutation counts (no downsampling, no simulation).

    num_muts is fixed to the literal "none" for this data_type — real samples don't
    share one target mutation count, so nothing here depends on it (see wildcard_constraints).
    rep is always 0. There is no ground truth for real samples: SigID (whose auto-penalty
    search logs a diagnostic against true weights) and evaluate_tool_performance require a
    weights_*.dat file and so cannot be used with data_type="realsim".
    """
    input:
        real_catalog=lambda wildcards: get_param("real_samples", "input_real_catalog"),
        script=os.path.join(workflow.basedir, "..", "scripts", "prepare_real_samples.py")
    output:
        data_files=expand(
            os.path.join(config["output_dir"], "generated_data", "realsim",
                         "counts_{{cancer_type}}-w_0-{{num_muts}}_{format}.dat"),
            format=get_required_output_formats()
        )
    params:
        output_formats=get_required_output_formats()
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    log:
        stderr=os.path.join(config["log_dir"], "create_realsim_catalog.{cancer_type}.{num_muts}.stderr.log"),
        stdout=os.path.join(config["log_dir"], "create_realsim_catalog.{cancer_type}.{num_muts}.stdout.log")
    shell:
        "(python '{input.script}' \
        --input_catalog {input.real_catalog} \
        --cancer_type {wildcards.cancer_type} \
        --num_muts_label {wildcards.num_muts} \
        --output_dir {config[output_dir]}/generated_data/realsim \
        --output_formats {params.output_formats} \
        ) > {log.stdout} 2> {log.stderr}"


rule prepare_external_catalog:
    """Prepare external mutational catalog for fitting tools.

    Only a mutation-channel catalog is required (contexts × samples) — no ground
    truth, just like realsim. external_datasets[{cancer_type}]["catalog"] is the
    only mandatory key. rep is always 0. Evaluated via cross-tool agreement, same
    as realsim, since there's no ground truth to compare against.
    """
    input:
        catalog=lambda wildcards: config["external_datasets"][wildcards.cancer_type]["catalog"],
        script=os.path.join(workflow.basedir, "..", "scripts", "fit_external_catalog.py")
    output:
        data_files=expand(
            os.path.join(config["output_dir"], "generated_data", "external",
                         "counts_{{cancer_type}}-w_0-{{num_muts}}_{format}.dat"),
            format=get_required_output_formats()
        )
    params:
        output_formats=get_required_output_formats(),
        matrix_order=config["matrix_order"],
        num_muts_label=lambda wildcards: wildcards.num_muts
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    log:
        stderr=os.path.join(config["log_dir"], "prepare_external_catalog.{cancer_type}.{num_muts}.stderr.log"),
        stdout=os.path.join(config["log_dir"], "prepare_external_catalog.{cancer_type}.{num_muts}.stdout.log")
    shell:
        "(python '{input.script}' \
        --input-catalog '{input.catalog}' \
        --dataset-name {wildcards.cancer_type} \
        --output-dir {config[output_dir]}/generated_data/external \
        --output-formats {params.output_formats} \
        --matrix-order {params.matrix_order} \
        --num-muts-label {params.num_muts_label}) > {log.stdout} 2> {log.stderr}"


rule create_custom_catalog:
    """Generate synthetic catalog from a user-specified custom signature configuration.

    Each dataset must give a 'combinations_file' (TSV: rows=combinations, cols=signatures);
    every sample randomly draws one row from the table, same logic as tissue data.

    Per-sample variation in the count matrix comes from multinomial mutation sampling.
    Per-dataset out_of_reference_sigs/weights override the top-level config values.
    """
    input:
        script=os.path.join(workflow.basedir, "..", "scripts", "MS_create_synthetic.py"),
        input_sigs=config["cosmic_signatures"]["base_file"]
    output:
        weights=os.path.join(config["output_dir"], "generated_data", "custom", "weights_{cancer_type}-w_{rep}-{num_muts}.dat"),
        data_files=expand(
            os.path.join(config["output_dir"], "generated_data", "custom", "counts_{{cancer_type}}-w_{{rep}}-{{num_muts}}_{format}.dat"),
            format=get_required_output_formats()
        )
    params:
        num_samples=config.get("num_samples", 100),
        random_seed=config.get("random_seed", 0),
        output_formats=get_required_output_formats(),
        mode_args=_custom_mode_args,
        out_of_reference_weights=_custom_out_of_reference_weights,
        out_sigs=_custom_out_sigs
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    log:
        stderr=os.path.join(config["log_dir"], "create_custom_catalog.{cancer_type}.{rep}.{num_muts}.stderr.log"),
        stdout=os.path.join(config["log_dir"], "create_custom_catalog.{cancer_type}.{rep}.{num_muts}.stdout.log")
    shell:
        "(python '{input.script}' \
        --num_samples {params.num_samples} \
        --num_muts {wildcards.num_muts} \
        --output_dir $(dirname {output.weights}) \
        --cancer_type {wildcards.cancer_type} \
        --realization {wildcards.rep} \
        --input_sigs {input.input_sigs} \
        --random_seed {params.random_seed} \
        --output_formats {params.output_formats} \
        {params.mode_args} \
        --out_of_reference_weights {params.out_of_reference_weights} \
        --out_sigs {params.out_sigs}) > {log.stdout} 2> {log.stderr}"
