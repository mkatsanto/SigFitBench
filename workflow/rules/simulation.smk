
"""Simulation rules for generating synthetic mutational catalogs"""

import os


def _custom_weights_file(wildcards):
    """The custom dataset's table of signature weights."""
    ds = config["custom_datasets"][wildcards.cancer_type]
    if "combinations_file" not in ds:
        raise ValueError(
            f"custom_datasets[{wildcards.cancer_type}] must have 'combinations_file'"
        )
    return ds["combinations_file"]

# Data generation rules
rule create_realsim_catalog:
    """Simulate samples for one tissue from signature activity in real tumours."""
    input:
        # Signature contributions file for specific cancer type
        signature_contributions_file=lambda wildcards: os.path.join(
            config["signature_contributions_dir"],
            f"signature_contributions_{config['wgs_or_wes']}_{wildcards.cancer_type}.dat"
        ),
        script=os.path.join(workflow.basedir, "..", "scripts", "MS_create_synthetic.py")
    output:
        weights=os.path.join(config["output_dir"], "generated_data", "realsim", "weights_{cancer_type}-w_{rep}-{num_muts}.dat"),
        # Dynamic outputs based on selected tools - at least one file must exist
        data_files=expand(
            os.path.join(config["output_dir"], "generated_data", "realsim", "counts_{{cancer_type}}-w_{{rep}}-{{num_muts}}_{format}.dat"),
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
        stderr=os.path.join(config["log_dir"], "create_realsim_catalog.{cancer_type}.{rep}.{num_muts}.stderr.log"),
        stdout=os.path.join(config["log_dir"], "create_realsim_catalog.{cancer_type}.{rep}.{num_muts}.stdout.log")
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


rule prepare_real_catalog:
    """Prepare one group of real samples for the fitting tools.

    Only a mutation-channel catalog is required (contexts × samples) — there is no
    ground truth. The {cancer_type} wildcard is resolved back to its real_datasets
    entry (and, if that catalog uses the "Tissue::SampleID" column convention, the
    tissue to filter to) via resolve_real_dataset() — see get_real_cancer_types() in
    common.smk. Catalogs without that convention are fit as one group named after the
    dataset key. Samples keep their own mutation counts (no downsampling).
    num_muts and rep are always "none"/0 in the output path: real samples have no
    target mutation count or repeated draws, so both are fixed placeholders in the
    shared output layout every tool rule uses (RUN_DIR in common.smk). Evaluated via
    cross-tool agreement, since there's no ground truth to compare against.
    """
    input:
        catalog=lambda wildcards: config["real_datasets"][resolve_real_dataset(wildcards.cancer_type)[0]]["catalog"],
        script=os.path.join(workflow.basedir, "..", "scripts", "prepare_real_catalog.py")
    output:
        data_files=expand(
            os.path.join(config["output_dir"], "generated_data", "real",
                         "counts_{{cancer_type}}-w_0-{{num_muts}}_{format}.dat"),
            format=get_required_output_formats()
        )
    params:
        output_formats=get_required_output_formats(),
        matrix_order=config["matrix_order"],
        num_muts_label=lambda wildcards: wildcards.num_muts,
        tissue_arg=lambda wildcards: (
            f"--tissue-prefix '{resolve_real_dataset(wildcards.cancer_type)[1]}'"
            if resolve_real_dataset(wildcards.cancer_type)[1] else ""
        )
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    log:
        stderr=os.path.join(config["log_dir"], "prepare_real_catalog.{cancer_type}.{num_muts}.stderr.log"),
        stdout=os.path.join(config["log_dir"], "prepare_real_catalog.{cancer_type}.{num_muts}.stdout.log")
    shell:
        "(python '{input.script}' \
        --input-catalog '{input.catalog}' \
        --dataset-name {wildcards.cancer_type} \
        --output-dir {config[output_dir]}/generated_data/real \
        --output-formats {params.output_formats} \
        --matrix-order {params.matrix_order} \
        --num-muts-label {params.num_muts_label} \
        {params.tissue_arg}) > {log.stdout} 2> {log.stderr}"


rule create_custom_catalog:
    """Simulate a custom dataset: each row of its table becomes one sample with num_muts mutations."""
    input:
        script=os.path.join(workflow.basedir, "..", "scripts", "MS_create_custom.py"),
        helpers=os.path.join(workflow.basedir, "..", "scripts", "MS_create_synthetic.py"),
        weights_file=_custom_weights_file,
        input_sigs=lambda wildcards: signature_catalogue("custom", wildcards.cancer_type)
    output:
        weights=os.path.join(config["output_dir"], "generated_data", "custom", "weights_{cancer_type}-w_{rep}-{num_muts}.dat"),
        data_files=expand(
            os.path.join(config["output_dir"], "generated_data", "custom", "counts_{{cancer_type}}-w_{{rep}}-{{num_muts}}_{format}.dat"),
            format=get_required_output_formats()
        )
    params:
        random_seed=config.get("random_seed", 0),
        output_formats=get_required_output_formats()
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
        --weights_file {input.weights_file} \
        --num_muts {wildcards.num_muts} \
        --output_dir $(dirname {output.weights}) \
        --dataset {wildcards.cancer_type} \
        --realization {wildcards.rep} \
        --input_sigs {input.input_sigs} \
        --random_seed {params.random_seed} \
        --output_formats {params.output_formats}) > {log.stdout} 2> {log.stderr}"
