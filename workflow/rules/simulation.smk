
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
rule create_synth_realistic_catalog:
    """Simulate samples for one tissue from signature activity in real tumours."""
    input:
        # Signature contributions file for specific cancer type
        signature_contributions_file=lambda wildcards: os.path.join(
            config["signature_contributions_dir"],
            f"signature_contributions_{config['wgs_or_wes']}_{wildcards.cancer_type}.dat"
        ),
        script=os.path.join(workflow.basedir, "..", "scripts", "MS_create_synthetic.py")
    output:
        known_answer=samples_file("known_answer_rep{rep}.dat", "synth_realistic"),
        counts=samples_file("counts_rep{rep}.dat", "synth_realistic")
    params:
        num_samples=get_param("ms_create_synthetic", "num_samples"),
        input_sigs=get_param("ms_create_synthetic", "input_sigs"),
        epsilon=get_param("ms_create_synthetic", "epsilon"),
        random_seed=get_param("ms_create_synthetic", "random_seed", "rng_seed"),
        wgs_or_wes=get_param("ms_create_synthetic", "wgs_or_wes"),
        out_of_reference_weights=get_param("ms_create_synthetic", "out_of_reference_weights"),
        out_sigs=config["out_of_reference_sigs"]
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    log:
        stderr=samples_log("synth_realistic", "stderr"),
        stdout=samples_log("synth_realistic", "stdout")
    shell:
        "(python {input.script:q}  \
        --num_samples {params.num_samples} \
        --num_muts {wildcards.num_muts} \
        --signature_contributions_file {input.signature_contributions_file:q} \
        --counts_file {output.counts:q} \
        --known_answer_file {output.known_answer:q} \
        --cancer_type {wildcards.cancer_type} \
        --realization {wildcards.rep} \
        --input_sigs {params.input_sigs:q} \
        --epsilon {params.epsilon} \
        --random_seed {params.random_seed} \
        --wgs_or_wes {params.wgs_or_wes} \
        --out_of_reference_weights {params.out_of_reference_weights} \
        --out_sigs {params.out_sigs}) > {log.stdout:q} 2> {log.stderr:q}"


rule prepare_real_catalog:
    """Prepare one group of real samples for the fitting tools.

    Only a mutation-channel catalog is required (contexts × samples) — there is no
    ground truth. The {cancer_type} wildcard is resolved back to its real_datasets
    entry (and, if that catalog uses the "Tissue::SampleID" column convention, the
    tissue to filter to) via resolve_real_dataset() — see get_real_cancer_types() in
    common.smk. Catalogs without that convention are fit as one group named after the
    dataset key. Samples keep their own mutation counts (no downsampling).
    Real samples have no target mutation count and are not drawn again, so their run
    folder is all_mutations and their repeat is 0 (REAL_NUM_MUTS in common.smk).
    Evaluated by comparing the tools with each other, since there is no known answer.
    """
    input:
        catalog=lambda wildcards: real_catalog_path(resolve_real_dataset(wildcards.cancer_type)[0]),
        # mutation order to use: the same as the simulated samples and the reference signatures
        order_like=config["cosmic_signatures"]["base_file"],
        script=os.path.join(workflow.basedir, "..", "scripts", "prepare_real_catalog.py")
    output:
        counts=samples_file("counts_rep{rep}.dat", "real")
    params:
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
        stderr=samples_log("real", "stderr"),
        stdout=samples_log("real", "stdout")
    shell:
        "(python {input.script:q} \
        --input-catalog {input.catalog:q} \
        --dataset-name {wildcards.cancer_type} \
        --counts-file {output.counts:q} \
        --order-like {input.order_like:q} \
        {params.tissue_arg}) > {log.stdout:q} 2> {log.stderr:q}"


rule create_synth_custom_catalog:
    """Simulate a custom dataset: each row of its table becomes one sample with num_muts mutations."""
    input:
        script=os.path.join(workflow.basedir, "..", "scripts", "MS_create_custom.py"),
        helpers=os.path.join(workflow.basedir, "..", "scripts", "MS_create_synthetic.py"),
        weights_file=_custom_weights_file,
        input_sigs=lambda wildcards: signature_catalogue("synth_custom", wildcards.cancer_type)
    output:
        known_answer=samples_file("known_answer_rep{rep}.dat", "synth_custom"),
        counts=samples_file("counts_rep{rep}.dat", "synth_custom")
    params:
        random_seed=config.get("random_seed", 0)
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    log:
        stderr=samples_log("synth_custom", "stderr"),
        stdout=samples_log("synth_custom", "stdout")
    shell:
        "(python {input.script:q} \
        --weights_file {input.weights_file:q} \
        --num_muts {wildcards.num_muts} \
        --counts_file {output.counts:q} \
        --known_answer_file {output.known_answer:q} \
        --dataset {wildcards.cancer_type} \
        --realization {wildcards.rep} \
        --input_sigs {input.input_sigs:q} \
        --random_seed {params.random_seed}) > {log.stdout:q} 2> {log.stderr:q}"
