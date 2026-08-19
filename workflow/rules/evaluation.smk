"""Evaluation rules for signature analysis results"""

import os

# ==============================================================================
# WILDCARD CONSTRAINTS
# ==============================================================================

wildcard_constraints:
    data_type=r"(synthetic|realsim|samples|external|custom)",
    cancer_type=r"[A-Za-z0-9_-]+",
    num_muts=r"(\d+|none)",
    identifier=r"[A-Za-z0-9_-]+",
    relevance=r"(standard|relevant)",
    tool=r"[A-Za-z0-9_]+",
    rep=r"\d+"


# ==============================================================================
# PER-RUN EVALUATION
# ==============================================================================

rule evaluate_tool_performance:
    """Evaluate signature fitting performance for a specific tool and dataset"""
    input:
        tool_results=os.path.join(
            config['output_dir'],
            "{data_type}", "{cancer_type}", "{relevance}", "{num_muts}", "rep_{rep}",
            "{tool}_contribution.dat"
        ),
        true_weights=lambda wildcards: os.path.join(
            config['output_dir'], "generated_data", wildcards.data_type,
            f"weights_{wildcards.cancer_type}-w_{'0' if wildcards.data_type == 'realsim' else wildcards.rep}-{wildcards.num_muts}.dat"
        ),
        reference_signatures=lambda wildcards: get_cosmic_signatures(wildcards.tool, wildcards.relevance),
        cosmic_file=config['cosmic_signatures']['base_file'],
        script=os.path.join(workflow.basedir, "..", "scripts", "evaluate_tool_performance.py")
    output:
        per_sample_metrics=os.path.join(
            config['output_dir'], "evaluation",
            "{data_type}_{cancer_type}_{num_muts}_{relevance}_{tool}_w{rep}_per_sample.tsv"
        )
    params:
        min_contribution=config["evaluation"]["min_contribution"],
        out_of_reference_sigs=OUT_OF_REF_SIGS_STR,
        penalty_arg=""
    log:
        stderr=os.path.join(config["log_dir"], "evaluation_{data_type}_{cancer_type}_{num_muts}_{relevance}_{tool}_w{rep}.stderr.log"),
        stdout=os.path.join(config["log_dir"], "evaluation_{data_type}_{cancer_type}_{num_muts}_{relevance}_{tool}_w{rep}.stdout.log")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    shell:
        "(python '{input.script}' \
        --tool-results {input.tool_results} \
        --true-weights {input.true_weights} \
        --reference-signatures {input.reference_signatures} \
        --cosmic-file {input.cosmic_file} \
        --output {output.per_sample_metrics} \
        --tool {wildcards.tool} \
        --num-muts {wildcards.num_muts} \
        --min-contribution {params.min_contribution} \
        --realization {wildcards.rep} \
        --cancer-type {wildcards.cancer_type} \
        --data-type {wildcards.data_type} \
        --out-of-reference-sigs {params.out_of_reference_sigs} \
        {params.penalty_arg}) \
        > {log.stdout} 2> {log.stderr}"


use rule evaluate_tool_performance as evaluate_penalty_method_performance with:
    input:
        tool_results=os.path.join(
            config['output_dir'],
            "{data_type}", "{cancer_type}", "{relevance}", "{num_muts}", "rep_{rep}",
            "{tool}_contribution_penalty_{penalty}.dat"
        ),
        true_weights=lambda wildcards: os.path.join(
            config['output_dir'], "generated_data", wildcards.data_type,
            f"weights_{wildcards.cancer_type}-w_{'0' if wildcards.data_type == 'realsim' else wildcards.rep}-{wildcards.num_muts}.dat"
        ),
        reference_signatures=lambda wildcards: get_cosmic_signatures(wildcards.tool, wildcards.relevance),
        cosmic_file=config['cosmic_signatures']['base_file'],
        script=os.path.join(workflow.basedir, "..", "scripts", "evaluate_tool_performance.py")
    output:
        per_sample_metrics=os.path.join(
            config['output_dir'], "evaluation",
            "{data_type}_{cancer_type}_{num_muts}_{relevance}_{tool}_penalty_{penalty}_w{rep}_per_sample.tsv"
        )
    params:
        min_contribution=config["evaluation"]["min_contribution"],
        out_of_reference_sigs=OUT_OF_REF_SIGS_STR,
        penalty_arg=lambda wildcards: f"--penalty-value {wildcards.penalty}"
    log:
        stderr=os.path.join(config["log_dir"], "evaluation_{data_type}_{cancer_type}_{num_muts}_{relevance}_{tool}_penalty_{penalty}_w{rep}.stderr.log"),
        stdout=os.path.join(config["log_dir"], "evaluation_{data_type}_{cancer_type}_{num_muts}_{relevance}_{tool}_penalty_{penalty}_w{rep}.stdout.log")


rule assess_fit_quality_bootstrap:
    """Assess signature fitting quality using bootstrap resampling."""
    input:
        tool_results=os.path.join(
            config['output_dir'],
            "{data_type}", "{cancer_type}", "{relevance}", "{num_muts}", "rep_{rep}",
            "{tool}_contribution.dat"
        ),
        input_catalog=lambda wildcards: get_input_data(wildcards.tool.lower(), wildcards),
        reference_signatures=lambda wildcards: get_cosmic_signatures(wildcards.tool, wildcards.relevance),
        script=os.path.join(workflow.basedir, "..", "scripts", "assess_fit_quality_bootstrap.py")
    output:
        assessment=os.path.join(
            config['output_dir'], "evaluation",
            "{data_type}_{cancer_type}_{num_muts}_{relevance}_{tool}_w{rep}_bootstrap_quality.tsv"
        )
    params:
        num_bootstrap=config["evaluation"]["num_bootstrap"],
        random_seed=lambda wildcards: config["random_seed"] + int(wildcards.rep)
    log:
        stderr=os.path.join(config["log_dir"], "bootstrap_quality_{data_type}_{cancer_type}_{num_muts}_{relevance}_{tool}_w{rep}.stderr.log"),
        stdout=os.path.join(config["log_dir"], "bootstrap_quality_{data_type}_{cancer_type}_{num_muts}_{relevance}_{tool}_w{rep}.stdout.log")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    shell:
        "(python '{input.script}' \
        --tool-results {input.tool_results} \
        --input-catalog {input.input_catalog} \
        --reference-signatures {input.reference_signatures} \
        --output {output.assessment} \
        --tool {wildcards.tool} \
        --num-bootstrap {params.num_bootstrap} \
        --random-seed {params.random_seed}) \
        > {log.stdout} 2> {log.stderr}"


rule evaluate_real_samples:
    """Evaluate signature fitting for real patient samples (without ground truth)"""
    input:
        tool_results=os.path.join(
            config['output_dir'],
            "samples", "{identifier}", "{relevance}",
            "{tool}_contribution.dat"
        ),
        script=os.path.join(workflow.basedir, "..", "scripts", "evaluate_real_samples.py")
    output:
        metrics=os.path.join(
            config['output_dir'], "evaluation",
            "samples_{identifier}_{relevance}_{tool}_metrics.dat"
        )
    log:
        stderr=os.path.join(config["log_dir"], "evaluation_samples_{identifier}_{relevance}_{tool}.stderr.log"),
        stdout=os.path.join(config["log_dir"], "evaluation_samples_{identifier}_{relevance}_{tool}.stdout.log")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    shell:
        "(python '{input.script}' \
        --tool-results {input.tool_results} \
        --output {output.metrics} \
        --tool {wildcards.tool} \
        --sample-id {wildcards.identifier}) \
        > {log.stdout} 2> {log.stderr}"


# ==============================================================================
# AGGREGATION
# ==============================================================================

rule compute_summary_metrics:
    """Concatenate per-sample TSVs and derive per-run summary statistics."""
    input:
        eval_files=get_eval_files(),
        script=os.path.join(workflow.basedir, "..", "scripts", "compute_summary_metrics.py")
    output:
        per_sample_metrics=os.path.join(config['output_dir'], "all_per_sample_metrics.tsv"),
        summary=os.path.join(config['output_dir'], "all_evaluation_metrics.tsv")
    log:
        stderr=os.path.join(config["log_dir"], "compute_summary_metrics.stderr.log"),
        stdout=os.path.join(config["log_dir"], "compute_summary_metrics.stdout.log")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    shell:
        "(python '{input.script}' \
        --input-files {input.eval_files} \
        --output-per-sample {output.per_sample_metrics} \
        --output {output.summary}) \
        > {log.stdout} 2> {log.stderr}"


rule aggregate_bootstrap_quality:
    """Aggregate bootstrap quality assessment results across all tools and conditions"""
    input:
        assessments=expand(
            os.path.join(
                config['output_dir'], "evaluation",
                "{data_type}_{cancer_type}_{num_muts}_{relevance}_{tool}_w{rep}_bootstrap_quality.tsv"
            ),
            data_type=active_data_types,
            cancer_type=config["cancer_types"],
            num_muts=config["num_muts"],
            relevance=config["relevance_modes"],
            tool=[t for t in get_selected_tools() if not has_penalty_sweep(t)],
            rep=range(config["num_realizations"])
        ) if config.get("evaluation", {}).get("run_bootstrap_quality", False) else [],
        script=os.path.join(workflow.basedir, "..", "scripts", "aggregate_bootstrap_quality.py")
    output:
        all_bootstrap=os.path.join(config['output_dir'], "all_bootstrap_quality.tsv")
    log:
        stderr=os.path.join(config["log_dir"], "aggregate_bootstrap_quality.stderr.log"),
        stdout=os.path.join(config["log_dir"], "aggregate_bootstrap_quality.stdout.log")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    shell:
        "(python '{input.script}' \
        --input-files {input.assessments} \
        --output {output.all_bootstrap}) \
        > {log.stdout} 2> {log.stderr}"


# ==============================================================================
# CROSS-TOOL AGREEMENT (no-ground-truth data types: realsim, external)
# ==============================================================================

rule evaluate_tool_agreement:
    """Compare every selected tool's fit pairwise for one no-ground-truth cell
    (realsim or external without ground_truth)."""
    input:
        tool_results=lambda wildcards: expand(
            os.path.join(
                config['output_dir'], "{data_type}", "{cancer_type}", "{relevance}",
                "{num_muts}", "rep_0", "{tool}_contribution.dat"
            ),
            data_type=wildcards.data_type,
            cancer_type=wildcards.cancer_type,
            relevance=wildcards.relevance,
            num_muts=wildcards.num_muts,
            tool=get_agreement_tools()
        ),
        script=os.path.join(workflow.basedir, "..", "scripts", "evaluate_tool_agreement.py")
    output:
        os.path.join(
            config['output_dir'], "evaluation",
            "{data_type}_{cancer_type}_{num_muts}_{relevance}_tool_agreement.tsv"
        )
    params:
        tools=get_agreement_tools(),
        min_contribution=config["evaluation"]["min_contribution"]
    log:
        stderr=os.path.join(config["log_dir"], "tool_agreement_{data_type}_{cancer_type}_{num_muts}_{relevance}.stderr.log"),
        stdout=os.path.join(config["log_dir"], "tool_agreement_{data_type}_{cancer_type}_{num_muts}_{relevance}.stdout.log")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    shell:
        "(python '{input.script}' \
        --tool-results {input.tool_results} \
        --tools {params.tools} \
        --min-contribution {params.min_contribution} \
        --data-type {wildcards.data_type} \
        --cancer-type {wildcards.cancer_type} \
        --relevance {wildcards.relevance} \
        --num-muts {wildcards.num_muts} \
        --output {output}) > {log.stdout} 2> {log.stderr}"


rule aggregate_tool_agreement:
    """Combine per-cell agreement TSVs from realsim and external datasets into one table."""
    input:
        agreement_files=get_agreement_files(),
        script=os.path.join(workflow.basedir, "..", "scripts", "aggregate_tool_agreement.py")
    output:
        os.path.join(config['output_dir'], "all_tool_agreement.tsv")
    log:
        stderr=os.path.join(config["log_dir"], "aggregate_tool_agreement.stderr.log"),
        stdout=os.path.join(config["log_dir"], "aggregate_tool_agreement.stdout.log")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    shell:
        "(python '{input.script}' \
        --input-files {input.agreement_files} \
        --output {output}) > {log.stdout} 2> {log.stderr}"


# ==============================================================================
# REPORTING
# ==============================================================================

rule create_html_report:
    """Create comprehensive HTML report with interactive visualizations"""
    input:
        summaries=os.path.join(config['output_dir'], "all_evaluation_metrics.tsv"),
        per_sample_metrics=os.path.join(config['output_dir'], "all_per_sample_metrics.tsv"),
        config_file=workflow.configfiles[0] if workflow.configfiles else []
    output:
        report=os.path.join(config['output_dir'], "sigfittest_report.html")
    params:
        title="SigFitTest Analysis Report",
        subtitle="Mutational Signature Fitting Tool Comparison and Evaluation",
        results_dir=config['output_dir']
    log:
        stderr=os.path.join(config["log_dir"], "html_report.stderr.log"),
        stdout=os.path.join(config["log_dir"], "html_report.stdout.log")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    script:
        os.path.join(workflow.basedir, "..", "scripts", "create_html_report.py")
