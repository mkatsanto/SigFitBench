"""Evaluation rules for signature analysis results"""

import os

# ==============================================================================
# WILDCARD CONSTRAINTS
# ==============================================================================

wildcard_constraints:
    data_type=r"(synth_realistic|real|samples|synth_custom)",
    cancer_type=r"[A-Za-z0-9_-]+",
    num_muts=r"(\d+|all)",
    identifier=r"[A-Za-z0-9_-]+",
    relevance=r"(|_relevant)",
    tool=r"[A-Za-z0-9_]+",
    rep=r"\d+"


# ==============================================================================
# PER-RUN EVALUATION
# ==============================================================================

rule evaluate_tool_performance:
    """Evaluate signature fitting performance for a specific tool and dataset"""
    input:
        tool_results=tool_result("{tool}"),
        true_weights=samples_file("known_answer_rep{rep}.dat"),
        reference_signatures=lambda wildcards: get_cosmic_signatures(
            wildcards.tool, wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type
        ),
        cosmic_file=lambda wildcards: signature_catalogue(wildcards.data_type, wildcards.cancer_type),
        script=os.path.join(workflow.basedir, "..", "scripts", "evaluate_tool_performance.py")
    output:
        per_sample_metrics=run_file("scores", "{tool}_rep{rep}.tsv")
    params:
        min_contribution=config["evaluation"]["min_contribution"],
        out_of_reference_sigs=OUT_OF_REF_SIGS_STR
    log:
        **run_log("scores_{tool}_rep{rep}")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    shell:
        "(python {input.script:q} \
        --tool-results {input.tool_results:q} \
        --true-weights {input.true_weights:q} \
        --reference-signatures {input.reference_signatures:q} \
        --cosmic-file {input.cosmic_file:q} \
        --output {output.per_sample_metrics:q} \
        --tool {wildcards.tool} \
        --num-muts {wildcards.num_muts} \
        --min-contribution {params.min_contribution} \
        --realization {wildcards.rep} \
        --cancer-type {wildcards.cancer_type} \
        --data-type {wildcards.data_type} \
        --out-of-reference-sigs {params.out_of_reference_sigs}) \
        > {log.stdout:q} 2> {log.stderr:q}"


rule assess_fit_quality_bootstrap:
    """Assess signature fitting quality using bootstrap resampling."""
    input:
        tool_results=tool_result("{tool}"),
        input_catalog=lambda wildcards: get_input_data(wildcards.tool.lower(), wildcards),
        reference_signatures=lambda wildcards: get_cosmic_signatures(
            wildcards.tool, wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type
        ),
        script=os.path.join(workflow.basedir, "..", "scripts", "assess_fit_quality_bootstrap.py")
    output:
        assessment=run_file("scores", "{tool}_bootstrap_rep{rep}.tsv")
    params:
        num_bootstrap=config["evaluation"]["num_bootstrap"],
        random_seed=lambda wildcards: config["random_seed"] + int(wildcards.rep)
    log:
        **run_log("scores_{tool}_bootstrap_rep{rep}")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    shell:
        "(python {input.script:q} \
        --tool-results {input.tool_results:q} \
        --input-catalog {input.input_catalog:q} \
        --reference-signatures {input.reference_signatures:q} \
        --output {output.assessment:q} \
        --tool {wildcards.tool} \
        --num-bootstrap {params.num_bootstrap} \
        --random-seed {params.random_seed}) \
        > {log.stdout:q} 2> {log.stderr:q}"


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
        "(python {input.script:q} \
        --tool-results {input.tool_results:q} \
        --output {output.metrics:q} \
        --tool {wildcards.tool} \
        --sample-id {wildcards.identifier}) \
        > {log.stdout:q} 2> {log.stderr:q}"


# ==============================================================================
# AGGREGATION
# ==============================================================================

rule compute_summary_metrics:
    """Concatenate per-sample TSVs and derive per-run summary statistics."""
    input:
        eval_files=get_eval_files(),
        script=os.path.join(workflow.basedir, "..", "scripts", "compute_summary_metrics.py")
    output:
        per_sample_metrics=summary_file("per_sample_metrics"),
        summary=summary_file("evaluation_metrics")
    log:
        **summary_log("evaluation_metrics")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    shell:
        "(python {input.script:q} \
        --input-files {input.eval_files:q} \
        --output-per-sample {output.per_sample_metrics:q} \
        --output {output.summary:q}) \
        > {log.stdout:q} 2> {log.stderr:q}"


# Run time and memory files of every tool run, for the realistic simulations
_benchmark_files = expand(
    tool_benchmark("{tool}"),
    data_type=active_data_types,
    cancer_type=config.get("cancer_types", ["Head-SCC"]),
    num_muts=config["num_muts"],
    relevance=relevance_endings(),
    tool=get_selected_tools(),
    rep=range(config.get("num_realizations", 1)),
)
# ...and for the custom simulations, whose runs are named after the dataset instead of a tissue
if "synth_custom" in config.get("data_types", []):
    _benchmark_files += expand(
        tool_benchmark("{tool}"),
        data_type="synth_custom",
        cancer_type=list(config.get("custom_datasets", {})),
        relevance=relevance_endings(),
        num_muts=config["num_muts"],
        rep=range(config.get("num_realizations", 1)),
        tool=get_selected_tools(),
    )


rule aggregate_tool_benchmarks:
    """Concatenate every tool rule's benchmark TSV (wall time + peak memory, written by
    Snakemake's `benchmark:` directive in tools.smk) into one summary table."""
    input:
        benchmark_files=_benchmark_files,
        script=os.path.join(workflow.basedir, "..", "scripts", "aggregate_tool_benchmarks.py")
    output:
        summary_file("tool_benchmarks")
    log:
        **summary_log("tool_benchmarks")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    shell:
        "(python {input.script:q} \
        --input-files {input.benchmark_files:q} \
        --output {output:q}) \
        > {log.stdout:q} 2> {log.stderr:q}"


rule aggregate_bootstrap_quality:
    """Aggregate bootstrap quality assessment results across all tools and conditions"""
    input:
        assessments=expand(
            run_file("scores", "{tool}_bootstrap_rep{rep}.tsv"),
            data_type=active_data_types,
            cancer_type=config["cancer_types"],
            num_muts=config["num_muts"],
            relevance=relevance_endings(),
            tool=get_selected_tools(),
            rep=range(config["num_realizations"])
        ) if config.get("evaluation", {}).get("run_bootstrap_quality", False) else [],
        script=os.path.join(workflow.basedir, "..", "scripts", "aggregate_bootstrap_quality.py")
    output:
        all_bootstrap=summary_file("bootstrap_quality")
    log:
        **summary_log("bootstrap_quality")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    shell:
        "(python {input.script:q} \
        --input-files {input.assessments:q} \
        --output {output.all_bootstrap:q}) \
        > {log.stdout:q} 2> {log.stderr:q}"


# ==============================================================================
# CROSS-TOOL AGREEMENT (real samples: no known answer, so the tools are compared with each other)
# ==============================================================================

rule evaluate_tool_agreement:
    """Compare every selected tool's fit pairwise for one no-ground-truth cell
    (real data)."""
    input:
        tool_results=lambda wildcards: expand(
            tool_result("{tool}"),
            rep=0,
            data_type=wildcards.data_type,
            cancer_type=wildcards.cancer_type,
            relevance=wildcards.relevance,
            num_muts=wildcards.num_muts,
            tool=get_agreement_tools()
        ),
        original_catalog=samples_file("counts_rep0.dat"),
        script=os.path.join(workflow.basedir, "..", "scripts", "evaluate_tool_agreement.py")
    output:
        pairwise=run_file("agreement", "tool_agreement.tsv"),
        signature_calls=run_file("agreement", "signature_calls.tsv"),
        sample_activity=run_file("agreement", "sample_activity.tsv"),
        signature_weights=run_file("agreement", "signature_weights.tsv")
    params:
        tools=get_agreement_tools(),
        min_contribution=config["evaluation"]["min_contribution"],
        relevance=lambda wildcards: relevance_mode(wildcards.relevance)
    log:
        **run_log("agreement")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    shell:
        "(python {input.script:q} \
        --tool-results {input.tool_results:q} \
        --tools {params.tools} \
        --min-contribution {params.min_contribution} \
        --data-type {wildcards.data_type} \
        --cancer-type {wildcards.cancer_type} \
        --relevance {params.relevance} \
        --num-muts {wildcards.num_muts} \
        --output {output.pairwise:q} \
        --output-signature-calls {output.signature_calls:q} \
        --output-sample-activity {output.sample_activity:q} \
        --output-signature-weights {output.signature_weights:q} \
        --original-catalog {input.original_catalog:q}) > {log.stdout:q} 2> {log.stderr:q}"


rule aggregate_tool_agreement:
    """Combine per-cell agreement TSVs from all real datasets into one table."""
    input:
        agreement_files=get_agreement_files(),
        script=os.path.join(workflow.basedir, "..", "scripts", "aggregate_tool_agreement.py")
    output:
        summary_file("tool_agreement")
    log:
        **summary_log("tool_agreement")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    shell:
        "(python {input.script:q} \
        --input-files {input.agreement_files:q} \
        --output {output:q}) > {log.stdout:q} 2> {log.stderr:q}"


rule aggregate_signature_calls:
    """Combine per-cell signature-call-consensus TSVs (which tools call which
    signatures active, and how often) into one table — same generic concat
    script as aggregate_tool_agreement, just a different input/output pair."""
    input:
        agreement_files=get_agreement_files("signature_calls"),
        script=os.path.join(workflow.basedir, "..", "scripts", "aggregate_tool_agreement.py")
    output:
        summary_file("signature_calls")
    log:
        **summary_log("signature_calls")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    shell:
        "(python {input.script:q} \
        --input-files {input.agreement_files:q} \
        --output {output:q}) > {log.stdout:q} 2> {log.stderr:q}"


rule aggregate_sample_activity:
    """Combine per-cell (sample, tool) mutation-burden / active-signature-count
    TSVs into one table — same generic concat script as aggregate_tool_agreement."""
    input:
        agreement_files=get_agreement_files("sample_activity"),
        script=os.path.join(workflow.basedir, "..", "scripts", "aggregate_tool_agreement.py")
    output:
        summary_file("sample_activity")
    log:
        **summary_log("sample_activity")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    shell:
        "(python {input.script:q} \
        --input-files {input.agreement_files:q} \
        --output {output:q}) > {log.stdout:q} 2> {log.stderr:q}"


rule aggregate_signature_weights:
    """Combine per-cell (tool, signature, sample) -> weight TSVs into one table —
    same generic concat script as aggregate_tool_agreement. This is the raw
    per-sample data behind the aggregated signature_calls consensus, used to
    draw per-signature, per-tool boxplots of predicted activity."""
    input:
        agreement_files=get_agreement_files("signature_weights"),
        script=os.path.join(workflow.basedir, "..", "scripts", "aggregate_tool_agreement.py")
    output:
        summary_file("signature_weights")
    log:
        **summary_log("signature_weights")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    shell:
        "(python {input.script:q} \
        --input-files {input.agreement_files:q} \
        --output {output:q}) > {log.stdout:q} 2> {log.stderr:q}"


rule create_agreement_report:
    """HTML report for real samples (no known answer): compares the tools with each other."""
    input:
        tool_agreement=summary_file("tool_agreement"),
        signature_calls=summary_file("signature_calls"),
        sample_activity=summary_file("sample_activity"),
        signature_weights=summary_file("signature_weights"),
        script=os.path.join(workflow.basedir, "..", "scripts", "create_agreement_report.py")
    output:
        report_file("real")
    params:
        title="SigFitBench: Real samples",
        subtitle="Agreement between tools (no known answer)"
    log:
        **summary_log("report_real")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    shell:
        "(python {input.script:q} \
        --tool-agreement {input.tool_agreement:q} \
        --signature-calls {input.signature_calls:q} \
        --sample-activity {input.sample_activity:q} \
        --signature-weights {input.signature_weights:q} \
        --output {output:q} \
        --title {params.title:q} \
        --subtitle {params.subtitle:q}) > {log.stdout:q} 2> {log.stderr:q}"


# ==============================================================================
# REPORTING
# ==============================================================================

rule create_html_report:
    """HTML report for one data type with known answers (realistic or custom simulations)."""
    input:
        summaries=summary_file("evaluation_metrics"),
        per_sample_metrics=summary_file("per_sample_metrics"),
        tool_benchmarks=summary_file("tool_benchmarks"),
        config_file=workflow.configfiles[0] if workflow.configfiles else []
    output:
        report=report_file("{report}")
    wildcard_constraints:
        report="|".join(name for name, _ in REPORTS.values())
    params:
        data_type=lambda wildcards: report_data_type[wildcards.report],
        title=lambda wildcards: "SigFitBench: " + REPORTS[report_data_type[wildcards.report]][1],
        subtitle="Mutational signature fitting tools compared with the known answer",
        results_dir=config['output_dir']
    log:
        **summary_log("report_{report}")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 8000 * attempt,
    script:
        os.path.join(workflow.basedir, "..", "scripts", "create_html_report.py")
