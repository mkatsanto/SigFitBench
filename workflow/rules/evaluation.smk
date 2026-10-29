"""Evaluation rules for signature analysis results"""

import os

# ==============================================================================
# WILDCARD CONSTRAINTS
# ==============================================================================

wildcard_constraints:
    data_type=r"(realsim|real|samples|custom)",
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
        tool_results=tool_result("{tool}"),
        true_weights=lambda wildcards: os.path.join(
            config['output_dir'], "generated_data", wildcards.data_type,
            f"weights_{wildcards.cancer_type}-w_{wildcards.rep}-{wildcards.num_muts}.dat"
        ),
        reference_signatures=lambda wildcards: get_cosmic_signatures(
            wildcards.tool, wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type
        ),
        cosmic_file=lambda wildcards: signature_catalogue(wildcards.data_type, wildcards.cancer_type),
        script=os.path.join(workflow.basedir, "..", "scripts", "evaluate_tool_performance.py")
    output:
        per_sample_metrics=os.path.join(
            config['output_dir'], "evaluation",
            "{data_type}_{cancer_type}_{num_muts}_{relevance}_{tool}_w{rep}_per_sample.tsv"
        )
    params:
        min_contribution=config["evaluation"]["min_contribution"],
        out_of_reference_sigs=OUT_OF_REF_SIGS_STR
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
        --out-of-reference-sigs {params.out_of_reference_sigs}) \
        > {log.stdout} 2> {log.stderr}"


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


# Benchmark TSVs for the uniform (non-custom) data types, built the same way get_eval_files()
# and rule all build their equivalent lists.
_benchmark_files = expand(
    tool_benchmark("{tool}"),
    data_type=active_data_types,
    cancer_type=config.get("cancer_types", ["Head-SCC"]),
    num_muts=config["num_muts"],
    relevance=config.get("relevance_modes", ["standard"]),
    tool=get_selected_tools(),
    rep=range(config.get("num_realizations", 1)),
)
# "custom" isn't part of active_data_types (dataset name stands in for cancer_type, no
# num_muts/rep sweep shared with the others), so it's appended separately -- same reason
# get_eval_files() and rule all's custom_tool_results handle it as its own list
# comprehension. Without this, a custom-only run (data_types: ["custom"]) leaves
# active_data_types empty and this rule would get zero input files.
if "custom" in config.get("data_types", []):
    _benchmark_files += expand(
        tool_benchmark("{tool}"),
        data_type="custom",
        cancer_type=list(config.get("custom_datasets", {})),
        relevance=config.get("relevance_modes", ["standard"]),
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
        os.path.join(config['output_dir'], "all_tool_benchmarks.tsv")
    log:
        stderr=os.path.join(config["log_dir"], "aggregate_tool_benchmarks.stderr.log"),
        stdout=os.path.join(config["log_dir"], "aggregate_tool_benchmarks.stdout.log")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    shell:
        "(python '{input.script}' \
        --input-files {input.benchmark_files} \
        --output {output}) \
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
            tool=get_selected_tools(),
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
# CROSS-TOOL AGREEMENT (no-ground-truth data type: real)
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
        original_catalog=lambda wildcards: get_input_data(get_agreement_tools()[0], wildcards),
        script=os.path.join(workflow.basedir, "..", "scripts", "evaluate_tool_agreement.py")
    output:
        pairwise=os.path.join(
            config['output_dir'], "evaluation",
            "{data_type}_{cancer_type}_{num_muts}_{relevance}_tool_agreement.tsv"
        ),
        signature_calls=os.path.join(
            config['output_dir'], "evaluation",
            "{data_type}_{cancer_type}_{num_muts}_{relevance}_signature_calls.tsv"
        ),
        sample_activity=os.path.join(
            config['output_dir'], "evaluation",
            "{data_type}_{cancer_type}_{num_muts}_{relevance}_sample_activity.tsv"
        ),
        signature_weights=os.path.join(
            config['output_dir'], "evaluation",
            "{data_type}_{cancer_type}_{num_muts}_{relevance}_signature_weights.tsv"
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
        --output {output.pairwise} \
        --output-signature-calls {output.signature_calls} \
        --output-sample-activity {output.sample_activity} \
        --output-signature-weights {output.signature_weights} \
        --original-catalog {input.original_catalog}) > {log.stdout} 2> {log.stderr}"


rule aggregate_tool_agreement:
    """Combine per-cell agreement TSVs from all real datasets into one table."""
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


rule aggregate_signature_calls:
    """Combine per-cell signature-call-consensus TSVs (which tools call which
    signatures active, and how often) into one table — same generic concat
    script as aggregate_tool_agreement, just a different input/output pair."""
    input:
        agreement_files=get_agreement_files("signature_calls"),
        script=os.path.join(workflow.basedir, "..", "scripts", "aggregate_tool_agreement.py")
    output:
        os.path.join(config['output_dir'], "all_signature_calls.tsv")
    log:
        stderr=os.path.join(config["log_dir"], "aggregate_signature_calls.stderr.log"),
        stdout=os.path.join(config["log_dir"], "aggregate_signature_calls.stdout.log")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    shell:
        "(python '{input.script}' \
        --input-files {input.agreement_files} \
        --output {output}) > {log.stdout} 2> {log.stderr}"


rule aggregate_sample_activity:
    """Combine per-cell (sample, tool) mutation-burden / active-signature-count
    TSVs into one table — same generic concat script as aggregate_tool_agreement."""
    input:
        agreement_files=get_agreement_files("sample_activity"),
        script=os.path.join(workflow.basedir, "..", "scripts", "aggregate_tool_agreement.py")
    output:
        os.path.join(config['output_dir'], "all_sample_activity.tsv")
    log:
        stderr=os.path.join(config["log_dir"], "aggregate_sample_activity.stderr.log"),
        stdout=os.path.join(config["log_dir"], "aggregate_sample_activity.stdout.log")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    shell:
        "(python '{input.script}' \
        --input-files {input.agreement_files} \
        --output {output}) > {log.stdout} 2> {log.stderr}"


rule aggregate_signature_weights:
    """Combine per-cell (tool, signature, sample) -> weight TSVs into one table —
    same generic concat script as aggregate_tool_agreement. This is the raw
    per-sample data behind the aggregated signature_calls consensus, used to
    draw per-signature, per-tool boxplots of predicted activity."""
    input:
        agreement_files=get_agreement_files("signature_weights"),
        script=os.path.join(workflow.basedir, "..", "scripts", "aggregate_tool_agreement.py")
    output:
        os.path.join(config['output_dir'], "all_signature_weights.tsv")
    log:
        stderr=os.path.join(config["log_dir"], "aggregate_signature_weights.stderr.log"),
        stdout=os.path.join(config["log_dir"], "aggregate_signature_weights.stdout.log")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    shell:
        "(python '{input.script}' \
        --input-files {input.agreement_files} \
        --output {output}) > {log.stdout} 2> {log.stderr}"


rule create_agreement_report:
    """HTML report for real samples (no known answer): compares the tools with each other."""
    input:
        tool_agreement=os.path.join(config['output_dir'], "all_tool_agreement.tsv"),
        signature_calls=os.path.join(config['output_dir'], "all_signature_calls.tsv"),
        sample_activity=os.path.join(config['output_dir'], "all_sample_activity.tsv"),
        signature_weights=os.path.join(config['output_dir'], "all_signature_weights.tsv"),
        script=os.path.join(workflow.basedir, "..", "scripts", "create_agreement_report.py")
    output:
        os.path.join(config['output_dir'], "sigfitbench_report_real_samples.html")
    params:
        title="SigFitBench: Real samples",
        subtitle="Agreement between tools (no known answer)"
    log:
        stderr=os.path.join(config["log_dir"], "agreement_report.stderr.log"),
        stdout=os.path.join(config["log_dir"], "agreement_report.stdout.log")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 2000 * attempt,
    shell:
        "(python '{input.script}' \
        --tool-agreement {input.tool_agreement} \
        --signature-calls {input.signature_calls} \
        --sample-activity {input.sample_activity} \
        --signature-weights {input.signature_weights} \
        --output {output} \
        --title '{params.title}' \
        --subtitle '{params.subtitle}') > {log.stdout} 2> {log.stderr}"


# ==============================================================================
# REPORTING
# ==============================================================================

rule create_html_report:
    """HTML report for one data type with known answers (realistic simulations or synthetic)."""
    input:
        summaries=os.path.join(config['output_dir'], "all_evaluation_metrics.tsv"),
        per_sample_metrics=os.path.join(config['output_dir'], "all_per_sample_metrics.tsv"),
        tool_benchmarks=os.path.join(config['output_dir'], "all_tool_benchmarks.tsv"),
        config_file=workflow.configfiles[0] if workflow.configfiles else []
    output:
        report=os.path.join(config['output_dir'], "sigfitbench_report_{report}.html")
    wildcard_constraints:
        report="|".join(name for name, _ in REPORTS.values())
    params:
        data_type=lambda wildcards: report_data_type[wildcards.report],
        title=lambda wildcards: "SigFitBench: " + REPORTS[report_data_type[wildcards.report]][1],
        subtitle="Mutational signature fitting tools compared with the known answer",
        results_dir=config['output_dir']
    log:
        stderr=os.path.join(config["log_dir"], "html_report_{report}.stderr.log"),
        stdout=os.path.join(config["log_dir"], "html_report_{report}.stdout.log")
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    script:
        os.path.join(workflow.basedir, "..", "scripts", "create_html_report.py")
