"""Tool rules for signature analysis methods"""

import os

# ==============================================================================
# WILDCARD CONSTRAINTS
# ==============================================================================

wildcard_constraints:
    data_type=r"(synthetic|realsim|samples|external|custom)",
    cancer_type=r"[A-Za-z0-9_-]+",
    identifier=r"[A-Za-z0-9_-]+",
    num_muts=r"(\d+|none)",
    relevance=r"(standard|relevant)",
    rep=r"\d+"


rule mutationalpatterns:
    input:
        data=lambda wildcards: get_input_data('mutationalpatterns', wildcards),
        cosmic=lambda wildcards: get_cosmic_signatures("mutationalpatterns", wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type),
    output:
        result=tool_result("MutationalPatterns", input_format="csv_no_index")
    params:
        strict=lambda wildcards: get_tool_config("mutationalpatterns")["strict"]
    log:
        **tool_logs("MutationalPatterns")
    benchmark:
        tool_benchmark("MutationalPatterns")
    conda:
        os.path.join(workflow.basedir, "envs", "mutationalpatterns.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    script:
        os.path.join(workflow.basedir, "..", "scripts", "MutationalPatterns.R")


rule deconstructsigs:
    input:
        data=lambda wildcards: get_input_data('deconstructsigs', wildcards),
        cosmic=lambda wildcards: get_cosmic_signatures("deconstructsigs", wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type),
    output:
        result=tool_result("deconstructSigs", input_format="tab_with_index")
    log:
        **tool_logs("deconstructSigs")
    benchmark:
        tool_benchmark("deconstructSigs")
    conda:
        os.path.join(workflow.basedir, "envs", "deconstructsigs.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    script:
        os.path.join(workflow.basedir, "..", "scripts", "deconstructSigs.R")

rule mmsig:
    input:
        data=lambda wildcards: get_input_data('mmsig', wildcards),
        cosmic=lambda wildcards: get_cosmic_signatures("mmsig", wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type),
    output:
        result=tool_result("mmsig", input_format="csv_with_index")
    log:
        **tool_logs("mmsig")
    benchmark:
        tool_benchmark("mmsig")
    conda:
        os.path.join(workflow.basedir, "envs", "mmsig.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    script:
        os.path.join(workflow.basedir, "..", "scripts", "mmsig.R")


rule sigfit:
    input:
        data=lambda wildcards: get_input_data('sigfit', wildcards),
        cosmic=lambda wildcards: get_cosmic_signatures("sigfit", wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type),
    output:
        result=tool_result("sigfit", input_format="csv_with_index"),
        lower90=tool_result("sigfit", "contribution_lower90")
    log:
        **tool_logs("sigfit")
    benchmark:
        tool_benchmark("sigfit")
    conda:
        os.path.join(workflow.basedir, "envs", "sigfit.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    script:
        os.path.join(workflow.basedir, "..", "scripts", "sigfit.R")

rule siglasso:
    input:
        data=lambda wildcards: get_input_data('siglasso', wildcards),
        cosmic=lambda wildcards: get_cosmic_signatures("siglasso", wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type),
    output:
        result=tool_result("siglasso", input_format="csv_quoted")
    log:
        **tool_logs("siglasso")
    benchmark:
        tool_benchmark("siglasso")
    conda:
        os.path.join(workflow.basedir, "envs", "siglasso.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    script:
        os.path.join(workflow.basedir, "..", "scripts", "sigLASSO.R")

rule sigminer:
    input:
        data=lambda wildcards: get_input_data('sigminer', wildcards),
        cosmic=lambda wildcards: get_cosmic_signatures("sigminer", wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type),
    output:
        result=tool_result("sigminer", input_format="csv_with_index")
    params:
        method=lambda wildcards: get_tool_config("sigminer")["method"]
    log:
        **tool_logs("sigminer")
    benchmark:
        tool_benchmark("sigminer")
    conda:
        os.path.join(workflow.basedir, "envs", "sigminer.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    script:
        os.path.join(workflow.basedir, "..", "scripts", "sigminer.R")

rule signature_tools:
    input:
        data=lambda wildcards: get_input_data('signature_tools', wildcards),
        cosmic=lambda wildcards: get_cosmic_signatures("signature_tools", wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type),
        script=os.path.join(workflow.basedir, "..", "scripts", "signature_tools.R"),
    output:
        result=tool_result("signature_tools", input_format="signature_tools")
    params:
        method=get_tool_config("signature_tools").get("method", "KLD"),
        exposureFilterType=get_tool_config("signature_tools").get("exposureFilterType", "giniScaledThreshold"),
        nboot=get_tool_config("signature_tools").get("nboot", 100),
    log:
        **tool_logs("signature_tools")
    benchmark:
        tool_benchmark("signature_tools")
    conda:
        os.path.join(workflow.basedir, "envs", "signature_tools.yaml")
    threads: lambda wildcards: get_tool_config("signature_tools").get("nparallel", 1)
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    shell:
        """
        Rscript --vanilla "{input.script}" \
            --data "{input.data}" \
            --cosmic "{input.cosmic}" \
            --output "{output.result}" \
            --method "{params.method}" \
            --filter-type "{params.exposureFilterType}" \
            --nboot {params.nboot} \
            --nparallel {threads} \
            > {log.stdout} 2> {log.stderr}""" + ON_TOOL_FAILURE

rule sigspack:
    input:
        data=lambda wildcards: get_input_data('sigspack', wildcards),
        cosmic=lambda wildcards: get_cosmic_signatures("sigspack", wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type),
    output:
        result=tool_result("sigspack", input_format="csv_with_index")
    log:
        **tool_logs("sigspack")
    benchmark:
        tool_benchmark("sigspack")
    conda:
        os.path.join(workflow.basedir, "envs", "sigspack.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    script:
        os.path.join(workflow.basedir, "..", "scripts", "SigsPack.R")

rule yapsa:
    input:
        data=lambda wildcards: get_input_data('yapsa', wildcards),
        cosmic=lambda wildcards: get_cosmic_signatures("yapsa", wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type),
    output:
        result=tool_result("yapsa", input_format="csv_with_index")
    log:
        **tool_logs("yapsa")
    benchmark:
        tool_benchmark("yapsa")
    conda:
        os.path.join(workflow.basedir, "envs", "yapsa.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    script:
        os.path.join(workflow.basedir, "..", "scripts", "YAPSA.R")


# ==============================================================================
# PYTHON-BASED TOOLS
# ==============================================================================

rule musical:
    input:
        data=lambda wildcards: get_input_data('musical', wildcards),
        weights=lambda wildcards: get_cosmic_signatures("musical", wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type),
        matrix_order=config["matrix_order"],
    output:
        result=tool_result("MuSiCal", input_format="csv_no_index")
    log:
        **tool_logs("MuSiCal")
    benchmark:
        tool_benchmark("MuSiCal")
    conda:
        os.path.join(workflow.basedir, "envs", "musical.yaml")
    script:
        os.path.join(workflow.basedir, "..", "scripts", "MuSiCal.py")


rule spa:
    input:
        data=lambda wildcards: get_input_data('spa', wildcards),
        signature_db=lambda wildcards: get_cosmic_signatures("spa", wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type),
    output:
        result=tool_result("SPA", input_format="tab_with_index")
    params:
        initial_remove_penalty=0.05
    log:
        **tool_logs("SPA")
    benchmark:
        tool_benchmark("SPA")
    conda:
        os.path.join(workflow.basedir, "envs", "spa.yaml")
    shadow:
        "minimal"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 8000 * attempt,
    script:
        os.path.join(workflow.basedir, "..", "scripts", "SPA.py")

rule spa_low:
    input:
        data=lambda wildcards: get_input_data('spa', wildcards),
        signature_db=lambda wildcards: get_cosmic_signatures("spa", wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type),
    output:
        result=tool_result("SPA_low", input_format="tab_with_index")
    params:
        initial_remove_penalty=0.02
    log:
        **tool_logs("SPA_low")
    benchmark:
        tool_benchmark("SPA_low")
    conda:
        os.path.join(workflow.basedir, "envs", "spa.yaml")
    shadow:
        "minimal"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    script:
        os.path.join(workflow.basedir, "..", "scripts", "SPA.py")

rule fspa:
    """GitHub-packaged FSPA (sirenfrappe/FSPA, installed via envs/fspa.yaml --
    see that file)"""
    input:
        data=lambda wildcards: get_input_data('fspa', wildcards),
        signature_db=lambda wildcards: get_cosmic_signatures("fspa", wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type),
    output:
        result=tool_result("FSPA", input_format="tab_with_index")
    params:
        initial_remove_penalty=0.05
    log:
        **tool_logs("FSPA")
    benchmark:
        tool_benchmark("FSPA")
    conda:
        os.path.join(workflow.basedir, "envs", "fspa.yaml")
    shadow:
        "minimal"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    script:
        os.path.join(workflow.basedir, "..", "scripts", "SPA.py")

rule beam_search:
    """GitHub-packaged beam_search FSPA (sirenfrappe/SigProfilerAssignment-beam-search, installed via envs/beam_search.yaml --
    see that file)"""
    input:
        data=lambda wildcards: get_input_data('beam_search', wildcards),
        signature_db=lambda wildcards: get_cosmic_signatures("beam_search", wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type),
    output:
        result=tool_result("beam_search", input_format="tab_with_index")
    params:
        initial_remove_penalty=0.05
    log:
        **tool_logs("beam_search")
    benchmark:
        tool_benchmark("beam_search")
    conda:
        os.path.join(workflow.basedir, "envs", "beam_search.yaml")
    shadow:
        "minimal"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    script:
        os.path.join(workflow.basedir, "..", "scripts", "beam_search.py")


rule msigact:
    """mSigAct (github.com/steverozen/mSigAct): maximum-likelihood fit (negative binomial by
    default), then removal of signatures while a likelihood-ratio test allows it."""
    input:
        data=lambda wildcards: get_input_data('msigact', wildcards),
        signature_db=lambda wildcards: get_cosmic_signatures("msigact", wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type),
    output:
        result=tool_result("mSigAct", input_format="tab_with_index")
    params:
        max_level=get_tool_config("msigact").get("max_level", 5),
        p_thresh=get_tool_config("msigact").get("p_thresh", 0.05),
        likelihood=get_tool_config("msigact").get("likelihood", "neg.binom"),
        # any other SparseAssignActivity arguments, by their mSigAct names (e.g. max.subsets)
        extra=get_tool_config("msigact").get("extra", {})
    log:
        **tool_logs("mSigAct")
    benchmark:
        tool_benchmark("mSigAct")
    conda:
        os.path.join(workflow.basedir, "envs", "msigact.yaml")
    shadow:
        "minimal"
    threads: get_tool_config("msigact").get("threads", 4)
    resources:
        mem_mb=lambda wildcards, attempt: 8000 * attempt,
    script:
        os.path.join(workflow.basedir, "..", "scripts", "mSigAct.R")


rule starsigndna:
    """StarSignDNA (github.com/uio-bmi/StarSignDNA): Poisson likelihood with an L1 penalty;
    with several samples, signatures are kept by their average exposure across the samples."""
    input:
        data=lambda wildcards: get_input_data('starsigndna', wildcards),
        signature_db=lambda wildcards: get_cosmic_signatures("starsigndna", wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type),
    output:
        result=tool_result("StarSignDNA", input_format="tab_with_index")
    params:
        n_iterations=get_tool_config("starsigndna").get("n_iterations", 1000),
        data_type="exome" if config.get("wgs_or_wes", "WGS") == "WES" else "genome"
    log:
        **tool_logs("StarSignDNA")
    benchmark:
        tool_benchmark("StarSignDNA")
    conda:
        os.path.join(workflow.basedir, "envs", "starsigndna.yaml")
    shadow:
        "minimal"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    script:
        os.path.join(workflow.basedir, "..", "scripts", "StarSignDNA.py")


rule sigtrim:
    """
    SigTrim: iterative z-score based unassigned-mutation trimming, followed by an
    SPA-style coarse pre-filter and greedy add/remove of known signatures.
    Runs the installed package (github.com/mkatsanto/SigTrim, pip-installed into
    sigtrim.yaml) via its `sigtrim` CLI entry point.
    """
    input:
        data=lambda wildcards: get_input_data('sigtrim', wildcards),
        weights=lambda wildcards: get_cosmic_signatures("sigtrim", wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type),
        matrix_order=config["matrix_order"],
    output:
        result=tool_result("SigTrim", input_format="csv_no_index")
    params:
        zscore_threshold=get_tool_config("sigtrim").get("zscore_threshold", 1.0),
        thresh_backward=get_tool_config("sigtrim").get("thresh_backward", 0.01),
        thresh_forward=get_tool_config("sigtrim").get("thresh_forward", 0.01),
        initial_remove_penalty=get_tool_config("sigtrim").get("initial_remove_penalty", 0.05),
        max_iter=get_tool_config("sigtrim").get("max_iter", 1000)
    log:
        **tool_logs("SigTrim")
    benchmark:
        tool_benchmark("SigTrim")
    conda:
        os.path.join(workflow.basedir, "envs", "sigtrim.yaml")
    shadow:
        "minimal"   # SigTrim also writes convergence_<sample>.csv files; keep them out of the results
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    shell:
        "(sigtrim \
        -i {input.data} \
        -w {input.weights} \
        -m {input.matrix_order} \
        -o {output.result} \
        --zscore-threshold {params.zscore_threshold} \
        --thresh-backward {params.thresh_backward} \
        --thresh-forward {params.thresh_forward} \
        --initial-remove-penalty {params.initial_remove_penalty} \
        --max-iter {params.max_iter}) > {log.stdout} 2> {log.stderr}" + ON_TOOL_FAILURE


rule sigid:
    """
    SigID: SPA-style greedy add/remove of known signatures with automatic per-sample
    L1 penalty selection for unknown-signature discovery. Produces a single
    contribution file (no penalty sweep required).
    Runs the installed package (github.com/mkatsanto/SigID, pip-installed into
    sigid.yaml) via its `sigid` CLI entry point.
    """
    input:
        data=lambda wildcards: get_input_data('sigid', wildcards),
        weights=lambda wildcards: get_cosmic_signatures("sigid", wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type),
        matrix_order=config["matrix_order"],
        # true_weights is optional for sigid itself (used only to log a diagnostic TAE,
        # never influences the fit — see SigID/sigid/cli.py's --true-weights help text);
        # no weights_*.dat is ever produced for realsim/external (no ground truth), so
        # this resolves to no input at all for those, rather than a missing-file error.
        true_weights=lambda wildcards: (
            os.path.join(
                config['output_dir'], "generated_data", wildcards.data_type,
                f"weights_{wildcards.cancer_type}-w_{'0' if wildcards.data_type == 'realsim' else wildcards.rep}-{wildcards.num_muts}.dat"
            ) if wildcards.data_type not in ("realsim", "external") else []
        ),
    output:
        result=tool_result("SigID", input_format="csv_no_index")
    params:
        thresh_backward=get_tool_config("sigid").get("thresh_backward", 0.001),
        max_iter=get_tool_config("sigid").get("max_iter", 1000),
        min_penalty=get_tool_config("sigid").get("min_penalty", 0.001),
        max_penalty=get_tool_config("sigid").get("max_penalty", 10000.0),
        initial_remove_penalty=get_tool_config("sigid").get("initial_remove_penalty", 0.05),
        true_weights_arg=lambda wildcards, input: (
            f"--true-weights {input.true_weights}" if input.true_weights else ""
        ),
    log:
        **tool_logs("SigID")
    benchmark:
        tool_benchmark("SigID")
    conda:
        os.path.join(workflow.basedir, "envs", "sigid.yaml")
    threads: 4
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    shell:
        "(sigid \
        -i {input.data} \
        -w {input.weights} \
        -m {input.matrix_order} \
        -o {output.result} \
        --thresh-backward {params.thresh_backward} \
        --max-iter {params.max_iter} \
        --min-penalty {params.min_penalty} \
        --max-penalty {params.max_penalty} \
        --initial-remove-penalty {params.initial_remove_penalty} \
        {params.true_weights_arg} \
        --n-jobs {threads}) > {log.stdout} 2> {log.stderr}" + ON_TOOL_FAILURE
