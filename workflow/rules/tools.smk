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
    rep=r"\d+",
    penalty=r"[0-9.]+"


rule mutationalpatterns:
    input:
        data=lambda wildcards: get_input_data('mutationalpatterns', wildcards),
        cosmic=lambda wildcards: get_cosmic_signatures("mutationalpatterns", wildcards.relevance),
    output:
        result=os.path.join(
            config['output_dir'],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "MutationalPatterns_contribution.dat"
        )
    params:
        strict=lambda wildcards: get_tool_config("mutationalpatterns")["strict"]
    log:
        stderr=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "mutationalpatterns_stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "mutationalpatterns_stdout.log"
        )
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
        cosmic=lambda wildcards: get_cosmic_signatures("deconstructsigs", wildcards.relevance),
    output:
        result=os.path.join(
            config['output_dir'],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "deconstructSigs_contribution.dat"
        )
    log:
        stderr=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "deconstructsigs_stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "deconstructsigs_stdout.log"
        )
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
        cosmic=lambda wildcards: get_cosmic_signatures("mmsig", wildcards.relevance),
    output:
        result=os.path.join(
            config['output_dir'],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "mmsig_contribution.dat"
        )
    log:
        stderr=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "mmsig_stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "mmsig_stdout.log"
        )
    conda:
        os.path.join(workflow.basedir, "envs", "r-base.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    script:
        os.path.join(workflow.basedir, "..", "scripts", "mmsig.R")


rule sigfit:
    input:
        data=lambda wildcards: get_input_data('sigfit', wildcards),
        cosmic=lambda wildcards: get_cosmic_signatures("sigfit", wildcards.relevance),
    output:
        result=os.path.join(
            config['output_dir'],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "sigfit_contribution.dat"
        ),
        lower90=os.path.join(
            config['output_dir'],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "sigfit_contribution_lower90.dat"
        )
    log:
        stderr=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "sigfit_stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "sigfit_stdout.log"
        )
    conda:
        os.path.join(workflow.basedir, "envs", "r-base.yaml")
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    script:
        os.path.join(workflow.basedir, "..", "scripts", "sigfit.R")

rule siglasso:
    input:
        data=lambda wildcards: get_input_data('siglasso', wildcards),
        cosmic=lambda wildcards: get_cosmic_signatures("siglasso", wildcards.relevance),
    output:
        result=os.path.join(
            config['output_dir'],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "siglasso_contribution.dat"
        )
    log:
        stderr=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "siglasso_stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "siglasso_stdout.log"
        )
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
        cosmic=lambda wildcards: get_cosmic_signatures("sigminer", wildcards.relevance),
    output:
        result=os.path.join(
            config['output_dir'],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "sigminer_contribution.dat"
        )
    params:
        method=lambda wildcards: get_tool_config("sigminer")["method"]
    log:
        stderr=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "sigminer_stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "sigminer_stdout.log"
        )
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
        cosmic=lambda wildcards: get_cosmic_signatures("signature_tools", wildcards.relevance),
        script=os.path.join(workflow.basedir, "..", "scripts", "signature_tools.R"),
    output:
        result=os.path.join(
            config['output_dir'],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "signature_tools_contribution.dat"
        )
    params:
        method=get_tool_config("signature_tools").get("method", "KLD"),
        exposureFilterType=get_tool_config("signature_tools").get("exposureFilterType", "giniScaledThreshold"),
        nboot=get_tool_config("signature_tools").get("nboot", 100),
    log:
        stderr=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "signature_tools_stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "signature_tools_stdout.log"
        )
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
            > {log.stdout} 2> {log.stderr}
        """

rule sigspack:
    input:
        data=lambda wildcards: get_input_data('sigspack', wildcards),
        cosmic=lambda wildcards: get_cosmic_signatures("sigspack", wildcards.relevance),
    output:
        result=os.path.join(
            config['output_dir'],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "sigspack_contribution.dat"
        )
    log:
        stderr=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "sigspack_stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "sigspack_stdout.log"
        )
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
        cosmic=lambda wildcards: get_cosmic_signatures("yapsa", wildcards.relevance),
    output:
        result=os.path.join(
            config['output_dir'],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "yapsa_contribution.dat"
        )
    log:
        stderr=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "yapsa_stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "yapsa_stdout.log"
        )
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
        weights=lambda wildcards: get_cosmic_signatures("musical", wildcards.relevance),
        matrix_order=config["matrix_order"],
    output:
        result=os.path.join(
            config['output_dir'],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "MuSiCal_contribution.dat"
        )
    log:
        stderr=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "musical_stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "musical_stdout.log"
        )
    conda:
        os.path.join(workflow.basedir, "envs", "musical.yaml")
    script:
        os.path.join(workflow.basedir, "..", "scripts", "MuSiCal.py")


rule spa:
    input:
        data=lambda wildcards: get_input_data('spa', wildcards),
        signature_db=lambda wildcards: get_cosmic_signatures("spa", wildcards.relevance),
    output:
        result=os.path.join(
            config['output_dir'],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "SPA_contribution.dat"
        )
    params:
        initial_remove_penalty=0.05
    log:
        stderr=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "spa_stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "spa_stdout.log"
        )
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
        signature_db=lambda wildcards: get_cosmic_signatures("spa", wildcards.relevance),
    output:
        result=os.path.join(
            config['output_dir'],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "SPA_low_contribution.dat"
        )
    params:
        initial_remove_penalty=0.01
    log:
        stderr=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "spa_low_stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "spa_low_stdout.log"
        )
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
    input:
        data=lambda wildcards: get_input_data('fspa', wildcards),
        signature_db=lambda wildcards: get_cosmic_signatures("fspa", wildcards.relevance),
    output:
        result=os.path.join(
            config['output_dir'],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "FSPA_contribution.dat"
        )
    params:
        initial_remove_penalty=0.05
    log:
        stderr=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "fspa_stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "fspa_stdout.log"
        )
    conda:
        os.path.join(workflow.basedir, "envs", "fspa.yaml")
    shadow:
        "minimal"
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt,
    script:
        os.path.join(workflow.basedir, "..", "scripts", "SPA.py")


rule sigtrim:
    """
    SigTrim: iterative z-score based unassigned-mutation trimming, followed by an
    SPA-style coarse pre-filter and greedy add/remove of known signatures.
    Runs the installed package (github.com/mkatsanto/SigTrim, pip-installed into
    sigtrim.yaml) via its `sigtrim` CLI entry point.
    """
    input:
        data=lambda wildcards: get_input_data('sigtrim', wildcards),
        weights=lambda wildcards: get_cosmic_signatures("sigtrim", wildcards.relevance),
        matrix_order=config["matrix_order"],
    output:
        result=os.path.join(
            config['output_dir'],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "SigTrim_contribution.dat"
        )
    params:
        zscore_threshold=get_tool_config("sigtrim")["zscore_threshold"],
        thresh_backward=get_tool_config("sigtrim")["thresh_backward"],
        thresh_forward=get_tool_config("sigtrim")["thresh_forward"],
        initial_remove_penalty=get_tool_config("sigtrim")["initial_remove_penalty"],
        max_iter=get_tool_config("sigtrim")["max_iter"]
    log:
        stderr=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "sigtrim_stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "sigtrim_stdout.log"
        )
    conda:
        os.path.join(workflow.basedir, "envs", "sigtrim.yaml")
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
        --max-iter {params.max_iter}) > {log.stdout} 2> {log.stderr}"


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
        weights=lambda wildcards: get_cosmic_signatures("sigid", wildcards.relevance),
        matrix_order=config["matrix_order"],
        true_weights=lambda wildcards: os.path.join(
            config['output_dir'], "generated_data", wildcards.data_type,
            f"weights_{wildcards.cancer_type}-w_{'0' if wildcards.data_type == 'realsim' else wildcards.rep}-{wildcards.num_muts}.dat"
        ),
    output:
        result=os.path.join(
            config['output_dir'],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "SigID_contribution.dat"
        )
    params:
        thresh_backward=get_tool_config("sigid").get("thresh_backward", 0.001),
        max_iter=get_tool_config("sigid").get("max_iter", 1000),
        min_penalty=get_tool_config("sigid").get("min_penalty", 0.001),
        max_penalty=get_tool_config("sigid").get("max_penalty", 10000.0),
        initial_remove_penalty=get_tool_config("sigid").get("initial_remove_penalty", 0.05),
    log:
        stderr=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "sigid_stderr.log"
        ),
        stdout=os.path.join(
            config["log_dir"],
            "{data_type}",
            "{cancer_type}",
            "{relevance}",
            "{num_muts}",
            "rep_{rep}",
            "sigid_stdout.log"
        )
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
        --true-weights {input.true_weights} \
        --n-jobs {threads}) > {log.stdout} 2> {log.stderr}"
