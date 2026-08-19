"""
Rules for identifying and managing relevant (active) signatures.

This module implements two-pass signature analysis:
1. Standard run: Fit with full signature set
2. Identify consistently active signatures
3. Relevant run: Re-fit with reduced signature set
"""

import os


rule identify_active_signatures:
    """
    Identify consistently active signatures from standard run results.
    
    This rule analyzes tool outputs to find signatures that are:
    - Active in multiple samples
    - Detected by multiple tools (consensus)
    - Contributing significantly to mutation burden
    
    Used for AUTO mode relevant signature detection.
    """
    input:
        # Standard run results for all selected tools
        tool_results=expand(
            os.path.join(
                config['output_dir'],
                "{{data_type}}_{{num_muts}}_standard-{tool}-contribution.dat"),
            tool=get_selected_tools()
        ),
        script=os.path.join(workflow.basedir, "..", "scripts", "identify_active_signatures.py")
    output:
        active_list=os.path.join(
            config['output_dir'],
            "relevant_signatures",
            "{data_type}_{num_muts}_active_signatures.txt"
        )
    params:
        min_contribution=lambda wildcards: config["relevant_detection"]["min_contribution"],
        min_sample_fraction=lambda wildcards: config["relevant_detection"]["min_sample_fraction"],
        min_avg_contribution=lambda wildcards: config["relevant_detection"]["min_avg_contribution"],
        consensus_tools=lambda wildcards: config["relevant_detection"]["consensus_tools"]
    conda:
        os.path.join(workflow.basedir, "envs", "python.yaml")
    log:
        stderr=os.path.join(
            config["log_dir"],
            "identify_active_signatures.{data_type}.{num_muts}.stderr.log"),
        stdout=os.path.join(
            config["log_dir"],
            "identify_active_signatures.{data_type}.{num_muts}.stdout.log")
    shell:
        """
        (python '{input.script}'  \
            --tool-results {input.tool_results} \
            --output {output.active_list} \
            --min-contribution {params.min_contribution} \
            --min-sample-fraction {params.min_sample_fraction} \
            --min-avg-contribution {params.min_avg_contribution} \
            --consensus-tools {params.consensus_tools} \
            --verbose) > {log.stdout} 2> {log.stderr}
        """


