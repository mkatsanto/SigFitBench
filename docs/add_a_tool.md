# Add your own tool

If your tool is a Python package on GitHub, adding it takes two files: an environment
that installs it, and a rule that runs it. Your tool then gets the same samples and
scores as the built-in tools and appears in the same report.

This is how `SigID` and `SPA_tree` are added; see their rules in
[`workflow/rules/tools.smk`](https://github.com/mkatsanto/SigFitBench/blob/main/workflow/rules/tools.smk).

## 1. The rule

Add this to the end of `workflow/rules/tools.smk`. It runs your tool's command, like the
`sigid` rule:

```python
rule mytool:
    input:
        data=lambda wildcards: get_input_data("mytool", wildcards),
        cosmic=lambda wildcards: get_cosmic_signatures("mytool", wildcards.relevance, data_type=wildcards.data_type, dataset=wildcards.cancer_type),
    output:
        result=tool_result("MyTool", input_format="tab_with_index")
    log:
        **tool_logs("MyTool")
    benchmark:
        tool_benchmark("MyTool")
    conda:
        os.path.join(workflow.basedir, "envs", "mytool.yaml")
    shell:
        "(mytool --counts {input.data:q} --signatures {input.cosmic:q} --output {output.result:q}"
        + extra_flags("mytool") + ") > {log.stdout:q} 2> {log.stderr:q}" + ON_TOOL_FAILURE
```

Replace `mytool` with your tool's name: in lower case for the rule and in quotes after
`get_input_data`, `get_cosmic_signatures` and `extra_flags`; as you want it shown in the
report for `MyTool`. Replace the command and its options with your tool's own.

## 2. The environment

Save as `workflow/envs/mytool.yaml`. The `pip` line installs your tool from GitHub:

```yaml
channels:
  - conda-forge
dependencies:
  - python>=3.9
  - pandas
  - pip
  - pip:
    - git+https://github.com/me/mytool.git@v1.0
```

Give a release tag (`@v1.0`) or a commit, not a branch, so the same version runs every
time.

## 3. Run it

Add `"MyTool"` to `tools_to_run` in a config, for example a copy of
`config/scenarios/quick_test.yaml`, and run it as usual. Snakemake builds the
environment the first time.

---

## What your tool reads and writes

| | File | Layout |
|---|---|---|
| Input | `{input.data}` | Tab-separated. One row per mutation type (96, first column `Type`), one column per sample, mutation counts |
| Input | `{input.cosmic}` | Tab-separated. One row per mutation type, one column per reference signature |
| Output | `{output.result}` | Comma-separated. One row per signature (first column = signature name), one column per sample |

The output must give the **number of mutations** for each signature in each sample, not
shares, and use the signature names from the reference.

## If your tool has no command

Use a short script that imports it instead, like the `spa_tree` rule does. Replace
`shell:` with:

```python
    script:
        os.path.join(workflow.basedir, "..", "scripts", "mytool.py")
```

and add `**extra_params("mytool")` under `params:`. The script finds the files in
`snakemake.input.data`, `snakemake.input.cosmic` and `snakemake.output.result`; see
[`scripts/SPA_tree.py`](https://github.com/mkatsanto/SigFitBench/blob/main/scripts/SPA_tree.py).

## Settings

Options under `tool_overrides` → `mytool` → `extra` in the config are added to the
command (or reach the script as `snakemake.params.extra`):

```yaml
tool_overrides:
  mytool:
    extra:
      --max-iter: 500
```
