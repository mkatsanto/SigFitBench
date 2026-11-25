"""Common functions and utilities for SigFitTest workflow

This module provides helper functions for:
- Tool configuration and inheritance
- Input/output path resolution
- COSMIC signature file management
- Tool selection and filtering
- Parameter inheritance across config sections
"""

import os
import re
import shlex
import sys


# ==============================================================================
# DATA TYPES
#   synth_realistic - simulated from signature activity in real tumours (known answer)
#   synth_custom    - simulated from the user's own table of signature mixes (known answer)
#   real            - real samples from real_datasets (no known answer: tools are
#                     compared with each other)
# active_data_types holds the realistic simulations only: custom and real runs are named
# after their datasets, not after the cancer_types list, so they are handled separately.
# ==============================================================================

def _peek_catalog_header(catalog_path):
    """Return the first non-comment, non-blank line of a catalog file (its header row),
    or None if the file can't be read. Used only to detect column-naming convention —
    doesn't load the full catalog (that's prepare_real_catalog.py's job at run time)."""
    if not catalog_path or not os.path.exists(catalog_path):
        return None
    with open(catalog_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                return line
    return None


def real_catalog_path(ds_name):
    """The catalog of one real_datasets entry: its `catalog`, or, if it gives a `vcf_dir` instead, the catalog the workflow
    builds from those VCFs (rule catalog_from_vcfs in vcf_input.smk)."""
    ds = config.get("real_datasets", {}).get(ds_name, {}) or {}
    if ds.get("catalog"):
        return ds["catalog"]
    if ds.get("vcf_dir"):
        return os.path.join(config["output_dir"], "internal", "catalogs_from_vcfs", f"{ds_name}_sbs96.dat")
    raise ValueError(f"real_datasets entry {ds_name!r} needs a 'catalog' or a 'vcf_dir'")


_REAL_CANCER_TYPES_CACHE = {}


def get_real_cancer_types(ds_name):
    """Cancer types (tissues) for one real_datasets entry.

    If the entry lists cancer_types, those are used. Otherwise they are read from the
    catalog's column headers, which may follow the "Tissue::SampleID" convention (only
    the header row is read). If any column does, every tissue found is used (one group
    per tissue). If none do, returns [ds_name]: the whole dataset is fit as one group
    named after the dataset key.
    """
    if ds_name in _REAL_CANCER_TYPES_CACHE:
        return _REAL_CANCER_TYPES_CACHE[ds_name]

    ds = config.get("real_datasets", {}).get(ds_name, {}) or {}
    result = list(ds.get("cancer_types") or [])
    if not result:
        header = _peek_catalog_header(real_catalog_path(ds_name))      # not there yet if built from VCFs: one group
        result = [ds_name]
        if header:
            delim = "\t" if "\t" in header else ("," if "," in header else None)
            cols = header.split(delim)[1:] if delim else header.split()[1:]  # skip index column
            tissues = sorted({c.strip().strip('"') .split("::")[0] for c in cols if "::" in c})
            if tissues:
                result = tissues

    _REAL_CANCER_TYPES_CACHE[ds_name] = result
    return result


def get_all_real_cancer_types():
    """Flat list of "cancer_type" wildcard values across every real_datasets entry —
    each dataset contributes either its own key (one group) or one entry per tissue
    (get_real_cancer_types). Used as the cancer_type dimension in rule all."""
    return [
        ct for ds_name in config.get("real_datasets", {})
        for ct in get_real_cancer_types(ds_name)
    ]


def resolve_real_dataset(cancer_type):
    """Map a "cancer_type" wildcard back to (dataset_name, tissue_prefix_or_None).

    tissue_prefix is None when the matched dataset is fit as one group
    (cancer_type == the dataset key itself); otherwise it's the tissue to filter the
    catalog to (passed to prepare_real_catalog.py's --tissue-prefix).
    """
    for ds_name in config.get("real_datasets", {}):
        tissues = get_real_cancer_types(ds_name)
        if tissues == [ds_name] and cancer_type == ds_name:
            return ds_name, None
        if cancer_type in tissues:
            return ds_name, cancer_type
    raise ValueError(
        f"cancer_type {cancer_type!r} does not match any real_datasets entry "
        f"(available: {list(config.get('real_datasets', {}).keys())})"
    )

active_data_types = [
    dt for dt in config.get("data_types", ["synth_realistic"])
    if dt not in ("synth_custom", "real")
]

if "real" in config.get("data_types", []) and not config.get("real_datasets"):
    logger.warning('data_types includes "real" but real_datasets is empty')

# Whether the run has simulations (scored against the known answer) and/or real samples
# (tools compared with each other); rule all asks only for the summaries that apply.
_has_ground_truth_data = bool(active_data_types) or "synth_custom" in config.get("data_types", [])
_has_no_ground_truth_data = "real" in config.get("data_types", [])

# One HTML report per data type with known answers: file name and title for each
REPORTS = {
    "synth_realistic": ("synth_realistic", "Realistic simulations"),
    "synth_custom": ("synth_custom", "Custom simulations"),
}
report_names = [REPORTS[dt][0] for dt in REPORTS if dt in config.get("data_types", ["synth_realistic"])]


def report_file(name):
    """One of the HTML reports, e.g. report_synth_realistic.html."""
    return os.path.join(config["output_dir"], f"report_{name}.html")


def summary_file(table):
    """One of the tables in summary/ that collect every run, e.g. evaluation_metrics.tsv."""
    return os.path.join(config["output_dir"], "summary", f"{table}.tsv")


def summary_log(name):
    """Log of a step that builds a summary table or a report."""
    base = os.path.join(config["log_dir"], "summary", name)
    return {"stderr": f"{base}.stderr.log", "stdout": f"{base}.stdout.log"}
report_data_type = {name: dt for dt, (name, _) in REPORTS.items()}


# ==============================================================================
# CROSS-TOOL AGREEMENT (real samples: no known answer, so the tools are compared with each other)
# ==============================================================================

def get_agreement_tools():
    """Tools usable for cross-tool agreement. SigID runs fine here — its true_weights
    input is optional (diagnostic-only), so it fits no-ground-truth data like every
    other tool."""
    return get_selected_tools()


def get_agreement_files(suffix="tool_agreement"):
    """Per-cell TSVs (tool_agreement or signature_calls — both outputs of the same
    evaluate_tool_agreement rule) for every active no-ground-truth dataset."""
    if "real" not in config.get("data_types", []):
        return []
    return expand(
        run_file("agreement", f"{suffix}.tsv"),
        data_type="real",
        cancer_type=get_all_real_cancer_types(),
        num_muts=REAL_NUM_MUTS,
        relevance=relevance_endings(),
    )


# ==============================================================================
# CONFIGURATION INHERITANCE FUNCTIONS
# ==============================================================================

def get_param(section, param, default_key=None):
    """Get parameter with inheritance: simulation[section] → config top-level → special mappings"""
    default_key = default_key or param

    sim_section = config.get("simulation", {}).get(section)
    if sim_section and isinstance(sim_section, dict):
        val = sim_section.get(param)
        if val is not None:
            return val

    val = config.get(param)
    if val is not None:
        return val

    special_mappings = {
        "input_sigs": config.get("cosmic_signatures", {}).get("base_file"),
        "rng_seed": config.get("random_seed"),
    }

    if param in special_mappings and special_mappings[param] is not None:
        return special_mappings[param]

    if default_key in special_mappings and special_mappings[default_key] is not None:
        return special_mappings[default_key]

    return None


def get_tool_config(tool_name):
    """Settings of one tool: the config's defaults, then its tool_overrides entry."""
    tool_config = config.get("defaults", {}).copy()

    for key, value in config.get("tool_overrides", {}).get(tool_name, {}).items():
        if isinstance(value, dict) and isinstance(tool_config.get(key), dict):
            tool_config[key].update(value)
        else:
            tool_config[key] = value

    return tool_config


# ==============================================================================
# EXTRA TOOL ARGUMENTS (tool_overrides: <tool>: extra: in the config)
# ==============================================================================

def extra_values(tool):
    """Extra arguments for one tool, from tool_overrides.<tool>.extra in the config.

    They are passed to the tool as they are and can replace what the rule sets itself
    (input, output, cores, ...), so use them with care.
    """
    return get_tool_config(tool).get("extra") or {}


def extra_params(tool):
    """For a rule's params. Adds nothing when no values are set, so that the job's params stay
    the same and Snakemake has no reason to rerun it."""
    values = extra_values(tool)
    return {"extra": values} if values else {}


def _shell_text(words):
    """Words to add to a shell command; empty if there are none, which leaves the command's text
    the same, so Snakemake has no reason to rerun it."""
    return "".join(" " + shlex.quote(str(w)).replace("{", "{{").replace("}", "}}") for w in words)


def extra_flags(tool):
    """Extra arguments as command-line flags: `true` gives the flag alone, `false` leaves it out."""
    words = []
    for flag, value in extra_values(tool).items():
        if value is True:
            words.append(flag)
        elif value is not False:
            words += [flag, value]
    return _shell_text(words)


def extra_name_values(tool):
    """Extra arguments as name=value words, for R scripts run from a shell command."""
    return _shell_text(f"{name}={value}" for name, value in extra_values(tool).items())


# ==============================================================================
# INPUT/OUTPUT PATH RESOLUTION
# ==============================================================================

# Reference signatures the tools fit with (base_file minus the hidden signatures), one file
# per format, e.g. internal/reference_signatures/tab_with_index.txt
REFERENCE_DIR = os.path.join("internal", "reference_signatures")


# Reference signature file each input format reads (every tool reads the same counts)
REFERENCE_FORMAT = {
    "tab_with_index": "tab_with_index",
    "csv_no_index": "tab_with_index",
    "csv_with_index": "csv_with_index",
    "csv_quoted": "csv_quoted",
    "signature_tools": "tab_with_index",
}

# Tools the workflow can run (name as in its result file -> input format).
# Filled by tool_result() when tools.smk is read.
REGISTERED_TOOLS = {}


# ==============================================================================
# RESULTS LAYOUT
#   {output_dir}/report_<data type>.html          one report per data type
#   {output_dir}/summary/<table>.tsv              one table each, all runs together
#   {output_dir}/runs/<data type>/<tissue or dataset>/<N>_mutations/
#       samples/counts_rep<r>.dat                 mutation counts the tools fit
#       samples/known_answer_rep<r>.dat           true signature mix of each sample
#       fits/<tool>_rep<r>.dat                    what each tool found
#       fits/extra/<tool>_<part>_rep<r>.dat       other files some tools write
#       scores/<tool>_rep<r>.tsv                  scores of each sample
#       benchmarks/<tool>_rep<r>.tsv              run time and memory
#       agreement/<table>.tsv                     tools compared with each other (real)
#   {output_dir}/internal/                        files the workflow needs, rarely opened
# Logs use the same runs/... folders under log_dir.
# ==============================================================================

# Folder ending for each relevance mode: nothing for the normal mode
RELEVANCE_ENDING = {"standard": "", "relevant": "_relevant"}

# Stands in for the mutation count of real samples, which keep their own counts
REAL_NUM_MUTS = "all"


def relevance_endings():
    """Folder endings of the relevance modes this run uses."""
    return [RELEVANCE_ENDING[m] for m in config.get("relevance_modes", ["standard"])]


def relevance_mode(ending):
    """Relevance mode ("standard" or "relevant") of a folder ending."""
    return "relevant" if ending == RELEVANCE_ENDING["relevant"] else "standard"


# Folder of one run: data type / tissue or dataset / mutations per sample (+ relevance)
RUN_DIR = os.path.join("runs", "{data_type}", "{cancer_type}", "{num_muts}_mutations{relevance}")
# The samples are the same in every relevance mode, so they are kept in the normal folder
SAMPLES_DIR = os.path.join("runs", "{data_type}", "{cancer_type}", "{num_muts}_mutations", "samples")


def run_file(*parts):
    """A file in one run's folder, e.g. run_file("scores", "{tool}_rep{rep}.tsv")."""
    return os.path.join(config["output_dir"], RUN_DIR, *parts)


def run_log(name):
    """stderr and stdout logs of one step of one run."""
    base = os.path.join(config["log_dir"], RUN_DIR, name)
    return {"stderr": f"{base}.stderr.log", "stdout": f"{base}.stdout.log"}


def samples_file(name, data_type="{data_type}"):
    """Counts or known answer of one run, e.g. samples_file("counts_rep{rep}.dat")."""
    return os.path.join(config["output_dir"], SAMPLES_DIR.replace("{data_type}", data_type), name)


def samples_log(data_type, stream):
    """Log of the step that makes one run's samples."""
    return os.path.join(config["log_dir"], SAMPLES_DIR.replace("{data_type}", data_type),
                        f"rep{{rep}}.{stream}.log")


def tool_result(tool, name="contribution", input_format=None):
    """Result file of one tool for one run, e.g. fits/SPA_rep0.dat. Other files a tool
    writes (name other than "contribution") go to fits/extra/.

    Tool rules also pass the input format the tool reads, which registers the tool.
    """
    if input_format:
        REGISTERED_TOOLS[tool] = input_format
    if name == "contribution":
        return run_file("fits", f"{tool}_rep{{rep}}.dat")
    return run_file("fits", "extra", f"{tool}_{name}_rep{{rep}}.dat")


# Added to the shell command of command-line tools: if the tool fails, write
# "# TOOL FAILED: <last line of its error log>" to its result file and print it, instead of
# stopping the workflow, so the HTML report is still built.
ON_TOOL_FAILURE = (
    ' || {{ reason=$(tail -n 1 {log.stderr:q} | tr -s "[:space:]" " "); '
    'echo "# TOOL FAILED: ${{reason:-see log}}" > {output.result:q}; '
    'echo "TOOL FAILED ({rule}): ${{reason:-see log}}" >&2; }}'
)


def tool_logs(tool):
    """stderr and stdout log files of one tool for one run."""
    return run_log(f"{tool}_rep{{rep}}")


def tool_benchmark(tool):
    """Run time and memory file of one tool for one run."""
    return run_file("benchmarks", f"{tool}_rep{{rep}}.tsv")


def get_input_data(tool_key, wildcards):
    """Mutation counts a tool fits: the run's samples/counts_rep<r>.dat. Every tool reads
    the same file (real samples have rep 0)."""
    return samples_file(f"counts_rep{wildcards.rep}.dat").format(
        data_type=wildcards.data_type, cancer_type=wildcards.cancer_type, num_muts=wildcards.num_muts)


# ==============================================================================
# COSMIC SIGNATURE FILE MANAGEMENT
# ==============================================================================

def _get_format_key(tool_name):
    """Input format of a tool, e.g. tab_with_index."""
    for tool, fmt in REGISTERED_TOOLS.items():
        if tool.lower() == tool_name.lower():
            return fmt
    return None


def signature_catalogue(data_type=None, dataset=None):
    """Signature catalogue to use: a custom dataset's signature_file, otherwise base_file."""
    base_file = config["cosmic_signatures"]["base_file"]
    if data_type == "synth_custom" and dataset:
        return config["custom_datasets"][dataset].get("signature_file", base_file)
    return base_file


def get_cosmic_signatures(tool_name, relevance, data_type=None, num_muts=None, dataset=None):
    """Reference signature file a tool fits with, in internal/reference_signatures/:
    <format>.txt, or custom_<dataset>/<format>.txt for a custom dataset, which may use
    its own catalogue. relevance is the mode or its folder ending."""
    catalogue = REFERENCE_FORMAT.get(_get_format_key(tool_name), "tab_with_index") + ".txt"
    reference_dir = os.path.join(config["output_dir"], REFERENCE_DIR)

    if data_type == "synth_custom" and dataset:
        return os.path.join(reference_dir, f"custom_{dataset}", catalogue)
    if relevance in ("relevant", RELEVANCE_ENDING["relevant"]):
        return os.path.join(reference_dir, f"relevant_{data_type}_{num_muts}", catalogue)
    return os.path.join(reference_dir, catalogue)


# ==============================================================================
# TOOL SELECTION AND FILTERING
# ==============================================================================

def get_selected_tools():
    """Get list of tools to run: 'all' or specific list"""
    tools_config = config.get("tools_to_run", "all")
    selection = (tools_config if isinstance(tools_config, (str, list))
                 else tools_config.get("selection", "all") if isinstance(tools_config, dict)
                 else "all")

    all_tools = list(REGISTERED_TOOLS)
    all_tools_lower = {t.lower() for t in all_tools}

    if selection == "all":
        return all_tools
    if isinstance(selection, list):
        # Names are compared ignoring case, but kept as written in tools_to_run
        unknown = [t for t in selection if t.lower() not in all_tools_lower]
        if unknown and not _warned_unknown_tools:
            print(f"WARNING: tools_to_run entries that are not known tools are skipped: "
                  f"{', '.join(unknown)}. Known tools: {', '.join(all_tools)}", file=sys.stderr)
            _warned_unknown_tools.append(True)
        return [t for t in selection if t.lower() in all_tools_lower]
    return all_tools


_warned_unknown_tools = []  # so the warning above is printed once per run


# The signatures left out of the reference, as one space-separated text for shell commands
OUT_OF_REF_SIGS_STR = ' '.join(config["out_of_reference_sigs"])


def get_eval_files():
    """Return all per-sample evaluation TSV paths for aggregation.

    real has no ground truth and is excluded from active_data_types entirely — it's
    evaluated separately via the tool-agreement pipeline, not this function.
    """
    _n_reps = config["num_realizations"]

    files = []
    for dt in active_data_types:
        _common = dict(
            data_type=[dt],
            cancer_type=config["cancer_types"],
            num_muts=config["num_muts"],
            relevance=relevance_endings(),
            rep=range(_n_reps),
        )
        files += expand(
            run_file("scores", "{tool}_rep{rep}.tsv"),
            tool=get_selected_tools(),
            **_common
        )
    custom_datasets = config.get("custom_datasets", {})
    if custom_datasets and "synth_custom" in config.get("data_types", []):
        _common_custom = dict(
            data_type=["synth_custom"],
            cancer_type=list(custom_datasets.keys()),
            num_muts=config["num_muts"],
            relevance=relevance_endings(),
            rep=range(_n_reps),
        )
        files += expand(
            run_file("scores", "{tool}_rep{rep}.tsv"),
            tool=get_selected_tools(),
            **_common_custom
        )

    return files


# ==============================================================================
# RELEVANT SIGNATURE HELPERS
# ==============================================================================

def get_relevant_signature_file(wildcards):
    """Return path to the relevant signatures file (auto-detected or manually created)."""
    rel_sigs = config.get('relevant_signatures')
    if rel_sigs == "auto":
        return os.path.join(
            config['output_dir'], "internal", "relevant_signatures",
            f"{wildcards.data_type}_{wildcards.num_muts}_active_signatures.txt"
        )
    elif isinstance(rel_sigs, list) and rel_sigs:
        manual_file = os.path.join(
            config['output_dir'], "internal", "relevant_signatures",
            "manual_relevant_signatures.txt"
        )
        os.makedirs(os.path.dirname(manual_file), exist_ok=True)
        with open(manual_file, 'w') as f:
            f.write("# Manually specified relevant signatures from config\n")
            for sig in rel_sigs:
                f.write(f"{sig}\n")
        return manual_file
    return None


def get_relevant_signatures_for_cosmic(wildcards):
    """Return comma-separated relevant signatures string for COSMIC file generation."""
    rel_sigs = config.get('relevant_signatures')
    if rel_sigs == "auto":
        sig_file = get_relevant_signature_file(wildcards)
        if not os.path.exists(sig_file):
            raise FileNotFoundError(
                f"AUTO mode enabled but active signatures file not found: {sig_file}\n"
                "Run standard mode first to generate active signatures."
            )
        sigs = [l.strip() for l in open(sig_file) if l.strip() and not l.startswith('#')]
        if not sigs:
            raise ValueError(f"No signatures found in {sig_file}")
        return ",".join(sigs)
    elif isinstance(rel_sigs, list) and rel_sigs:
        return ",".join(rel_sigs)
    raise ValueError(
        f"relevant_signatures must be 'auto' or a non-empty list. Got: {rel_sigs}"
    )
