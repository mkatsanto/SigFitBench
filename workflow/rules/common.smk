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
import sys


# ==============================================================================
# ACTIVE DATA TYPES
# Computed once here so every included rule file can reference active_data_types.
# realsim (like custom/external) is handled by its own dedicated target block in
# rule all — it has no ground truth and doesn't use the shared cancer_types/num_muts
# lists, so it's excluded from the generic per-data-type expansion here.
# ==============================================================================

def _realsim_configured():
    real_samples_cfg = config.get("simulation", {}).get("real_samples", {}) or {}
    return bool(real_samples_cfg.get("input_real_catalog"))

def get_realsim_cancer_types():
    """Cancer types to run for realsim, from simulation.real_samples.cancer_types."""
    return config.get("simulation", {}).get("real_samples", {}).get("cancer_types", [])


def _peek_catalog_header(catalog_path):
    """Return the first non-comment, non-blank line of a catalog file (its header row),
    or None if the file can't be read. Used only to detect column-naming convention —
    doesn't load the full catalog (that's fit_external_catalog.py's job at run time)."""
    if not catalog_path or not os.path.exists(catalog_path):
        return None
    with open(catalog_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if line and not line.startswith("#"):
                return line
    return None


_EXTERNAL_CANCER_TYPES_CACHE = {}


def get_external_cancer_types(ds_name):
    """Cancer types (tissues) for one external_datasets entry.

    Auto-detected from the catalog's column headers using the "Tissue::SampleID"
    convention (same one realsim's PCAWG file uses) — peeks at just the header row,
    not the whole file. If any column follows that convention, returns the sorted
    list of distinct tissue prefixes found (per-tissue splitting, like realsim). If
    none do, returns [ds_name]: the whole dataset is fit as one flat bucket named
    after the dataset key itself (unsplit — the original external behavior, still
    correct for catalogs that don't use per-tissue column naming).
    """
    if ds_name in _EXTERNAL_CANCER_TYPES_CACHE:
        return _EXTERNAL_CANCER_TYPES_CACHE[ds_name]

    catalog_path = config.get("external_datasets", {}).get(ds_name, {}).get("catalog")
    header = _peek_catalog_header(catalog_path)
    result = [ds_name]
    if header:
        delim = "\t" if "\t" in header else ("," if "," in header else None)
        cols = header.split(delim)[1:] if delim else header.split()[1:]  # skip index column
        tissues = sorted({c.strip().strip('"') .split("::")[0] for c in cols if "::" in c})
        if tissues:
            result = tissues

    _EXTERNAL_CANCER_TYPES_CACHE[ds_name] = result
    return result


def get_all_external_cancer_types():
    """Flat list of "cancer_type" wildcard values across every external_datasets
    entry — each dataset contributes either its own key (unsplit) or one entry per
    auto-detected tissue (get_external_cancer_types). Used as the cancer_type
    dimension in rule all's expand()-based external_tool_results."""
    return [
        ct for ds_name in config.get("external_datasets", {})
        for ct in get_external_cancer_types(ds_name)
    ]


def resolve_external_dataset(cancer_type):
    """Map a "cancer_type" wildcard back to (dataset_name, tissue_prefix_or_None).

    tissue_prefix is None when the matched dataset isn't using per-tissue splitting
    (cancer_type == the dataset key itself); otherwise it's the tissue to filter the
    catalog to (passed to fit_external_catalog.py's --tissue-prefix).
    """
    for ds_name in config.get("external_datasets", {}):
        tissues = get_external_cancer_types(ds_name)
        if tissues == [ds_name] and cancer_type == ds_name:
            return ds_name, None
        if cancer_type in tissues:
            return ds_name, cancer_type
    raise ValueError(
        f"cancer_type {cancer_type!r} does not match any external_datasets entry "
        f"(available: {list(config.get('external_datasets', {}).keys())})"
    )

active_data_types = [
    dt for dt in config.get("data_types", ["synthetic"])
    if dt not in ("custom", "realsim", "external")
]

if "realsim" in config.get("data_types", []) and not _realsim_configured():
    logger.warning("realsim requested but simulation.real_samples.input_real_catalog is not set")

# Whether this run has any ground-truth data (synthetic/custom) vs. any no-ground-truth
# data (realsim/external) — rule all uses these to request the right summary/report
# target(s) only, instead of unconditionally requesting both (which would crash
# compute_summary_metrics for a no-ground-truth-only run: get_eval_files() returns []
# and its --input-files argument requires at least one file).
_has_ground_truth_data = bool(active_data_types) or "custom" in config.get("data_types", [])
_has_no_ground_truth_data = (
    "realsim" in config.get("data_types", []) or "external" in config.get("data_types", [])
)

# One HTML report per data type with known answers: file name and title for each
REPORTS = {
    "synthetic": ("realistic_simulations", "Realistic simulations"),
    "custom": ("synthetic", "Synthetic"),
}
report_names = [REPORTS[dt][0] for dt in REPORTS if dt in config.get("data_types", ["synthetic"])]
report_data_type = {name: dt for dt, (name, _) in REPORTS.items()}


# ==============================================================================
# CROSS-TOOL AGREEMENT (no-ground-truth data types: realsim, external)
# ==============================================================================

def get_agreement_tools():
    """Tools usable for cross-tool agreement. SigID runs fine here — its true_weights
    input is optional (diagnostic-only), so it fits no-ground-truth data like every
    other tool."""
    return get_selected_tools()


def get_agreement_files(suffix="tool_agreement"):
    """Per-cell TSVs (tool_agreement or signature_calls — both outputs of the same
    evaluate_tool_agreement rule) for every active no-ground-truth dataset."""
    eval_dir = os.path.join(config['output_dir'], "evaluation")
    relevance_modes = config.get("relevance_modes", ["standard"])
    files = []

    if "realsim" in config.get("data_types", []):
        files += expand(
            os.path.join(eval_dir, f"realsim_{{cancer_type}}_none_{{relevance}}_{suffix}.tsv"),
            cancer_type=get_realsim_cancer_types(),
            relevance=relevance_modes,
        )

    if "external" in config.get("data_types", []):
        for ds_name in config.get("external_datasets", {}):
            for cancer_type in get_external_cancer_types(ds_name):
                files += expand(
                    os.path.join(eval_dir, f"external_{cancer_type}_none_{{relevance}}_{suffix}.tsv"),
                    relevance=relevance_modes,
                )

    return files


# ==============================================================================
# CONFIGURATION INHERITANCE FUNCTIONS
# ==============================================================================

def _get_gen_name(tool_key):
    """Name of the counts file a tool reads, e.g. deconstructSigs for data_for_deconstructSigs.dat."""
    for format_data in _format_groups().values():
        if tool_key.lower() in [t.lower() for t in format_data.get("tools", [])]:
            base = os.path.basename(format_data.get("input_data_base", ""))
            m = re.match(r"data_for_(.+)\.dat$", base)
            return m.group(1) if m else tool_key
    return tool_key


def get_required_output_formats():
    """Get list of output formats needed for selected tools"""
    return sorted({_get_gen_name(t) for t in get_selected_tools()})


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
    """Get tool config with inheritance: defaults → format_group → tool_overrides"""
    tool_config = config.get("defaults", {}).copy()

    for format_data in _format_groups().values():
        if tool_name.lower() in [t.lower() for t in format_data.get("tools", [])]:
            cs = format_data.get("cosmic_signatures", {}) or {}
            tool_config.update({
                "cosmic_signatures": {
                    "standard": cs.get("standard", ""),
                    "relevant": cs.get("relevant", "")
                },
                "input_data": format_data.get("input_data_base", ""),
                **{k: v for k, v in format_data.items()
                   if k not in ["tools", "cosmic_signatures", "input_data_base"]}
            })
            break

    for key, value in config.get("tool_overrides", {}).get(tool_name, {}).items():
        if isinstance(value, dict) and isinstance(tool_config.get(key), dict):
            tool_config[key].update(value)
        else:
            tool_config[key] = value

    return tool_config


# ==============================================================================
# INPUT/OUTPUT PATH RESOLUTION
# ==============================================================================

# Counts file and signature catalogue each input format uses
FORMAT_FILES = {
    "tab_with_index": ("data/data_for_deconstructSigs.dat",
                       "{output_dir}/cosmic_signatures/COSMIC_v3_SBS_GRCh38_tab_with_index.txt"),
    "csv_no_index": ("data/data_for_MutationalPatterns.dat",
                     "{output_dir}/cosmic_signatures/COSMIC_v3_SBS_GRCh38_tab_with_index.txt"),
    "csv_with_index": ("data/data_for_sigfit.dat",
                       "{output_dir}/cosmic_signatures/COSMIC_v3_SBS_GRCh38_csv_with_index.txt"),
    "csv_quoted": ("data/data_for_sigLASSO.dat",
                   "{output_dir}/cosmic_signatures/COSMIC_v3_SBS_GRCh38_csv_quoted.txt"),
    "signature_tools": ("data/data_for_signature_tools.dat",
                        "{output_dir}/cosmic_signatures/COSMIC_v3_SBS_GRCh38_tab_with_index.txt"),
}

# Tools the workflow can run (name as in its result file -> input format).
# Filled by tool_result() when tools.smk is read.
REGISTERED_TOOLS = {}


def _format_groups():
    """Tools grouped by input format: the config's tool_formats if it has one, otherwise
    built from the tool rules."""
    if config.get("tool_formats"):
        return config["tool_formats"]
    groups = {}
    for tool, fmt in REGISTERED_TOOLS.items():
        data_file, catalogue = FORMAT_FILES[fmt]
        group = groups.setdefault(fmt, {"tools": [], "input_data_base": data_file,
                                        "cosmic_signatures": {"standard": catalogue}})
        group["tools"].append(tool)
    return groups


# Folder of one fitting run: data type / sample set / signature set / mutations per sample /
# replicate. Every tool result, log and benchmark path is built from this one template.
RUN_DIR = os.path.join("{data_type}", "{cancer_type}", "{relevance}", "{num_muts}", "rep_{rep}")


def tool_result(tool, name="contribution", input_format=None):
    """Result file of one tool for one run, e.g. SPA_contribution.dat.

    Tool rules also pass the input format the tool reads, which registers the tool.
    """
    if input_format:
        REGISTERED_TOOLS[tool] = input_format
    return os.path.join(config["output_dir"], RUN_DIR, f"{tool}_{name}.dat")


# Added to the shell command of command-line tools: if the tool fails, write
# "# TOOL FAILED: <last line of its error log>" to its result file and print it, instead of
# stopping the workflow, so the HTML report is still built.
ON_TOOL_FAILURE = (
    ' || {{ reason=$(tail -n 1 {log.stderr} | tr -s "[:space:]" " "); '
    'echo "# TOOL FAILED: ${{reason:-see log}}" > {output.result}; '
    'echo "TOOL FAILED ({rule}): ${{reason:-see log}}" >&2; }}'
)


def tool_logs(tool):
    """stderr and stdout log files of one tool for one run."""
    base = os.path.join(config["log_dir"], RUN_DIR, tool.lower())
    return {"stderr": f"{base}_stderr.log", "stdout": f"{base}_stdout.log"}


def tool_benchmark(tool):
    """Run time and memory file of one tool for one run."""
    return os.path.join(config["output_dir"], "benchmarks", RUN_DIR, f"{tool}_benchmark.tsv")


def get_input_data(tool_key, wildcards):
    """Return input data path based on data_type: synthetic/realsim/samples/external.

    Realsim uses the same {cancer_type}-w_0-{num_muts} naming as synthetic (rep fixed at 0),
    so downstream rules need no special-casing.
    """
    gen_name = _get_gen_name(tool_key)
    dt = getattr(wildcards, 'data_type', None)

    if dt == 'external':
        return os.path.join(config["output_dir"], "generated_data", "external",
                            f"counts_{wildcards.cancer_type}-w_0-{wildcards.num_muts}_{gen_name}.dat")

    if dt == 'custom':
        return os.path.join(config["output_dir"], "generated_data", "custom",
                            f"counts_{wildcards.cancer_type}-w_{wildcards.rep}-{wildcards.num_muts}_{gen_name}.dat")

    if dt in ('synthetic', 'realsim'):
        # realsim always uses rep 0; synthetic uses the actual rep wildcard
        rep = 0 if dt == 'realsim' else wildcards.rep
        return os.path.join(config["output_dir"], "generated_data", dt,
                            f"counts_{wildcards.cancer_type}-w_{rep}-{wildcards.num_muts}_{gen_name}.dat")

    if dt == 'samples' and hasattr(wildcards, 'identifier'):
        sample_dir = config.get("sample_input_dir", "input/samples")
        for path in [f"{wildcards.identifier}_mutations.dat", f"{wildcards.identifier}_{gen_name}.dat"]:
            full_path = os.path.join(sample_dir, path)
            if os.path.exists(full_path):
                return full_path
        if "input_catalog" in config.get("samples", {}):
            return config["samples"]["input_catalog"]

    override = config.get("tool_overrides", {}).get(tool_key, {}).get("input_data")
    if override:
        return override

    for format_info in _format_groups().values():
        if tool_key.lower() in [t.lower() for t in format_info.get("tools", [])]:
            return format_info.get("input_data_base") or []
    return []


# ==============================================================================
# COSMIC SIGNATURE FILE MANAGEMENT
# ==============================================================================

def _get_format_key(tool_name):
    """Input format of a tool, e.g. tab_with_index."""
    for fmt_key, format_data in _format_groups().items():
        if tool_name.lower() in [t.lower() for t in format_data.get("tools", [])]:
            return fmt_key
    return None


def signature_catalogue(data_type=None, dataset=None):
    """Signature catalogue to use: a custom dataset's signature_file, otherwise base_file."""
    base_file = config["cosmic_signatures"]["base_file"]
    if data_type == "custom" and dataset:
        return config["custom_datasets"][dataset].get("signature_file", base_file)
    return base_file


def get_cosmic_signatures(tool_name, relevance, data_type=None, num_muts=None, dataset=None):
    """Signature file a tool fits with. Custom datasets get their own file."""
    cs = get_tool_config(tool_name).get("cosmic_signatures", {}) or {}
    base_path = cs.get("standard") or ""

    if "{output_dir}" in base_path:
        base_path = base_path.format(output_dir=config["output_dir"])
    if base_path.startswith("results/"):
        base_path = base_path.replace("results/", config["output_dir"] + "/", 1)

    if data_type == "custom" and dataset:
        # Same file name as the tool's standard file, with _CUSTOM_{dataset} added
        if base_path.startswith(config["output_dir"]):
            return base_path.replace(".txt", f"_CUSTOM_{dataset}.txt")
        fmt_key = _get_format_key(tool_name)
        stem = os.path.splitext(os.path.basename(base_path))[0]
        return os.path.join(
            config["output_dir"], "cosmic_signatures",
            f"{stem}_{fmt_key}_CUSTOM_{dataset}.txt"
        )

    if relevance == "relevant":
        if base_path.startswith(config["output_dir"]):
            # Standard file is already in output dir: insert _RELEVANT before .txt
            base_path = base_path.replace(".txt", "_RELEVANT.txt")
        else:
            # Standard file is a raw input (e.g. csv_no_index uses input/ directly):
            # build the format-specific relevant path in output dir
            fmt_key = _get_format_key(tool_name)
            stem = os.path.splitext(os.path.basename(base_path))[0]
            base_path = os.path.join(
                config["output_dir"], "cosmic_signatures",
                f"{stem}_{fmt_key}_RELEVANT.txt"
            )

        if data_type and num_muts:
            base_path = base_path.replace("_RELEVANT.txt", f"_RELEVANT_{data_type}_{num_muts}.txt")

    return base_path


# ==============================================================================
# TOOL SELECTION AND FILTERING
# ==============================================================================

def get_selected_tools():
    """Get list of tools to run: 'all' or specific list"""
    tools_config = config.get("tools_to_run", "all")
    selection = (tools_config if isinstance(tools_config, (str, list))
                 else tools_config.get("selection", "all") if isinstance(tools_config, dict)
                 else "all")

    all_tools = [t for fmt in _format_groups().values()
                 for t in fmt.get("tools", [])]
    all_tools_lower = {t.lower() for t in all_tools}

    if selection == "all":
        return all_tools
    if isinstance(selection, list):
        # Case-insensitive membership check so "SPA"/"MuSiCal" in tools_to_run
        # match "spa"/"musical" in the known tools, while preserving the case from
        # tools_to_run (which must match the hardcoded output filenames in tools.smk).
        unknown = [t for t in selection if t.lower() not in all_tools_lower]
        if unknown and not _warned_unknown_tools:
            print(f"WARNING: tools_to_run entries that are not known tools are skipped: "
                  f"{', '.join(unknown)}. Known tools: {', '.join(all_tools)}", file=sys.stderr)
            _warned_unknown_tools.append(True)
        return [t for t in selection if t.lower() in all_tools_lower]
    return all_tools


_warned_unknown_tools = []  # so the warning above is printed once per run


# Pre-joined string for passing out-of-reference signatures to shell commands
OUT_OF_REF_SIGS_STR = ' '.join(config["out_of_reference_sigs"])


def get_eval_files():
    """Return all per-sample evaluation TSV paths for aggregation.

    realsim and external have no ground truth and are excluded from active_data_types
    entirely — they're evaluated separately via the tool-agreement pipeline, not this
    function.
    """
    _eval_dir = os.path.join(config['output_dir'], "evaluation")
    _n_reps = config["num_realizations"]

    files = []
    for dt in active_data_types:
        _common = dict(
            data_type=[dt],
            cancer_type=config["cancer_types"],
            num_muts=config["num_muts"],
            relevance=config["relevance_modes"],
            rep=range(_n_reps),
        )
        files += expand(
            os.path.join(_eval_dir, "{data_type}_{cancer_type}_{num_muts}_{relevance}_{tool}_w{rep}_per_sample.tsv"),
            tool=get_selected_tools(),
            **_common
        )
    custom_datasets = config.get("custom_datasets", {})
    if custom_datasets and "custom" in config.get("data_types", []):
        _common_custom = dict(
            data_type=["custom"],
            cancer_type=list(custom_datasets.keys()),
            num_muts=config["num_muts"],
            relevance=config["relevance_modes"],
            rep=range(_n_reps),
        )
        files += expand(
            os.path.join(_eval_dir, "{data_type}_{cancer_type}_{num_muts}_{relevance}_{tool}_w{rep}_per_sample.tsv"),
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
            config['output_dir'], "relevant_signatures",
            f"{wildcards.data_type}_{wildcards.num_muts}_active_signatures.txt"
        )
    elif isinstance(rel_sigs, list) and rel_sigs:
        manual_file = os.path.join(
            config['output_dir'], "relevant_signatures",
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
