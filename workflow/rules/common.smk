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

active_data_types = [
    dt for dt in config.get("data_types", ["synthetic"])
    if dt not in ("custom", "realsim")
]

if "realsim" in config.get("data_types", []) and not _realsim_configured():
    logger.warning("realsim requested but simulation.real_samples.input_real_catalog is not set")


# ==============================================================================
# CROSS-TOOL AGREEMENT (no-ground-truth data types: realsim, external)
# ==============================================================================

def get_agreement_tools():
    """Tools usable for cross-tool agreement: excludes penalty-sweep tools (no single
    contribution.dat) and SigID (its auto-penalty search requires a true_weights input
    we don't have for no-ground-truth data types)."""
    excluded = {"SigID"}
    return [t for t in get_selected_tools() if not has_penalty_sweep(t) and t not in excluded]


def get_agreement_files():
    """Per-cell tool-agreement TSVs for every active no-ground-truth dataset."""
    eval_dir = os.path.join(config['output_dir'], "evaluation")
    relevance_modes = config.get("relevance_modes", ["standard"])
    files = []

    if "realsim" in config.get("data_types", []):
        files += expand(
            os.path.join(eval_dir, "realsim_{cancer_type}_none_{relevance}_tool_agreement.tsv"),
            cancer_type=get_realsim_cancer_types(),
            relevance=relevance_modes,
        )

    if "external" in config.get("data_types", []):
        for ds_name, ds in config.get("external_datasets", {}).items():
            files += expand(
                os.path.join(eval_dir, f"external_{ds_name}_{ds['num_muts_label']}_{{relevance}}_tool_agreement.tsv"),
                relevance=relevance_modes,
            )

    return files


# ==============================================================================
# CONFIGURATION INHERITANCE FUNCTIONS
# ==============================================================================

def _get_gen_name(tool_key):
    """Derive input data format name from tool_formats config."""
    for format_data in config.get("tool_formats", {}).values():
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

    for format_data in config.get("tool_formats", {}).values():
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

    for format_info in config.get("tool_formats", {}).values():
        if tool_key in format_info.get("tools", []):
            return format_info.get("input_data_base") or []
    return []


# ==============================================================================
# COSMIC SIGNATURE FILE MANAGEMENT
# ==============================================================================

def _get_format_key(tool_name):
    """Return the tool_formats group key for a given tool."""
    for fmt_key, format_data in config.get("tool_formats", {}).items():
        if tool_name.lower() in [t.lower() for t in format_data.get("tools", [])]:
            return fmt_key
    return None


def get_cosmic_signatures(tool_name, relevance, data_type=None, num_muts=None):
    """Get COSMIC signatures path for tool and relevance (standard/relevant)."""
    cs = get_tool_config(tool_name).get("cosmic_signatures", {}) or {}
    base_path = cs.get("standard") or ""

    if "{output_dir}" in base_path:
        base_path = base_path.format(output_dir=config["output_dir"])
    if base_path.startswith("results/"):
        base_path = base_path.replace("results/", config["output_dir"] + "/", 1)

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

    all_tools = [t for fmt in config.get("tool_formats", {}).values()
                 for t in fmt.get("tools", [])]
    all_tools_lower = {t.lower() for t in all_tools}

    if selection == "all":
        return all_tools
    if isinstance(selection, list):
        # Case-insensitive membership check so "SPA"/"MuSiCal" in tools_to_run
        # match "spa"/"musical" in tool_formats, while preserving the case from
        # tools_to_run (which must match the hardcoded output filenames in tools.smk).
        return [t for t in selection if t.lower() in all_tools_lower]
    return all_tools


def has_penalty_sweep(tool_name):
    """Check if a tool is in the global penalty_tools list"""
    return tool_name in config.get("penalty_tools", [])


def get_penalty_values(tool_name):
    """Get the global penalty sweep values for a penalty tool"""
    if not has_penalty_sweep(tool_name):
        return []
    return [str(float(p)) for p in config["penalty_values"]]


# Pre-joined string for passing out-of-reference signatures to shell commands
OUT_OF_REF_SIGS_STR = ' '.join(config["out_of_reference_sigs"])


def get_eval_files():
    """Return all per-sample evaluation TSV paths (regular + penalty sweep) for aggregation.

    realsim has no ground truth and is excluded from active_data_types entirely —
    it's evaluated separately via the tool-agreement pipeline, not this function.
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
            tool=[t for t in get_selected_tools() if not has_penalty_sweep(t)],
            **_common
        )
        files += expand(
            os.path.join(_eval_dir, "{data_type}_{cancer_type}_{num_muts}_{relevance}_{tool}_penalty_{penalty}_w{rep}_per_sample.tsv"),
            tool=[t for t in get_selected_tools() if has_penalty_sweep(t)],
            penalty=[str(float(p)) for p in config.get("penalty_values", [])],
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
            tool=[t for t in get_selected_tools() if not has_penalty_sweep(t)],
            **_common_custom
        )
        files += expand(
            os.path.join(_eval_dir, "{data_type}_{cancer_type}_{num_muts}_{relevance}_{tool}_penalty_{penalty}_w{rep}_per_sample.tsv"),
            tool=[t for t in get_selected_tools() if has_penalty_sweep(t)],
            penalty=[str(float(p)) for p in config.get("penalty_values", [])],
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
