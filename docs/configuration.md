# Configuration reference

Every run needs a config file passed with `--configfile`. The configs in
[`config/scenarios/`](https://github.com/mkatsanto/SigFitBench/tree/main/config/scenarios)
are complete examples; the [Guide](guide.md) explains when to use each.

Settings marked **required** must be present even when the data type that uses them is
not active.

## What to run

| Setting | Example | Meaning |
|---|---|---|
| `tools_to_run` **required** | `["SPA", "MuSiCal"]` | Tools to compare. Names as in the [tool table](https://github.com/mkatsanto/SigFitBench#fitting-tools), spelled exactly. |
| `data_types` | `["synth_realistic"]` | Where catalogs come from: `synth_realistic`, `synth_custom`, `real`. Several allowed. |
| `relevance_modes` | `["standard"]` | Which list of reference signatures the tools fit against. `"standard"`: all signatures in the COSMIC file, minus `out_of_reference_sigs`. `"relevant"`: only the signatures listed in `relevant_signatures`; results go in folders ending in `_relevant`. `["standard", "relevant"]` runs both. Most runs use `["standard"]`. |
| `relevant_signatures` | `null` | Used only by the `"relevant"` mode: the signatures the tools may use, e.g. `["SBS1", "SBS5", "SBS18"]`. Leave `null` otherwise. |

## Simulated data (`synth_realistic`, `synth_custom`)

| Setting | Example | Meaning |
|---|---|---|
| `cancer_types` | `["Head-SCC", "Breast"]` | Tissues to simulate (`synth_realistic` only). Each needs `cosmic_tissue_data/signature_contributions_{WGS or WES}_{tissue}.dat`. |
| `num_muts` **required** | `[1000, 10000]` | Mutations per sample. One run per value. Use `[]` when only `real` is active. |
| `num_samples` | `30` | Different signature mixes per tissue (`synth_realistic` only). |
| `num_realizations` | `1` | Independent samples generated from each signature mix. With `3`, every mix gives 3 samples that have the same known answer but different random mutations. |
| `out_of_reference_sigs` **required** | `["SBS86", "SBS90"]` | Unknown signatures: taken out of the list the tools fit against, but still used to simulate samples (`synth_realistic`). Use `[]` for none. See [Unknown signatures](guide.md#unknown-signatures). |
| `out_of_reference_weights` | `[0.1]` | One value per unknown signature added to each sample: its share of mutations. `[0]` or `[]` adds none. |
| `random_seed` | `0` | Same seed, same simulated data. |
| `wgs_or_wes` | `"WGS"` | Tissue profiles for whole genomes (`WGS`) or exomes (`WES`). Fewer tissues exist for `WES`. |

Tissues available for `WGS`: Breast, CNS-GBM, ColoRect-AdenoCA, Eso-AdenoCA, Head-SCC,
Liver-HCC, Lung-AdenoCA, Lymph-BNHL, Ovary-AdenoCA, Panc-AdenoCA, Skin-Melanoma,
Stomach-AdenoCA.

### Custom tables

```yaml
custom_datasets:
  MyMix:                                          # name in result folders and the report
    combinations_file: "config/scenarios/custom_mixtures.tsv"   # your table: rows = samples, columns = signatures
    out_of_reference_sigs: ["SBS90"]              # optional: in the table, but not given to the tools
    signature_file: "input/COSMIC_v3.4_SBS_GRCh38.txt"   # optional: other COSMIC version
```

`num_samples` and the top-level `out_of_reference_*` settings are not used for custom
tables. See [Custom mixtures](guide.md#custom-mixtures).

## Real data (`real`)

```yaml
real_datasets:
  PCAWG:
    catalog: "real_mutational_catalogs/real_samples_WGS_PCAWG-146_input_profiles.dat"
    cancer_types: ["Breast-AdenoCA", "Skin-Melanoma"]   # optional: only these tissues
  MyCohort:
    catalog: "input/my_cohort.tsv"                      # 96 mutation types x samples
  MyVCFs:
    vcf_dir: "input/my_vcfs"                            # one VCF per sample instead of a catalog
    reference_fasta: "input/ref/genome.fa"              # same genome build as the VCFs, with a .fai index
```

An entry with `vcf_dir` instead of `catalog` has its catalog built by the workflow (SigVCF), one column per VCF. SigVCF runs
from its source folder: set `sigvcf_source` (default: the `SigVCF` folder next to SigFitBench).

Each entry is one catalog. Columns named `Tissue::Sample` are split into one group per
tissue; `cancer_types` keeps only the tissues listed. Catalogs without `::` in their
column names are fitted as one group, named after the entry (`MyCohort`).

See [Own catalog](guide.md#own-catalog) and [Real PCAWG](guide.md#real-pcawg).

## Files and folders

| Setting | Example | Meaning |
|---|---|---|
| `output_dir` **required** | `"results_my_run"` | All results. Use a new one per run. |
| `log_dir` **required** | `"logs_my_run"` | Logs of every job. |
| `input_dir` | `"input"` | Folder with reference files. |
| `signature_contributions_dir` **required** | `"cosmic_tissue_data"` | Tissue profiles used for simulation. |
| `cosmic_signatures.base_file` **required** | `"input/COSMIC_v3.3.1_SBS_GRCh38.txt"` | COSMIC signatures, used to simulate and (minus `out_of_reference_sigs`) to fit. |
| `matrix_order` **required** | `"input/mut_matrix_order_MutationalPatterns.dat"` | Order of the 96 mutation types. Leave as is. |
| `epsilon` | `0.00000001` | Small number to avoid dividing by zero. Leave as is. |

COSMIC files in `input/`, named `COSMIC_v<version>_SBS_<build>.txt`:

| Version | Genome builds | Notes |
|---|---|---|
| 3.6 | GRCh37, GRCh38 | Newest. Adds `SBS22c` and `SBS100`–`SBS113` to v3.4. |
| 3.4 | GRCh37, GRCh38 | Splits `SBS22` and `SBS40` into `SBS22a/b` and `SBS40a/b/c`; adds `SBS96`–`SBS99`. |
| 3.3.1 | GRCh37, GRCh38 | Default. Needed for `synth_realistic`, whose tissue profiles use its signature names. |
| 3, 3.1, 3.2, 2, 1 | GRCh37 and/or GRCh38 | Older releases. |

Use the build your mutations were called on. Newer versions work for `synth_custom`
(through `signature_file`) and `real`.

## Scoring

```yaml
evaluation:
  min_contribution: 10.0          # fewer mutations than this = signature counted as absent
  run_bootstrap_quality: false    # in development; leave false
  num_bootstrap: 1
```

All three keys are **required**.

## Tool settings

Under `tool_overrides`, keyed by the tool name in lower case. Two are **required** when
their tool runs:

```yaml
tool_overrides:
  mutationalpatterns:
    strict: true                 # true = drop signatures that barely improve the fit
  sigminer:
    method: "NNLS"               # "NNLS", "QP" or "SA"
```

Optional:

| Tool | Setting | Default | Meaning |
|---|---|---|---|
| `signature_tools` | `method` | `"KLD"` | `"KLD"` or `"NNLS"` |
| | `exposureFilterType` | `"giniScaledThreshold"` | or `"fixedThreshold"` |
| | `nboot` | `100` | Bootstrap repeats used to drop uncertain signatures |
| | `nparallel` | `1` | Cores per job |
| `msigact` | `max_level` | `5` | Remove up to this many signatures at once when testing. Above 2 is slow with all COSMIC signatures. |
| | `p_thresh` | `0.05` | Remove a signature only if the test p-value is at least this |
| | `likelihood` | `"neg.binom"` | or `"multinom"` |
| | `threads` | `4` | Cores per job |
| | `extra` | `{}` | Any other mSigAct `SparseAssignActivity` argument, e.g. `max.subsets: 2500` (any tool can have this, see below) |
| `starsigndna` | `n_iterations` | `1000` | Optimisation steps |

`SPA` and `SPA_low` have no settings; they differ only in the removal penalty (0.05 and
0.02).

## Extra tool arguments

To use a tool option that has no setting above, add it under `extra` for that tool, in the run's
config. Use the tool's own argument name:

```yaml
tool_overrides:
  sigminer:
    extra:
      auto_reduce: true      # argument of the R function sig_fit()
  sigtrim:
    extra:
      --verbose: true        # flag of the sigtrim command: true gives the flag alone, false leaves it out
```

- **Priority:** a value here replaces the workflow's own value for that argument.
- **Blocked arguments:** the workflow stops with an error if you set an argument it controls itself,
  such as the input data, the output or the number of cores. The full list is
  `PROTECTED_ARGUMENTS` in `workflow/rules/common.smk`.
- **Wrong arguments:** a tool that rejects an argument is reported like any other tool failure:
  its result file says `# TOOL FAILED`, and the reason is in its log.
- **Reruns:** only the tool whose values changed reruns. Tools with no values set are not affected.
- **Not included:** SigNoise and SigNoiseID, which take all their options from `tool_overrides`.
- **`signature_tools`:** single values only, no lists.

## Changing defaults without editing the config

Any setting can be given on the command line:

```bash
snakemake ... --configfile config/scenarios/quick_test.yaml \
  --config num_samples=10 output_dir=results_10_samples log_dir=logs_10_samples
```
