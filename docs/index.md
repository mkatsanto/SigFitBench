# SigFitBench

A Snakemake workflow that compares mutational signature fitting tools. It runs 14
published tools on the same mutation catalogs, scores every result against the same
known answer, and writes one interactive HTML report.

## Get started

```bash
git clone https://github.com/mkatsanto/SigFitBench.git
cd SigFitBench
conda env create -f install/environment.yml
conda activate sigfitbench

snakemake --profile profiles/local-conda \
  --snakefile workflow/Snakefile --configfile config/scenarios/quick_test.yaml
```

The first run builds one conda environment per tool (30–60 minutes).

## Pages

- **[Guide](guide.md)**: the questions SigFitBench answers, each with a ready-to-run config
- **[Configuration](configuration.md)**: every setting
- **[Results and metrics](results.md)**: output files, scores and the report

Installation, the list of tools, cluster use and citations are in the
[README](https://github.com/mkatsanto/SigFitBench#readme).
