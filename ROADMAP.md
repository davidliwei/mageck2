# MAGeCK2 Roadmap

MAGeCK2 is the next generation of [MAGeCK](https://sourceforge.net/projects/mageck/) (Model-based
Analysis of Genome-wide CRISPR-Cas9 Knockout), the widely adopted open-source toolkit for analyzing
pooled CRISPR screens. This document describes where the project is, where it is going, and how the
community can get involved. It is a living document and will be revised as priorities evolve.

## Vision

MAGeCK2 aims to be the de facto AI-ready, open-source tool for pooled CRISPR screen analysis —
callable from Python notebooks, agentic copilots, and standard bioinformatics pipelines, with data
structures that interoperate with the broader single-cell and perturbation-modeling ecosystem.

## Where we are today

MAGeCK2 builds on MAGeCK and currently adds:

- Paired-sample analysis
- UMI-aware counting
- Paired-guide (dual-sgRNA) screen counting

Core functionality carried forward from MAGeCK: sgRNA counting, quality control, RRA and MLE
essentiality testing with copy-number correction, and functional enrichment.

## Planned directions

### 1. An AI-ready Python core with a standardized data representation

- A documented, stable **Python API** for the core operations (counting, testing/ranking, QC), with
  the command line re-expressed as a thin layer over that API.
- **R access via a [reticulate](https://rstudio.github.io/reticulate/) wrapper** over the same Python
  codebase, rather than a separately maintained R implementation.
- Adoption of **[AnnData](https://anndata.readthedocs.io/)/[MuData](https://mudata.readthedocs.io/)**
  as a standard CRISPR-screen data representation, proposed to the
  [scverse](https://scverse.org/) community via a public RFC, so MAGeCK2 outputs flow into the
  single-cell and perturbation-modeling stack.
- Integration points for LLM agents and scientific copilots, including
  [MCP](https://modelcontextprotocol.io/)-compatible endpoints and vendor-neutral agent skill
  packages.
- Standardized benchmark datasets and an evaluation harness for screen-analysis tasks.

### 2. Scalability and ecosystem integration

- **CPU-parallel performance** improvements on the counting and ranking paths to support
  genome-scale and emerging large-library screens.
- Modernized first-class integrations with the wider ecosystem (e.g.
  [nf-core/crisprseq](https://nf-co.re/crisprseq), [Galaxy](https://galaxyproject.org/), and
  [Latch.bio](https://latch.bio/)) onto the new API and data schema.
- Refreshed tutorials and example datasets, and community engagement (including workshops).

### Engineering practices (ongoing)

Robust tests, continuous integration, versioned releases, and clear documentation are maintained
continuously as part of every change, not as a separate effort. Releases are distributed through
standard channels (PyPI, Bioconda, Bioconductor).

## Out of scope (for now)

To keep the project focused, the following are intentionally **not** planned at this stage:

- A separately maintained, parallel R reimplementation (R is supported via the reticulate wrapper).
- GPU acceleration (the performance work targets CPU parallelism).
- Bundling external multi-omics datasets (e.g. DepMap) as an integrated data layer.
- Development of new AI/ML models. MAGeCK2 is the *software* that prepares data, features, and
  benchmarks for such models.

## Getting involved

Suggestions and contributions are welcome. Please open an
[issue](https://github.com/davidliwei/mageck2/issues) to propose features or report problems, see
[CONTRIBUTING.md](CONTRIBUTING.md) for how to contribute code, and join the
[MAGeCK user forum](https://groups.google.com/g/mageck) for discussion.
