# vibrio-integron-pipeline — Development log

**Student:** Liam Nguyen (48596804)
**Course:** SCIE3220 research project
**Supervisor:** Jan Engelstädter (Engelstädter Lab)
**Compute:** UQ Bunya HPC cluster (micromamba environments, SLURM)

This repository holds the day-to-day development of the project: detection-tool
benchmarking, the genus-wide survey, and the first interactive visualiser prototype.
The visualiser was later packaged and released as a standalone tool in its own
repository, [NqynX/intervis](https://github.com/NqynX/intervis); see that repository's
development log for the release period (late September onward). This log records the
work in this repository by period, with the commits made in each.

## Repository layout

- `data/` — reference genomes and inputs
- `scripts/` — detection, reconstruction, annotation and analysis scripts
- `results/` — detection outputs, cassette tables and figures
- `intervis/` — the visualiser prototype (engine, CLI, R/Shiny app), later split into NqynX/intervis
- `db/` — DIAMOND Swiss-Prot database
- `tool-benchmark/` — synthetic-genome runs for the final benchmark

## Development by period

### Period 1 — Project setup and first reconstruction (7–10 August 2026)

Set the project up and got from a genome to a drawn cassette array. The detection
environment was pinned in the first week so every later result is reproducible, and the
reconstruction step was built first because every downstream analysis and the visualiser
depend on it.

- **7 Aug 2026 · `9d7d79f` · Initial commit** — created the repository and project structure.
- **7 Aug 2026 · `556352e` · README.md** — stated the aim: detect and visualise integron cassette arrays in Vibrio genomes.
- **7 Aug 2026 · `964984c` · Initial repository structure** — set up the data, scripts and results folders.
- **7 Aug 2026 · `9119b5e` · Add files via upload** — added the first input files and working notes.
- **7 Aug 2026 · `17de5c0` · Record IntegronFinder 2.0.2 environment** — fixed the main detection tool's version so runs are reproducible.
- **10 Aug 2026 · `f88a5f3` · Add cassette array reconstruction from IntegronFinder output** — the base step that all later analysis and the viewer depend on.
- **10 Aug 2026 · `50b82f5` · Add integron metadata and cassette length annotation** — added integron metadata and per-cassette lengths to the reconstruction.
- **10 Aug 2026 · `b7950fe` · Add first integron cassette array visualisation** — drew the first picture of a reconstructed array.
- **10 Aug 2026 · `8092880` · Improve integron cassette array visualisation layout** — reworked the layout so it reads more clearly.
- **10 Aug 2026 · `47af50e` · Add initial project meeting notes** — recorded notes from the first meeting with Jan (scope and direction).

### Period 2 — Detection tools and tool comparison (27–28 August 2026)

Installed the second detector and compared the detectors on a real superintegron on a
common basis — a shared genome with a swept HattCI score threshold, not a single
arbitrary cut-off. The I-VIP decision was recorded with its reason rather than silently
dropped.

- **27 Aug 2026 · `9700e80` · Phase 1: Install and record HattCI environment** — installed HattCI and recorded its environment for reproducibility.
- **27 Aug 2026 · `2eb00f5` · Update tool config with HattCI installation details** — added the exact install details to the tool config.
- **27 Aug 2026 · `e5de480` · Phase 1 complete: document I-VIP exclusion, finalise tool config and survey** — recorded why I-VIP is set aside (class-1 scope; no integrons on these genomes).
- **28 Aug 2026 · `9e53940` · Phase 2-3: N16961 genome, IntegronFinder and HattCI initial run on chrII** — ran both detectors on the N16961 chromosome II superintegron.
- **28 Aug 2026 · `e34bfda` · Phase 2-3: initial tool runs on N16961 chrII with Vscore analysis** — added HattCI Viterbi-score analysis to choose a sensible threshold.
- **28 Aug 2026 · `bef0894` · Phase 3.2: tool comparison script with Vscore threshold analysis** — compared the detectors across a swept threshold, not one fixed cut-off.
- **28 Aug 2026 · `fa96b07` · Phase 3: complete tool comparison stats and summary paragraph** — finished the comparison statistics and a written summary.
- **28 Aug 2026 · `431e181` · Phase 4: cassette table and sequence extraction (59 cassettes from N16961 chrII)** — extracted the cassette table and sequences for annotation.
- **28 Aug 2026 · `f0ae0cd` · Add zoomed integron visualisation script** — added a zoomed view to inspect single cassettes and their attC sites.

### Period 3 — Annotation and multi-genome pipeline (29 August 2026)

Gave cassettes real names and automated the whole path, then ran it across more genomes.
The low and highly variable cassette counts were the first clear sign that detection
settings strongly affect the result — the observation that later motivated the formal
benchmark. A wrong-genome input was corrected as soon as it was spotted, to stop it
propagating.

- **29 Aug 2026 · `da14f9a` · Phase 6: cassette annotation (BLASTx Swiss-Prot + CARD/RGI)** — first real gene names on the cassettes.
- **29 Aug 2026 · `10d160f` · Add automated pipeline script for integron detection, cassette reconstruction, and annotation** — combined the three steps into one pipeline script.
- **29 Aug 2026 · `1134471` · Phase 7: Run pipeline on 4 additional Vibrio genomes** — variable counts showed detection settings strongly change the result.
- **29 Aug 2026 · `5403638` · Add comparative BLASTx analysis, update .gitignore, remove old wrong-genome results** — removed earlier results that used the wrong genome, before the error spread.
- **29 Aug 2026 · `f8e84f3` · Add files via upload** — added supporting data files.

### Period 4 — Genus-wide survey and quantitative analyses (31 August 2026)

Put the single-genome observations on a population footing. The Cayley/permutation result
separates cassette-order rearrangement from overall genome relatedness, which whole-genome
ANI alone cannot do, and recording versions here fixed the reproducibility baseline for the
whole project.

- **31 Aug 2026 · `7282594` · Genus-wide survey: 988 complete Vibrionaceae genomes** — ran the detection pipeline across the genus.
- **31 Aug 2026 · `f89a7e2` · HattCI concordance (84%), per-genome DefenseFinder (3741 systems, 78% of genomes)** — measured HattCI–IntegronFinder concordance and ran DefenseFinder.
- **31 Aug 2026 · `e1bd318` · Cayley/permutation analysis: order conserved above 99% ANI, random below** — measured cassette-order rearrangement separately from genome relatedness.
- **31 Aug 2026 · `47d2836` · Ignore SLURM logs, regenerable intermediates, and unused I-VIP submodule** — kept only source and inputs tracked.
- **31 Aug 2026 · `40eaad7` · Record software and database versions** — reproducibility baseline for the project.

### Period 5 — Interactive visualiser and distance network (11–18 September 2026)

Turned the static plots into an interactive, attC-resolved viewer with a comparison mode.
Ribbons were made a similarity gradient instead of a flat shared-family colour so each link
shows how conserved it is, and kept neutral until the identity step runs rather than showing
fabricated values. The deliverable shape was reviewed with J. Engelstädter: an R/Shiny app
is the webpage tool, the input is a genome, attC belongs to its cassette and is revealed on
zoom, comparison is pairwise, and annotation comes from public databases.

- **11 Sep 2026 · `5815a89` · Add v0.1 integron visualiser: .integrons parser + plot_integron()** — first version of the interactive viewer.
- **11 Sep 2026 · `a51a5dc` · First N16961 visualiser render (v0.1)** — confirmed the viewer works end to end on real data.
- **11 Sep 2026 · `0c739da` · v0.2 multi-track region view: gene arrows + attC + length** — attC drawn as a feature, not just a coordinate.
- **11 Sep 2026 · `28497bc` · N16961 per-cassette function categories** — cassettes can be coloured by function.
- **11 Sep 2026 · `3830fa2` · Multi-genome comparison data (N16961, V. cholerae, V. parahaemolyticus)** — input for the pairwise comparison view.
- **11 Sep 2026 · `6a57cbf` · Add Cayley/MTF distance-network builder and survey network_data.json** — shows cassette-order change between genomes, which ANI does not.
- **11 Sep 2026 · `3b7cdb7` · comparison_v2: attC sites + per-family % identity** — ribbons shaded by per-pair amino-acid identity; each link shows how conserved it is.
- **18 Sep 2026 · `9376603` · Add example .integrons (2 Vibrio superintegrons) for intervis comparison test** — a reproducible comparison test case.
- **18 Sep 2026 · `14c1318` · Add intervis annotation tables (functions, families, identities)** — the annotation tables the views use; the tool was then split into the separate NqynX/intervis repository for release.

## Reproducibility

Detection used IntegronFinder 2.0.2 (Prodigal 2.6.3, HMMER 3.3.2, INFERNAL 1.1.4) and
HattCI 1.0b, run in micromamba environments on Bunya. Software and database versions were
recorded in-repo (`40eaad7`). The visualiser prototype's engine, CLI and R/Shiny app, and
the full released tool, are in [NqynX/intervis](https://github.com/NqynX/intervis).
