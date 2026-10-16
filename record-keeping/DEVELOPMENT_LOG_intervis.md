# intervis — Development log

**Student:** Liam Nguyen (48596804)
**Course:** SCIE3220 research project
**Supervisor:** Jan Engelstädter (Engelstädter Lab)

intervis is an integron-native viewer for Vibrionaceae cassette arrays. It runs
IntegronFinder, reconstructs cassette intervals between consecutive *attC* sites, names
cassette proteins by DIAMOND blastp against Swiss-Prot, and renders a zoomable,
*attC*-resolved view with built-in pairwise array comparison, behind both an R/Shiny
interface and a command-line interface.

The tool began as a prototype inside the
[vibrio-integron-pipeline](https://github.com/NqynX/vibrio-integron-pipeline) repository
and was split into this repository for release. See that repository's development log for
the earlier history (project setup, detection benchmarking, the genus-wide survey, and the
visualiser prototype — August to mid-September 2026). This log records the release period.

## Development by period

### Period 6 — Packaging and v1.0.0 release (28 September 2026)

Brought the prototype together into one engine behind an R/Shiny interface and a CLI, and
made the first public release. Detection was kept as an established method (IntegronFinder)
rather than writing a new *attC* classifier; the contribution is the integron-native,
end-to-end path with *attC*-resolved, zoomable cassettes and built-in pairwise comparison.

- One engine (run IntegronFinder `--local-max` → reconstruct cassette intervals → name cassette proteins by DIAMOND blastp against Swiss-Prot → emit viewer data), with the R/Shiny app and CLI on top.
- Tool split from the development repository into this repository; **v1.0.0** tagged and released.

### Period 7 — Out-of-the-box fix (9 October 2026, `6cb0bca`)

A clean clone by the supervisor had shown only the demo genomes and nothing on *Run*,
because the reference genomes were git-ignored and the environment file was incomplete.

- Committed the built-in reference genomes so a fresh clone ships with them.
- Completed `environment.yml` (IntegronFinder + DIAMOND + Python; R/Shiny for the GUI).
- Fixed backend path resolution so the app runs out of the box.
- Commit `6cb0bca` (`git describe`: v1.0.0-3-g6cb0bca) — the analysis version cited in the report.

Built-in examples are auto-discovered from `app/examples/`: any file stem with a `.fna`
(and optional `.integrons`/GenBank) becomes a selectable reference genome, with the label
taken from the FASTA header.

### Period 8 — Final benchmark and synthetic benchmark (early–mid October 2026)

Evaluated the detection step formally. Five genomes with experimental evidence were scored
(IntegronFinder default and `--local-max`, HattCI, I-VIP, BacAnt); F1 and threshold-independent
average precision were computed against the IntegronFinder local-search reference. A synthetic
benchmark then measured recall against a known 80-site *attC* catalogue under random and
compensatory divergence (seed 20260101). Run on Bunya; scripts and outputs added to this repo.

### Period 9 — Analysis, figures and records added (14 October 2026)

Published the analysis scripts, figures, the detection-benchmark workbook and the project
record to the repository (commits `49730c5`, `9aa648b`, `fd714da`), so the tool repository
also carries the evaluation that supports the report.

## Reproducibility

- **Detection:** IntegronFinder 2.0.2 (`--local-max`, `--circ`, `--keep-palindromes`, max inter-*attC* 10 kb).
- **Annotation:** DIAMOND 2.2.5 blastp of translated cassette proteins vs UniProtKB/Swiss-Prot (575,748 sequences, 11 Jun 2026); naming by top hit (`--max-target-seqs 1 --evalue 1e-5 --more-sensitive`).
- **Comparison:** reciprocal best hits (`-k 50 --evalue 1e-5 --more-sensitive`), retained at ≥30% amino-acid identity and ≥0.5 coverage of the shorter protein; homology families by single-linkage clustering of threshold-passing hits.
- **Environment:** `environment.yml` (micromamba/conda); R ≥ 4.1 with `shiny` for the GUI.
- **Case study:** N16961 chromosome II, NC_002506.1 (assembly GCF_000006745.1).
