# Development of a Bioinformatics Pipeline for Identifying and Visualising Integron Cassette Arrays in Vibrio Genomes

## Project overview

Integrons enable bacteria to capture, rearrange and express gene cassettes, with cassette position influencing gene expression. This project investigates large chromosomal integrons in publicly available *Vibrio* genomes, including complex superintegrons.

Existing bioinformatics tools will be evaluated for their ability to:

- identify integron regions;
- detect gene cassettes;
- reconstruct cassette order; and
- visualise cassette-array organisation.

Based on this evaluation, a reproducible bioinformatics pipeline will be developed to extract integron regions, organise the resulting information and generate clear visualisations of cassette organisation. The workflow will provide a computational foundation for future quantitative studies of integron dynamics and for testing predictions from evolutionary models using genomic data.

## Project details

- **Project period:** 27 July 2026 – 25 November 2026
- **Laboratory:** Engelstädter Lab
- **Primary supervisor:** Jan Engelstädter
- **Student:** Liam Nguyen
- **Institution:** The University of Queensland

## Repository purpose

This repository is used for project record keeping, pipeline development and version control. It will contain:

- literature and software evaluation notes;
- data acquisition and processing scripts;
- integron-detection workflows;
- cassette-array reconstruction scripts;
- visualisation scripts;
- configuration and environment files;
- analysis logs and progress records; and
- documentation required to reproduce the workflow.

## Repository structure

```text
vibrio-integron-pipeline/
├── README.md
├── docs/
│   ├── meeting_notes/
│   ├── literature_notes/
│   └── workflow_documentation/
├── config/
├── scripts/
│   ├── data_download/
│   ├── preprocessing/
│   ├── integron_detection/
│   ├── cassette_reconstruction/
│   └── visualisation/
├── notebooks/
├── results/
│   ├── tables/
│   └── figures/
├── logs/
├── environments/
└── tests/
