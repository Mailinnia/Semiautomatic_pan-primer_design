# Semiautomatic pan-primer design
Semiautomatic pan-primer design is a semi-automated, multi-step workflow for clustering-based reduction of viral sequence databases and downstream pan-primer design with varVAMP.
The workflow is intended to reduce the influence of heavily represented viral groups while retaining representative sequence diversity for downstream multiple sequence alignment and broad-range primer design. It
was developed and evaluated using _Orthoflavivirus_ and _Alphavirus_ as proof-of-concept datasets.
## Table of Contents
- [Workflow overview](#overview_)
- [Installation and configuration](#install_)
- [Input requirements](#input_)
- [Quick start](#quick_)
- [Documentation](#doc_)
  - [Preparing the input database](#prep_)
  - [Clustering and database reduction](#install_)
    - [Scoring system](#score_)
    - [Clustering report and outputs](#cluster_out_)
    - [Tools and performance notes](#tools_)
  - [Extracting the reduced database](#extract_)
  - [Alignment of reduced database](#align_)
  - [Primer design with varVAMP](#design_)
    - [Primer prediction](#pred_)
      - [Local BLAST databases](#blast_)
  - [Citation](#cite_)
  - [Disclaimer](#disc_)


# <a name="overview_"></a> Workflow overview
<p align="center">
<img src="https://github.com/user-attachments/assets/a2986f5c-70f1-4cc5-908d-d18100036fa9" width="800" height="" >
</p>

The Semiautomatic pan-primer design workflow consists of the following sequential stages:

0. Database preparation – A curated viral sequence database and corresponding metadata are prepared according to the target virus group and intended application. Database construction and filtering are not standardized within the Semiautomatic pan-primer design workflow, as appropriate inclusion criteria are database-specific.
1. Clustering and database reduction – Sequences are grouped into major clusters and subclusters, representative sequences are selected, and the clustering and database-reduction results are summarized in an HTML report.
2. Multiple sequence alignment – The reduced database is aligned. MAFFT is supported by the provided alignment utility, but an existing compatible alignment can also be supplied.
3. Primer design and prediction – The alignment is processed with varVAMP across combinations of consensus thresholds and maximum ambiguity settings. Candidate primer sets can be evaluated against on-target and off-target BLAST databases, and the results are collated in an HTML report. Primer prediction can also be run independently on existing primer sets.

# <a name="install_"></a> Installation and configuration

