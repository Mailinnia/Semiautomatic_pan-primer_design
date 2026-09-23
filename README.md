# Semiautomatic pan-primer design
Semiautomatic pan-primer design is a semi-automated, multi-step workflow for clustering-based reduction of viral sequence databases and downstream pan-primer design with varVAMP.
The workflow is intended to reduce the influence of heavily represented viral groups while retaining representative sequence diversity for downstream multiple sequence alignment and broad-range primer design. It
was developed and evaluated using _Orthoflavivirus_ and _Alphavirus_ as proof-of-concept datasets.

[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.22110656.svg)](https://doi.org/10.5281/zenodo.22110656)

## Table of Contents
- [Workflow overview](#overview_)
- [Installation and configuration](#install_)
- [Input requirements](#input_)
- [Quick start](#quick_)
- [Documentation](#doc_)
  - [Preparing the input database](#prep_)
  - [Clustering and database reduction](#cluster_)
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
### Installation
Install the required dependencies listed in `semiautomatic.yml` in an appropriate Conda environment. For example:

```
conda create -p <path/to/environment>
conda activate <path/to/environment>
conda env update -f semiautomatic.yml
```
The environment includes the Python and R dependencies and command-line
tools required by the workflow, including Vclust, VSEARCH, CD-HIT-EST,
Mash, MAFFT, varVAMP, NCBI BLAST+, SeqKit, and Krona.

### Configuration
Two scripts, `functions_primer_check.py` and `db_taxonomy.py` require a local Krona taxonomy database.
After installing Krona, navigate to the Krona installation directory and run `updateTaxonomy.sh`.
```
cd <path/to/krona>
./updateTaxonomy.sh
```

By default, the taxonomy database `taxonomy.tab`is created in `<path/to/krona/taxonomy`. Provide the full path to `taxonomy.tab` in both `functions_primer_check.py` and `db_taxonomy.py`.
```
KT_TAXONOMY_DB = "/path/to/krona/taxonomy/taxonomy.tab"
```

### Command-line access

For convenient command-line use, it is recommended to make the main workflow scripts executable:

```
chmod +x cluster_db.py cluster_extract.py mafft_entropy.py primer_design.py db_taxonomy.py primer_check.py
```
The scripts can then be linked to a directory included in the user's `PATH`, for example a personal `bin` directory, using symbolic links:
```
ln -s </full/path/to/semiautomatic>/cluster_db.py <path/to/bin>/cluster_db.py
...
ln -s </full/path/to/semiautomatic>/primer_check.py <path/to/bin>/primer_check.py
```
Ensure that the selected `bin` directory is included in `PATH`.


# <a name="input_"></a> Input requirements
The required input depends on the stage of the workflow:

### Clustering and database reduction

- A nucleotide FASTA file containing the input sequence database.
- A matching TSV metadata file containing accession identifiers and relevant metadata.

Database construction and filtering should be tailored to the target viral group and intended application.

### Multiple sequence alignment

- The extracted reduced database in FASTA format.

The provided alignment utility uses MAFFT. Alternatively, an existing compatible multiple sequence alignment can be supplied directly to the primer-design workflow.

### Primer design

- The aligned reduced database in FASTA format.

The alignment is used as input for varVAMP and evaluated across the specified primer-design parameter combinations.

# <a name="quick_"></a> Quick start
A typical analysis proceeds as follows.

**1. Cluster and reduce the input database**
```
cluster_db.py -f FASTA -m META
```
**2. Extract the selected reduced database**
```
cluster_extract.py -f FASTA -d REDUCED-DB
```
**3. Align the reduced database:**

```
mafft_entropy.py -f REDUCED-FASTA
```
**4. Perform primer design on the aligned reduced database**
```
primer_design.py  -a REDUCED-ALN
```

Run any script with `-h` for the complete set of options.

# <a name="doc_"></a> Documentation
The sections below describe each workflow stage, scoring system, output structure, primer prediction, and construction of local BLAST databases.

## <a name="prep_"></a> Preparing the input database
Semiautomatic takes a curated viral sequence database as input. Database construction and filtering are intentionally not standardized by the workflow, as appropriate inclusion criteria depend on the target viral group and study objective.

One approach is to identify the viral family, genus, or species of interest in NCBI, apply appropriate filtering criteria (e.g. sequence length), and download the corresponding GenBank files. The resulting dataset should be curated before use, including removal of sequences that are unsuitable for representing naturally occurring viral diversity, such as synthetic constructs, chimeric sequences, and other sequences that do not meet the study-specific inclusion criteria. The GenBank files can subsequently be processed to extract sequence metadata.

At minimum, the metadata file must contain the columns `accession` and `viral_species`, with accessions corresponding to those in the input FASTA file. The use of current ICTV species nomenclature for `viral_species` is recommended.

## <a name="cluster_"></a> Clustering and database reduction
The script performs hierarchical clustering of genomic sequences and reduction of the database in three main stages: major clustering followed by subclustering and sequence selection.

```
cluster_db.py [-h] -f FASTA -m META [-c MAJOR_CLUSTER] [-s SUBCLUSTER] [-n NUMSEQ] [-p PREFIX] [-o OUT_DIR] [-z THREADS] [--max-len MAX_LEN] [--all] [--overcluster-cutoff OVERCLUSTER_CUTOFF]
                     [--align-dir ALIGN_DIR] [--mash MASH] [--report-title REPORT_TITLE] [--skip-major SKIP_MAJOR]
                     [--cluster-method {vclust-single,vclust-complete,vclust-uclust,vclust-cd-hit,vclust-set-cover,vclust-leiden,vsearch,cd-hit-est}] [-a PREFERRED_ACCESSIONS]

Cluster database

options:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Fasta file containing the original database
  -m META, --meta META  TSV file containing meta-data for the fasta file, such as "species"
  -c MAJOR_CLUSTER, --major-cluster MAJOR_CLUSTER
                        Range for major clustering. (default: 0.70 0.75 0.80 0.85)
  -s SUBCLUSTER, --subcluster SUBCLUSTER
                        Range for sub-clustering. (default: 0.96 0.97 0.98 0.99)
  -n NUMSEQ, --numseq NUMSEQ
                        Number of sequences per major cluster (default: 100)
  -p PREFIX, --prefix PREFIX
                        Prefix. Default is input fasta
  -o OUT_DIR, --out-dir OUT_DIR
                        Output directory. Default is current working directory
  -z THREADS, --threads THREADS
                        Number of CPU threads to use. (default: 20)
  --max-len MAX_LEN     Maximum sequence length for vsearch. (default: 50000)
  --all                 Subcluster all major clusters regardless of silhouette score selection (default: False)
  --overcluster-cutoff OVERCLUSTER_CUTOFF
                        Sets the cut-off for overclustering penalty for calculating the max acceptable sub_id (default: 0.2)
  --align-dir ALIGN_DIR
                        Provide the directory of the vclust alignment to skip this step. This is useful if vclust align has already been performed on the input database
  --mash MASH           Provide a precomputed mash distance matrix file to skip the matrix generation step. This is useful if you've already run <mash dist -t> on the input database and want to reuse the
                        output. (default: None)
  --report-title REPORT_TITLE
                        Title of report. Default is input filename
  --skip-major SKIP_MAJOR
                        Skip major clustering by recreating the clusters from the provided cluster tsv
  --cluster-method {vclust-single,vclust-complete,vclust-uclust,vclust-cd-hit,vclust-set-cover,vclust-leiden,vsearch,cd-hit-est}
                        Method to use for clustering (default: vclust-single)
  -a PREFERRED_ACCESSIONS, --preferred-accessions PREFERRED_ACCESSIONS
                        Path to file containing accessions to preferentially retain during database reduction. For each subcluster listed accessions are selected first (if present) and remaining slots are
                        filled by random sampling
```
**1. Major clustering**<br>
Clustering is first performed using a range of major cluster identity thresholds. Each clustering result is evaluated using *silhouette scores*, which quantify how well each cluster is separated from the others based on pairwise distances. Silhouette scores are calculated from a precomputed total average nucleotide identity (tANI) or Mash distance matrix (an alignment-free method that estimates pairwise sequence dissimilarity through a
_k-mer_-based sketch approximation), depending on the chosen clustering method. 

As **default**, the major cluster identity producing the highest silhouette score is selected for further subclustering. To continue with all major cluster IDs (regardless of which achieved the best silhouette score), use the flag: `--all`

**2. Subclustering**<br>
Within each major cluster, subclustering is performed to select representative sequences while preserving diversity and avoiding overclustering. The weighted overclustering penalty score is then calculated for each subclustering threshold. This score quantifies how much a subclustering exceeds the user-defined target number of sequences per major cluster (`--numseq`). A penalty of 0 indicates that the number of subclusters is within or below the target, while increasing values represent progressively stronger overclustering. 

**3. Sequence selection**<br>
Approximately *n* sequences (`--numseq`) are randomly selected from each major cluster using Python's *random.sample* method. At least one sequence is selected from every subcluster to maintain diversity. As a result: More than *n* sequences may be selected if a major cluster contains more than *n* subclusters. If a major cluster contains fewer than *n* sequences, all sequences are retained.

The choice of _n_ can have a substantial effect on the composition of the reduced database and should be adjusted to suit the input dataset. Consider both the number of sequences in each major cluster and the number of subclusters generated at the selected subclustering threshold. Higher subclustering thresholds may generate many subclusters, in which case the final number of selected sequences can exceed _n_. The relative abundance of different viral groups should also be considered. If a group of particular interest contains far fewer sequences than heavily sampled groups, reducing _n_ may help prevent the larger groups from dominating the reduced database and provide a more balanced representation.

### <a name="score_"></a> Scoring system
**Silhouette scores**<br> 
Major clustering scores are calculated using silhouette scores based on a distance matrix. The silhouette score measures how well each sequence fits within its assigned cluster compared to neighboring clusters. Scores range from -1 to 1, where values approaching 1 indicate well-separated clusters, values near 0 indicate overlapping clusters, and negative values indicate that sequences may be more similar to a neighboring cluster than to their assigned cluster.

**Overclustering penalty scores**<br>
The overclustering penalty is used to identify subclustering thresholds that split major clusters into more subclusters than needed for the requested database reduction. It is calculated relative to the user-defined number of target sequences per major cluster (`--numseq`).

* A penalty of 0 means that the number of subclusters does not exceed `--numseq`.
* A higher penalty indicates increasing fragmentation into more subclusters than the target.
* The mean penalty gives equal weight to every major cluster.
* The weighted penalty gives greater influence to major clusters containing more sequences and is used by the workflow when evaluating subclustering thresholds.

Because the penalty depends directly on `--numseq`, it should always be interpreted in the context of the selected target number of sequences. The value of `--numseq` should therefore be chosen with consideration of the size and composition of the input database and the level of subclustering required.

### <a name="cluster_out_"></a> Clustering report and outputs
<p align="center">
<img src="https://github.com/user-attachments/assets/66c98638-a22e-46e2-b129-b5e11f20d4b7" width="600" height="" >
</p>

An HTML report summarizes the clustering and database-reduction results using interactive sunburst plots, clustering scores, and tables highlighting clusters that contain multiple viral_species as well as viral_species split across multiple clusters.

**Directories and files**

|Directory | Outputs|
|---|---|
|1_dist| Contains the vclust-align files or the mash distance matrix|
|2_clustering | Contains the output for the major clusters and the subclusters, respectively |
|3_scores | Contains the calculated scores for the major clusters and subclusters, respectively |
|4_reduced_db | Contains the tsv files for each of the suggested reduced databases |

### <a name="tools_"></a> Tools and performance notes
The workflow supports Vclust, VSEARCH, and CD-HIT-EST. CD-HIT-EST is limited to sequence identity thresholds \>=0.80. Vclust supports several types of clustering methods: vclust-single, vclust-complete,
vclust-uclust, vclust-cd-hit, vclust-set-cover, and vclust-leiden. See [Clustering algorithms in Vclust](https://www.nature.com/articles/s41592-025-02701-7/figures/3) for more information.

#### Vclust

When using Vclust, the most time-consuming step is usually generating the ANI matrix. However, if the matrix has been computed previously, it can be reused with `--align-dir`. The script will then automatically
search for files ending in `align.tsv` and `align.ids.tsv` in the provided directory.

#### VSEARCH and CD-HIT

When using VSEARCH or CD-HIT, the most time-consuming steps are usually computing the Mash distance matrix, and performing major clustering itself. In the benchmark datasets used for development, VSEARCH required
substantially more runtime than CD-HIT-EST, whereas CD-HIT-EST produced less coherent clusters, particularly for datasets with greater variation in genome length.

If a Mash distance matrix has already been computed previously, you can reuse it with: `--mash <path/to/mash_matrix.tsv>`<br>
In the benchmark datasets used for development, Mash-derived silhouette scores were more sensitive to variation in genome length than scores derived from converted tANI values.

#### Skip reclustering
If major clusters have already been generated previously, you can skip re-clustering by providing the existing file using `--skip-major <path/to/major_cluster.tsv`. <br>
This will recreate the corresponding clusters without rerunning clustering. It is important to use the same clustering method that was used for generating the original file.
 
## <a name="extract_"></a> Extracting the reduced database
Following evaluation of the clustering results and HTML report, select the reduced database corresponding to the preferred major clustering and subclustering thresholds. The selected reduced database is provided as a TSV file containing the representative sequences chosen during database reduction. `cluster_extract.py` uses this file together with the original FASTA database to extract the selected sequences and generate the reduced FASTA file for downstream alignment.

```
 cluster_extract.py [-h] -f FASTA -d DB [-p PREFIX] [-o OUT_DIR] [-z THREADS]

Extract sequences from reduced database

options:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Fasta file containing the original database
  -d DB, --db DB        TSV file containing the reduced database meta data
  -p PREFIX, --prefix PREFIX
                        Prefix. Default is input database
  -o OUT_DIR, --out-dir OUT_DIR
                        Output directory. Default is current working directory
  -z THREADS, --threads THREADS
                        Number of CPU threads to use. (default: 20)
```


## <a name="align_"></a> Alignment of reduced database
The extracted reduced database is aligned prior to primer design. The provided `mafft_entropy.py` utility can be used to generate a multiple sequence alignment with MAFFT and visualize sequence variability across the alignment using an entropy plot. The entropy plot can be used to assess the overall alignment and identify conserved and variable regions.
```
mafft_entropy.py [-h] -f FASTA [-m {auto,einsi,linsi,ginsi}] [-p PREFIX] [-d OUT_DIR] [--aa] [-z THREADS] [-t PLOT_TITLE] [--skip-aln] [-i] [--memsave]

Create MAFFT alignment with entropy plot

options:
  -h, --help            show this help message and exit
  -f FASTA, --fasta FASTA
                        Fasta file for alignment
  -m {auto,einsi,linsi,ginsi}, --method {auto,einsi,linsi,ginsi}
                        Mafft method.(default: auto)
  -p PREFIX, --prefix PREFIX
                        Prefix. Default is input fasta
  -d OUT_DIR, --out-dir OUT_DIR
                        Output directory. Default is current working directory
  --aa                  Amino acid mode
  -z THREADS, --threads THREADS
                        Number of CPU threads to use. (default: 20)
  -t PLOT_TITLE, --plot-title PLOT_TITLE
                        Title of entropy plot. Default is input filename
  --skip-aln            Skip alignmnet. -f is assumed to be aligned file.
  -i, --plot-interactive
                        Interactive plot (default: False)
  --memsave             Use memsave instead of nomemsave when running MAFFT (default: False)
```

## <a name="design_"></a> Primer design with varVAMP
Once the alignment is ready, it can be supplied to the primer-design script, which runs varVAMP across user-defined parameter combinations.

```
primer_design.py [-h] -a ALN [-t THRESHOLDS] [-n MAX_AMBIGS] [-m {qpcr,single}] [-p PREFIX] [-b PRIMER_PREFIX] [-x OFF_TARGET_DB] [--allowance-3-prime ALLOWANCE_3_PRIME] [-y TARGET_DB]
                        [--target-mismatch TARGET_MISMATCH] [--opt-amp OPT_AMP] [--max-amp MAX_AMP] [--min-len MIN_LEN] [--opt-len OPT_LEN] [--max-len MAX_LEN] [--min-qamp MIN_QAMP] [--max-qamp MAX_QAMP]
                        [--min-temp MIN_TEMP] [--max-temp MAX_TEMP] [--opt-temp OPT_TEMP] [-o OUT_DIR] [-z THREADS] [-r REPORT_TITLE] [--blast-remote]

Primer design using varVamp

options:
  -h, --help            show this help message and exit
  -a ALN, --aln ALN     Alignment file
  -t THRESHOLDS, --thresholds THRESHOLDS
                        Consensus thresholds (default: "0.75 0.80 0.85")
  -n MAX_AMBIGS, --max-ambigs MAX_AMBIGS
                        List of max number of ambiguous characters in a primer (default: "3 4 5 6")
  -m {qpcr,single}, --mode {qpcr,single}
                        varVamp mode (default: qpcr). Automatically switches to single mode if no primers are found for qpcr.
  -p PREFIX, --prefix PREFIX
                        Output prefix. Default uses alignment name
  -b PRIMER_PREFIX, --primer-prefix PRIMER_PREFIX
                        Primer prefix (default: Primer).
  -x OFF_TARGET_DB, --off-target-db OFF_TARGET_DB
                        Path to off-target blast database. If provided, primers will be checked for potential off-targets
  --allowance-3-prime ALLOWANCE_3_PRIME
                        Number of mismatches allowed in the first 3 bases of the 3' end for off-target (default: 0)
  -y TARGET_DB, --target-db TARGET_DB
                        Path to target blast database. If provided, primers will be checked whether they can potentially catch all targets
  --target-mismatch TARGET_MISMATCH
                        Max number of mismatches when comparing to target blast db (default: 3).
  --opt-amp OPT_AMP     Optimal amplicon length when single mode is used (default: 150).
  --max-amp MAX_AMP     Max amplicon length when single mode is used (default: 300).
  --min-len MIN_LEN     Min primer length (default: 18).
  --opt-len OPT_LEN     Optimal primer length (default: 21).
  --max-len MAX_LEN     Max primer length (default: 24).
  --min-qamp MIN_QAMP   Min q-amplicon length when qPCR mode is used (default: 70).
  --max-qamp MAX_QAMP   Max q-amplicon length when qPCR mode is used (default: 200).
  --min-temp MIN_TEMP   Min primer temp (default: 56).
  --max-temp MAX_TEMP   Max primer temp (default: 63).
  --opt-temp OPT_TEMP   Optimal primer temp (default: 60).
  -o OUT_DIR, --out-dir OUT_DIR
                        Output directory. Default is <current directory>/primer-design_<prefix>_<date>
  -z THREADS, --threads THREADS
                        Number of CPU threads to use. (default: 10)
  -r REPORT_TITLE, --report-title REPORT_TITLE
                        Title of report. Default is prefix

```
**Designing primers**<br> 
The primer design workflow takes an alignment as input which is fed into varVAMP with different combinations of consensus thresholds (`--thresholds`) and different maximum ambiguous characters (`--max-ambig`) for the primers, in increasing order. The probe ambiguity is automatically calculated by varVAMP as _n–1_, where _n_ is the max primer ambiguity.

The workflow can be used in both qPCR mode and single mode. When the qPCR mode is used, if it fails for a given threshold and max ambiguity, it will automatically switch to single mode for those parameters. If the
qPCR mode fails for all ambiguities for a given threshold, the workflow will permanently switch to single mode as the consensus threshold increases. Additionally, if a given ambiguity fails in both qPCR mode
and single mode for a given threshold, that ambiguity will not be used as the consensus threshold increases.

It is possible to adjust the size ranges of the amplicons, size of primers, and primer temperatures. The qPCR probe temperature will automatically be adjusted when changing the primer temperature.

The outputs of each threshold are collated and primer sets are compared and duplicates are filtered out, e.g. a primer set found at threshold 0.85 could also be found at threshold 0.80, if so, only the primer set
at 0.85 is kept (but both will be shown in the figures). Duplicate primer sets may have different probes associated. The probes in that case will be given a suffix (a, b, c, etc), and shown as one primer set
with two alternative probes. A given primer _sequence_ may have several primer names, as names are given based on primer pairs, e.g. primer_1_F and primer_2_F are identical, however primer_1_R and primer_2_R are not, therefore the forward primer sequence has two different names, as it is part of two different pairs.

If `--off-target-db` and/or `--on-target-db` are provided, primer prediction will be performed on the designed primers.

The primer-design workflow produces a TSV file containing the candidate primer sets and an interactive HTML report summarizing the design results.

The HTML report includes:
* Entropy plots for the evaluated consensus thresholds, showing sequence variability across the alignment together with potential primer-binding regions and the locations of the designed primers and probes.
* Primer tables summarizing the resulting primer pairs and, where applicable, probes. If primer prediction is performed using on-target and/or off-target databases, the corresponding predicted hits are included in these tables.
* Run status for each combination of consensus threshold and maximum ambiguity setting, indicating whether primer design was successful and whether qPCR or single mode was used.

Together, these outputs allow results from the different varVAMP parameter combinations to be reviewed and compared in a single report.

### <a name="pred_"></a> Primer prediction

**Off-target prediction**<br>
Provide an off-target database using `--off-target-db` to identify sequences that may be unintentionally amplified by a candidate primer set. BLAST hits are filtered according to primer-binding criteria, including:
* Matching of the first three bases at the 3′ end;
* A minimum alignment length of 50% of the primer;
* A mismatch allowance of one mismatch per five bases outside the first three bases at the 3′ end.
* The predicted amplicon length (_X_) should fall within an acceptable range (0.5 × amplicon length ≤ _X_ ≤ max(1000 bp, 1.5 × amplicon length))

Primer pairs passing these criteria are evaluated according to their orientation and the distance between predicted binding sites to identify potential off-target amplicons. Probe hits, when present, are additionally required to occur within the predicted amplicon.

**On-target prediction**<br>
Provide an on-target database using `--target-db` to estimate how many target sequences may be amplified by each candidate primer set. Primer variants are queried against the database and filtered according to alignment length, mismatch allowance, primer orientation, and expected amplicon length. The maximum number of allowed mismatches can be adjusted using `--target-mismatch`.

* Matching of the first three bases at the 3′ end;
* A minimum alignment length of 75% of the primer;
* A mismatch allowance of _n_ mismatch per five bases outside the first three bases at the 3′ end.
* Directionality of the primers, ergo forward primers must be forward, etc.
* The predicted amplicon length falls within a dynamically calculated range of the expected product size

**Standalone primer prediction**<br>
Primer prediction can also be run independently of the primer-design workflow, allowing existing or published primer sets to be evaluated using the same prediction criteria.
The primer prediction can be run separately from the workflow above.<br>
It requires a TSV-file with `primer_name`, `seq`, `amp_name`,`oligo_type`, `amp_len` columns.
```
primer_check.py [-h] -f FILE [-p PREFIX] [-x OFF_TARGET_DB] [-y TARGET_DB] [--target-mismatch TARGET_MISMATCH] [-o OUT_DIR] [-z THREADS] [-c ON_TARGET_COLUMN] [-g OFF_TARGET_COLUMN]
                       [--allowance-3-prime ALLOWANCE_3_PRIME] [--blast-remote]

Primer check using BLASTn

options:
  -h, --help            show this help message and exit
  -f FILE, --file FILE  TSV file containing primers. Must have columns "[primer_name, seq, amp_name, amp_len, oligo_type (left, right, probe)]"
  -p PREFIX, --prefix PREFIX
                        Output prefix. Default uses tsv file name
  -x OFF_TARGET_DB, --off-target-db OFF_TARGET_DB
                        Path to off-target blast database. If provided, primers will be checked for potential off-targets
  -y TARGET_DB, --target-db TARGET_DB
                        Path to target blast database. If provided, primers will be checked whether they can potentially catch all targets
  --target-mismatch TARGET_MISMATCH
                        Max number of mismatches when comparing to target blast db (default: 3).
  -o OUT_DIR, --out-dir OUT_DIR
                        Output directory. Default is <current directory>/primer-check_<prefix>_<date>
  -z THREADS, --threads THREADS
                        Number of CPU threads to use. (default: 10)
  -c ON_TARGET_COLUMN, --on-target-column ON_TARGET_COLUMN
                        Prefix for on-target column names. (default: on_target)
  -g OFF_TARGET_COLUMN, --off-target-column OFF_TARGET_COLUMN
                        Prefix for off-target column names. (default: off_target)
  --allowance-3-prime ALLOWANCE_3_PRIME
                        Number of mismatches allowed in the first 3 bases of the 3' end for off-target (default: 0)

```


### <a name="blast_"></a> Local BLAST databases
To run the workflow with `--target-db` and `--off-target-db` it is necessary to create local blast databases.<br> For the **off-target database** it is necessary to also provide taxids for the accessions, at minimum at species level.

If `[ERROR] BLAST Database error: Database memory map file error` assign more memory to the instance.

#### Off-target database

To create the taxonomy file, see [DB taxonomy file](#db_taxonomy_).

Create a blast database with taxids:
```
makeblastdb -dbtype nucl -hash_index -parse_seqids -in <FASTA> -out <DB-NAME> -taxid_map <TAXONOMY.TSV>
```

The taxonomy file must consist of accession and taxid:
```
KY325469.1      3048459
KY586889.1      3052464
MT261963.1      3052464
PQ155013.1      3052464
MN244562.1      3052464
```

#### On-target database

The on-target database does not need taxids, however, it is no issue if they are present.
```
makeblastdb -dbtype nucl -hash_index -in <FASTA> -out <DB-NAME> 
```

### <a name="db_taxonomy_"></a> Database taxonomy file

To create the database taxonomy you must have a meta file in tsv format, which contains columns `accession` and `species` and a fasta file containing the same accessions as in the meta file:
```
db_taxonomy -f <fasta_file> -m <meta-file.tsv>
```

# <a name="cite_"></a> Citation
If you use the Semiautomatic pan-primer design workflow in your research, please cite the associated manuscript. Citation details and DOI will be added upon publication.

The software is archived on Zenodo:<br>
**Semiautomatic pan-primer design, v1.0.0.**;  https://doi.org/10.5281/zenodo.22110655

The Semiautomatic pan-primer design workflow builds upon several existing tools. Please also cite the software used in your analysis, including:

* varVAMP – primer design: Fuchs et al., 2025; https://doi.org/10.1038/s41467-025-60175-9
* Vclust – sequence clustering and tANI calculation: Zielezinski et al., 2025; https://doi.org/10.1038/s41592-025-02701-7
* VSEARCH – sequence clustering: Rognes et al., 2016; https://doi.org/10.1093/bioinformatics/btl158
* CD-HIT – sequence clustering: Li and Godzik, 2006; https://doi.org/10.1093/bioinformatics/btl158
* MAFFT – multiple sequence alignment: Katoh and Standley, 2013; https://doi.org/10.1093/molbev/mst010
* Mash – genomic distance estimation: Ondov et al., 2016; https://doi.org/10.1186/s13059-016-0997-x

Additional software dependencies and their respective references are described in the associated manuscript.

# <a name="disc_"></a> Disclaimer
The Semiautomatic pan-primer design workflow generates and evaluates candidate primer designs for downstream assay development. _In silico_ predictions do not replace experimental validation and assay optimization.



