# Sanger Sequencing Data Analysis Results

This directory contains the results of processing Sanger sequencing data (.ab1 files) using Niangao.

## Overview

The application automates the process of cleaning, assembling (when paired-end reads are provided), and optionally performing BLAST searches on Sanger sequencing data. This README provides a description of the files produced by the analysis pipeline.

## File List and Descriptions

### Core Output Files

*   **`seq_trimed_all.fasta`**: All quality-controlled (trimmed) sequences from the input .ab1 files, stored in FASTA format. This file contains all the sequences *before* any assembly step.
*   **`trimming_report.csv`**: A report detailing the trimming process for each .ab1 file. See the [Trimming Report Details](#trimming-report-details) section for column descriptions.
*   **`seq_trimed_single.zip`**: A compressed archive (.zip) containing individual FASTA files for each of the trimmed sequences (one file per .ab1 input). These are the same sequences as in `seq_trimmed_all.fasta`, but separated into individual files.

### Sequence Assembly Results (Paired-End Reads Only)

*These files are generated only if the input data includes paired-end reads (identified by `_Fwd` and `_Rev` suffixes in the filename).*

*   **`seq_trimed_assembled_all.fasta`**: This FASTA file contains the final sequences after the assembly process. It includes:
    *   Successfully assembled contigs (merged sequences from paired reads). These are named with the `_assemble` suffix.
    *   For pairs where assembly *failed*, the *longer* of the trimmed forward (Fwd) or reverse (Rev) sequence is included. These are named with `_Fwd` or `_Rev` suffix.
    *   Sequences from .ab1 files that had *only* a forward or *only* a reverse read (but not both). These are named with the original sample name and the `_Fwd` or `_Rev` suffix.
    *   Sequences from .ab1 files that did *not* have `_Fwd` or `_Rev` suffixes, *and* were not part of any paired-end set. These are named by the original sample.

*   **`cap3_aseembly_details.ace`**: The complete output file from the CAP3 assembly program, in ACE format. This file contains detailed information about how the contigs were assembled, including alignments and quality scores.

*   **`seq_trimed_assembled_single.zip`**: A compressed archive (.zip) containing individual FASTA files for each of the assembled sequences *and* the sequences that failed to assemble (as described for `seq_trimed_assembled_all.fasta`).

### BLAST Alignment Results (Optional)

*These files are generated only if the "Do BLAST?" option is set to "YES".*

*   **`blast_result_final.tsv`**: A tab-separated value (TSV) file containing the results of the BLAST search against the selected database. See the [BLAST Results Details](#blast-results-details) section for column descriptions.

## Trimming Report Details

The `trimming_report.csv` file provides a summary of the sequence trimming process for each input .ab1 file.

| Column Name             | Description                                                                                         |
| :---------------------- | :-------------------------------------------------------------------------------------------------- |
| `file_name`             | The base name of the original .ab1 file (without the .ab1 extension).                                |
| `len._raw`              | The length (in bases) of the raw sequence extracted from the .ab1 file.                              |
| `len._trimmed`           | The length (in bases) of the sequence *after* trimming.                                              |
| `ave._quality_raw`       | The average Phred quality score of the *raw* sequence.                                                 |
| `ave._quality_trimmed`   | The average Phred quality score of the *trimmed* sequence.                                             |
| `len._deleted`           | The total number of bases removed during trimming (`len._raw` - `len._trimmed`).                  |
| `ratio_deleted`         | The proportion of bases removed during trimming (`len._deleted` / `len._raw`).                         |
| `deleted_count1`        | The number of bases removed during the first trimming step (sliding window average quality filtering). |
| `deleted_count2`        | The number of bases removed during the second trimming step (longest continuous high-quality region).  |
| `deleted_count3`        | The number of bases removed during the third trimming step (low-quality base filtering).             |

## BLAST Results Details

The `blast_result_final.tsv` file provides the results of the BLAST search, with lineage annotations.

| Column Name        | Description                                                                                                |
|--------------------|------------------------------------------------------------------------------------------------------------|
| `qseqid`           | Query sequence ID (matches sequence names in FASTA files).                                                 |
| `sseqid`           | Subject sequence ID (from the BLAST database).                                                               |
| `staxids`          | NCBI Taxonomy ID(s) of the subject sequence.                                                                   |
| `pident`           | Percentage of identical matches in the alignment.                                                             |
| `length`           | Length of the alignment.                                                                                       |
| `qcovs`            | Query coverage per subject (percentage of the query sequence that aligns to the subject).                   |
| `mismatch`         | Number of mismatches in the alignment.                                                                       |
| `gapopen`          | Number of gap openings in the alignment.                                                                     |
| `qstart`           | Start position of the alignment in the query sequence.                                                       |
| `qend`             | End position of the alignment in the query sequence.                                                         |
| `sstart`           | Start position of the alignment in the subject sequence.                                                       |
| `send`             | End position of the alignment in the subject sequence.                                                         |
| `sstrand`          | Strand of the subject sequence in the alignment (plus or minus).                                             |
| `evalue`           | Expect value (E-value) of the alignment. Lower values indicate higher significance.                          |
| `bitscore`         | Bit score of the alignment (measure of the quality of the alignment).                                       |
| `taxids`           | NCBI Taxonomy ID(s).                                                                                         |

## Authors

*   Wenhao Zhou ([zhouwenhao@westlake.edu.cn](mailto:zhouwenhao@westlake.edu.cn))
*   Xinyu Wang ([wangxinyu30@westlake.edu.cn](mailto:wangxinyu30@westlake.edu.cn))

## Citation list

This project utilizes several R packages and external tools.

**R Packages:**

*   **shiny:** [https://github.com/rstudio/shiny](https://github.com/rstudio/shiny).

*   **shinyjs:** [https://github.com/daattali/shinyjs](https://github.com/daattali/shinyjs).

*   **DT:** [https://github.com/rstudio/DT](https://github.com/rstudio/DT).

*   **zip:** [https://github.com/r-lib/zip](https://github.com/r-lib/zip).

*   **sangerseqR:** [https://github.com/jonathonthill/sangerseqR](https://github.com/jonathonthill/sangerseqR).

*   **RcppRoll:** [https://github.com/kevinushey/RcppRoll](https://github.com/kevinushey/RcppRoll).

*   **seqinr:** [https://github.com/lbbe-software/seqinr](https://github.com/lbbe-software/seqinr).

*   **sys:** [https://github.com/jeroen/sys](https://github.com/jeroen/sys).

**External Tools:**

*   **CAP3:** Huang X, Madan A. CAP3: A DNA Sequence Assembly Program. Genome Res. 1999;9(9):868-877. doi:[10.1101/gr.9.9.868](https://doi.org/10.1101/gr.9.9.868)
*   **BLAST:** Altschul SF, Gish W, Miller W, Myers EW, Lipman DJ. Basic local alignment search tool. J Mol Biol. 1990;215(3):403-410. doi:[10.1016/S0022-2836(05)80360-2](https://doi.org/10.1016/S0022-2836(05)80360-2)
*   **Taxonkit:** Shen W, Ren H. TaxonKit: A practical and efficient NCBI taxonomy toolkit. J Genet Genomics. 2021;48(9):844-850. doi:[10.1016/j.jgg.2021.03.006](https://doi.org/10.1016/j.jgg.2021.03.006)

**Databases:**

*   [https://ftp.ncbi.nlm.nih.gov/blast/db/ITS_RefSeq_Fungi.tar.gz](https://ftp.ncbi.nlm.nih.gov/blast/db/ITS_RefSeq_Fungi.tar.gz)
*   [https://ftp.ncbi.nlm.nih.gov/blast/db/ITS_eukaryote_sequences.tar.gz](https://ftp.ncbi.nlm.nih.gov/blast/db/ITS_eukaryote_sequences.tar.gz)
*   [https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz](https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz)
*   [https://ftp.ncbi.nlm.nih.gov/blast/db/28S_fungal_sequences.tar.gz](https://ftp.ncbi.nlm.nih.gov/blast/db/28S_fungal_sequences.tar.gz)
*   [https://ftp.ncbi.nlm.nih.gov/blast/db/18S_fungal_sequences.tar.gz](https://ftp.ncbi.nlm.nih.gov/blast/db/18S_fungal_sequences.tar.gz)