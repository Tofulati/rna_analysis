# RNA Modification Pipeline

[![Web Application](https://img.shields.io/badge/Website-Live_Site-success)](https://rna-analysis.vercel.app/)
[![GitHub](https://img.shields.io/badge/GitHub-Repository-blue)](https://github.com/Tofulati/rna_analysis)
[![Jupyter Notebook](finalized.ipynb)]

Welcome to the **RNA Modification Pipeline**! This project processed and visualized complex genetic data—specifically, chemical modifications made to RNA molecules—originating from Oxford Nanopore sequencing technology. 

---

## Project Overview & Biological Background

**What are RNA Modifications?**
Just like DNA, RNA is made up of molecular "letters" (A, C, G, U). Sometimes, the body chemically alters these letters after they are created, which can change how genes function. Two critical modifications we tracked in this project were:
* **A-to-I (Adenosine to Inosine):** Often plays a role in how the immune system and brain function.
* **m6A (Adenosine methylation):** One of the most common modifications, influencing how RNA is processed and degraded.

**The Samples**
We analyzed two specific cellular samples for this project:
* **MR01-1:** T-Cell Leukemia Samples
* **MR01-2:** iPS Generated T-Cell Samples

---

## Tech Stack

This project was divided into a **Data Processing Pipeline** (run on a high-performance computing cluster) and a **Web Application**.

**Data Processing:**
* **Python** (Pandas, NumPy) - For data manipulation and calculations.
* **modkit & samtools** - Industry-standard tools for reading raw genomic sequencing files.
* **SLURM** - Job scheduler used to manage the computing power required to process thousands of genes.

**Web Dashboard:**
* **React** - Interactive user interface.
* **Supabase** - Cloud database storing our processed genetic data.
* **Vercel** - Web hosting platform.

---

## Setup & Installation

If you would like to run the data processing pipeline locally, you will need to set up a Python virtual environment and install the required tools.

**1. Create and Activate a Virtual Environment:**
```bash
# Create the virtual environment
python3 -m venv rnaenv

# Activate the virtual environment (Linux/macOS)
source rnaenv/bin/activate

# Activate the virtual environment (Windows)
rnaenv\Scripts\activate

```

**2. Install Python Dependencies:**

```bash
pip install pandas numpy intervaltree supabase

```

**3. Install Genomic Tools:**
You will also need `samtools` and `modkit` installed on your system. These are best installed via `conda` or downloaded directly from their respective repositories.

> **Note for SDSC Expanse Users:** If you are running this pipeline on the Expanse cluster (or a similar SLURM environment), you can use the provided `.sh` and python submission scripts (e.g., `submit_extract.py`, `submit_output.py`) to orchestrate the jobs across multiple nodes.

---

## Code Architecture & Pipeline Mechanics

The backend of this project is a high-throughput computational pipeline written primarily in Python. It transforms raw, compressed sequencing data into aggregated, web-ready statistics. Here is a step-by-step breakdown of how the codebase processed the data.

### Phase 1: Data Ingestion & Conversion

Nanopore sequencing machines output data in `.bam` (Binary Alignment Map) files. These files are highly compressed and contain the millions of short RNA reads aligned to the human genome.

* **The Code:** The pipeline uses the `pysam` library to open and read these `.bam` files.
* **The Process:** Because genomic files are too massive to analyze all at once, the code splits the data chromosome by chromosome. It extracts only the essential information—such as the genetic sequence and the modification tags—and saves them into Pandas DataFrames (stored as `.pkl` or "pickle" files). This converts the data from a raw bioinformatics format into a standard, tabular format that is much faster to query.

### Phase 2: Gene Identification

Before looking for modifications, the pipeline needs to know exactly where to look.

* **The Code:** Using `pymysql`, the script connects directly to the UCSC Genome Browser's public database.
* **The Process:** It fetches the exact start and end coordinates for all human genes (using the standard hg38 genome build). The pipeline then cross-references these locations with our sample data, filtering out genes that lack enough sequencing coverage (requiring at least 10 reads in both samples). Genes that passed this test are categorized as "Primary" targets for deep analysis.

### Phase 3: Decoding the Modifications (The Core Engine)

This is the most computationally intensive part of the pipeline, where the code determines if a specific RNA "letter" was modified.

* **The Code:** This step relies heavily on string parsing and regular expressions (`re`).
* **The Process:** 1. **Alignment:** The code reads the "CIGAR string"—a compact code provided by the sequencer that explains exactly how a specific RNA read matches up against the standard human genome (e.g., noting if a letter was deleted or inserted).
2. **Probability Extraction:** It reads the `MM` (Modification) and `ML` (Probability) tags. These tags tell us where the sequencer *thinks* a modification occurred and how confident it is (on a scale of 0 to 255).
3. **Classification:** The script calculates the probabilities for the sequence and assigns the most likely state to every single position: Unmodified, A-to-I (Inosine), or m6A.

### Phase 4: Functional Annotation & Aggregation

Once the modifications are pinpointed, the code must figure out what part of the gene they belong to.

* **The Code:** The pipeline uses a data structure called an `IntervalTree`.
* **The Process:** An Interval Tree is highly efficient at answering the question, "Does position X fall within region Y?" The code maps every discovered modification to specific transcript features defined by a BED file: Exons, Introns, 5' UTRs, and 3' UTRs. Finally, the script aggregates these millions of mapped points into clean summary statistics, calculating total counts, averages, standard deviations, and Counts Per Kilobase (CPK) to normalize the data across genes of different lengths.

### Phase 5: Cloud Database Integration

The final step bridges the gap between the heavy bioinformatics processing and the interactive web frontend.

* **The Code:** A dedicated upload script uses the `supabase` Python client.
* **The Process:** It reads the final aggregated data tables, cleans up any missing values or incompatible data types (ensuring they are JSON-safe), and uploads them in batches of 100 to a Supabase PostgreSQL database. Once the data hits the database, it immediately becomes searchable on the live React website.

---

## Using the Web Dashboard

You can explore the processed data interactively at our [Live Webpage](https://rna-analysis.vercel.app/).

**Instructions:**

1. **Filter your view:** At the top of the page, select your sample (`MR01-1` or `MR01-2`), the region of interest (e.g., `UTR 3'`), and the Modification Type (`A-to-I`, `m6A`, or `Both`).
2. **Explore the Graph:** A scatter plot will populate showing genes that match your criteria. Hover over the dots to see quick details.
3. **Deep Dive into a Gene:** Scroll down to the data table or use the search bar to find a specific gene (e.g., `AZIN1`).
4. **Compare Data:** Clicking on a gene reveals a detailed breakdown of its modification probabilities, averages, and counts. You can swap samples at the top of the page to compare how the modifications changed between the samples.

---

## Contributors & Acknowledgements

This project was developed within the **Hojun Li Lab**.

* **Albert Ho** - Main Developer
* **Hojun Li** - Principal Investigator (PI)
* **Fay Jiang** - Collaborator

---

## References

* [Oxford Nanopore Modkit Documentation](https://nanoporetech.github.io/modkit/quick_start.html)
* [SAM File Specifications](https://samtools.github.io/hts-specs/SAMv1.pdf)
* [UCSC Genome Browser](https://genome.ucsc.edu/)