<h1 align="center">Niangao: Sanger Sequencing Results Processing Tool</h1>
<p align="center">
  <strong>🔬 A Fully Automated Tool for Microbial Identification Using Sanger Sequencing</strong>
</p>

## Introduction

> If you work with Sanger sequencing data and need a fast, efficient way to process `.ab1` files, **Niangao** is the perfect tool for you! 🚀 

Niangao provides a streamlined workflow for:

✅ **Sequence trimming**  
✅ **Paired-end assembly (by CAP3)**  
✅ **BLAST alignment (against NCBI Databases)** 

Use it as a **Web tool**, an **R package**, or a **Docker container** based on your needs! 👇  
## Platforms
### Web
An [online version](https://niangao.lab.westlake.edu.cn) of Niangao is available for non-coders. Please check the [tips and features](#using-niangao-web---tips-and-features) for Niangao before using.
### R package
If you are familiar with R and want to build your own pipeline with R code, please try our [R package version](https://github.com/WeStrainGroup/Niangao_R_package) of Niangao. **Note** that this version ***ONLY*** have the trimming functionality! See [using R package](#using-r-package) for detailed code.
### Docker image
If you need to build Niangao at your local server for data security and faster uploads, we recommend you to use Dockerized Niangao, which have exactly the same functions and interface as Niangao Web. See [using Docker image](#using-docker-image) for the detailed code to build it.
## User manual
### Using Niangao Web - Tips and Features
- **BLAST Options:** Select "YES" to run BLAST. Choose your desired database and the number of aligned subject sequence to display.
- **File Upload:** You can upload `.ab1` files by dragging and dropping them onto the file input area or by clicking the "Browse..." button.
- **Paired-End Assembly:** For paired-end sequence assembly, name your forward and reverse read files with a consistent pattern: `SampleName_Fwd.ab1` and `SampleName_Rev.ab1`.  Only files with matching names and the correct suffixes (`_Fwd`, `_Rev`) will be assembled.
- **Failed Assemblies:** If CAP3 fails to assemble a pair of sequences (often due to low quality), the longer of the two sequences will be included in the output.
- **Downloading Results:**  Download your processed sequences and reports using the download button located next to the file input. **Important:** Results are ***NOT*** saved after you refresh or close the page.
### Using R package 
#### Step 1: Install tool package for Niangao installation
Install R package '*devtools*' and '*BiocManager*' on your device if you haven't installed yet:
``` R
if (!requireNamespace("devtools", quietly = TRUE)) {
  install.packages("devtools")
}

if (!requireNamespace("BiocManager", quietly = TRUE)) {
  install.packages("BiocManager")
}

```
#### Step 2: Install Niangao package
``` R
devtools::install_github("WeStrainGroup/Niangao_R_package")
```
- The relevant dependency package will be automatically installed. 
- If you failed at this step, please check if it is a network error and then try again.
#### Step 3: Load the package
``` R
library(Niangao)
``` 
#### Step 4: Call the main function
``` R
results <- niangao("path/to/your/ab1/files") # Specify the directory containing your .ab1 files
```
#### Step 5: Access the results
``` R
print(results$result_table)  # View the summary table
print(results$dna_string_set)  # View the DNAStringSet object
	
# You can then create your own pipeline with the results of Niangao!
```
### Using Docker Image
#### Prerequisites
Before you begin, ensure you have the following:
*   **Docker:**
    *   Windows/macOS: [Docker Desktop](https://www.docker.com/products/docker-desktop)
    *   Linux: [Docker Engine](https://docs.docker.com/engine/install/)
*   **Hardware Resources:** 
	-   CPU: 4+ cores if you need to run BLAST for hundreds of files
#### Step 1: Pull the image from Docker Hub
```bash
docker pull westraingroup/niangao
```
#### Step 2: Run the Docker container
```bash
docker run -d -p 3838:3838 westraingroup/niangao
 ```
*   `-d`: Runs the container in detached mode (in the background).
*   `-p 3838:3838`: Maps port 3838 on your host machine to port 3838 inside the container.  If port 3838 is already in use, change the *first* 3838 to a different available port (e.g., `-p 8080:3838`).
#### Step 3: Access the application
Open your web browser and go to:  `http://localhost:3838` (or `http://localhost:8080` if you used a different host port in the previous step).
#### If you need
If you what to check the container id to stop it (which is **highly recommended** after you finished your work) or check detailed running log for troubleshooting, you should use these code:
1.  **Find the Container ID:**
```bash
docker ps
```
This command lists all running containers.  Look for the container with the image name `niangao:latest` and note its `CONTAINER ID`.

2.  **Check the running log (optional, useful for trouble shooting)**
```bash
docker logs <container_id>
```
Replace `<container_id>` with the actual container ID from the previous step.

3.  **Stop the Container:**
```bash
docker stop <container_id>
```
Replace `<container_id>` with the actual container ID from the previous step.
## Citation
This project utilizes several R packages and external tools.
**R Packages:** [shiny](https://github.com/rstudio/shiny), [shinyjs](https://github.com/daattali/shinyjs), [DT](https://github.com/rstudio/DT), [zip](https://github.com/r-lib/zip), [sangerseqR](https://github.com/jonathonthill/sangerseqR),  [RcppRoll](https://github.com/kevinushey/RcppRoll), [seqinr](https://github.com/lbbe-software/seqinr), [sys](https://github.com/jeroen/sys).

**External Tools:** [CAP3](https://faculty.sites.iastate.edu/xqhuang/cap3-and-pcap-sequence-and-genome-assembly-programs), [BLAST](https://blast.ncbi.nlm.nih.gov/Blast.cgi), [Taxonkit](https://github.com/shenwei356/taxonkit/).

**Databases:**
*   [ITS_RefSeq_Fungi](https://ftp.ncbi.nlm.nih.gov/blast/db/ITS_RefSeq_Fungi.tar.gz)
*   [ITS_eukaryote_sequences](https://ftp.ncbi.nlm.nih.gov/blast/db/ITS_eukaryote_sequences.tar.gz)
*   [16S_ribosomal_RNA](https://ftp.ncbi.nlm.nih.gov/blast/db/16S_ribosomal_RNA.tar.gz)
*   [28S_fungal_sequences](https://ftp.ncbi.nlm.nih.gov/blast/db/28S_fungal_sequences.tar.gz)
*   [18S_fungal_sequences](https://ftp.ncbi.nlm.nih.gov/blast/db/18S_fungal_sequences.tar.gz)
## License
This project is licensed under the MIT License - see the [LICENSE](LICENSE) file for details.
## Citation
If you use Niangao in your research or work, please cite this tool - see the [CITATION](CITATION) file for details, or click 'Cite this repository' on this page.
## Contact us
For questions, suggestions, or bug reports, please contact us via email at [zhouwenhao@westlake.edu.cn](mailto:zhouwenhao@westlake.edu.cn) and [wangxinyu30@westlake.edu.cn](mailto:wangxinyu30@westlake.edu.cn).
We greatly appreciate your contributions and feedback! ❤