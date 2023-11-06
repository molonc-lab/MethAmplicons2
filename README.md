# MethAmplicons2
Command line tool written in Python for generation of lollipop and ridgleline plots from targeted bisulfite sequencing.

Generate lollipop plots that show the most frequent epialleles for a region of interest as well as a "combined" epiallele that is an average of less frequent (<5%) epialleles e.g. the RAD51C promoter region: 

<img width="515" alt="34059-Tumour_S9_L001__001_RAD51C (RAD51C)_5perc_barplot" src="https://github.com/molonc-lab/MethAmplicons2/assets/128128145/bbf1f840-7d72-4c71-9014-538330ca23c9">

Can also plot the combined/average epialleles for multiple samples of interest together:
<img width="390" alt="RAD51C methylation combined" src="https://github.com/molonc-lab/MethAmplicons2/assets/128128145/f8ab39a3-0a96-4b2a-9797-0f8a1279e664">

Lastly, you can plot ridgeline plots that allow one to view the distribution of the number of methylated CpGs per epiallele by sample. This provides another visual tool to assess the proportion of non-methylated, partially methylated, and fully methylated epialleles in a sample:  
<img width="1346" alt="ridgeline_combined_homo" src="https://github.com/molonc-lab/MethAmplicons2/assets/128128145/06a1efd0-4077-46a3-84ff-35efebe99279">

This repo provides a generalised "CLI tool" version of the code from [MethAmplicons](https://github.com/okon/MethAmplicons) which is a collection of scripts used for the analysis of RAD51C methylation data in:
[Nesic, K, Kondrashova, O et al. "Acquired *RAD51C* promoter methylation loss causes PARP inhibitor resistance in high grade serous ovarian carcinoma." Cancer Research (2021)](https://cancerres.aacrjournals.org/content/early/2021/07/27/0008-5472.CAN-21-0774)


This tool uses [FLASH](https://ccb.jhu.edu/software/FLASH/) paired-end read merging software for merging reads: [FLASH: Fast length adjustment of short reads to improve genome assemblies. T. Magoc and S. Salzberg. Bioinformatics 27:21 (2011), 2957-63](https://doi.org/10.1093/bioinformatics/btr507)


# Getting Started
- To get started with the tool, follow the steps INSTALLATION and USE below.  

## INSTALLATION: 

Currently, to install and use methamplicons: 
1. (Recommended) Create a virtual environment for installation 
2. Clone or download the MethAmplicons2 repository (put the files in your working directory) then install the methamplicons tool using flit.
3. Run methamplicons and approve moving the flash read merging tool to your system's bin folder or the bin folder of your virtual environment. 

### Step 1. Create and activate a virtual environment using a Python version between 3.6.8 and 3.9.6 (recommended):
- It is recommended to create a virtual environment where the tool can be installed
  ```bash
  #Command for creating virtual environment named 'ma_env'
  python3 -m venv ma_env
  #Command for activating virtual environment 
  source ma_env/bin/activate
  ``` 

### Step 2.1 - Getting the files from the GitHub repo:
- Clone the methamplicons repository in the directory where you want the code files to go (alternatively download the repo folder and move it to this directory): 

  ```bash
  # use ssh link if on an HPC
  git clone https://github.com/molonc-lab/methamplicons.git
  
  # cd into the MethAmplicons2 directory
  cd MethAmplicons2
  ``` 

### Step 2.2 - Install the methamplicons tool using flit: 

  ```bash
  # Install flit if you have not already 
  pip install flit
  # FOR VIRTUAL ENVIRONMENT: in the repo folder (MethAmplicons2) install methamplicons
  # and all its requirements with flit by specifying the path to the virtual environment's python 
  flit install --python ma_env/bin/python

  # ALTERNATIVELY: if you are using a non-organisation computer and have decided not to use a virtual environment
  flit install --symlink
  ``` 

## USE 

### Step 1: Activate the virtual environment from before if not activated already (using the path to the activate file)
```bash
  source ma_env/bin/activate
```

### Step 2: Run methamplicons, either by specifying absolute paths to files or giving relative paths 

#### Example command
```bash
#specify the location of methamplicons - not required if installed to usr/bin of Mac (non-organisation computers)
  ma_env/bin/methamplicons --PE_read_dir test \
    --amplicon_info test/BS_primers_amplicons_CDS_RC.tsv \
    --sample_labels test/SampleID_labels_amplicon_meth.csv \
    --output_dir output \
    --min_seq_freq 0.01 \
    --verbose true

```

### Step 2. (First time use) Approve moving flash to your system or virtual environment's bin: 
- methamplicons WILL NOT run without a binary for FLASH in the bin folder containing the version of Python being used to run the command. 
For Macs:
- When running methamplicons for the first time, respond 'y' when prompted with "Flash binary not found. Would you like to move a copy of the flash binary for Unix to your bin (y/n)?". A precompiled binary of Flash for Macs will be moved from src/methamplicons to the usr/bin folder or the bin of your active virtual environment (e.g. /Users/brettl/ma_env/bin/). If you happen to already have a binary for Flash in either of these folders you will not be prompted. If you respond with 'n' the program will exit.

Credit to the creators of FLASH: 
FLASH: Fast length adjustment of short reads to improve genome assemblies. T. Magoc and S. Salzberg. Bioinformatics 27:21 (2011), 2957-63.

#### Requirements for directories and files provided as arguments: 
- Example tsv and csv files are provided under tests

##### --PE_read_dir - directory containing paired end read files:
- An assumption of the program is that the last instance of R1 or R2 before the file extension (.fastq, .fastq.gz) indicates that a file contains read 1s of a pair or read 2s of a pair. 
- The tested files had read 1s at the same line (position) as the read 2s in the other file, however order shouldn't be important as each fastq files reads are placed in dictionaries and so a read's counterpart can be searched. 

##### --amplicon_info - tsv file containing the information about the amplified regions: 
- The tab-delimited (tsv) file should have data organised into 'columns':
    - AmpliconName, Primer1, Primer2, Sequence, and CDS
      
- Columns should contain: 
    - AmpliconName is the name given to the amplicon.
    - Primer1 and Primer2 are the primers that will match the reads (they ARE bisulfite converted). Primer2 is assumed to be the reverse primer and therefore its reverse complement is used by the program (the user can provide the reverse primer sequence as is).
    - Sequence - reference sequence for the amplified region that is NOT bisulfite converted. The reference sequence should start with primer 1's sequence and ends with primer 2's (reverse complement) sequence.
    - CDS is the distance of the first base in the reference sequence relative to the first base in the CDS. For genes on the reverse strand the direction of the amplicon and primer 1 and 2 need to be considered. 0 may be put as a stand in value. 

- If multiple regions are targeted in a single analysis, multiple amplicon entries can be provided (including overlapping regions/amplicons). These will be extracted from the reads and analysed separately.
  
Example tsv file:  
| Amplicon_Name |	Primer1  | 	Primer2  | 	Sequence | CDS |
| ------------- | ------------- | ------------- | ------------- | ------------- |
| RAD51C	| GAAAATTTATAAGATTGCGTAAAGTTGTAAGG |	CTAACCCCGAAACAACCAAACTCC | GAAAATTTACAAGACTGCGCAAAGCTGCAAGGCCCGGAGCCCCGTGCGGCCAGGCCGCAGAGCCGGCCCCTTCCGCTTTACGTCTGACGTCACGCCGCACGCCCCAGCGAGGGCGTGCGGAGTTTGGCTGCTCCGGGGTTAG	| -156 |
| BRCA1_l	| TTGTTGTTTAGCGGTAGTTTTTTGGTT	| AACCTATCCCCCGTCCAAAAA |	CTGCTGCTTAGCGGTAGCCCCTTGGTTTCCGTGGCAACGGAAAAGCGCGGGAATTACAGATAAATTAAAACTGCGACTGCGCGGCGTGAGCTCGCTGAGACTTCCTGGACGGGGGACAGGCT |	-1361 |
| BRCA1_s	| TTGTTGTTTAGCGGTAGTTTTTTGGTT	| CAATCGCAATTTTAATTTATCTATAATTCCC |	CTGCTGCTTAGCGGTAGCCCCTTGGTTTCCGTGGCAACGGAAAAGCGCGGGAATTACAGATAAATTAAAACTGCGACTG	| -1361 |

##### --sample_labels - csv file containing sample label information (optional):
- This file is not required, however it can be used to map the Sample Id (name used in the fastq files) to the SampleLabel or ShortLabel if the CSV includes the following columns:
    - SampleID, SampleLabel, ShortLabel

#### Example output files and directories: 
5_perc_sample1 (BRCA1_l) <br />
5_perc_sample1 (BRCA1_l)_no_legend <br />
5_perc_sample1 (BRCA1_s) <br />
5_perc_sample1 (BRCA1_s)_no_legend <br />
5_perc_sample2 (BRCA1_s) <br />
5_perc_sample2 (BRCA1_s)_no_legend <br />
5_perc_sample2 (RAD51C) <br />
5_perc_sample2 (RAD51C)_no_legend <br />
All_samples_combined_colour_meth_BRCA1_l.pdf <br />
All_samples_combined_colour_meth_BRCA1_s.pdf <br />
All_samples_combined_colour_meth_RAD51C.pdf <br />
All_samples_combined_colour_unmeth_BRCA1_l.pdf <br />
All_samples_combined_colour_unmeth_BRCA1_s.pdf <br />
All_samples_combined_colour_unmeth_RAD51C.pdf <br />
demultiplexed <br />
df_alleles_sort_all.csv <br />
df_exclude_unmeth-alleles_freq.csv <br />
df_meth_freq.csv <br />
merged <br />
ridgeline_plot_BRCA1_l.pdf <br />
ridgeline_plot_BRCA1_s.pdf <br />
ridgeline_plot_RAD51C.pdf <br />

## Argument info 

```
usage: methamplicons [-h] [--PE_read_dir PE_READ_DIR]
                     [--amplicon_info AMPLICON_INFO]
                     [--sample_labels SAMPLE_LABELS] [--output_dir OUTPUT_DIR]
                     [--min_seq_freq MIN_SEQ_FREQ] [--verbose {true,false}]
                     [--save_data {true,false}] [--lollipop {true,false}]
                     [--ridge {true,false}]

CLI tool for plotting targeted bisulfite sequencing

optional arguments:
  -h, --help            show this help message and exit
  --PE_read_dir PE_READ_DIR
                        Desired input directory with fastq files, gzipped or
                        not
  --amplicon_info AMPLICON_INFO
                        Path to the amplicon info file in tsv format, e.g.:
                        AmpliconName Primer1 Primer2 ReferenceSequence
  --sample_labels SAMPLE_LABELS
                        Path to sample labels file in csv format
  --output_dir OUTPUT_DIR
                        Desired output directory
  --min_seq_freq MIN_SEQ_FREQ
                        Threshold frequency an extracted epiallele sequence
                        must have to be included in analysis
  --verbose {true,false}
                        Print all output after file parsing (default: true).
  --save_data {true,false}
                        Save processed data in csv format (default: true).
  --lollipop {true,false}
                        Save a lollipop graph (default: true).
  --ridge {true,false}  Save a ridgeline graph (default: true).
```
## Alternative OS support:
Linux:
- Take the flash binary from linux_flash and move it to the appropriate bin on your system or the src/methamplicons folder before running the program (and then respond 'y' as above).
- If the attached binary does not work try other binary files from: https://ccb.jhu.edu/software/FLASH/
- For a Linux system, extract the binary from the FLASH-1.2.11-Linux-x86_64.tar.gz

Windows:
- methamplicons attempts to move the flash binary to the bin folder which does not exist on Windows, however such a bin would exist in a Python or Anaconda environment: a virtual environment is therefore required to run methamplicons on Windows. 

