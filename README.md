# methamplicons
Command line tool written in Python for generation of lollipop and ridgleline plots from targeted bisulfite sequencing.
- To get started with the tool, follow the steps under INSTALLATION and USE below.  

## Example command
  ```bash
methamplicons --PE_read_dir test --amplicon_info test/BS_primers_amplicons_CDS_RC.tsv --sample_labels test/SampleID_labels_amplicon_meth.csv --output_dir output
  ``` 

## Requirements for directories and files provided as arguments: 
- Example tsv and csv files are provided under tests

### --PE_read_dir - directory containing paired end read files:
- An assumption of the program is that the last instance of R1 or R2 before the file extension (.fastq, .fastq.gz) indicates that a file contains read 1s of a pair or read 2s of a pair. 
- The tested files had read 1s at the same line (position) as the read 2s in the other file, however order shouldn't be important as each fastq files reads are placed in dictionaries and so a read's counterpart can be searched. 

### --amplicon_info - tsv file containing the information about the amplified regions: 
- The tsv file should have data organised into 'columns' AmpliconName, Primer1, Primer2, Sequence, and CDS where: 
    - Primer1 and Primer2 are the primers that will match the reads (they ARE bisulfite converted) and the reference sequence is NOT bisulfite converted.
        - Primer2 is assumed to be the reverse primer and therefore its reverse complement is used by the program (the user can provide the reverse primer sequence as is).
    - CDS is the distance of the first base in the reference sequence relative to the first base in the CDS (this can be found by using BLAT with the reference sequence in IGV). 0 may be put as a stand in value. 
    - Reference sequence for the amplified region: an assumption is that the reference sequence starts with primer 1's sequence and ends with primer 2's (reverse complement) sequence.

Example tsv file: 
Amplicon_Name 	Primer1   	Primer2   	Sequence	CDS
RAD51C	GAAAATTTATAAGATTGCGTAAAGTTGTAAGG	CTAACCCCGAAACAACCAAACTCC	GAAAATTTACAAGACTGCGCAAAGCTGCAAGGCCCGGAGCCCCGTGCGGCCAGGCCGCAGAGCCGGCCCCTTCCGCTTTACGTCTGACGTCACGCCGCACGCCCCAGCGAGGGCGTGCGGAGTTTGGCTGCTCCGGGGTTAG	-156
BRCA1_l	TTGTTGTTTAGCGGTAGTTTTTTGGTT	AACCTATCCCCCGTCCAAAAA	CTGCTGCTTAGCGGTAGCCCCTTGGTTTCCGTGGCAACGGAAAAGCGCGGGAATTACAGATAAATTAAAACTGCGACTGCGCGGCGTGAGCTCGCTGAGACTTCCTGGACGGGGGACAGGCT	-1361
BRCA1_s	TTGTTGTTTAGCGGTAGTTTTTTGGTT	CAATCGCAATTTTAATTTATCTATAATTCCC	CTGCTGCTTAGCGGTAGCCCCTTGGTTTCCGTGGCAACGGAAAAGCGCGGGAATTACAGATAAATTAAAACTGCGACTG	-1361

### csv file:
This file is not required, however it can be used to map the SampleId (which may be in a file name) to the SampleLabel or ShortLabel if the CSV includes the following columns: 
SampleID,SampleLabel,ShortLabel

# INSTALLATION: 

Currently, to install and use methamplicons:
1. (Recommended) Create a virtual environment for installation 
2. Download and extract the binary for the flash read merging tool and move it to 
the user's bin folder or the bin in the virtual environment. 
3A. Clone or download the MethAmplicons2 repository (put the files in your working directory) then install the tool using flit.
3B. Install the methamplicons tool from testpypi: pip install -i https://test.pypi.org/simple/ methamplicons

## Step 1. Create and activate a virtual environment (recommended):
- It is recommended to create a virtual environment where the tool can be installed
  ```bash
  #Command for creating virtual environment named 'ma_env'
  python3 -m venv ma_env
  #Command for activating virtual environment 
  source ma_env/bin/activate
  ``` 

## Step 2. (Optional) Download the flash binary for your system and move it to the appropriate bin folder: 
- If you already have a binary of Johns Hopkins University's FLASH in your bin or the bin of your active virtual environment, methamplicons will use this to run. 
- Otherwise: when running methamplicons for the first time, respond 'y' when prompted with "Flash binary not found. Would you like to move a copy of the flash binary for Unix to your bin (y/n)?" If you respond with 'n' the program will exit and you can either try again or manually obtain the binary. 
- methamplicons WILL NOT run without a binary for FLASH in the user's bin so if you prefer not to automatically do this: 
        - Download and extract the binary for FLASH yourself and move it to the appropriate bin: 
            - Get your system's flash binary file from: https://ccb.jhu.edu/software/FLASH/
            - For a Linux system, extract the binary from the FLASH-1.2.11-Linux-x86_64.tar.gz
 ```bash
 #extract the binary 
  tar -xzvf FLASH-1.2.11-Linux-x86_64.tar.gz

  ``` 

## Step 3 Option A: Installing from GitHub

### Step 3A.1 - Getting the files from the GitHub repo:
- Clone the methamplicons repository in the directory where you want the code files to go (alternatively download the repo folder and move it to this directory): 

  ```bash
  # use ssh link if on an HPC
  git clone https://github.com/molonc-lab/methamplicons.git
  
  # cd into the MethAmplicons2 directory
  cd MethAmplicons2

  ``` 

## Step 3A.2 - Install the methamplicons tool using flit: 
  ```bash
  # Install flit if you have not already 
  pip install flit
  # In the repo folder (MethAmplicons2) install methamplicons and all its requirements with:
  flit install --symlink
  ``` 

## Step 3 Option B: Installing from testpypi
- The command below should install methamplicons and all its requirements (but the version may be less recent than GitHub)
  ```bash
  pip3 install -i https://test.pypi.org/simple/methamplicons

  ``` 

# USE 

## Step 1: Activate the virtual environment from before if not activated already (using the path to the activate file)
source mba_env/bin/activate

## Step 2: Run methamplicons, either by specifying absolute paths to files or giving relative paths 
- Example with relative paths, using the aforementioned test folder in labdata/Lab_OlgaK/brettL/test: 
- The output folder does not need to already exist
methamplicons --PE_read_dir test --amplicon_info test/BS_primers_amplicons_CDS_RC.tsv --sample_labels test/SampleID_labels_amplicon_meth.csv --output_dir output

Example output files and directories: 
5_perc_183--T4_S11 (BRCA1_l) 
5_perc_183--T4_S11 (BRCA1_l)_no_legend 
5_perc_183--T4_S11 (BRCA1_s) 
5_perc_183--T4_S11 (BRCA1_s)_no_legend 
5_perc_183-R2-T1_S14 (BRCA1_s) 
5_perc_183-R2-T1_S14 (BRCA1_s)_no_legend 
5_perc_13375_S32 (RAD51C) 
5_perc_13375_S32 (RAD51C)_no_legend 
5_perc_NTC--T1_S3 (BRCA1_l) 
5_perc_NTC--T1_S3 (BRCA1_s) 
5_perc_NTC--T1_S3 (BRCA1_s)_no_legend 
5_perc_NTC--T1_S3 (RAD51C) 
5_perc_OV--T4_S10 (BRCA1_l) 
5_perc_OV--T4_S10 (BRCA1_l)_no_legend 
5_perc_OV--T4_S10 (BRCA1_s) 
5_perc_OV--T4_S10 (BRCA1_s)_no_legend 
5_perc_OV--T4_S10 (RAD51C) 
5_perc_OV-R2-T1_S13 (BRCA1_s) 
5_perc_OV-R2-T1_S13 (BRCA1_s)_no_legend 
All_samples_combined_colour_meth_BRCA1_l.pdf 
All_samples_combined_colour_meth_BRCA1_s.pdf 
All_samples_combined_colour_meth_RAD51C.pdf 
All_samples_combined_colour_unmeth_BRCA1_l.pdf 
All_samples_combined_colour_unmeth_BRCA1_s.pdf 
All_samples_combined_colour_unmeth_RAD51C.pdf 
demultiplexed 
df_alleles_sort_all.csv 
df_exclude_unmeth-alleles_freq.csv 
df_meth_freq.csv 
merged
ridgeline_plot_BRCA1_l.pdf 
ridgeline_plot_BRCA1_s.pdf 
ridgeline_plot_RAD51C.pdf

## Argument info 

usage: methamplicons [-h] [--PE_read_dir PE_READ_DIR]
                        [--amplicon_info AMPLICON_INFO]
                        [--sample_labels SAMPLE_LABELS]
                        [--output_dir OUTPUT_DIR] [--save_data {true,false}]
    

optional arguments:
  -h, --help            show this help message and exit
  --PE_read_dir PE_READ_DIR
                        Desired input directory with fastq files, gzipped or
                        not
  --amplicon_info AMPLICON_INFO
                        Path to the amplicon info file in tsv format, e.g.:
                        AmpliconName Primer1 Primer2 ReferenceSequence
  --sample_labels SAMPLE_LABELS
                        Path to sample labels file in csv format - not currently used by the tool
  --output_dir OUTPUT_DIR
                        Desired output directory 
  --save_data {true,false}
                        Save processed data in csv format (default: true).
