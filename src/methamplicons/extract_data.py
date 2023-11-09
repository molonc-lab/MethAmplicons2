import gzip
import subprocess
import os
import csv
from pkg_resources import resource_filename


class DataExtractionError(Exception):
    pass

class ExtractData:
    """
    This class contains attributes and methods for extracting data from 
    amplicon sequencing (not necessarily targeted bisulfite sequencing)
    """
    COMPLEMENT = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}

    def __init__(self):
        self.read1_seqs = {}
        self.read2_seqs = {}

    def open_file(self, sequences_file: str):
        # Open the file depending on whether it's gzipped or not
        if sequences_file.endswith(".gz"):
            file_handle = gzip.open(sequences_file, 'rt')  # Open gzipped file in text mode
        else:
            file_handle = open(sequences_file, 'r')
        return file_handle
    
    def valid_characters(self, sequence, seq_type): 
        if any(base not in "ACTGN" for base in sequence):
            raise DataExtractionError(f"Invalid characters in {seq_type} file")
        
    def reverse_complement(self, seq): 
        rev_comp = "".join(self.COMPLEMENT.get(base, base) for base in reversed(seq))
        return rev_comp
        
    def read_fastq(self, fastq_file):
        seqs = {}
        seq_info = None
        file_handle = self.open_file(fastq_file)

        with file_handle as rs_file:
            lines = rs_file.readlines()
            for i in range(0, len(lines), 4): 
                # Remove leading @ and whitespace
                seq_info = lines[i].strip()
                #THIS HANDLES A SPECIFIC PAIRED END READ DATA format
                seq_info = seq_info.split()[0]
                sequence = lines[i + 1].strip()
                phred = lines[i + 3].strip()
                ##print(f"Sequence info {seq_info}: sequence: {sequence}")
                self.valid_characters(sequence, "Sequences")
                seqs[seq_info] = [sequence, phred]

        return seqs
    
    def sort_fqs_by_primer(self, read1_seqs_file, read2_seqs_file, refseqs, amplicon_info, base_name, out_dir, runflash = True):
        """
        place the reads into as many pairs of files as there are unique sets of primers
        when a new pair of files is done being created, run flash on them 
        """

        read1_seqs, read2_seqs = self.load_paired_read_files(read1_seqs_file, read2_seqs_file)

        #create a directory "demultiplexed" in the specified output directory
        parsed_out_dir = os.path.join(out_dir, "demultiplexed")
        if not os.path.exists(parsed_out_dir):
            os.makedirs(parsed_out_dir)

        # will have read 1 and read 2 files for a given sample and will later append the gene name
        new_r1_base_name = os.path.join(parsed_out_dir, base_name + "_R1")
        new_r2_base_name = os.path.join(parsed_out_dir, base_name + "_R2")

        #print(f"The new base names from the demulitplexed function (point 1) are {new_r1_base_name} and {new_r2_base_name}")

        # get average read length
        if runflash: 
            read_lengths = [len(lst[0]) for lst in read1_seqs.values()]
            avg_read_len = sum(read_lengths) / len(read_lengths)

        #print(f"amplicon_info at sort_fqs: {amplicon_info}")

        for amplicon_name, primer_info in amplicon_info.items(): 
            r1s_for_region = new_r1_base_name + "_" + amplicon_name + ".fastq"
            r2s_for_region = new_r2_base_name + "_" + amplicon_name + ".fastq"

            #print(f"The new base names from the demulitplexed function (point 2) are {r1s_for_region} and {r2s_for_region}")

            new_fastq_r1 = open(r1s_for_region, "w")
            new_fastq_r2 = open(r2s_for_region, "w")
            new_fastq_r1.write("")
            new_fastq_r2.write("")
            new_fastq_r1.close()
            new_fastq_r2.close()
            
            new_fastq_r1 = open(r1s_for_region, "a")
            new_fastq_r2 = open(r2s_for_region, "a")
            
            #the amplicon_info dictionary contains the primers to search for in the PE read files
            #primers_to_search = primer_info[0:2]
            primers_to_search = [primer[:10] for primer in primer_info[0:2]]

            for read_id, read1_seq_qual in read1_seqs.items(): 
                read2_seq_qual = read2_seqs[read_id]
                if any(primer in read1_seq_qual[0] for primer in primers_to_search):
                    # this read corresponds to the region of interest
                    # put the read1 and read2 into the two region-specific fq files
                    ##print("match found")

                    r1_fq_entry = read_id + "1:N:0:32\n" + read1_seq_qual[0] + "\n+\n" + read1_seq_qual[1] + "\n"
                    r2_fq_entry = read_id + "2:N:0:32\n" + read2_seq_qual[0] + "\n+\n" + read2_seq_qual[1] + "\n"
                    ##print(r1_fq_entry)
                    new_fastq_r1.write(r1_fq_entry)
                    new_fastq_r2.write(r2_fq_entry)

            new_fastq_r1.close()
            new_fastq_r2.close()

            if runflash: 
                base_name_reg = base_name + "_" + amplicon_name
                #print(f"The base name for the merged file is {base_name_reg}")
                refseq_len = len(refseqs[amplicon_name])
                # pass required aruments to run_flash, converting required arguments to strings
                self.run_flash(r1s_for_region, r2s_for_region, base_name_reg, out_dir, str(int(avg_read_len)), str(refseq_len))
    def get_flash_binary_path(self):
        #may need to specify specific version of tool
        return resource_filename('methamplicons', 'flash')
    
    def set_verbose(self, verbose):
        self.verbose = verbose

    def run_flash(self, r1s_for_region, r2s_for_region, base_name_reg, out_dir, avg_read_len, refseq_len):

        output_dir = os.path.join(out_dir, "merged")
        # make output directory if it doesn't exist
        if not os.path.exists(output_dir):
            os.makedirs(output_dir)

        flash_binary = self.get_flash_binary_path()
        # ALL ARGUMENTS MUST BE STRINGS
        cmd = [flash_binary, "-m", "10", "-M", refseq_len, "-x", "0.25", "-O", "-r", avg_read_len, \
                "-f", refseq_len, r1s_for_region, r2s_for_region, "-d", output_dir, "-o", base_name_reg]
        
        #subprocess.run(cmd, universal_newlines=True, check=True)

        if (self.verbose == "true"):
            subprocess.run(cmd, universal_newlines=True, check=True)
        else:
            output_file = os.path.join(out_dir, "flash_stdout.txt")
            with open(output_file, "ab") as outf:
                subprocess.run(cmd, stdout=outf, stderr=subprocess.PIPE, universal_newlines=True, check=True)

    def merge_reads(self, r1_seqs_file, r2_seqs_file, refseqs, amplicon_info, output_dir):
        #print(f"The R1 and R2 files at merge_reads are {r1_seqs_file} and {r2_seqs_file}")
        """
        By default, place the reads in a directory called merged in the preferred output directory,
        which is the current working directory by default, with the base name for the 
        paired end read files as a subdirectory
        """
        base_name = ""
        #get the base name of the fastq files
        #can change this code to remove surrounding underscores/better stylise the sample name
        if r1_seqs_file.endswith("fastq.gz"):
            base_name = os.path.basename(r1_seqs_file).replace("R1", "").replace(".fastq.gz", "")
        elif r1_seqs_file.endswith("fq.gz"):
            base_name = os.path.basename(r1_seqs_file).replace("R1", "").replace(".fq.gz", "")
        elif r1_seqs_file.endswith(".fastq"):
            base_name = os.path.basename(r1_seqs_file).replace("R1", "").replace(".fastq", "")
        elif r1_seqs_file.endswith(".fq"):
            base_name = os.path.basename(r1_seqs_file).replace("R1", "").replace(".fq", "")

        #print(f"The file base name at merge_reads is {base_name}")
        # sort/demultiplex the fastq files
        # run run_flash in this function as all the required info is there
        # read1_seqs, read2_seqs, refseqs, primer info, base_name, out_dir,
        self.sort_fqs_by_primer(r1_seqs_file, r2_seqs_file, refseqs, amplicon_info, base_name, output_dir)

    def check_tsv_entry(self, entry_values): 

        is_valid_entry = True

        if any(value.isspace() for value in entry_values):
            is_valid_entry = False
        
        if any(len(value) == 0 for value in entry_values): 
            is_valid_entry = False

        if not sum(all(base in 'ACTG' for base in seq) for seq in entry_values[1:]):
            is_valid_entry = False

        # THIS data validation statement needs to be fixed
        #if not entry_values[-1].strip("-").isdigit(): 
            #is_valid_entry = False

        return is_valid_entry
    
    def parse_tsv(self, amplicon_info_tsv):

        primer_dict = {}
        refseqs = {} 

        for line in amplicon_info_tsv:
            #for a tsv file
            line_split = line.split("\t")

            if line.startswith("Amplicon_Name"): 
                continue 

            try:
                amplicon_name, primer1, primer2, refseq, pos_relative_CDS = line.split("\t")
                if not self.check_tsv_entry((amplicon_name, primer1, primer2, refseq, pos_relative_CDS)): 
                    continue
                #also check that primer1 primer2 and refseq contain valid characters
            except: 
                raise DataExtractionError("Incomplete entry in amplicon info file")
            
            fwd_pos = len(primer1) - 1
            #remove newline
            refseq = refseq.strip()
            rev_pos = len(refseq) - len(primer2) - 1

            primer2 = self.reverse_complement(primer2)

            #should check for duplicate amplicon names!!! 
            if not (amplicon_name in primer_dict.keys() or amplicon_name in refseqs.keys()):
                primer_dict[amplicon_name] = [primer1, primer2, fwd_pos, rev_pos, int(pos_relative_CDS)]
                #remove newline character if applicable
                refseqs[amplicon_name] = refseq

        #print(f"primer_dict:\n {primer_dict}\nRefseqs:\n{refseqs}")

        return primer_dict, refseqs

    def read_primer_seq_file(self, amplicon_info_file): 

        try:
            with open(amplicon_info_file) as amplicon_info:
                primer_dict, refseqs = self.parse_tsv(amplicon_info)
        except: 
            with open(amplicon_info_file, encoding = 'latin-1') as amplicon_info:
                primer_dict, refseqs = self.parse_tsv(amplicon_info)
        
        if len(primer_dict) == 0: 
            raise DataExtractionError("The amplicon info file was empty")

        return primer_dict, refseqs

    def load_paired_read_files(self, tb_seqs_file1: str, tb_seqs_file2: str):
        
        print("\n")
        read1_seqs = self.read_fastq(tb_seqs_file1)
        if len(read1_seqs) == 0: 
            raise DataExtractionError("The read 1 sequences file was empty")
        else:
            print(f"{len(read1_seqs)} read 1 sequences in {tb_seqs_file1}")

            read2_seqs = self.read_fastq(tb_seqs_file2)
            if len(read2_seqs) == 0: 
                raise DataExtractionError("The read 2 sequences file was empty")
            else:
                print(f"{len(read2_seqs)} read 2 sequences in {tb_seqs_file2}")
        
        return read1_seqs, read2_seqs

                

