from methamplicons.plotter import Plotter
from methamplicons.extract_meth import ExtractMeth
import argparse
import os
import shutil
import sys
import pandas as pd
from functools import reduce
import re
import numpy as np
from collections import defaultdict

class MethAmplicon:
    
    def __init__(self):
        self.extract_meth = ExtractMeth()
        self.plotter = Plotter()

        # use argparse for parsing of command line input 
        self.parser = argparse.ArgumentParser(description='CLI tool for plotting targeted bisulfite sequencing')
        self.setup_parser()

        self.refseq_cpgs = {}
        self.amplicon_info = {}

        self.sample_efficiencies = {}

    def valid_file(self, file_path, seq_type, extensions):
        """
        This function checks if a file at a provided path has a given extension and does exist
        This way, an obvious error in the file paths is pointed out right away
        And the correctness of the actual contents of the file will be checked by downstream processing
        """
        if not os.path.exists(file_path):
            raise argparse.ArgumentTypeError(f'{seq_type} file: {file_path} does not exist')
        if not any(file_path.endswith(ext) for ext in extensions):
            raise argparse.ArgumentTypeError(f'{seq_type} file must end with one of the \
                                             following extensions: {", ".join(extensions)}')
        return file_path
    
    def valid_PE_read_dir(self, dir_path): 
        if not os.path.isdir(dir_path): 
            raise argparse.ArgumentTypeError(f'Input directory for PE reads {dir_path} is not a directory')
        return dir_path

    def valid_amplicon_info_file(self, file_path):
        return self.valid_file(file_path, "Primers and Reference Sequences", ['.tsv'])
    
    def valid_labels_file(self, file_path):
        return self.valid_file(file_path, "Labels", ['.csv'])
    
    def valid_out_dir(self, dir_path): 
        if not os.path.isdir(dir_path): 
            os.makedirs(dir_path, exist_ok=True)

        return dir_path
    
    def valid_thresh(self, freq_thresh):
        freq_thresh = float(freq_thresh)
        if not (freq_thresh >= 0 and freq_thresh < 1): 
            raise argparse.ArgumentTypeError(f'{freq_thresh} is an invalid frequency - frequency threshold must be >= 0 and <1')
        
        return freq_thresh

    def setup_parser(self):
        """
        Create the arguments that the command line tool will use
        Set some of the arguments to default values
        """
        self.parser.add_argument('--PE_read_dir', type=self.valid_PE_read_dir, \
                                 default = os.getcwd(), help="Desired input directory with fastq files, gzipped or not")
    
        self.parser.add_argument('--amplicon_info', type=self.valid_amplicon_info_file, \
                                 help="Path to the amplicon info file in tsv format, e.g.: \
                                 AmpliconName   Primer1   Primer2  ReferenceSequence")
        
        self.parser.add_argument('--sample_labels', type=self.valid_labels_file, nargs='?', \
                                 help="Path to sample labels file in csv format", default=None)
        
        #desired directory for all output from merged files to dataframes and graphs
        #will be current working directory by default. 
        self.parser.add_argument('--output_dir', type=self.valid_out_dir, \
                                 default = os.getcwd(), help="Desired output directory")
        
        self.parser.add_argument('--min_seq_freq', type=self.valid_thresh, \
                                 default = 0.01, help="Threshold frequency an extracted epiallele sequence must have to be included in analysis (default:0.01)")
        
        self.parser.add_argument('--verbose', type=str, choices=['true', 'false'], \
                                 default='true', help="Print all output after file parsing (default: true).")
        
        # the save_data argument is true by default, and the user can also set it to false with --save_data false
        self.parser.add_argument('--save_data', type=str, choices=['true', 'false'], \
                                 default='true', help="Save processed data in csv format (default: true).")
        
        # add an option to save or delete intermediate files -default will be delete
        self.parser.add_argument('--save_intermediates', type=str, choices=['true', 'false'], \
                                 default='true', help="Save 'demultiplexed' and merged read files for all combinations of samples and amplicons (default: true).")
        
        # add an option to save or delete intermediate files -default will be delete
        self.parser.add_argument('--combine_lanes', type=str, choices=['true', 'false'], \
                                 default='true', help="Combine fastq files from different sequencing lanes (L001, L002) into single R1 and R2 files (default: true).")

    @staticmethod
    def replace_last(source_string, replace_what, replace_with):
        head, _sep, tail = source_string.rpartition(replace_what)
        return head + replace_with + tail
    
    def combine_lanes(self):
        """
        Combine all lanes for R1 and R2 files with the same basename.
        """
        files = os.listdir(self.args.PE_read_dir)
        grouped_files = defaultdict(lambda: {'R1': [], 'R2': []})

        # group R1 R2 files by basename 
        for file in files:
            if file.endswith('.fastq') or file.endswith('.fq') \
               or file.endswith('.fastq.gz') or file.endswith('.fq.gz'):
                basename_match = re.search(r'^(.*?)_L00[0-9]', file)
                if basename_match:
                    basename = basename_match.group(1)
                    if "R1" in file:
                        r2_file = self.replace_last(file, "R1", "R2")
                        if r2_file in files:    
                            grouped_files[basename]["R1"].append(file)
                            grouped_files[basename]["R2"].append(r2_file)
                    #else continue

        print(f"grouped files:\n{grouped_files}")

        # combine files for each group (lanes for a sample)
        for basename, reads in grouped_files.items():

            # will do R1 files then R2 files
            for read_type, files in reads.items():
                #{R1: [L001_R1.fq, L002_R1.fq], R2: ..}
                combined_seqs = {}
                for file in files:
                    # will get dict with {@VH01336:71:AACC3JLHV:1:1101:28570:1000:TGTTGTT+ACCTATC : ATGTG...}
                    # the seq info will match two reads and will only correspond to one lane
                    seqs = self.extract_meth.read_fastq(os.path.join(self.args.PE_read_dir, file))
                    # all R1 reads for multiple lanes now together 
                    combined_seqs.update(seqs)

                # Construct the filename for the combined file
                # this will get rid of whatever came after R1/R2 but that's probably fine
                extension = '.fastq'  # Assuming standard extension; adjust if necessary
                combined_file_name = f"{basename}_all_lanes_{read_type}{extension}"
                
                # Write combined sequences to a new file
                with open(os.path.join(self.read_dir, combined_file_name), 'w') as f:
                    for seq_info, (sequence, phred) in combined_seqs.items():
                        f.write(f"@{seq_info}\n{sequence}\n+\n{phred}\n")

        
    def get_paired_files(self):
        """
        Get all pairs of fastq files in a directory
        """
        fastq_files = [f for f in os.listdir(self.read_dir) \
                       if f.endswith('.fastq') or f.endswith('.fq')\
                        or f.endswith('.fastq.gz') or f.endswith('.fq.gz')]
        ##print(f"Found fastq files in input directory: {fastq_files}")
        paired_files = []
        for f in fastq_files:
            if "R1" in f:
                # Replace the last occurrence of R1 with R2 to find the pair
                r2_file = self.replace_last(f, "R1", "R2")
                
                if r2_file in fastq_files:
                    paired_files.append((os.path.join(self.read_dir, f), os.path.join(self.read_dir, r2_file)))

        print(f"paired files: {paired_files}")

        return paired_files
        
    def merge_loop(self):
        """
        Iterate over all files in the given input directory, find paired files, and created merged read files
        """ 
        paired_files = self.get_paired_files()
        
        for read_file_pair in paired_files: 
            # merge the reads
            self.extract_meth.merge_reads(read_file_pair[0], read_file_pair[1], self.refseqs, self.amplicon_info, self.args.output_dir)

    def get_amp_name(self,file_name): 
        """
        Return the name of the amplicon based on the merged read file name
        Use this to group amplicon data into directories with data for multiple samples
        """
        amplicon_name = ""
        #avoid the error caused by parsing file with delimiter
        for amp_name in self.refseqs.keys(): 
            if amp_name in file_name:
                # if the name of a region is in the file name it is most likely the name
                # but will also check that the only thing that comes after the name is .extendedFrags.fastq
                if file_name.split(amp_name)[-1:] == [".extendedFrags.fastq"]:
                    amplicon_name = amp_name
                    #print(f"amplicon name: {amplicon_name} in file name")

        return amplicon_name
    
    @staticmethod
    def use_basename(file_name):
        pattern = r".extendedFrags.fastq"
        replacement = ""
        sid = re.sub(pattern, replacement, os.path.basename(file_name))
        return sid

    def get_sname(self, file_name, amplicon_name): 
        """
        Get the sample name from the file name - for merged read files
        """

        if self.args.combine_lanes:
            match = re.search(r'^(.*?)_all_lanes', file_name)
            sid = match.group(1)
        else:
            match = re.search(r'^(.*?)_L00[0-9]', file_name)
            sid = match.group(1)

        # with the sid, try to see if there is a corresponding sample name in the 
        # sample name csv
        try: 
            sname=self.labels_df.loc[sid]['ShortLabel']
            if pd.isnull(sname):
                sname=self.labels_df.loc[sid]['SampleLabel']
        except:
            sname = None

        return sname
    
    def plot_per_sample_lollipop(self, alleles_sort,refseq, fwd_pos, rev_pos, filtered_reads, pos_rel_CDS, sname, amplicon_name, output_dir):
        """
        Plots lollipop plots containing the specific epiallales for a given sample by calling methods in Plotter
        """
        # Prepare individual sample plots grouping alleles <2% and 5% - actually just 5% but will leave in a loop
        # in the event multiple minimum frequencies are desired
        for freq_min in [5]:
            df=self.extract_meth.convert_to_df(alleles_sort,refseq, fwd_pos, rev_pos, filtered_reads,freq_min)
            ##print(f"Dataframe for first plot {df}")
            df_below_freq=self.extract_meth.calculate_meth_fraction_min(alleles_sort, refseq, fwd_pos, rev_pos, filtered_reads,freq_min)

            sum_variable = df.variable + pos_rel_CDS
            df.variable = np.where(sum_variable >= 0, sum_variable + 1, sum_variable)
            sum_variable_below = df_below_freq.variable + pos_rel_CDS
            df_below_freq.variable = np.where(sum_variable_below >= 0, sum_variable_below + 1, sum_variable_below)
            
            if df_below_freq.freq.sum() > 0:  
                # if you have epialleles with frequency below 5%      
                self.plotter.plot_lollipop_combined(df,df_below_freq,sname,output_dir,freq_min, amplicon_name)
            else:     
                self.plotter.plot_lollipop(df,sname,output_dir,freq_min, amplicon_name)

    
    def do_combined_lollipop(self, df_alt, amplicon_names): 
        df_alts_by_region = {}

        df_alt= df_alt.rename_axis('position').reset_index()
        #print(f"df_alt columns {df_alt.columns}")
        
        for amplicon_name in amplicon_names: 
            for col_name in df_alt.columns: 
                #print(f"1. col_name is {col_name}")
                if amplicon_name in col_name: 
                    # need to create a dataframe with alleles specific to that amplicon (remove NAN for that column)
                    """
                    Filter out NAN values for a column corresponding to a given amplicon-sample combination when creating a new dataframe for a given amplicon, e.g. RAD51C
                    , and also before merging to an existing dataframe for a given amplicon so that only alleles corresponding to a given amplicon are included in its ridgeline plot
                    """
                    filtered_df_alt = df_alt[df_alt[col_name].notna()]
                    #reset indices to get rid of mismatched indices warning message - index does not provide information 
                    # in this case so it should not be important
                    df_alt = df_alt.reset_index(drop=True)

                    if amplicon_name not in df_alts_by_region: 
                        df_alts_by_region[amplicon_name] = pd.DataFrame()
                        df_alts_by_region[amplicon_name]["position"] = filtered_df_alt["position"]
                    #could change this to use only the sample name
                    #need to remove all NAs, then merge, then convert NAs to zeros
                    df_alts_by_region[amplicon_name][col_name] = filtered_df_alt[col_name]

        for amplicon_name, df_alt_for_region in df_alts_by_region.items():
            pos_rel_CDS = self.amplicon_info[amplicon_name][-1]
            pos_promoter = list(map(lambda x: x + pos_rel_CDS + 1 if (x + pos_rel_CDS) >= 0 else x + pos_rel_CDS, self.extract_meth.get_cpg_positions(self.refseqs[amplicon_name], self.amplicon_info[amplicon_name][2],self.amplicon_info[amplicon_name][3] )))
            df_alt_for_region['pos'] = pos_promoter
            df_alt_for_region.drop(columns=["position"], inplace=True)
            amp_out_dir = os.path.join(self.args.output_dir, amplicon_name)

            if not os.path.exists(amp_out_dir):
                os.makedirs(amp_out_dir)
            
            if self.args.save_data:
                # want to save this df_alt_for_region in the corresponding amplicon folder
                df_alt_for_region.to_csv(os.path.join(amp_out_dir,f"average_{amplicon_name}_meth_by_sample.csv"))

            self.plotter.plot_lollipop_colour(df=df_alt_for_region, outpath=amp_out_dir,
                             outname=f"All_samples_combined_avgd_meth_{amplicon_name}.pdf")


    def meth_amplicon_loop(self):
        
        allele_sort_dfs = []

        df_alt = pd.DataFrame() 
        merged_path = os.path.join(self.args.output_dir, "merged")

        #After merging of reads is completed, program flow is essentially the same as in MethAmplicon
        for i,file in enumerate([f for f in os.listdir(merged_path) \
                       if f.endswith("extendedFrags.fastq")]):

            #print(f"Identified a merged file {file}")
            amplicon_name = self.get_amp_name(file)
            
            # create a directory for the specific amplicon if it does not already exist
            amp_out_dir = os.path.join(self.args.output_dir, amplicon_name)
            if not os.path.exists(amp_out_dir):
                os.makedirs(amp_out_dir)

            sname = self.get_sname(file, amplicon_name)

            if sname == None:
                # a short name for the file could not be found
                sname = self.use_basename(file)
                #print(f"full file sname:{sname}")
            else:
                sname = sname + "_parse_" + amplicon_name
                #print(f"short sname:{sname}")

            file_path = os.path.join(merged_path, file)
            fwd_primer, rev_primer, fwd_pos, rev_pos, pos_rel_CDS = tuple(self.amplicon_info[amplicon_name])

            print(f"Getting epiallele counts for sample '{self.extract_meth.parse_name(sname)}' for amplicon '{amplicon_name}'")

            # Generate dictionary with all reads for that amplicon in the merged file 
            d, below_thresh =self.extract_meth.get_all_reads(file_path, fwd_primer, rev_primer)
            #print(f'read counts dictionary\n{d}')

            refseq = self.refseqs[amplicon_name]

            # Count only CpG sites in alleles
            alleles_sort,filtered_reads, filt_for_length, filt_for_CpG_AG, reads_n =self.extract_meth.count_alleles(d, refseq, fwd_pos, rev_pos)

            ### make this a function - save_info() in the future
                
            self.extract_meth.set_threshold(0)
            all_reads_df, zero_reads =self.extract_meth.get_all_reads(file_path, fwd_primer, rev_primer,printout=False)
            self.extract_meth.set_threshold(self.args.min_seq_freq)
            
            num_ts_obs, exp_ts, num_reads_used, num_non_cpg_cs = self.extract_meth.get_efficiency_vals(all_reads_df, refseq, fwd_pos, rev_pos)

            # split sname into sample and amplicon
            if '_parse_' in sname:
                sample, amplicon = sname.split('_parse_', 1)
            elif '_all_lanes_' in sname:
                sample, amplicon = sname.split('_all_lanes_', 1)
            else:
                sample = sname
                amplicon = amplicon  

            key = (sample, amplicon)

            if key not in self.sample_efficiencies:
                if exp_ts == "Empty":
                    self.sample_efficiencies[key] = ["No_reads", None, None, None, None, filtered_reads, filt_for_length, filt_for_CpG_AG, reads_n, below_thresh]
                elif exp_ts == "Badseqs":
                    self.sample_efficiencies[key] = ["None_w_length_refseq", None, None, None, None, filtered_reads, filt_for_length, filt_for_CpG_AG, reads_n, below_thresh]
                elif not exp_ts == 0:
                    self.sample_efficiencies[key] = [num_ts_obs / exp_ts, num_ts_obs, exp_ts, num_reads_used, num_non_cpg_cs, filtered_reads, filt_for_length, filt_for_CpG_AG, reads_n, below_thresh]
                    # print(f"for sample {sample} and amplicon {amplicon}: num_ts_obs={num_ts_obs}, exp_ts={exp_ts}")
                else:
                    self.sample_efficiencies[key] = ["No_non_CpG_cs", None, None, None, None]
            else:
                self.sample_efficiencies[key] = [f"Sample amplicon pair name is not unique but efficiency is {num_ts_obs / exp_ts}", None, None, None, None, filtered_reads, filt_for_length, filt_for_CpG_AG, reads_n, below_thresh]
                print("Attempted to record bisulfite conversion efficiency for a \
                    sample amplicon pair twice, there should only be one of each \
                    sample amplicon pair, if two samples have the same name, please rename one sample")
            ### end of save_info function 
                    
            if alleles_sort == []: 
                print(f"No epialleles were found for amplified region: {amplicon_name} for {sname}, trying next region")
                print("\n")
                #should not have to use continue 
                continue
                
            df_alleles_sort = pd.DataFrame(alleles_sort, columns=['allele',sname])
            #print(f"Alleles sorted as dataframe for {sname} and {amplicon_name}: \n {df_alleles_sort.to_string()}")

            allele_sort_dfs.append(df_alleles_sort.copy())

            #Calculate methylation fraction per CpG site
            df_sample=self.extract_meth.calculate_meth_fraction(alleles_sort, refseq, fwd_pos, rev_pos)

            amp_out_dir = os.path.join(self.args.output_dir, amplicon_name)
            if not os.path.exists(amp_out_dir):
                os.makedirs(amp_out_dir)
            
            if self.args.save_data:
                sname_parsed = self.extract_meth.parse_name(sname)
                # want to save this df_alt_for_region in the corresponding amplicon folder
                df_alleles_sort.to_csv(os.path.join(amp_out_dir,f"{sname_parsed}_df_alleles_sort.csv")) 
            
            df_sample.columns=[sname]
            if i == 0:
                df_alt=df_sample
            elif i > 0:
                df_alt=df_alt.join(df_sample, how='outer')

            self.plot_per_sample_lollipop(alleles_sort,refseq, fwd_pos, rev_pos, filtered_reads, pos_rel_CDS, sname, amplicon_name, amp_out_dir)
            
            print("\n")

        dfs = [df.set_index('allele') for df in allele_sort_dfs]

        df_alleles_sort_all2 = reduce(lambda left, right: pd.merge(left, right, on='allele', how='outer'), dfs)

        #get the names of the different amplicons
        amplicon_names = self.refseqs.keys()
        #plot a ridgeline plot based on the accumulated data from multiple samples
        self.plotter.ridgeline(df_alleles_sort_all2, self.refseqs, self.args.output_dir, self.args.save_data, self.amplicon_info)

        self.do_combined_lollipop(df_alt, amplicon_names)

    def choose_read_dir(self):

        if self.args.combine_lanes:
            # make the combine lanes directory name 
            self.all_lanes_dir = os.path.join(self.args.output_dir, "combined_lanes")
            # set the read directory to that value
            self.read_dir = self.all_lanes_dir
            if not os.path.exists(self.all_lanes_dir):
                # if the folder doesn't exist yet, the combined fastq files haven't been made
                os.makedirs(self.all_lanes_dir)
                # so run the function to make them - possibly time consuming so only want to run 
                # if this folder doesn't exist yet
                self.combine_lanes()
            else:
                if os.listdir(self.all_lanes_dir) == []:
                    self.combine_lanes()
                else:
                    print("Using fastq files in existing combined lanes directory.\n \
                     If you would like to remake these combined fastq files, either:\n \
                      1) Specify a new output directory to make a new combined lanes directory \
                      2) Delete this combined lanes directory in the given output directory")

                self.read_dir = self.all_lanes_dir

        else:
            self.read_dir = self.args.PE_read_dir

    def save_info_csv(self):
        keys = pd.DataFrame.from_records(list(self.sample_efficiencies.keys()), columns=['Sample', 'Amplicon'])
        values = pd.DataFrame.from_records(list(self.sample_efficiencies.values()), \
                                            columns=['BS_Conv_Eff', 'Num_Ts_Obs', \
                                                    'Num_Exp_Ts_Total', 'Num_Reads_Used_Non_CpG', \
                                                    'Num_Non_CpG_Cs', "Retained_for_CpG_Total", \
                                                    "Excl_for_CpG_length", "Excl_for_CpG_AG", \
                                                    "Reads_above_thresh", "Reads_below_thresh_ct"])
        efficiency_df = pd.concat([keys, values], axis=1)
        efficiency_df.to_csv(os.path.join(self.args.output_dir, "bisulfite_seq_info.csv"), index=False)

        
    def run(self):
        # when the app is run from the command line, parse the arguments
        self.args = self.parser.parse_args()

        if not self.args.PE_read_dir:
            raise argparse.ArgumentTypeError("--PE_read_dir (Paired end reads directory) is missing!")
        if not self.args.amplicon_info:
            raise argparse.ArgumentTypeError("--amplicon_info (amplicon info tsv file) is missing!")
        
        seq_freq_threshold = self.args.min_seq_freq
        self.extract_meth.set_threshold(seq_freq_threshold)

        self.choose_read_dir()
        
        #create df for short/sample name lookup from file name from the provided csv
        try:
            self.labels_df = pd.read_csv(self.args.sample_labels, index_col = 0)
        except:
            self.labels_df = pd.DataFrame()

        # want to get the column for shortlabels in order
        self.plotter.set_labels(self.labels_df)
        
        # process tsv file
        self.amplicon_info, self.refseqs = self.extract_meth.read_primer_seq_file(self.args.amplicon_info)
        
        # disable print
        if self.args.verbose == "false":
            sys.stdout = open(os.devnull, 'w')
        # pass the verbose value to extract_meth to optionally redirect flash subprocess output
        self.extract_meth.set_verbose(self.args.verbose)

        # iterate over the paired end read files and process data 
        self.merge_loop()
        self.meth_amplicon_loop()

        self.save_info_csv()

        if self.args.save_intermediates == "false":
            # delete the merged and demultiplexed directories and all files they contain
            dirs_to_delete = ["merged", "demultiplexed", "combined_lanes"]
            print("Deleting combined lane, merged and demultiplexed read directories - to run without deleting, run with --save_intermediates true")

            for dir_to_delete in dirs_to_delete:
                full_path = os.path.join(self.args.output_dir, dir_to_delete)
                if os.path.exists(full_path):
                    shutil.rmtree(full_path)


def main():   

    app = MethAmplicon()
    app.run()

if __name__ == "__main__":
    main()