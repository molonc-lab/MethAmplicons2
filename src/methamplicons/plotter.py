import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1.inset_locator import inset_axes
import seaborn as sns
#from joypy import joyplot
import pandas as pd
pd.options.mode.chained_assignment = None  # default='warn'

import numpy as np
import os
from methamplicons.extract_meth import ExtractMeth

class Plotter:

    def ridgeline(self, df_alleles_sort_all, refseqs, outpath,  save_data, amplicon_info, outname = "ridgeline_plot"): 
        # Show relative frequencies for the different numbers of methylated CpGs/epiallele by sample
        # Also we are only interested in the same region - 1 facet grid per amplicon with 1 plot per sample 

        #print(f"df_alleles_sort_all at ridgeline: \n{df_alleles_sort_all}")
        df_alleles_sort_all= df_alleles_sort_all.rename_axis('allele').reset_index()
        #print(f"df_alleles_sort_all: \n{df_alleles_sort_all}")

        data_by_amplicon = {}

        ext_meth = ExtractMeth()
        
        for amplicon_name in refseqs.keys(): 
            for col_name in df_alleles_sort_all.columns: 
                if col_name.endswith(amplicon_name):  # check if amplicon_name is at the end of col_name
                    #filtered_df = filtered_df.fillna(0) 

                    if amplicon_name not in data_by_amplicon: 
                        #initialize the allele column which is like rownames
                        data_by_amplicon[amplicon_name] = pd.DataFrame()
                        data_by_amplicon[amplicon_name]["allele"] = df_alleles_sort_all["allele"]
                    #for the corresponding amplicon dataframe, add in that sample's column
                    data_by_amplicon[amplicon_name][col_name] = df_alleles_sort_all[col_name]

        for amplicon_name, allele_data_by_sample in data_by_amplicon.items():
            fwd_pos, rev_pos = tuple(amplicon_info[amplicon_name])[2:4]
      
            num_cpg = len(ext_meth.get_cpg_positions(refseqs[amplicon_name], fwd_pos, rev_pos))

            #print(f"the number of cpgs for {amplicon_name} is {num_cpg}")

            allele_data_by_sample['allele_length'] = allele_data_by_sample['allele'].str.len()
            allele_data_by_sample = allele_data_by_sample[allele_data_by_sample['allele_length'] == num_cpg]
            allele_data_by_sample.drop(columns=['allele_length'], inplace=True)

            # number of Cs in each allele
            allele_data_by_sample['cpg'] = allele_data_by_sample['allele'].str.count('C')

            #print(f"\nallele_data_by_sample:\n{allele_data_by_sample}")

            #print(f"Allele data by sample {allele_data_by_sample.to_string()}")
            #max_cpg = allele_data_by_sample['allele'].apply(lambda x: len(x)).max()
            #print(f"The value of max cpg is {max_cpg}")


            melted_df = allele_data_by_sample.melt(id_vars=["allele", "cpg"], 
                                                var_name="sample", 
                                                value_name="count")

            melted_df = melted_df.fillna(0)

            # counts for a given number of CpGs/epiallele for each sample
            total_counts = melted_df.groupby('sample')['count'].sum()

            # relative frequencies
            melted_df['rel_freq'] = melted_df.apply(lambda row: row['count'] / total_counts[row['sample']], axis=1)

            # group by 'sample' and 'cpg', summing counts and relative frequencies
            grouped_df = melted_df.groupby(['sample', 'cpg']).agg({
                'count': 'sum',
                'rel_freq': 'sum'
            }).reset_index()

            # all CpG counts for each sample
            all_cpgs = np.arange(0, num_cpg + 1)

            # Create a dataframe with all combinations of sample and cpg count
            all_samples = melted_df['sample'].unique()
            all_combinations = pd.MultiIndex.from_product([all_samples, all_cpgs], names=['sample', 'cpg']).to_frame(index=False)

            # Merge this with the grouped_df to ensure all combinations exist
            grouped_df = all_combinations.merge(grouped_df, on=['sample', 'cpg'], how='left').fillna(0)

            # sort by sample and cpg for plotting
            sorted_df = grouped_df.sort_values(by=['sample', 'cpg'])

            pal = sns.color_palette(palette='Set2', n_colors=len(all_samples)) 
            g = sns.FacetGrid(sorted_df, row="sample", hue="sample", height=2, aspect=15, palette=pal)

            # plot relative frequencies
            g.map_dataframe(sns.lineplot, x='cpg', y='rel_freq')

            # fill the area underneath the lineplot
            g.map_dataframe(plt.fill_between, x='cpg', y1=0, y2='rel_freq', alpha=0.5)

            # add a horizontal line for each plot
            g.map(plt.axhline, y=0, lw=2, clip_on=False)

            # Setting x-ticks to integers and adjusting x-axis limit
            for ax in g.axes.flat:
                ax.set_xticks(np.arange(0, num_cpg + 1))
                ax.set_xlim(0, num_cpg)

                # Add sample name text to plots
                min_x_value = min(ax.get_xlim())
                label_x_position = min_x_value + (2 * min_x_value)
                ax.text(label_x_position, 0.02, ax.get_title(), fontweight='bold', fontsize=15, color='k')
                ax.set_title('')  

            g.set(facecolor="None")

            # space between plots
            g.fig.subplots_adjust(hspace=-0.5)

            # remove yticks
            g.set(yticks=[])
            g.despine(bottom=True, left=True)

            plt.setp(ax.get_xticklabels(), fontsize=15, fontweight='bold')
            plt.xlabel('Number of meCpGs/epiallele', fontweight='bold', fontsize=15)

            amp_out_dir = os.path.join(outpath, amplicon_name)
            if not os.path.exists(amp_out_dir):
                os.makedirs(amp_out_dir)
            
            if save_data:
                # want to save this df_alt_for_region in the corresponding amplicon folder
                sorted_df.to_csv(os.path.join(amp_out_dir,f"Ridgeline_data_{amplicon_name}.csv"))
                #melted_df.to_csv(os.path.join(amp_out_dir,f"melted_df_{amplicon_name}.csv"))
                allele_data_by_sample.to_csv(os.path.join(amp_out_dir,f"Specific_allele_data_{amplicon_name}.csv"))

            filename = f"{outname}_{amplicon_name}.pdf"
            fullpath = os.path.join(amp_out_dir, filename)
            g.savefig(fullpath)

    def plot_lollipop_colour (self, df, outpath, outname="All_samples_combined_colour.pdf"):  
        
        #print(f"Dataframe for combined samples pre-melt {df}")
        # Changing default font to Arial
        plt.rcParams['font.sans-serif'] = "Arial"
        plt.rcParams['font.family'] = "sans-serif"
        
        df_melt = df.melt(id_vars="pos")
        df_melt['variable']= df_melt["variable"].str.split('_').str[0]
        df_melt = df_melt.sort_index(ascending=False)
        #df_melt = df_melt.sort_values(by=['variable'])

        #print(f"Dataframe for combined samples melted {df_melt}")

        plt.set_cmap('coolwarm')
        plt.figure()
        fig, ax = plt.subplots(figsize=(5,4))
        ax.hlines(df_melt['variable'],min(df_melt['pos']),max(df_melt['pos']), label='_nolegend_', zorder=1)
        im = ax.scatter(df_melt['pos'],df_melt['variable'],label="meC", c=df_melt['value'],edgecolors ="black", s=50, zorder=2)
        im.set_clim(0,1)
        cbar = fig.colorbar(im, ax=ax,ticks=[0.1, 0.5, 0.9])
        cbar.ax.tick_params(labelsize=6.5, # label size
                            length=1.5, # length of ticks
                            pad = 0.4) # distance between ticks and labels

        ax.axes.set_xticks(list(df_melt['pos'].unique()))
        ax.tick_params(axis='x', which='major', labelsize=6.5, rotation=45)
        ax.tick_params(axis='y', which='major', labelsize=6.5)
        
        plt.tight_layout()
        fig.savefig(outpath + "/" + outname)
        plt.close()
    
    def plot_lollipop (self, df,sname,outpath,freq_min, amplicon_name):
    
        # Changing default font to Arial
        plt.rcParams['font.sans-serif'] = "Arial"
        plt.rcParams['font.family'] = "sans-serif"
        
        plt.set_cmap('coolwarm')

        df_C = df[df['value']=="C"]
        df_T = df[df['value']=="T"]
        
        fig, (ax1, ax2) = plt.subplots(1, 2, sharey=True, figsize=(10,5))

        ax1.hlines(df['seq'],min(df['variable']),max(df['variable']), label='_nolegend_', zorder=1)
        ax1.scatter(df_C['variable'],df_C['seq'],label="meC", color="#B30326",edgecolors ="black", s=50, zorder=2)
        ax1.scatter(df_T['variable'],df_T['seq'],label="C", color="#3A4CC0",edgecolors ="black", s=50, zorder=3)
        ax1.legend(loc='lower center', bbox_to_anchor=(0.5, -0.3), ncol=2)
        ax1.axes.set_xticks(list(df['variable'].unique()))
        ax1.tick_params(axis='x', which='major', labelsize=6.5, rotation=45)
        xlabel = amplicon_name + " CpG site"
        ax1.axes.set_xlabel(xlabel)

        ax2.barh(df['seq'], df['freq'], align='center', color='grey')
        ax2.axes.set_yticks([])
        ax2.axes.set_xlabel("Frequency (%)")
        ax2.set_xlim([0, 100])

        plt.suptitle(sname)# + f"\nMethylation alleles detected at >{freq_min}% frequency") 

        fig.tight_layout(rect=[0, 0.03, 1, 0.9])

        fig.savefig(f"{outpath}/{sname}_alleles_above_{freq_min}_perc_freq_barplot.pdf")
        
        plt.close()

    def plot_lollipop_combined (self, df,df_below_freq,sname,outpath,freq_min, amplicon_name, colbar=True):

        # Changing default font to Arial
        plt.rcParams['font.sans-serif'] = "Arial"
        plt.rcParams['font.family'] = "sans-serif"

        plt.set_cmap('coolwarm')

        df_C = df[df['value']=="C"]
        df_T = df[df['value']=="T"]

        hb=df['seq'].drop_duplicates().count()

        fig, ((ax3, ax4),(ax1, ax2)) = plt.subplots(2, 2, sharex='col',sharey='row',
                                                    gridspec_kw=dict(width_ratios=[10,5], height_ratios=[1,hb],hspace=0))    

        # Merged epialleles
        ax3.hlines(df_below_freq['seq'],min(df_below_freq['variable']),max(df_below_freq['variable']), label='_nolegend_', zorder=1)
        im = ax3.scatter(df_below_freq['variable'],df_below_freq['seq'],label='_nolegend_', c=df_below_freq['shade'],edgecolors ="black", s=50, zorder=2)
        im.set_clim(0,1)
        ax3.tick_params(axis='x', which='major', labelsize=8, rotation=45)
        ax3.axes.set_xticks([])

        # Merged epialleles - frequency
        ax4.barh(df_below_freq['seq'], df_below_freq['freq'], align='center', color='grey')
        ax4.axes.set_yticks([])
        ax4.set_xlim(0,100)

        # Individual epialleles
        ax1.hlines(df['seq'],min(df['variable']),max(df['variable']), label='_nolegend_', zorder=1)
        ax1.scatter(df_C['variable'],df_C['seq'],label="meC", color="#B30326",edgecolors ="black", s=50, zorder=2)
        ax1.scatter(df_T['variable'],df_T['seq'],label="C", color="#3A4CC0",edgecolors ="black", s=50, zorder=3)
        ax1.axes.set_xticks(list(df['variable'].unique()))
        ax1.tick_params(axis='x', which='major', labelsize=6.5, rotation=45)
        xlabel = amplicon_name + " CpG site"
        ax1.axes.set_xlabel(xlabel,size=8,labelpad=4)
        ax1.axes.set_yticks([])

        # Individual epialleles - frequency
        ax2.barh(df['seq'], df['freq'], align='center', color='grey')
        ax2.axes.set_yticks([])
        ax2.tick_params(axis='x', labelsize=6.5)
        ax2.axes.set_xlabel("Frequency (%)",size=8,labelpad=10)
        ax2.set_xlim(0,100)

        #Text labels
        ax3.text(-0.135,0.4, f"<{freq_min}% frequency \n(merged)", size=6.5, ha="center", 
                transform=ax3.transAxes)
        ax1.text(-0.14,0.9, f"≥{freq_min}% frequency", size=6.5, ha="center", 
                transform=ax1.transAxes)

        # Figure title
        fig.suptitle(sname, size=8, weight='bold')

        fig.tight_layout(rect=[0, 0.03, 1, 0.9])
        
        if colbar:
            # After tight layout, otherwise issues warnings
            # Controlling the placement of the colour bar
            adj_par = len(df.seq.unique())
            adj_height=(2.5+adj_par*2.5)
            adj_bbox_coord=(1.0436+0.0432*adj_par)

            axins = inset_axes(ax3,
                            width="60%",  
                            height=f"{adj_height}%",
                            #height="7.5%",  
                            loc='lower left',
                            bbox_to_anchor=(0, adj_bbox_coord, 1, 1), 
                            #bbox_to_anchor=(0, 1.25963685, 1, 1), 
                            bbox_transform=ax3.transAxes,
                            borderpad=0)
            # Colour bar
            cbar = fig.colorbar(im, 
                                #label='CpG methylation',
                                ax=ax3,
                                cax=axins, 
                                ticks=[0.1, 0.5, 0.9],
                                orientation="horizontal")

            cbar.ax.tick_params(labelsize=6.5, # label size
                                length=1.5, # length of ticks
                                pad = 0.4) # distance between ticks and labels

            cbar.set_label(label='CpG methylation',
                        labelpad=-25,
                        size=8)
            
            #fig.savefig(f"{outpath}/{sname}_{freq_min}perc_barplot.pdf")
            fig.savefig(f"{outpath}/{sname}_allele_freq_w_below_{freq_min}_perc_avgd_barplot.pdf")
        else:
            fig.savefig(f"{outpath}/{sname}_{freq_min}perc_barplot_nolegend.pdf")
            

        plt.close()
        