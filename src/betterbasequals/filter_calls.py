from collections import defaultdict
from betterbasequals.utils import eprint
import numpy as np
import pandas as pd


class BBQFilter:
    def __init__(
        self,
        vcf_file,
        outfile,
        cutoff_lower, 
        cutoff_upper
    ):

        self.calls =  pd.read_csv(vcf_file, sep = "\t", header = None, names = ['CHROM', 'POS', 'ID', 'REF', 'ALT', 'QUAL', 'FILTER', 'INFO'])
        self.cutoff_lower = cutoff_lower
        self.cutoff_upper = cutoff_upper
        self.outfile = outfile


    def filter_BBQ(self):
        n_var = self.calls.shape[0]
        eprint(f'Total number of calls before filtering: {n_var}')
        eprint(f'Number of PASS calls before filtering: {len(self.calls[self.calls["FILTER"]=="PASS"])}')

        # separate N and N_total info columns
        info_strings = '{"' + self.calls.INFO.str.split(';').str.join('","').str.replace('=','":"').str.replace("\"\",", "") + '"}' 
        info_df = pd.json_normalize(info_strings.apply(eval))
        calls_info = pd.concat([self.calls.reset_index(drop=True), info_df[["N", "N_total"]]], axis=1)
        calls_info['N'] = calls_info['N'].astype('int')
        calls_info['N_total'] = calls_info['N_total'].astype('int')

        # filtered coverage quantiles of PASS variants
        N = calls_info["N"][calls_info['FILTER'] == "PASS"]
        N_lower, N_upper = self.quantiles(N)
        
        # total coverage quantiles of PASS variants
        N_total = calls_info["N_total"][calls_info['FILTER'] == "PASS"]
        N_t_lower, N_t_upper = self.quantiles(N_total)

        # filter PASS variants based on coverage quantiles
        filtered_calls = calls_info[["CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO"]][(calls_info['FILTER'] != "PASS") | ((calls_info['N'].between(N_lower, N_upper)) & (calls_info['N_total'].between(N_t_lower, N_t_upper)) & (calls_info['FILTER'] == "PASS"))]
        filtered_calls = filtered_calls.sort_values(['CHROM', 'POS'], ascending=[True, True])

        # write output
        filtered_calls.to_csv(self.outfile,  sep='\t', header = False, index = False)
        
        n_filtered = n_var - filtered_calls.shape[0]

        eprint(f'Total number of calls after filtering: {filtered_calls.shape[0]}')
        eprint(f'Number of PASS calls after filtering: {len(filtered_calls[filtered_calls["FILTER"]=="PASS"])}')
    
        return n_filtered

    def quantiles(self, cov_list):
        # calculate and return quantiles of pandas series
        quantiles = cov_list.quantile([self.cutoff_lower, self.cutoff_upper])
        lower = quantiles[self.cutoff_lower]
        upper = quantiles[self.cutoff_upper]
        
        return lower, upper
                


    