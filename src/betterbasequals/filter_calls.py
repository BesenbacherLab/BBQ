from betterbasequals.utils import eprint
from operator import itemgetter
import numpy as np


class BBQFilter:
    def __init__(
        self,
        vcf_file,
        outfile,
        cutoff_lower, 
        cutoff_upper
    ):

        self.vcf = vcf_file
        self.cutoff_lower = cutoff_lower
        self.cutoff_upper = cutoff_upper
        self.outfile = outfile


    def filter_BBQ(self):
        cov_list = []
        cov_total_list = []
        calls = open(self.vcf, 'r')
        n_pass = 0

        # create a list of PASS variant coverages
        for line in calls:
            if line.split("\t")[6] == "PASS":
                n_pass += 1
                cov = itemgetter(1, 4)(line.split("\t")[7].split(";"))
                cov_list.append(int(cov[0].split("=")[1]))
                cov_total_list.append(int(cov[1].split("=")[1]))
        calls.close()
        eprint(f'Number of PASS calls before filtering: {n_pass}')
        
        cov_arr = np.array(cov_list)
        cov_total_arr = np.array(cov_total_list)

        # filtered coverage quantiles of PASS variants
        cov_q_l = np.quantile(cov_arr, self.cutoff_lower)
        cov_q_u = np.quantile(cov_arr, self.cutoff_upper)

        # total coverage quantiles of PASS variants
        cov_total_q_l = np.quantile(cov_total_arr, self.cutoff_lower)
        cov_total_q_u = np.quantile(cov_total_arr, self.cutoff_upper)

        # filter PASS variants based on coverage quantiles 
        pass_count = 0
        n_filtered = 0
        n_pass = 0
        with open(self.vcf, 'r') as calls:
            for line in calls:
                if line.split("\t")[6] == "PASS":
                    cov = cov_list[pass_count]
                    cov_total = cov_total_list[pass_count]
                    pass_count += 1
                    if cov_q_l <= cov <= cov_q_u:
                        if cov_total_q_l <= cov_total <= cov_total_q_u:
                            n_pass += 1
                            self.outfile.write(line)
                        else: 
                            n_filtered += 1
                    else: 
                        n_filtered += 1
                else: 
                    self.outfile.write(line)

        eprint(f'Number of PASS calls after filtering: {n_pass}')
        
        return n_filtered
                


    
