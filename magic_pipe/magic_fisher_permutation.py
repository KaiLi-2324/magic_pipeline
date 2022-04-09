import time
import random
import numpy as np
from tqdm import tqdm
import multiprocessing
import scipy.stats as stats
from collections import defaultdict
import statsmodels.stats.multitest as multi
from .local_exceptions import *

"""
@Author: Kai Li
@email: likai@wiucas.ac.cn
@Description: This script is designed to perform fisher exact test and permutation for each gene
"""


class ReadSnpMatrix:

    def __init__(self, path, path_case_control, path_pathway):
        self.path = path
        self.path_case_control = path_case_control
        self.path_pathway = path_pathway
        self.pathway_genes = defaultdict(set)
        self.sample_types = dict()
        self.case = []
        self.control = []
        self.genes_variant_samples_num = defaultdict(dict)

    def get_case_control_samples(self):
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Reading from case and control file {self.path_case_control}")
        with open(self.path_case_control) as file:
            for line in file:
                items = line.strip().split("\t")
                self.sample_types[items[0]] = items[1]
                if items[1] == "case":
                    self.case.append(items[0])
                elif items[1] == "control":
                    self.control.append(items[0])
                else:
                    raise InvalidSampleTypeError("Expected only case/control sample types, got {}".format(items[1]))

    def get_pathway_genes(self):
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Reading from pathway {self.path_pathway}")
        with open(self.path_pathway) as file:
            for line in file:
                items = line.strip().split("\t")
                self.pathway_genes[items[0]].update(set(items[2:]))

    def read_in_snp_matrix(self):
        self.get_case_control_samples()
        count = 0
        total_genes = set()
        genes_samples_with_variant = defaultdict(lambda: defaultdict(set))
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: parsing snp matrix {self.path}")
        with open(self.path) as file:
            line = file.readline()
            samples = line.strip().split("\t")[5:]
            line = file.readline()
            while line:
                items = line.strip().split("\t")
                gene = items[1]
                total_genes.add(gene)
                genotypes = items[5:]
                for each_item in genotypes:
                    if each_item == "1" or each_item == "2":
                        if self.sample_types[samples[count]] == "case":
                            genes_samples_with_variant[gene]["case"].add(samples[count])
                        elif self.sample_types[samples[count]] == "control":
                            genes_samples_with_variant[gene]["control"].add(samples[count])
                        else:
                            raise ValueError("Expected only case/control types, got {}".format(self.sample_types[samples[count]]))
                    count += 1
                count = 0
                line = file.readline()

        if not self.path_pathway:
            for each_gene in total_genes:
                self.genes_variant_samples_num[each_gene]["case"] = len(genes_samples_with_variant[each_gene]["case"])
                self.genes_variant_samples_num[each_gene]["control"] = len(genes_samples_with_variant[each_gene]["control"])
        else:
            self.get_pathway_genes()
            for each_pathway in self.pathway_genes:
                this_pathway_genes = self.pathway_genes[each_pathway]
                this_pathway_case_samples = {genes_samples_with_variant[each_gene]["case"] for each_gene in this_pathway_genes}
                this_pathway_control_samples = {genes_samples_with_variant[each_gene]["control"] for each_gene in this_pathway_genes}
                self.genes_variant_samples_num[each_pathway]["case"] = len(this_pathway_case_samples)
                self.genes_variant_samples_num[each_pathway]["control"] = len(this_pathway_control_samples)


def fisher_perm_each_gene(each_gene, this_gene_case_sample_num,
                          this_gene_control_sample_num, case_samples, control_samples, hypothesis, path_pathway):
    # in this function, each_gene stands for a gene and is a bytes object, not str
    try:
        variant_sample_num = this_gene_case_sample_num + this_gene_control_sample_num

        def permutation():
            # perform permutation for each gene for 1000 times
            random_samples = random.sample(case_samples + control_samples, variant_sample_num)
            case_variant_intersect_samples = list(frozenset(random_samples).intersection(case_samples))
            control_variant_intersect_samples = list(frozenset(random_samples).intersection(control_samples))
            # yield the p value
            return stats.fisher_exact([[len(case_variant_intersect_samples),
                                       len(case_samples) - len(case_variant_intersect_samples)],
                                      [len(control_variant_intersect_samples),
                                       len(control_samples) - len(control_variant_intersect_samples)]],
                                      alternative=hypothesis)[1]

        odds_ratio, this_gene_p_value = stats.fisher_exact([[this_gene_case_sample_num,
                                                             len(case_samples) - this_gene_case_sample_num],
                                                            [this_gene_control_sample_num,
                                                             len(control_samples) - this_gene_control_sample_num]],
                                                           alternative=hypothesis)
        # if pathway is given, we return the median of the 1000 permutation values
        # if no pathway is specified, we return the 1000 permutation values
        if path_pathway:
            this_gene_p_permutation = np.median([permutation() for _ in range(1000)])
        else:
            this_gene_p_permutation = [permutation() for _ in range(1000)]
        return (each_gene, this_gene_p_value, this_gene_p_permutation, odds_ratio)
    except Exception as e:
        return (each_gene, str(e), "error")


def parse_snp_matrix(path_matrix, path_case_control, path_pathway, method, path_out, num_process,
                     progress_bar, hypothesis):
    snp_matrix = ReadSnpMatrix(path=path_matrix, path_case_control=path_case_control, path_pathway=path_pathway)
    snp_matrix.read_in_snp_matrix()
    genes = list(snp_matrix.genes_variant_samples_num.keys())
    gene_case_control_samples = defaultdict(dict)
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start performing fisher test and permutation")

    mp_list = []
    pbar = ""
    if progress_bar:
        pbar = tqdm(desc=f"  Fisher test/Permutation", total=len(genes))

    def call_back(result):
        if progress_bar:
            pbar.update()
        mp_list.append(result)

    pool = multiprocessing.Pool(processes=int(num_process))
    for each_gene in genes:
        this_gene_variant_case_samples_num = snp_matrix.genes_variant_samples_num[each_gene]["case"]
        this_gene_variant_control_samples_num = snp_matrix.genes_variant_samples_num[each_gene]["control"]
        gene_case_control_samples[each_gene]["case"] = this_gene_variant_case_samples_num
        gene_case_control_samples[each_gene]["control"] = this_gene_variant_control_samples_num
        pool.apply_async(func=fisher_perm_each_gene, args=(each_gene, this_gene_variant_case_samples_num,
                                                           this_gene_variant_control_samples_num,
                                                           snp_matrix.case, snp_matrix.control, hypothesis,
                                                           path_pathway,), callback=call_back)

    pool.close()
    pool.join()

    if progress_bar:
        pbar.close()

    genes_pvals_pairs = dict()
    genes_with_error = dict()
    total_genes_permutation = []
    if not mp_list:
        raise EmptyResultError("No results got!")
    for each_item in mp_list:
        each_gene = each_item[0]
        if each_item[-1] == "error":
            genes_with_error[each_gene] = each_item[1]
        else:
            gene_case_control_samples[each_gene]["case_wild"] = len(snp_matrix.case) - gene_case_control_samples[each_gene]["case"]
            gene_case_control_samples[each_gene]["control_wild"] = len(snp_matrix.control) - gene_case_control_samples[each_gene]["control"]
            # each_item[3] is odds ratio
            gene_case_control_samples[each_gene]["odds_ratio"] = each_item[3]
            # each_item[1] is p value
            genes_pvals_pairs[each_gene] = each_item[1]
            if path_pathway:
                # if pathway is input, the each_item[2] is a number, if not, it's a python list
                gene_case_control_samples[each_gene]["pperm"] = each_item[2]
            else:
                total_genes_permutation.extend(sorted(each_item[2]))

    if genes_with_error:
        path_error = f"{path_out}.err"
        print(f"[ERROR]: Error occurred when performing fisher/permutation for certain genes, Writing failed parsing "
              f"genes to {path_error}")
        with open(path_error, "w+") as err:
            for each_gene in genes_with_error:
                err.write(f"{each_gene}\t{genes_with_error[each_gene]}\n")
        raise GenesParsingError(genes_with_error)

    # Because genes in the Manager.dict is out of order after the multiprocessing method, we have to reorder the genes
    # Calculate adjusted p value for each gene, and sort the genes according to it's p value by ascending order
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start sorting genes and adjusting p values")
    # sort genes_pvals_pairs(dict) by value(p value), and return sorted keys of genes_pvals_pairs
    genes_sorted = sorted(genes_pvals_pairs, key=lambda k: genes_pvals_pairs[k])
    pvals_sorted = [genes_pvals_pairs[each_gene] for each_gene in genes_sorted]
    pvals_adjusted = dict(zip(genes_sorted, multi.multipletests(pvals_sorted, method=method)[1].tolist()))

    # if pathway is not input, merge all the p perm into a sorted list, and calculate the permutation for each gene by
    # each 1000 p perms in the list, if it's input, just assign total_genes_permutation to total_genes_permutation_mean
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start calculating p permutation")
    if path_pathway:
        total_genes_permutation_mean = [gene_case_control_samples[each_gene]["pperm"] for each_gene in genes_sorted]
    else:
        total_genes_permutation_mean = list(map(np.mean, np.array_split(sorted(total_genes_permutation), len(genes_sorted))))
    if len(total_genes_permutation_mean) != len(genes_sorted):
        raise GenesPermNumberNotEqualError("Gene numbers not equal to the calculated p permutation numbers")
    genes_pperm_pairs = dict(zip(genes_sorted, total_genes_permutation_mean))

    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start writing to {path_out}")
    with open(path_out, "w+") as out:
        out.write("gene\tcase\tcase_wild\tcontrol\tcontrol_wild\tp_value\todds_ratio\tp_adj\tp_permutation\n")
        for each_gene in genes_sorted:
            out.write("{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\t{}\n".format(each_gene,
                                                                    gene_case_control_samples[each_gene]["case"],
                                                                    gene_case_control_samples[each_gene]["case_wild"],
                                                                    gene_case_control_samples[each_gene]["control"],
                                                                    gene_case_control_samples[each_gene]["control_wild"],
                                                                    genes_pvals_pairs[each_gene],
                                                                    gene_case_control_samples[each_gene]["odds_ratio"],
                                                                    pvals_adjusted[each_gene],
                                                                    genes_pperm_pairs[each_gene]))
