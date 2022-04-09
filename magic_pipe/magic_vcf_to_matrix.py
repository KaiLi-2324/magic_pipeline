import re
import time
import portion
from collections import Counter
from .local_exceptions import *


"""
@Author: Kai Li
@Email: likai@wiucas.ac.cn
@Description: This script is designed to analyze vcf from GATK
"""


def replace(snp):
    nums = snp.split("/")
    if re.search(r"\./\.", snp):
        return "."
    else:
        return str(int(nums[0]) + int(nums[1]))


def parse_vcf(path_vcf, i, I, trans_type, variant_type, e, path_out):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start reading from {path_vcf}")
    lof_types = ["stop_gained", "frameshift_variant", "start_lost", "splice_acceptor_variant", "splice_donor_variant"]
    missense_types = ["inframe_deletion", "inframe_insertion", "missense_variant", "stop_lost"]
    output_file = open(path_out, "w+")
    formats = []
    sample_names = []
    snp_out_format = ["0/0", "0/1", "1/1", "1/0", "./."]
    pattern_AF = re.compile(r";AF=(.*?);AN=")
    pattern_CSQ = re.compile(r";CSQ=(.*?)$")
    with open(path_vcf) as file:
        for line in file:
            if line.startswith("#"):
                if re.search(r'Format: (.*?)">', line):
                    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: start parsing format line")
                    formats = re.findall(r'Format: (.*?)">', line)[0].split("|")
                elif line.startswith("#CHROM"):
                    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: start collecting sample names")
                    sample_names = line.strip().split("\t")[9:]
                    output_file.write("variant\tsymbol\tensemble\tstat\tinfo\t{}\n".format("\t".join(sample_names)))
            elif not formats:
                raise FormatLineNotFoundError(f"No format line found in {path_vcf}! 'Format: Allele|Consequence|IMPACT|SYMBOL|Gene|Feature_type|Feature...'")
            elif not sample_names:
                raise SamplesNotFoundError(f"No sample names found in {path_vcf}!")
            else:
                items = line.strip().split("\t")
                AF = float(pattern_AF.findall(items[7])[0])
                # generate an interval [i, I], if AF not in interval, then pass this line
                section = portion.closed(i, I)
                if AF not in section:
                    continue
                transcripts = pattern_CSQ.findall(items[7])[0].split(",")
                this_line_gene_trans_out = dict()
                # iterate over each transcript
                for each_transcript in transcripts:
                    this_transcript_items = each_transcript.split("|")
                    # replace the empty string in this_transcript_items with 0
                    this_transcript_format_dict_rep0 = dict(zip(formats, list(map(lambda k: k if k else 0, this_transcript_items))))
                    if not this_transcript_format_dict_rep0[trans_type]:
                        continue
                    # if not max(float(this_transcript_format_dict_rep0["gnomADg_AF"]), float(this_transcript_format_dict_rep0["gnomAD_AF"])) <= e:
                    #   continue
                    if not float(this_transcript_format_dict_rep0["MAX_AF"]) <= e:
                        continue
                    if max(float(this_transcript_format_dict_rep0["SpliceAI_pred_DS_AG"]),
                           float(this_transcript_format_dict_rep0["SpliceAI_pred_DS_AL"]),
                           float(this_transcript_format_dict_rep0["SpliceAI_pred_DS_DG"]),
                           float(this_transcript_format_dict_rep0["SpliceAI_pred_DS_DL"])) < 0.5:
                        splice_ai = True
                    else:
                        splice_ai = False
                    if variant_type == "lof":
                        if any(list(map(lambda k: 1 if re.search(r"(^{}$|^{}&|&{}&|&{}$)".format(k, k, k, k), this_transcript_format_dict_rep0["Consequence"]) else 0, lof_types))):
                            this_line_gene_trans_out[f"{this_transcript_format_dict_rep0['SYMBOL']}_{this_transcript_format_dict_rep0['Gene']}"] = each_transcript
                        # any of the splice ai must be greater than or equal to 0.5
                        elif max(float(this_transcript_format_dict_rep0["SpliceAI_pred_DS_AG"]),
                                 float(this_transcript_format_dict_rep0["SpliceAI_pred_DS_AL"]),
                                 float(this_transcript_format_dict_rep0["SpliceAI_pred_DS_DG"]),
                                 float(this_transcript_format_dict_rep0["SpliceAI_pred_DS_DL"])) >= 0.5:
                            if this_transcript_format_dict_rep0["SpliceAI_pred_SYMBOL"] == this_transcript_format_dict_rep0["SYMBOL"] and this_transcript_format_dict_rep0["LoF"] != "LC":
                                this_line_gene_trans_out[f"{this_transcript_format_dict_rep0['SYMBOL']}_{this_transcript_format_dict_rep0['Gene']}"] = each_transcript
                    elif variant_type == "missense_benign_1":
                        if any(list(map(lambda k: 1 if re.search(r"(^{}$|^{}&|&{}&|&{}$)".format(k, k, k, k), this_transcript_format_dict_rep0["Consequence"]) else 0, missense_types))) and float(this_transcript_format_dict_rep0["CADD_PHRED"]) < 15:
                            if splice_ai:
                                this_line_gene_trans_out[f"{this_transcript_format_dict_rep0['SYMBOL']}_{this_transcript_format_dict_rep0['Gene']}"] = each_transcript
                    elif variant_type == "missense_benign_2":
                        if re.search(r"^tolerated\(.*\)$", str(this_transcript_format_dict_rep0["SIFT"])) and re.search(r"^benign\(.*?\)$", str(this_transcript_format_dict_rep0["PolyPhen"])):
                            if splice_ai:
                                this_line_gene_trans_out[f"{this_transcript_format_dict_rep0['SYMBOL']}_{this_transcript_format_dict_rep0['Gene']}"] = each_transcript
                    elif variant_type == "missense_damage_1":
                        if any(list(map(lambda k: 1 if re.search(r"(^{}$|^{}&|&{}&|&{}$)".format(k, k, k, k), this_transcript_format_dict_rep0["Consequence"]) else 0, missense_types))) and float(this_transcript_format_dict_rep0["CADD_PHRED"]) >= 15:
                            if splice_ai:
                                this_line_gene_trans_out[f"{this_transcript_format_dict_rep0['SYMBOL']}_{this_transcript_format_dict_rep0['Gene']}"] = each_transcript
                    elif variant_type == "missense_damage_2":
                        if re.search(r"^deleterious\(.*\)$", str(this_transcript_format_dict_rep0["SIFT"])) and re.search(r"^probably_damaging\(.*?\)$", str(this_transcript_format_dict_rep0["PolyPhen"])):
                            if splice_ai:
                                this_line_gene_trans_out[f"{this_transcript_format_dict_rep0['SYMBOL']}_{this_transcript_format_dict_rep0['Gene']}"] = each_transcript
                    elif variant_type == "synonymous":
                        if this_transcript_format_dict_rep0["Consequence"] == "synonymous_variant":
                            if splice_ai:
                                this_line_gene_trans_out[f"{this_transcript_format_dict_rep0['SYMBOL']}_{this_transcript_format_dict_rep0['Gene']}"] = each_transcript
                if this_line_gene_trans_out:
                    samples_info_items = list(map(lambda k: k.split(":")[0], items[9:]))
                    # count the frequency of "0/0", "1/1", "1/0", "0/1" and "./."
                    snp_type_counts = Counter(samples_info_items)
                    snp_type_out = []
                    for each_snp in snp_out_format:
                        if each_snp in snp_type_counts:
                            snp_type_out.append(str(snp_type_counts[each_snp]))
                        else:
                            snp_type_out.append("0")
                    # replace "0/0" with "0", "0/1" with "1", "1/1" with "2", "1/0" with "1", "./." with "."
                    samples_info_replace = list(map(replace, samples_info_items))
                    for each_gene in this_line_gene_trans_out:
                        gene, ens = each_gene.split("_")
                        if all([gene, ens]):
                            output_file.write("{}\t{}\t{}\t{}\t{}\t{}\n".format(f"{items[0].split(':')[-1]}_{items[1]}_{items[3]}_{items[4]}", gene, ens, "/".join(snp_type_out), this_line_gene_trans_out[each_gene], "\t".join(samples_info_replace)))
    output_file.close()
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Finish writing to {path_out}")
