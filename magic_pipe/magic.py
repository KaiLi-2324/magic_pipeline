import re
import os
import time
import argparse
import subprocess
import multiprocessing
from .local_exceptions import VersionError


def require_args(args_list):
    # warn the user to input the required args
    class RequireArgs(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if values not in args_list:
                msg = "Args of {} must be chosen from {}".format(f"-{self.dest}", ",".join(args_list))
                raise argparse.ArgumentTypeError(msg)
            setattr(args, self.dest, values)
    return RequireArgs


def require_argl():
    class RequireArgs(argparse.Action):
        def __call__(self, parser, args, values, option_string=None):
            if type(values) != int:
                msg = "Args of {} must be integer".format(f"-{self.dest}")
                raise argparse.ArgumentTypeError(msg)
            elif values < 1:
                msg = "Args of {} must be greater than or equal to 1".format(f"-{self.dest}")
                raise argparse.ArgumentTypeError(msg)
            else:
                cpu_count = multiprocessing.cpu_count()
                if values > cpu_count:
                    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Local cpu count: {cpu_count}, "
                          f"processes user provide: {values} Adjusting to {cpu_count}")
                    # values = cpu_count
            setattr(args, self.dest, values)

    return RequireArgs


def _fastq_to_vcf_(args):
    from .magic_fastq_to_vcf import fastq_to_vcf_parallel
    from magic_pipe import tools, data
    cpu_counts = multiprocessing.cpu_count()
    cpu_to_be_use = cpu_counts // 2
    if args.process * args.threads > cpu_counts:
        print(f"[WARNING]: Allocated threads plus processes have exceeded the cpu count of this machine, adjusting "
              f"threads to {cpu_to_be_use} and processes to {cpu_to_be_use}")
        args.process = cpu_to_be_use
        args.threads = cpu_to_be_use
    path_tools = tools.__path__[0]
    if not args.gatk:
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: GATK not provided, using the in-built one "
              f"instead")
        args.gatk = f"{path_tools}/GenomeAnalysisTK.jar"
        os.chmod(args.gatk, int("755", base=8))
    if not args.picard:
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Picard not provided, using the in-built one "
              f"instead")
        args.picard = f"{path_tools}/picard.jar"
        os.chmod(args.picard, int("755", base=8))
    if not args.trim:
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Trimmomatic not provided, using the in-built "
              f"one instead")
        args.trim = f"{path_tools}/trimmomatic-0.39.jar"
        os.chmod(args.trim, int("755", base=8))
    if not args.sambamba:
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Sambamba not provided, using the in-built one"
              f" instead")
        args.sambamba = f"{path_tools}/sambamba-0.8.0"
        os.chmod(args.sambamba, int("755", base=8))
    if not args.path_adapter:
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Adapters not provided, using the in-built one"
              f" instead")
        path_data = data.__path__[0]
        args.path_adapter = f"{path_data}/adapters4qc.fasta"
    for each_item in [args.trim, args.java, args.bwa, args.gatk, args.samtools, args.sambamba, args.path_adapter]:
        if not os.path.isfile(each_item):
            raise FileNotFoundError(f"Could not find {each_item}. Please provide the absolute path")
    pattern_gatk_version = re.compile(r"^3\.\d-\d-.*")
    gatk_version = subprocess.check_output(f"{args.java} -jar {args.gatk} --version", shell=True).decode()
    if not pattern_gatk_version.search(gatk_version):
        raise VersionError(f"Only GATK 3 is supported! Got {gatk_version.strip()}")
    fastq_to_vcf_parallel(path_sample_list=args.path_sample_list,
                          num_process=args.process,
                          java=args.java,
                          threads=args.threads,
                          path_out=args.path_out,
                          path_adapter=args.path_adapter,
                          path_trim=args.trim,
                          path_bwa=args.bwa,
                          path_ref=args.path_ref,
                          path_sambamba=args.sambamba,
                          path_probe=args.path_probe,
                          path_gatk=args.gatk,
                          path_mills=args.path_mills,
                          path_dbsnp=args.path_dbsnp,
                          path_samtools=args.samtools,
                          path_picard=args.picard,
                          path_bgzip=args.bgzip,
                          path_tabix=args.tabix,
                          force=args.force)


def _vcf_to_matrix_(args):
    import portion
    from .magic_vcf_to_matrix import parse_vcf
    trans_type = args.t
    variant_type = args.v
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: trans type {trans_type}  variant type "
          f"{variant_type}")
    if args.i == "-inf":
        i = -portion.inf
    else:
        i = float(args.i)
    if args.I == "inf":
        I = portion.inf
    else:
        I = float(args.I)
    parse_vcf(path_vcf=args.path_vcf,
              i=i,
              I=I,
              trans_type=trans_type,
              variant_type=variant_type,
              e=float(args.e),
              path_out=args.path_out)


def _fisher_perm_(args):
    from .magic_fisher_permutation import parse_snp_matrix
    path_pathway = args.w if args.w else ""
    if not os.path.exists(os.path.dirname(args.o)):
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Creating directory {os.path.dirname(args.o)}")
        os.makedirs(os.path.dirname(args.o))
    parse_snp_matrix(path_matrix=args.p,
                     path_case_control=args.c,
                     path_pathway=path_pathway,
                     method=args.m,
                     path_out=args.o,
                     num_process=args.l,
                     progress_bar=args.g,
                     hypothesis=args.s)


def _skat_acat_(args):
    from .magic_run_skat_acat import get_case_control_samples, parse_covariate, read_in_snp_matrix, calculate_skat
    case_samples_num, control_samples_num = get_case_control_samples(args.c)
    if args.X:
        covariate = parse_covariate(args.X, case_samples_num, control_samples_num, args.c)
    else:
        covariate = None
    genes_snp_genotypes = read_in_snp_matrix(args.p)
    calculate_skat(genes_snp_genotypes, case_samples_num, control_samples_num, args.o, covariate)


def _count_variant_(args):
    from .cy_vcf_utils import CCountVariantsNum
    if not os.path.exists(os.path.dirname(args.o)):
        os.makedirs(os.path.dirname(args.o))
    count_variant_num = CCountVariantsNum(args.p, args.o, args.w)
    if args.w:
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Counting snp numbers for all samples with pathway {args.w}")
        count_variant_num.read_in_pathway()
        count_variant_num.count_snp_matrix_with_pathway()
    else:
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Counting snp numbers for all samples without pathway table")
        count_variant_num.count_snp_matrix()


def main():
    parser = argparse.ArgumentParser(prog='magic',
                                     description="A scalable pipeline for WES/WGS data analysis\n"
                                                 "magic contains five modules: \n"
                                                 "\t%(prog)s fastq_to_vcf: Analyzing the raw fastq format files to "
                                                 "Variant Call Format files. This step includes trimming, mapping and"
                                                 "haplotype calling\n"
                                                 "\t%(prog)s vcf_to_matrix: Extract biological features and information"
                                                 " from the vcf file generated by GATK\n"
                                                 "\t%(prog)s fisher_perm: Perform fisher test and permutation test for"
                                                 "genes extracted from the vcf_to_matrix command\n"
                                                 "\t%(prog)s skat_acat: Perform SKAT and ACAT analysis for genes "
                                                 "extracted from the vcf_to_matrix command\n"
                                                 "\t%(prog)s count_variant: Count genes with snp for each sample\n",
                                     formatter_class=argparse.RawTextHelpFormatter)
    subparsers = parser.add_subparsers(title="modules", help='magic modules, use -h/--help for help')
    sub_fastq_to_vcf = subparsers.add_parser("fastq_to_vcf", description="Perform trimming, mapping and haplotype "
                                                                         "calling for the raw fastq files")
    sub_vcf_to_matrix = subparsers.add_parser("vcf_to_matrix", description="Convert the vcf generated by GATK to a"
                                                                           "matrix, each row stands for a snp")
    sub_fisher_perm = subparsers.add_parser("fisher_perm", description="Perform fisher test and permutation test for"
                                                                       "each gene generated by vcf_to_matrix")
    sub_skat_acat = subparsers.add_parser("skat_acat", description="Perform SKAT and ACAT analysis for each gene "
                                                                   "gneraed by vcf_to_matrix")
    sub_count_variant = subparsers.add_parser("count_variant", description="Count genes with snp for each sample")

    # fastq to vcf
    sub_fastq_to_vcf.add_argument("-p", help="A file which contains path to fastq files for each sample, each row "
                                             "represents a sample and MUST have three columns with tab delimited, the "
                                             "first column is sample name, the second is path to fastq read1 and the "
                                             "third is path to fastq read2. eg.\nsample\tsample_read1\tsample_read2\n"
                                             "sample1\tsample1_read1\tsample1_read2\n...",
                                  action="store", dest="path_sample_list", type=str, required=True)
    sub_fastq_to_vcf.add_argument("-r", help="Path to reference fasta file, such as hg38.fa/GRCh38.fa", action="store",
                                  dest="path_ref", type=str, required=True)
    sub_fastq_to_vcf.add_argument("-b", help="Path to probe file", action="store", dest="path_probe", type=str,
                                  required=True)
    sub_fastq_to_vcf.add_argument("-m", help="Path to mills vcf, such as Mills_and_1000G_gold_standard.indels.hg38.vcf."
                                             "gz. Must be provided in the base recalibrating step", action="store",
                                  dest="path_mills", type=str, default=None)
    sub_fastq_to_vcf.add_argument("-d", help="Path to dbsnp vcf, such as dbsnp_146.hg38.vcf.gz", action="store",
                                  dest="path_dbsnp", type=str, default=None)
    sub_fastq_to_vcf.add_argument("-a", help="Path to fasta file which contains adapter sequences", action="store",
                                  dest="path_adapter", type=str, default=None)
    sub_fastq_to_vcf.add_argument("-l", help="Processes to use in the program, default is 5", action="store",
                                  dest="process", type=int, default=5)
    sub_fastq_to_vcf.add_argument("-o", help="Directory to put the output results", action="store", dest="path_out",
                                  type=str, default=5)
    sub_fastq_to_vcf.add_argument("-threads", help="Threads to use in trimmomatic/bwa/GATK/sambamba, default is 8",
                                  action="store", dest="threads", type=int, default=8)
    sub_fastq_to_vcf.add_argument("-trim", help="Path to trimmomatic executable, eg. /usr/bin/trimmomatic-0.38.jar",
                                  action="store", dest="trim", type=str, default=None)
    sub_fastq_to_vcf.add_argument("-java", help="Path to java executable, eg. /usr/bin/java", action="store",
                                  dest="java", type=str, default="/usr/bin/java")
    sub_fastq_to_vcf.add_argument("-bwa", help="Path to bwa executable, eg. /usr/bin/bwa", action="store", dest="bwa",
                                  type=str, required=True)
    sub_fastq_to_vcf.add_argument("-gatk", help="Path to GATK java package, only GATK 3 is supported", action="store",
                                  dest="gatk", type=str, default=None)
    sub_fastq_to_vcf.add_argument("-samtools", help="Path to samtools executable, eg. /usr/bin/samtools",
                                  action="store", dest="samtools", type=str, required=True)
    sub_fastq_to_vcf.add_argument("-sambamba", help="Path to sambamba executable, eg. /usr/bin/sambamba",
                                  action="store", dest="sambamba", type=str, default=None)
    sub_fastq_to_vcf.add_argument("-picard", help="Path to picard java package", action="store", dest="picard",
                                  type=str, default=None)
    sub_fastq_to_vcf.add_argument("-tabix", help="Path to tabix executable, eg. /usr/bin/tabix", action="store",
                                  dest="tabix", type=str, required=True)
    sub_fastq_to_vcf.add_argument("-bgzip", help="Path to bgzip executable, eg. /usr/bin/bgzip", action="store",
                                  dest="bgzip", type=str, required=True)
    sub_fastq_to_vcf.add_argument("-force", help="Force to rerun the steps even if they have been run before, not "
                                                 "skipping them. Must be chosen from trim/map/sort/markdup/recal/bqsr"
                                                 "/index/haplotype/all/none. If all is given, all the steps will be "
                                                 "rerun. If none is given, all steps will  be skipped if they have"
                                                 "been run. Multiple args are accepted, -force trim mapping index "
                                                 "means force rerun the trim mapping and index step if they have been"
                                                 "run before",
                                  choices=["trim", "map", "sort", "markdup", "recal", "bqsr", "index", "haplotype",
                                           "all", "none"], nargs='+', type=str, default="all")
    sub_fastq_to_vcf.set_defaults(func=_fastq_to_vcf_)

    # vcf to matrix
    sub_vcf_to_matrix.add_argument("-p", help="Path to reference file", action="store", dest="path_vcf", type=str,
                                   required=True)
    sub_vcf_to_matrix.add_argument("-i", help="Lower limit to AF value, default is -inf, e.g. 0.001", action="store",
                                   dest="i", type=str, default="-inf")
    sub_vcf_to_matrix.add_argument("-I", help="Upper limit to AF value, default is inf, e.g. 0.01", action="store",
                                   dest="I", type=str, default="inf")
    sub_vcf_to_matrix.add_argument("-t", help="Transcript types to filter the vcf, must be chosen from CANONICAL/CCDS"
                                              "/RefSeq", action=require_args(["CANONICAL", "CCDS", "RefSeq"]), type=str,
                                   required=True)
    sub_vcf_to_matrix.add_argument("-v", help="Variant types to filter the vcf, must be chosen from lof"
                                              "/missense_benign_1/missense_benign_2/missense_damage_1/missense_damage_2"
                                              "/synonymous", action=require_args(["lof", "missense_benign_1",
                                                                                  "missense_benign_2",
                                                                                  "missense_damage_1",
                                                                                  "missense_damage_2",
                                                                                  "synonymous"]), type=str,
                                   required=True)
    sub_vcf_to_matrix.add_argument("-o", help="Path to put the output", action="store", dest="path_out", type=str,
                                   required=True)
    sub_vcf_to_matrix.add_argument("-e", help="Filter for MAX_AF, default is inf", action="store", dest="e", type=str,
                                   default="inf")
    sub_vcf_to_matrix.set_defaults(func=_vcf_to_matrix_)

    # fisher test and permutation test
    sub_fisher_perm.add_argument("-p", help="Path to input matrix with columns representing samples and rows "
                                            "representing snps", action="store", dest="p", type=str, required=True)
    sub_fisher_perm.add_argument("-c", help="Path to input file with case and control samples, must be two columns, "
                                            "first col for samples, the other for case/control", action="store",
                                 dest="c", type=str, required=True)
    sub_fisher_perm.add_argument("-o", help="Path to put output matrix with the following information: gene\tcase\t"
                                            "case_wild\tcontrol\tcontrol_wild\tp_value\tp_adj\tp_permutation",
                                 action="store", dest="o", type=str, required=True)
    sub_fisher_perm.add_argument("-m", help="Method to use in the adjustment of p values, must be chosen from "
                                            "bonferroni/sidak/holm-sidak/holm/simes-hochberg/hommel/fdr_bh/fdr_by"
                                            "/fdr_tsbh/fdr_tsbky, default is fdr_bh",
                                 action=require_args(["bonferroni", "sidak", "holm-sidak", "holm", "simes-hochberg",
                                                     "hommel", "fdr_bh", "fdr_by", "fdr_tsbh", "fdr_tsbky"]),
                                 type=str, default="fdr_bh")
    sub_fisher_perm.add_argument("-l", help="Processes to use in the fisher test and permutation step, default is 10",
                                 action=require_argl(), type=int, default=10)
    sub_fisher_perm.add_argument("-g", help="Whether to show the progress bar, default is off", action="store_true",
                                 default=False)
    sub_fisher_perm.add_argument("-s", help="Alternative hypothesis for fisher exact test, must be chosen from two-"
                                            "sided/less/greater, less and greater is one-sided, default is two-sided",
                                 action=require_args(["two-sided", "less", "greater"]), type=str, default="two-sided")
    sub_fisher_perm.add_argument("-w", help="Path to pathway table, each row stands for a pathway, the 1st and 2nd col "
                                            "are pathway name while the other cols are genes in this pathway. e.g."
                                            "pathway\tpathway\tgene1\tgene2...", action="store", default=None)
    sub_fisher_perm.set_defaults(func=_fisher_perm_)

    # SKAT and ACAT analysis
    sub_skat_acat.add_argument("-p", help="Path to input matrix with columns representing samples and rows representing"
                                          " snps", action="store", dest="p", type=str, required=True)
    sub_skat_acat.add_argument("-c", help="Path to input file with case and control samples, must be two columns, first"
                                          " col for samples, the other for case/control", action="store", dest="c",
                               type=str, required=True)
    sub_skat_acat.add_argument("-X", help="Path to covariates file, default id not provided", action="store", dest="X",
                               type=str, default=None)
    sub_skat_acat.add_argument("-o", help="Path to put output matrix with the following information: "
                                          "gene\tburden_p1\tburden_p2\tskat_p1\tskat_p2\tacatv_p1\tacatv_p2\tskato_p1"
                                          "\tskato_p2\tacato_p", action="store", dest="o", type=str, required=True)
    sub_skat_acat.set_defaults(func=_skat_acat_)

    # count variant numbers for each sample
    sub_count_variant.add_argument("-p", help="Path to input matrix with columns representing samples and rows "
                                              "representing snps", action="store", dest="p", type=str, required=True)
    sub_count_variant.add_argument("-o", help="Path to put output file with two columns, one for sample and the other "
                                              "for variant nums", action="store", dest="o", type=str, required=True)
    sub_count_variant.add_argument("-w", help="Path to a gene table of a certain pathway, must be one column and each "
                                              "row stands for a gene, default is not provided", action="store",
                                   dest="w", type=str, default="")
    sub_count_variant.set_defaults(func=_count_variant_)

    args = parser.parse_args()
    if hasattr(args, "func"):
        args.func(args)
    else:
        parser.print_help()


if __name__ == '__main__':
    main()
