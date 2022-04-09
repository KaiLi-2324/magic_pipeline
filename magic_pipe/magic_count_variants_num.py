import os
import time
import argparse
from .cy_vcf_utils import CCountVariantsNum


def main():
    parser = argparse.ArgumentParser(description="This script is designed to stat samples with variant numbers")
    parser.add_argument("-p", help="Path to input matrix with columns representing samples and rows representing snps",
                        action="store", dest="p", type=str, required=True)
    parser.add_argument("-o", help="Path to put output file with two columns, one for sample and the other for "
                                   "variant nums", action="store", dest="o", type=str, required=True)
    parser.add_argument("-w", help="Path to a gene table of a certain pathway, must be one column and each row stands "
                                   "for a gene, default is not provided", action="store", dest="w", type=str,
                        default="")
    args = parser.parse_args()
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


if __name__ == '__main__':
    main()
