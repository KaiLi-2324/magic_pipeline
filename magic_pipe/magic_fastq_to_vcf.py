import os
import time
import subprocess
import multiprocessing
from .local_exceptions import *


def create_dir(path):
    if not os.path.exists(path):
        os.makedirs(path)
    return path


def create_mark_file(temp_dir, file_name):
    if not os.path.isdir(temp_dir):
        os.makedirs(temp_dir)
    mark = f"{temp_dir}/{file_name}"
    open(mark, "w+").close()


def create_sym_link(src_file, dst_dir):
    file_name = os.path.basename(src_file)
    dst_file = f"{dst_dir}/{file_name}"
    if not os.path.exists(dst_file):
        os.symlink(src=src_file, dst=dst_file)
    else:
        os.remove(dst_file)
        os.symlink(src=src_file, dst=dst_file)
    return dst_file


def parse_sample_list(path_sample_list):
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Start reading from {path_sample_list}")
    sample_fastqs = []
    with open(path_sample_list) as file:
        for line in file:
            items = line.strip().split("\t")
            sample_fastqs.append(items)
    return sample_fastqs


def index_ref(path_ref, path_mills, path_dbsnp, path_out, path_samtools, java, path_picard, path_bgzip, path_tabix):
    ref_dir = create_dir(f"{path_out}/ref")
    path_ref_link = create_sym_link(src_file=path_ref, dst_dir=ref_dir)
    for each_type in [".amb", ".ann", ".bwt", ".pac", ".sa"]:
        if not os.path.exists(f"{path_ref}{each_type}"):
            raise FileNotFoundError(f"BWA index {path_ref}{each_type} for {path_ref} does not exist! Please run command"
                                    f": bwa index {path_ref}")
        else:
            create_sym_link(src_file=f"{path_ref}{each_type}", dst_dir=ref_dir)
    path_log = create_dir(path=f"{path_out}/log")
    path_temp = create_dir(path=f"{path_out}/temp")
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Indexing fasta {path_ref_link}")
    cmd_index_fasta = f"{path_samtools} faidx {path_ref_link} > {path_log}/index.fasta.log 2>&1"
    subprocess.check_output(cmd_index_fasta, shell=True)
    if not os.path.exists(f"{path_ref_link}.dict"):
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Creating dict for {path_ref_link}")
        cmd_dict_fasta = f"{java} -Djava.io.tmpdir={path_temp} -jar {path_picard} CreateSequenceDictionary R={path_ref_link} O={path_ref_link}.dict > {path_log}/dict.fasta.log 2>&1"
        subprocess.check_output(cmd_dict_fasta, shell=True)
    path_mills_link = create_sym_link(src_file=path_mills, dst_dir=ref_dir)
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: Indexing vcf {path_mills_link}")
    if not path_mills_link.endswith(".gz"):
        path_mills_return = f"{path_mills_link}.gz"
        subprocess.check_output(f"{path_bgzip} -f {path_mills_link}", shell=True)
        subprocess.check_output(f"{path_tabix} -p vcf -f {path_mills_link}.gz", shell=True)
    else:
        path_mills_return = path_mills_link
        subprocess.check_output(f"{path_tabix} -p vcf -f {path_mills_link}", shell=True)
    path_dbsnp_link = create_sym_link(src_file=path_dbsnp, dst_dir=ref_dir)
    print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]:: Indexing vcf {path_dbsnp_link}")
    if not path_dbsnp_link.endswith(".gz"):
        path_dbsnp_return = f"{path_dbsnp_link}.gz"
        if not os.path.exists(f"{path_dbsnp_link}.gz"):
            subprocess.check_output(f"{path_bgzip} -f {path_dbsnp_link}", shell=True)
        if not os.path.exists(f"{path_dbsnp_link}.gz.tbi"):
            subprocess.check_output(f"{path_tabix} -p vcf -f {path_dbsnp_link}.gz", shell=True)
    else:
        path_dbsnp_return = path_dbsnp_link
        if not os.path.exists(f"{path_dbsnp_link}.tbi"):
            subprocess.check_output(f"{path_tabix} -p vcf -f {path_dbsnp_link}", shell=True)
    return path_ref_link, path_mills_return, path_dbsnp_return


def get_trim_out(path_out, sample):
    this_sample_trim_out = f"{path_out}/clean/{sample}"
    if not os.path.exists(this_sample_trim_out):
        os.makedirs(this_sample_trim_out)
    read1_clean_out = f"{this_sample_trim_out}/{sample}_R1.fastq.gz"
    read2_clean_out = f"{this_sample_trim_out}/{sample}_R2.fastq.gz"
    read1_unpaired_out = f"{this_sample_trim_out}/{sample}_R1.unpaired.fastq.gz"
    read2_unpaired_out = f"{this_sample_trim_out}/{sample}_R2.unpaired.fastq.gz"
    return read1_clean_out, read2_clean_out, read1_unpaired_out, read2_unpaired_out


class SampleReRun:
    trim = True
    map = True
    sort = True
    markdup = True
    recal = True
    bqsr = True
    index = True
    haplotype = True

    def run_all(self):
        self.trim = False
        self.map = False
        self.sort = False
        self.markdup = False
        self.recal = False
        self.bqsr = False
        self.index = False
        self.haplotype = False


class SampleFastqToVcf(SampleReRun):
    def __init__(self, sample, read1, read2, java, threads, path_out, path_adapter, path_trim, path_bwa, path_ref,
                 path_sambamba, path_probe, path_gatk, path_mills, path_dbsnp, path_samtools, force):
        self.sample = sample
        self.read1 = read1
        self.read2 = read2
        self.java = java
        self.path_trim = path_trim  # path_trim is path to trimmomatic
        self.path_bwa = path_bwa  # path_bwa is path to bwa
        self.path_sambamba = path_sambamba  # path_picard is path to picard tools
        self.threads = int(threads)
        self.adapter = path_adapter
        self.path_ref = path_ref
        self.path_probe = path_probe
        self.path_gatk = path_gatk
        self.path_mills = path_mills
        self.path_dbsnp = path_dbsnp
        self.path_samtools = path_samtools
        self.path_temp = create_dir(path=f"{path_out}/temp/{sample}")
        self.path_log = create_dir(path=f"{path_out}/log/{sample}")
        self.path_mapping = create_dir(path=f"{path_out}/mapping/{sample}")
        self.path_vcf = create_dir(path=f"{path_out}/vcf/{sample}")
        self.read1_clean_out, self.read2_clean_out, self.read1_unpaired_out, self.read2_unpaired_out = get_trim_out(path_out=path_out, sample=sample)
        if force != "none":
            if force == "all":
                self.run_all()
            else:
                setattr(self, force, False)
        self.analyze()

    def analyze(self):
        if os.path.isfile(f"{self.path_temp}/trim.done") and self.trim:
            print(f" Trimming for {self.sample} has been done, skipping it")
        else:
            print(f" Trimming: {self.sample}")
            self.trimming()
        if os.path.isfile(f"{self.path_temp}/mapping.done") and self.map:
            print(f" Mapping for {self.sample} has been done, skipping it")
        else:
            print(f" Aligning: {self.sample}")
            self.mapping()
        if os.path.isfile(f"{self.path_temp}/sorting.done") and self.sort:
            print(f" Sorting for {self.sample} has been done, skipping it")
        else:
            print(f" Sorting bam: {self.sample}")
            self.sort_bam()
        if os.path.isfile(f"{self.path_temp}/markdup.done") and self.markdup:
            print(f" Marking dup for {self.sample} has been done, skipping it")
        else:
            print(f" Marking duplicates: {self.sample}")
            self.mark_dup()
        if os.path.isfile(f"{self.path_temp}/baserecal.done") and self.recal:
            print(f" Base Recalibrating for {self.sample} has been done, skipping it")
        else:
            print(f" Base Recalibrating: {self.sample}")
            self.base_recalibrator()
        if os.path.isfile(f"{self.path_temp}/applybqsr.done") and self.bqsr:
            print(f" Applying BQSR for {self.sample} has been done, skipping it")
        else:
            print(f" Applying BQSR: {self.sample}")
            self.apply_bqsr()
        if os.path.isfile(f"{self.path_temp}/index.done") and self.index:
            print(f" Indexing bam for {self.sample} has been done, skipping it")
        else:
            print(f" Indexing bam: {self.sample}")
            self.index_bam()
        if os.path.isfile(f"{self.path_temp}/haplotype.done") and self.haplotype:
            print(f" Calling Haplotype for {self.sample} has been done, skipping it")
        else:
            print(f" Calling Haplotype: {self.sample}")
            self.haplotype_caller()

    def trimming(self):
        cmd = f"{self.java} -Djava.io.tmpdir={self.path_temp} -jar {self.path_trim} PE -threads {self.threads} {self.read1} {self.read2} {self.read1_clean_out} {self.read1_unpaired_out} {self.read2_clean_out} {self.read2_unpaired_out} LEADING:15 TRAILING:15 SLIDINGWINDOW:5:20 AVGQUAL:20 MINLEN:36 ILLUMINACLIP:{self.adapter}:2:20:10:8:keepBothReads > {self.path_log}/{self.sample}.trim.log 2>&1"
        subprocess.check_output(cmd, shell=True)
        create_mark_file(temp_dir=self.path_temp, file_name="trim.done")

    def mapping(self):
        cmd = f'{self.path_bwa} mem -M -t {self.threads} -R "@RG\\tID:{self.sample}\\tSM:{self.sample}\\tLB:{self.sample}\\tPL:illumina" {self.path_ref} {self.read1_clean_out} {self.read2_clean_out} | {self.path_samtools} view -bS -F 4 -q 20 - -o {self.path_mapping}/{self.sample}.bam'
        subprocess.check_output(cmd, shell=True)
        create_mark_file(temp_dir=self.path_temp, file_name="mapping.done")

    def sort_bam(self):
        cmd = f"{self.path_sambamba} sort -t {self.threads} {self.path_mapping}/{self.sample}.bam {self.sample}.sort"
        subprocess.check_output(cmd, shell=True)
        create_mark_file(temp_dir=self.path_temp, file_name="sort.done")

    def mark_dup(self):
        cmd = f"{self.path_sambamba} markdup -t {self.threads} {self.path_mapping}/{self.sample}.sorted.bam {self.path_mapping}/{self.sample}.sort.rmdup.bam > {self.path_log}/{self.sample}.mkdup.log 2>&1"
        subprocess.check_output(cmd, shell=True)
        create_mark_file(temp_dir=self.path_temp, file_name="markdup.done")

    def base_recalibrator(self):
        cmd = f"{self.java} -Djava.io.tmpdir={self.path_temp} -jar {self.path_gatk} -T BaseRecalibrator -R {self.path_ref} -I {self.path_mapping}/{self.sample}.sort.rmdup.bam -L {self.path_probe} --known-sites {self.path_dbsnp} --known-sites {self.path_mills} -o {self.path_mapping}/{self.sample}.table > {self.path_log}/{self.sample}.recalibrator.log 2>&1"
        subprocess.check_output(cmd, shell=True)
        create_mark_file(temp_dir=self.path_temp, file_name="baserecal.done")

    def apply_bqsr(self):
        if os.path.exists(f"{self.path_mapping}/{self.sample}.bam"):
            os.remove(f"{self.path_mapping}/{self.sample}.bam")
        cmd = f"{self.java} -Djava.io.tmpdir={self.path_temp} -jar {self.path_gatk} -T ApplyBQSR -R {self.path_ref} -I {self.path_mapping}/{self.sample}.sort.rmdup.bam -bqsr {self.path_mapping}/{self.sample}.table -o {self.path_mapping}/{self.sample}.bam > {self.path_log}/{self.sample}.bqsr.log 2>&1"
        subprocess.check_output(cmd, shell=True)
        create_mark_file(temp_dir=self.path_temp, file_name="applybqsr.done")

    def index_bam(self):
        cmd = f"{self.path_sambamba} index -t {self.threads} {self.path_mapping}/{self.sample}.bam"
        subprocess.check_output(cmd, shell=True)
        create_mark_file(temp_dir=self.path_temp, file_name="index.done")

    def haplotype_caller(self):
        cmd = f"{self.java} -Djava.io.tmpdir={self.path_temp} -jar {self.path_gatk} -T HaplotypeCaller -R {self.path_ref} -I {self.path_mapping}/{self.sample}.bam -o {self.path_vcf}/{self.sample}.g.vcf > {self.path_log}/{self.sample}.hapcall.log 2>&1"
        subprocess.check_output(cmd, shell=True)
        create_mark_file(temp_dir=self.path_temp, file_name="haplotype.done")


def fastq_to_vcf_each_sample(sample, read1, read2, java, threads, path_out, path_adapter, path_trim,
                             path_bwa, path_ref, path_sambamba, path_probe, path_gatk, path_mills,
                             path_dbsnp, path_samtools, mp_list, force):
    try:
        each_sample_fastq_to_vcf = SampleFastqToVcf(sample=sample,
                                                    read1=read1,
                                                    read2=read2,
                                                    java=java,
                                                    threads=threads,
                                                    path_out=path_out,
                                                    path_adapter=path_adapter,
                                                    path_trim=path_trim,
                                                    path_bwa=path_bwa,
                                                    path_ref=path_ref,
                                                    path_sambamba=path_sambamba,
                                                    path_probe=path_probe,
                                                    path_gatk=path_gatk,
                                                    path_mills=path_mills,
                                                    path_dbsnp=path_dbsnp,
                                                    path_samtools=path_samtools,
                                                    force=force)
        mp_list.append((sample, "done"))
    except Exception as e:
        mp_list.append((sample, "error", str(e)))


def fastq_to_vcf_parallel(**kwargs):
    path_ref, path_mills, path_dbsnp = index_ref(path_ref=kwargs["path_ref"],
                                                 path_mills=kwargs["path_mills"],
                                                 path_dbsnp=kwargs["path_dbsnp"],
                                                 path_out=kwargs["path_out"],
                                                 path_samtools=kwargs["path_samtools"],
                                                 java=kwargs["java"],
                                                 path_picard=kwargs["path_picard"],
                                                 path_bgzip=kwargs["path_bgzip"],
                                                 path_tabix=kwargs["path_tabix"])
    sample_fastqs = parse_sample_list(path_sample_list=kwargs["path_sample_list"])
    pool = multiprocessing.Pool(processes=int(kwargs["num_process"]))
    mp_list = multiprocessing.Manager().list()
    for each_item in sample_fastqs:
        sample, read1, read2 = each_item
        pool.apply_async(func=fastq_to_vcf_each_sample, args=(sample,
                                                              read1,
                                                              read2,
                                                              kwargs["java"],
                                                              kwargs["threads"],
                                                              kwargs["path_out"],
                                                              kwargs["path_adapter"],
                                                              kwargs["path_trim"],
                                                              kwargs["path_bwa"],
                                                              path_ref,
                                                              kwargs["path_sambamba"],
                                                              kwargs["path_probe"],
                                                              kwargs["path_gatk"],
                                                              path_mills,
                                                              path_dbsnp,
                                                              kwargs["path_samtools"],
                                                              mp_list,
                                                              kwargs["force"],))
    pool.close()
    pool.join()

    error_info = []
    if not mp_list:
        raise EmptyResultError("No results got!")
    for each_item in mp_list:
        if each_item[1] == "error":
            error_info.append(f"{each_item[0]}\t{each_item[2]}")
    if error_info:
        path_error = f'{kwargs["path_out"]}/error.txt'
        with open(path_error, "w+") as error:
            error.write("\n".join(error_info) + "\n")
        raise UnExpectedError(f"Error occurred while parsing the input data, error information has been written "
                              f"to {path_error}")
    else:
        print(f"[{time.strftime('%Y-%m-%d %H:%M:%S', time.localtime())}]: All jobs have been finished")
