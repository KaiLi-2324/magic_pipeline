# distutils: language=c++
cimport cython
from .cy_vcf_utils cimport *
from .local_exceptions import SamplesNotFoundError
# from cython.operator cimport dereference as deref, preincrement as preinc


cdef class CReadMatrix:

    @cython.embedsignature(True)
    def __cinit__(self, str path):
        self.path = path
        self.fp = NULL
        self.closed = False
        self.open()

    def __enter__(self):
        return self

    def __iter__(self):
        if self.fp == NULL:
            self.open()
        return self._iterate_()

    def __exit__(self, exc_type, exc_val, exc_tb):
        if not self.closed:
            self.close()

    def _iterate_(self):
        cdef:
            char *line = NULL
            size_t length = 0
            ssize_t read = 0

        while read != -1:
            read = getline(&line, &length, self.fp)
            if read != -1:
                yield line

    cpdef CReadMatrix open(self):
        cdef char *file = NULL
        file_name = self.path.encode("UTF-8")
        file = file_name
        self.fp = fopen(file, "rb")
        if self.fp == NULL:
            raise IOError("Failed to open file {}".format(self.path))
        return self


    cpdef bytes readline(self):
        # the getline method of libc will return a C/C++ string (Python bytes) not Python string
        cdef:
            char *line = NULL
            size_t length = 0
            ssize_t read

        read = getline(&line, &length, self.fp)
        if read == -1:
            return None
        return line

    cpdef void close(self):
        cdef int errno
        if not self.closed:
            errno = fclose(self.fp)
            if errno != 0:
                raise IOError("Failed to close file {}".format(self.path))
            else:
                self.closed = True


cdef class CCountVariantsNum:

    @cython.embedsignature(True)
    def __cinit__(self, str path, str path_out, str path_pathway):
        self.path = path
        self.path_out = path_out
        self.path_pathway = path_pathway

    cpdef int read_in_pathway(self) except +:
        cdef:
            bytes line, gene
        if not self.path_pathway:
            raise ValueError("pathway table is not provided!")
        with CReadMatrix(self.path_pathway) as pathway:
            for line in pathway:
                gene = line.strip()
                if gene:
                    self.pathway_genes[gene] = 1

    cpdef int count_snp_matrix(self) except +:
        cdef:
            int count = 0
            list samples, items
            bytes line, each_item

        with CReadMatrix(self.path) as file:
            line = file.readline()
            samples = line.strip().split(b"\t")[5:]
            if not samples:
                raise SamplesNotFoundError("No samples found in {}".format(self.path))
            line = file.readline()
            while line:
                items = line.rstrip().split(b"\t")[5:]
                for each_item in items:
                    if strcmp(each_item, "1") == 0:
                        self.samples_variants_count[samples[count]] += 1
                    elif strcmp(each_item, "2") == 0:
                        self.samples_variants_count[samples[count]] += 1
                    count += 1
                line = file.readline()
                count = 0
        return 0

    cpdef int count_snp_matrix_with_pathway(self) except +:
        cdef:
            int count = 0
            list samples, items, total_items
            bytes line, each_item, gene
        self.read_in_pathway()
        with CReadMatrix(self.path) as file:
            line = file.readline()
            samples = line.strip().split(b"\t")[5:]
            if not samples:
                raise SamplesNotFoundError("No samples found in {}".format(self.path))
            line = file.readline()
            while line:
                total_items = line.rstrip().split(b"\t")
                items = total_items[5:]
                gene = total_items[1]
                if self.pathway_genes[gene] == 1:
                    for each_item in items:
                        if strcmp(each_item, "1") == 0:
                            self.samples_variants_count[samples[count]] += 1
                        elif strcmp(each_item, "2") == 0:
                            self.samples_variants_count[samples[count]] += 1
                        count += 1
                line = file.readline()
                count = 0
        return 0
