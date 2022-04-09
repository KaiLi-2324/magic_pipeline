from libc.stdlib cimport *
from libc.stdio cimport *
from cpython cimport bool
from libcpp.string cimport string
from libcpp.unordered_map cimport unordered_map
from libc.string cimport strcmp, strtok, strcspn, strcmp


cdef class CReadMatrix:
    cdef:
        str path
        FILE *fp
        bool closed
    cpdef CReadMatrix open(self)
    cpdef bytes readline(self)
    cpdef void close(self)


cdef class CCountVariantsNum:
    cdef:
        str path
        str path_out
        str path_pathway
        public unordered_map[string, int] samples_variants_count
        public unordered_map[string, int] pathway_genes
    cpdef int read_in_pathway(self) except +
    cpdef int count_snp_matrix(self) except +
    cpdef int count_snp_matrix_with_pathway(self) except +


