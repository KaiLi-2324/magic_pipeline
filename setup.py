# from distutils.core import setup
import os
from setuptools import setup, find_packages, Extension

try:
    from Cython.Build import cythonize
except ModuleNotFoundError:
    print("[WARNING]: Cython is required to build this package. Tyr to install Cython by: pip install cython")
    exit(-1)


rlib_path_dir = ""
if os.environ['LD_LIBRARY_PATH']:
    lib_path = os.environ['LD_LIBRARY_PATH'].split(":")
else:
    lib_path = []
if lib_path:
    for each_dir in lib_path:
        if "libRblas.so" in os.listdir(each_dir):
            rlib_path_dir = each_dir
if not rlib_path_dir:
    print(f"[WARNING]: libRblas.so not found in LD_LIBRARY_PATH, there might be some problem in installing rpy2")

setup(
    name="magic",
    version="1.0.0",
    packages=find_packages(),
    author="Kai Li",
    python_requires='>=3.7',
    author_email="likai@ucas.ac.cn",
    ext_modules=cythonize([Extension("magic_pipe.cy_vcf_utils",
                                     ["magic_pipe/cy_vcf_utils.pyx"],
                                     extra_compile_args=["-std=c++11"])],
                          compiler_directives={"language_level": 3,
                                               "embedsignature": True}),
    requires=['Cython'],
    install_requires=['portion==2.1.4',
                      'tqdm==4.42.1',
                      'scipy==1.3.3',
                      'statsmodels==0.12.1',
                      'pandas==0.25.3',
                      'numpy==1.19.4',
                      'rpy2==3.4.5'],
    entry_points={'console_scripts': ['magic=magic_pipe.magic:main']},
    zip_safe=False,
    include_package_data=True,
    package_data={'magic_pipe': ['data/adapters4qc.fasta',
                                 "tools/sambamba-0.8.0",
                                 "tools/trimmomatic-0.39.jar",
                                 "tools/GenomeAnalysisTK.jar",
                                 "tools/picard.jar"]}
)
