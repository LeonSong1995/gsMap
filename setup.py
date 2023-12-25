from setuptools import setup, find_packages

setup(
    name='GPS',
    version='0.1a',
    packages= find_packages('src'),
    package_dir={'': 'src'},
    url='https://github.com/LeonSong1995/GPS/tree/main',
    license='MIT',
    author='liyang,wenhao',
    author_email='songliyang@westlake.edu.cn',
    description='Genetics-Informed Pathogenic Spatial Mapping for Spatial Transcriptomics',
    long_description='''
         GPS (Genetics-Informed Pathogenic Spatial Mapping) is a Python package 
         for integrating GWAS data with spatial transcriptomics to identify traits-related 
         spatial regions.
     ''',
)
