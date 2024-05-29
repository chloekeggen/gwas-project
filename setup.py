from setuptools import setup, find_packages
from linear_gwas.__version__ import __version__

setup(
    name='linear_gwas',
    version=__version__,
    author='Ieva Sereiva',
    author_email='isereiva@ucsd.edu',
    description='Simple GWAS on continuous phenotypes performing association analysis via linear regression',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'matplotlib',
        'statsmodels',
        'cyvcf2',
        'scikit-learn'
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Linux',
    ],
    entry_points={
        "console_scripts": [
            "gwas-tools-cli=linear_gwas.gwas_tools_cli:main",
            "gwas-plots-cli=linear_gwas.gwas_plots_cli:main"
        ],
    },
)

