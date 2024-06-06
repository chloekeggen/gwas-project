from setuptools import setup, find_packages
from linear_gwas.__version__ import __version__

setup(
    name='linear_gwas',
    version=__version__,
    author='Ieva Sereiva, Chloe Keggen, Manaswini Kanaparthy',
    author_email='isereiva@ucsd.edu',
    description='Simple GWAS on continuous phenotypes performing association analysis via linear regression',
    packages=find_packages(),
    install_requires=[
    'pandas~=2.0.1',
    'numpy~=1.24.1',
    'matplotlib~=3.7.1',
    'statsmodels~=0.14.0',
    'pyvcf3~=1.0.0'
    ],
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: Linux',
    ],
    entry_points={
        "console_scripts": [
            "gwas-tools-cli=linear_gwas.gwas_tools_cli:main"
        ],
    },
)


