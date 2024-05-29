from setuptools import setup, find_packages

setup(
    name='gwas_package',
    version='0.1.0',
    packages=find_packages(),
    install_requires=[
        'pandas',
        'numpy',
        'cyvcf2'
    ],
    author='Chloe Keggen A17081400, Ieva Sereiva A16764544, Manaswini Kanaparthy A17169602
',
    author_email='ckeggen@ucsd.edu',
    description='A package for preprocessing GWAS data',
    long_description=open('README.md').read(),
    long_description_content_type='text/markdown',
    url='https://github.com/chloekeggen/gwas-project',
    classifiers=[
        'Programming Language :: Python :: 3',
        'License :: OSI Approved :: MIT License',
        'Operating System :: OS Independent',
    ],
    python_requires='>=3.6',
)
