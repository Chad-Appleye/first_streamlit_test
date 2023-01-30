import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="codx-biotools",
    version="0.0.2",
    description="General CoDx bioinformatics tools",
    long_description=long_description,
    long_description_content_type="text/markdown",
    packages=setuptools.find_packages(),
    include_package_data=True,
    scripts=[
        'bin/combine_fastqataloguer'
    ],
    install_requires=[
        # 'pysam',
        'pandas',
        'httplib2',
        'bs4',
        'numpy',
        'biopython', 
        'plotly', 
        'tqdm',
        'Bio', 
        'xlrd',
        'scipy',
        'sklearn',
        'matplotlib',
        'openpyxl',
        'scikit-allel'
    ]
)