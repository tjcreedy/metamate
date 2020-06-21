import setuptools

with open("README.md", "r") as fh:
    long_description = fh.read()

setuptools.setup(
    name="NUMTdumper-tjcreedy", # Replace with your own username
    version="0.1b",
    author="Thomas J. Creedy",
    author_email="thomas.creedy@gmail.com",
    description="NUMTdumper: validated NUMT removal for mitochondrial "
                "metabarcoding",
    long_description=long_description,
    long_description_content_type="text/markdown",
    url="https://github.com/tjcreedy/numtdumper",
    packages=setuptools.find_packages(),
    classifiers=[
        "Programming Language :: Python :: 3",
        "Programming Language :: R",
        "License :: OSI Approved :: GNU General Public License (GPL)",
        "Development Status :: 4 - Beta",
        "Environment :: Console",
        "Intended Audience :: Science/Research",
        "Operating System :: OS Independent",
        "Topic :: Scientific/Engineering :: Bio-Informatics"
    ],
    python_requires='>=3.6',
    package_data={
        'Rscripts':[
            'numtdumper/getclades.R',
            'numtdumper/maketree.R'
            ],
        'specifications':[
            'specifications.txt'
            ]
    },
    entry_points={
        "console_scripts":[
            "numtdumper = numtdumper.numtdumper:main",
            "filtertranslate = numtdumper.filtertranslate:main"
            #"filterlength = numtdumper.filterlength:main"
            #"filterreference = numtdumper.filterreference:main"
            ]}
)

#TODO add url to published paper