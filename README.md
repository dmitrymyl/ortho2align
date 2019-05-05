# ortho2align
Sequence alignment tool based on syntenic protein neighbourhood derived from OrthoDB.

# What is it?
`ortho2align` is the package for alignmnent of nucleotide sequences from one species against another species' genome. It uses OrthoDB data of orthologous proteins to construct regions of quasi-synteny to utilise them as guides for alignment.

# For what tasks this package can be used?
* general alignment of any sequences from the query genome to the target genome;
* alignment of low-conservative sequences to the target genome (ncRNAs etc.).

# Dependencies
The package was tested under following dependencies:
* python 3.6
* pandas 0.24?
* blat v. 35
* bedtools 2.25?
* OrthoDB v10

# How to run
TODO: add sample code to run:
1. Getting orthodb map file.
2. Running ortho2align.py

# How does it work?
Input files are:
* genome of query species
* genome of target species
* protein annotation of query species
* genome annotation of target species
* OrthoDB map file
* coordinates of query sequences in the genome of query species

First, the packages inferes protein neighbourhood of query sequences in the genome of query species at the given radius with `bedtools window -w radius`. Next, it retrieves orthologs of neighbouring proteins in the target genome and construct quasi-syntenic regions. In case orthologous protein in the target genome are placed within merge distance, then they will be merged. After there are no possibilities to merge, derived syntenic ranges can be flanked to some extent. This might be helpful in case protein neihbourhood of one query sequence contain only one protein. There might be more than one quasi-syntenic regions for one query sequence due to paralogues.

Then query sequences are aligned against their syntenic regions with BLAT. User can define tile size and minimal identity of sequences to report. BLAT was chosen for its convenient psl3 format of alignments, that provides exon-intron-like structure of aligned regions. Alignment of many sequences can take a lot of time so user can specify how many cores can be used for alignment process.

Main output file is a `json` containing array of dictionary records, one query sequence per record. Each record contains information about query sequence,
protein neighbourhood, syntenic regions and found alignments. Intermediate files are created in the working directory. THe main working format is json due to its convenient interoperability with python data structures.

# Structure of the package
All scripts produce json files.
`get_neighbourhood.py` retrieves protein neighbourhood of query sequences in query genome based on supplied protein annotation in gtf|gff format.
`extract_mapping.py` extracts orthodb mapping data from bulk file for one species.
`chromsizes_fasta.py` gets chromosome sizes of given genome in fasta format and returns them in json format.
`annotation2json.py`translates given genome annotation of target species into json format as a list of dictionaries.
`map_synteny.py`maps protein neighbourhood to orthologs in the target species and compose syntenic ranges.
`get_fasta.py` retrieves query and syntenic target sequences one per file from given query and target genomes.
`grid_alignment.py`performs alignment of query sequences to target syntenies.
`ortho2align.py` master script to ~~rule them all~~ run listed above scripts in sequential manner. All output files produced by the scripts have fixed names, so user can run each step separately as long as one follows naming conventions.
TODO: add image of data flow within scripts.

# TODO
* complete README.md
* add synteny map example
* add orthodb file processing
* add examples folder