# Download and unpack human hg38 genome.
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz
gunzip hg38.fa.gz
# Download and unpack mouse mm10 genome.
wget https://hgdownload.soe.ucsc.edu/goldenPath/mm10/bigZips/mm10.fa.gz
gunzip mm10.fa.gz
# Download liftOver chain file.
wget https://hgdownload.soe.ucsc.edu/goldenPath/hg38/liftOver/hg38ToMm10.over.chain.gz
# Download mouse gencode vM25 annotation file.
wget http://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_mouse/release_M25/gencode.vM25.annotation.gff3.gz
gunzip gencode.vM25.annotation.gff3.gz
grep -P $"\tgene\t" gencode.vM25.annotation.gff3 > gencode.vM25.annotation.genes.gff3
# Run ortho2align
ortho2align @test.config
