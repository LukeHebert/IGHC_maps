# Rationale
This repository records the steps taken & tools used to generate a map of the immunoglobulin C (constant) gene loci among four species. This map helps depict what might have been a gene deletion event in an ancestor of the ferret and mink that is not an ancestor of the dog or human.

Fortunately, IMGT provides C gene sequences for all these species. Unfortunately, their sequence record descriptions include start & stop positions that cannot trivially be mapped back to an organism's genome. Therefore, some other mapping tool must be employed to re-annotate their location within a genome. This repository uses a BLAST wrapper.

## I. Acquire IMGT C genes & Genomes
### A. Use NCBI Datasets tool to download genome FASTAs: https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/command-line/datasets/
- Mustel putorius furo
    - `datasets download genome accession GCF_011764305.1 --include genome`
    - (assembly details: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_011764305.1/)
- Canis lupus familiaris
    - `datasets download genome accession GCF_000002285.5 --include genome`
    - (assembly details: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000002285.5/)
- Homo sapiens
    - `datasets download genome accession GCF_000001405.40 --include genome`
    - (assembly details: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_000001405.40/
- Neogale vison
    - `datasets download genome accession GCF_020171115.1 --include genome`
    - (assembly details: https://www.ncbi.nlm.nih.gov/datasets/genome/GCF_020171115.1/)
### B. Use IMGT's gene DB download tool to get previously-annotated entire genomic C gene sequences
- https://www.imgt.org/genedb/
- after narrowing selection by species, gene type C genes, immunoglobulin heavy chain, etc., select nucleotide FASTA output, "C-GENE" download option to get the full genomic C gene locus incuding exons & introns
- For dog IGHG1, IGHG3, and IGHG4, the "C-GENE" is not downloadable via IMGT. It appears these three genes have not been mapped within the immunoglobulin locus in the dog genome assembly, and instead their exons actually map best to a curated mRNA database. Therefore, I have used IGHG2's top 4 hits, hoping to catch the three IGHG paralogs in the dog genome assembly.

## II. Run scripts in this order to create C gene exon map
### A. Per Species
1. keep_one_allele.py
2. download [makeblastdb](https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/) or use [conda install](https://anaconda.org/bioconda/blast) and run `makeblastdb -in species_genome.fasta -parse_seqids -dbtype nucl`
3. cgenes_BLASTn_genome.py
4. besthits_perquery.py
### B. Now using all Species Data
5. combine_hits.py
6. cgene_map.py