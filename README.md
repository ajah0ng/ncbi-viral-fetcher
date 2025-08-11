This Python script involves building a bioinformatics pipeline script to retrieve, annotate, and analyse viral genomic sequences from the NCBI Virus Genome database.

In simpler words, it's a Google search for virus DNA/RNA, but for scientists.

- You give it a virus ID (from the NCBI database).
- It looks up that virus’s complete genome sequence from NCBI.
- It creates a file with all the gene annotations — basically a list of important parts in the virus genome and what they do.
- It produces a FASTA file which is the raw genome sequence, GTF file (gene and coding sequence annotations) and Codon usage table, counts of every codon in the genome.

Short Python script that makes retrieving genomic data easier <3


