# 🧬 This document records how to manipulate some files of bioinformatics
🚨 If you don't believe in me, copy and paste this page into the chat box of an AI you trust to see whether they are correct. I ⚠️**highly recommend**⚠️ you to do this, because I maintain this repo myself, so it is easy for me to make some mistakes.

## 1. Environment Setup
This practice requires the following environment (only words after $ are the commands, words after # are comments for clarity, and you don't need to type them in the terminal):
```bash
$ conda create -n yeast-practice -y
	# -y means skipping confirmation
$ conda activate yeast-practice
$ conda install wget seqkit python
```

## 2. Get Data
In this practice, we will be using data from SGD (Saccharomyces Genome Database). You can type this command in your terminal emulator:
```bash
$ wget http://sgd-archive.yeastgenome.org/sequence/S288C_reference/S288C_reference_genome_Current_Release.tgz
```

## 3. Unpack and Inspect
You have a `.tgz` file in your directory now, let's gain some basic notions about it. Type this command:
```bash
$ tar -xzvf S288C_reference_genome_Current_Release.tgz
	# x means extract, z means input file is in .gz format
	# v means verbose, f designates the input file name
```
This will create a folder called `S288C_reference_genome_R64-5-1_20240529`.
```bash
$ ls S288C_reference_genome_R64-5-1_20240529
```
It will return:
```bash
gene_association_R64-5-1_20240529.sgd
NotFeature_R64-5-1_20240529.fasta.gz
orf_coding_all_R64-5-1_20240529.fasta.gz
orf_trans_all_R64-5-1_20240529.fasta.gz
other_features_genomic_R64-5-1_20240529.fasta.gz
rna_coding_R64-5-1_20240529.fasta.gz
S288C_reference_sequence_R64-5-1_20240529.fsa.gz
saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz
```
So what are these files? Let's see:

| File | Description |
|:--|:--|
| `S288C_reference_sequence_*.fsa.gz` | The reference genome in FASTA format. |
| `saccharomyces_cerevisiae_*.gff.gz` | Genome annotation — genes, transcripts, CDS, etc. |
| `orf_coding_all_*.fasta.gz` | All ORF (open reading frame) sequences. |
| `orf_trans_all_*.fasta.gz` | All transcripts (mRNA sequences). |
| `rna_coding_*.fasta.gz` | All non-coding RNAs. |
| `other_features_genomic_*.fasta.gz` | Miscellaneous genomic features. |
| `NotFeature_*.fasta.gz` | Sequences not annotated as features (intergenic or unknown). |
| `gene_association_*.sgd` | Gene ontology (GO) annotation mapping file. |

Now that we know some basics, let's start exploring the details.

## 4. Genome Overview
Go to this new folder.
```bash
$ cd S288C_reference_genome_R64-5-1_20240529
```
Get a basic notion of yeast genome.
```bash
$ seqkit stat S288C_reference_sequence_R64-5-1_20240529.fsa.gz
file                                              format  type  num_seqs     sum_len  min_len    avg_len    max_len
S288C_reference_sequence_R64-5-1_20240529.fsa.gz  FASTA   DNA         17  12,157,105   85,779  715,123.8  1,531,933
```
It turns out that their are 17 sequences. The shortest one is of 85,779 bp, while the longest one is of 1,531,933 bp. Now let's dive deeper into the genome.
```bash
$ seqkit stat -a S288C_reference_sequence_R64-5-1_20240529.fsa.gz	# -a means --all
file                                              format  type  num_seqs     sum_len  min_len    avg_len    max_len       Q1       Q2       Q3  sum_gap      N50  N50_num  Q20(%)  Q30(%)  AvgQual  GC(%)  sum_n
S288C_reference_sequence_R64-5-1_20240529.fsa.gz  FASTA   DNA         17  12,157,105   85,779  715,123.8  1,531,933  439,888  745,751  948,066        0  924,431        6       0       0        0  38.15      0
```
We get more information.

1️⃣ **Assembly Quality Metrics**
| Metric | Value | Interpretation |
|--------|-------|----------------|
| **N50** | 924,431 bp | **Definition**: Length of the shortest sequence at 50% of the total genome length when all sequences are ordered from longest to shortest.<br> Indicates 50% of the genome is contained in sequences ≥924 kb.<br> High value suggests excellent assembly continuity. |
| **N50_num** | 6 | Number of sequences needed to reach N50 threshold.<br> Only 6 longest chromosomes contain half of the genome. |
| **sum_gap** | 0 | Zero gaps in assembly (perfect continuity). |
| **sum_n** | 0 | No ambiguous bases (N's) in reference. |

2️⃣ **Sequence Composition**
| Metric | Value | Interpretation |
|--------|-------|----------------|
| **GC(%)** | 38.15% | Typical GC content for S. cerevisiae.<br> Slightly lower than gene-rich regions (which average ~40%). |
| **Q20(%)/Q30(%)** | 0% | Not applicable for reference genomes (quality scores only relevant for raw sequencing data). |

3️⃣ **Size Distribution (Quartiles)**
| Metric | Value | Interpretation |
|--------|-------|----------------|
| **Q1** | 439,888 bp | 25% of chromosomes are ≤439 kb. |
| **Q2 (Median)** | 745,751 bp | Half of chromosomes are ≤746 kb. |
| **Q3** | 948,066 bp | 75% of chromosomes are ≤948 kb. |

We have obtained some general information of yeast genome. It's time to focus on each individual sequence. Let's check out what does this FASTA file look like first.
```bash
$ gzcat S288C_reference_sequence_R64-5-1_20240529.fsa.gz | head
>ref|NC_001133| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=I]
CCACACCACACCCACACACCCACACACCACACCACACACCACACCACACCCACACACACA
CATCCTAACACTACCCTAACACAGCCCTAATCTAACCCTGGCCAACCTGTCTCTCAACTT
ACCCTCCATTACCCTGCCTCCACTCGTTACCCTGTCCCATTCAACCATACCACTCCGAAC
CACCATCCATCCCTCTACTTACTACCACTCACCCACCGTTACCCTCCAATTACCCATATC
CAACCCACTGCCACTTACCCTACCATTACCCTACCATCCACCATGACCTACTCACCATAC
TGTTCTTCTACCCACCATATTGAAACGCTAACAAATGATCGTAAATAACACACACGTGCT
TACCCTACCACTTTATACCACCACCACATGCCATACTCACCCTCACTTGTATACTGATTT
TACGTACGCACACGGATGCTACAGTATATACCATCTCAAACTTACCCTACTCTCAGATTC
CACTTCACTCCATGGCCCATCTCTCACTGAATCAGTACCAAATGCACTCACATCATTATG
```
Recall how did we introduce FASTA file format. `>ref|NC_001133| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=I]
` is the ID of the sequence followed by some descriptions. In this case, we find that this is chromosome I of yeast, and is from S288C strain. You can also inspect the whole file.
```bash
$ zless -S S288C_reference_sequence_R64-5-1_20240529.fsa.gz
```
It will open a text viewer in your terminal. But FASTA files are always very long, it is not very useful to look at the whole file directly to get information. For example if we want ot know how many chromosomes does yeast have, rather than counting the sequences manually in `zless`, we can do the following. Of course we have obtained this stat by the previous `seqkit stat`, I just want to show more about the command line.
```bash
$ gzcat S288C_reference_sequence_R64-5-1_20240529.fsa.gz | grep '>' | wc -l
      17
$ seqkit seq -n S288C_reference_sequence_R64-5-1_20240529.fsa.gz | wc -l
      17
```
Single quotes are use to prevent shell from parsing the redirection operator (`>`). It gives us that yeast has 17 chromosomes. You can remove the `wc -l` command (pipe `|` before `wc` also removed) to get the details, or we may apply some professional bioinformatics command line tools.
```bash
$ seqkit fx2tab -n -l -g S288C_reference_sequence_R64-5-1_20240529.fsa.gz
ref|NC_001133| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=I]	230218	39.27
ref|NC_001134| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=II]	813184	38.34
ref|NC_001135| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=III]	316620	38.53
ref|NC_001136| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=IV]	1531933	37.91
ref|NC_001137| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=V]	576874	38.51
ref|NC_001138| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=VI]	270161	38.73
ref|NC_001139| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=VII]	1090940	38.06
ref|NC_001140| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=VIII]	562643	38.50
ref|NC_001141| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=IX]	439888	38.90
ref|NC_001142| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=X]	745751	38.37
ref|NC_001143| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XI]	666816	38.07
ref|NC_001144| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XII]	1078177	38.48
ref|NC_001145| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XIII]	924431	38.20
ref|NC_001146| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XIV]	784333	38.64
ref|NC_001147| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XV]	1091291	38.16
ref|NC_001148| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XVI]	948066	38.06
ref|NC_001224| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [location=mitochondrion] [top=circular]	85779	17.11
```
`fx2tab` means FASTA to table, `-n` means output names, `-l` means output the length, `-g` means output GC content. Recall what did we get previously. The longest chromosome turns out to be chromosome IV, and the shortest one is mitochondrial DNA. The last column is the GC content. We also notice that in the last row, there is a description `[top=circular]`. `top` means topology, so this means mitocondrial DNA is circular, which is a common sense. Although we have seen that there is no ambiguous base, we still want to know how to check for each individual sequence.
```bash
$ seqkit fx2tab -n -l -B N S288C_reference_sequence_R64-5-1_20240529.fsa.gz | awk -F'\t' '$3 == 0 {print}'
```
`-B` means `--base-content` and `N` represents any base. `awk` is a command line tool that processes text row by row. `-F'\t'` means designate tab as the delimiter. `$3` means the third field (column) in the current row. Since we know that there is no ambiguous base in our reference, this command will output everything. You can verify this by removing the command after pipe `|` (also removed), and then do it again.

## 5. Understand the Annotations
Note that there is a GFF file in our directory.
```bash
$ ls saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz
saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz
```
We will start with looking at the raw file.
```bash
$ gzcat saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz | head -n 25
##gff-version 3
#!date-produced 2024-05-29 16:45:20
#!data-source SGD
#!assembly R64-4-1
#!refseq-version GCF_000146045.2
#
# Saccharomyces cerevisiae S288C genome (version=R64-4-1)
#
# Features from the 16 nuclear chromosomes labeled chrI to chrXVI,
# plus the mitochondrial genome labeled chrmt.
#
# Created by Saccharomyces Genome Database (http://www.yeastgenome.org/)
#
# Weekly updates of this file are available for download from:
# https://downloads.yeastgenome.org/latest/saccharomyces_cerevisiae.gff.gz
#
# Please send comments and suggestions to sgd-helpdesk@lists.stanford.edu
#
# SGD is funded as a National Human Genome Research Institute Biomedical Informatics Resource from
# the U. S. National Institutes of Health to Stanford University.
#
chrI	SGD	chromosome	1	230218	.	.	.	ID=chrI;dbxref=NCBI:BK006935.2;Name=chrI
chrI	SGD	telomere	1	801	.	-	.	ID=TEL01L;Name=TEL01L;Note=Telomeric%20region%20on%20the%20left%20arm%20of%20Chromosome%20I%3B%20composed%20of%20an%20X%20element%20core%20sequence%2C%20X%20element%20combinatorial%20repeats%2C%20and%20a%20short%20terminal%20stretch%20of%20telomeric%20repeats;display=Telomeric%20region%20on%20the%20left%20arm%20of%20Chromosome%20I;dbxref=SGD:S000028862;curie=SGD:S000028862
chrI	SGD	X_element	337	801	.	-	.	Parent=TEL01L;Name=TEL01L_X_element
chrI	SGD	X_element_combinatorial_repeat	63	336	.	-	.	Parent=TEL01L;Name=TEL01L_X_element_combinatorial_repeat
```
Try to recall what did we discuss about GFF file format in File-Format directory and identify what does each field represent. First let's see how many lines do we have. So what identifier can we use to determine which line is the actual annotation?
```bash
$ gzcat saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz | grep '#'
##gff-version 3
#!date-produced 2024-05-29 16:45:20
#!data-source SGD
#!assembly R64-4-1
#!refseq-version GCF_000146045.2
#
# Saccharomyces cerevisiae S288C genome (version=R64-4-1)
#
# Features from the 16 nuclear chromosomes labeled chrI to chrXVI,
# plus the mitochondrial genome labeled chrmt.
#
# Created by Saccharomyces Genome Database (http://www.yeastgenome.org/)
#
# Weekly updates of this file are available for download from:
# https://downloads.yeastgenome.org/latest/saccharomyces_cerevisiae.gff.gz
#
# Please send comments and suggestions to sgd-helpdesk@lists.stanford.edu
#
# SGD is funded as a National Human Genome Research Institute Biomedical Informatics Resource from
# the U. S. National Institutes of Health to Stanford University.
#
chrVII	SGD	gene	899527	902329	.	-	.	ID=YGR200C;Name=YGR200C;gene=ELP2;so_term_name=protein_coding_gene;Alias=ELP2,Elongator%20subunit%20ELP2,KTI3,TOT2;Ontology_term=GO:0002098,GO:0005634,GO:0005737,GO:0006357,GO:0008017,GO:0032447,GO:0033588,SO:0000236;Note=Subunit%20of%20Elongator%20complex%3B%20binds%20to%20microtubules%20via%20conserved%20alkaline%20residues%3B%20has%20two%20seven-bladed%20WD40%20&#946%3B%20propellers%3B%20Elongator%20complex%20is%20required%20for%20modification%20of%20wobble%20nucleosides%20in%20tRNA%3B%20target%20of%20Kluyveromyces%20lactis%20zymocin;display=Subunit%20of%20Elongator%20complex;dbxref=SGD:S000003432;orf_classification=Verified;curie=SGD:S000003432
chrXII	SGD	gene	406844	408240	.	-	.	ID=YLR132C;Name=YLR132C;gene=USB1;so_term_name=protein_coding_gene;Alias=USB1,phosphoric%20diester%20hydrolase;Ontology_term=GO:0005634,GO:0005634,GO:0005634,GO:0005737,GO:0005739,GO:0005739,GO:0005739,GO:0034477,GO:0034477,GO:1990838,SO:0000236;Note=Putative%20poly%28U%29-specific%203'-to-5'%20RNA%20exonuclease%3B%20involved%20in%203'-end%20processing%20of%20U6%20snRNA%20removing%20uridines%20and%20generating%20a%20terminal%202&#8242%3B%2C3&#8242%3B%20cyclic%20phosphate%3B%20essential%20protein%20that%20localizes%20to%20the%20nucleus%20and%20mitochondria%3B%20overexpression%20suppresses%20the%20respiratory%20defects%20of%20oxa1%20and%20mtf2%20mutants%3B%20homolog%20of%20S.pombe%20gene%2C%20mpn1%20and%20human%20gene%2C%20hUSB1%3B%20mutations%20in%20hUSB1%20are%20associated%20with%20a%20rare%20genodermatosis%2C%20poikiloderma%20with%20neutropenia%20%28OMIM%20604173%29;display=Putative%20poly%28U%29-specific%203'-to-5'%20RNA%20exonuclease;dbxref=SGD:S000004122;orf_classification=Verified;curie=SGD:S000004122
###
##FASTA
```
It turns out that some annotation lines also contains hashtag (`#`) and this GFF also contains FASTA sequences. We want to extract solely GFF part first.
```bash
$ gzcat saccharomyces_cerevisiae_R64-5-1_20240529.gff.gz | sed -n '/^###/q; p' | gzip > s_cerevisiae_annotation_only.gff.gz
```
`-n` means not to output each line automatically, `^###` means match the line starting with `###`, `q` means quit, `p` means print. Now we get the GFF file excluding FASTA and we can count how many annotation lines do we have.
```bash
$ gzcat s_cerevisiae_annotation_only.gff.gz | grep -v '^#' | wc -l
   28347
```
`-v` means inverse match. There are a total of 28347 annotation lines. Not too many. Now let's see what type the features are.
```bash
$ gzcat s_cerevisiae_annotation_only.gff.gz | grep -v '^#' | cut -d $'\t' -f 3 | sort | uniq -c | sort -nr
11125 mRNA
7072 CDS
6613 gene
 543 ARS
 497 noncoding_exon
 384 long_terminal_repeat
 378 intron
 299 tRNA_gene
 299 tRNA
 196 ARS_consensus_sequence
  91 transposable_element_gene
  78 region
  77 snoRNA_gene
  77 snoRNA
  64 telomere
  50 LTR_retrotransposon
  47 plus_1_translational_frameshift
  32 X_element
  32 centromere
  31 telomeric_repeat
  31 ncRNA_gene
  31 ncRNA
  28 X_element_combinatorial_repeat
  24 rRNA_gene
  24 rRNA
  24 five_prime_UTR_intron
  21 uORF
  19 Y_prime_element
  17 chromosome
  16 centromere_DNA_Element_III
  16 centromere_DNA_Element_II
  16 centromere_DNA_Element_I
  12 pseudogene
  12 blocked_reading_frame
   8 origin_of_replication
   8 matrix_attachment_site
   8 internal_transcribed_spacer_region
   8 external_transcribed_spacer_region
   6 snRNA_gene
   6 snRNA
   4 silent_mating_type_cassette_array
   3 Z1_region
   3 Y_region
   3 X_region
   3 non_transcribed_region
   2 Z2_region
   2 W_region
   2 mating_type_region
   2 intein_encoding_region
   1 telomerase_RNA_gene
   1 telomerase_RNA
   1 recombination_enhancer
```
`cut` used to extract the fields, `-d` means designate delimiter, `$'\t'` means tab, `-f 3` means the third field, `uniq -c` means count how many times does each word appear, `sort -nr` means sort by number in descending order. Amazing. There are many types. You may want to explore what does each feature mean. I guess these are SO (Sequence Ontology) terms, so you can visit their official website. Simply google it. Very easy. But according to my experience, those terms are used ambiguously in practice. You can see that the number of mRNA > gene, which reveals that there exists alternative splicing. Now let's see how many features does each chromosome have.
```bash
$ gzcat s_cerevisiae_annotation_only.gff.gz | grep -v '^#' | cut -d $'\t' -f 1 | sort | uniq -c | sort -nr
3507 chrIV
2563 chrVII
2495 chrXII
2472 chrXV
2194 chrXVI
2184 chrXIII
1912 chrII
1826 chrXIV
1697 chrX
1453 chrXI
1426 chrV
1349 chrVIII
1020 chrIX
 855 chrIII
 649 chrVI
 504 chrI
 241 chrmt
```
Obviously the longest chromosome has the most features, but how about the density? Here we will only count the number of genes to calculate the density. First, we get the number of genes on each chromosome.
```bash
$ gzcat s_cerevisiae_annotation_only.gff.gz | awk -F'\t' '$3=="gene" {print $1}' | sort | uniq -c | awk '{print $2, $1}' OFS='\t'
chrI	117
chrII	456
chrIII	184
chrIV	837
chrIX	241
chrmt	28
chrV	323
chrVI	140
chrVII	585
chrVIII	323
chrX	399
chrXI	349
chrXII	579
chrXIII	507
chrXIV	437
chrXV	597
chrXVI	511
```
This command is very complicated, but don't be scare. Try to ask AI to let it explain this command. Now we want to get the length of each chromosome.
```bash
$ seqkit fx2tab -n -l S288C_reference_sequence_R64-5-1_20240529.fsa.gz
ref|NC_001133| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=I]	230218
ref|NC_001134| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=II]	813184
ref|NC_001135| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=III]	316620
ref|NC_001136| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=IV]	1531933
ref|NC_001137| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=V]	576874
ref|NC_001138| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=VI]	270161
ref|NC_001139| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=VII]	1090940
ref|NC_001140| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=VIII]	562643
ref|NC_001141| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=IX]	439888
ref|NC_001142| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=X]	745751
ref|NC_001143| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XI]	666816
ref|NC_001144| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XII]	1078177
ref|NC_001145| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XIII]	924431
ref|NC_001146| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XIV]	784333
ref|NC_001147| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XV]	1091291
ref|NC_001148| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [chromosome=XVI]	948066
ref|NC_001224| [org=Saccharomyces cerevisiae] [strain=S288C] [moltype=genomic] [location=mitochondrion] [top=circular]	85779
```
But the name provided by `seqkit` does not match FASTA, so we need to make a little change.
```bash
$ seqkit fx2tab -n -l S288C_reference_sequence_R64-5-1_20240529.fsa.gz | awk -F'\t' '{chr=""; if($0 ~ /chromosome=/) {split($0, a, /chromosome=/); split(a[2], b, /]/); chr="chr"b[1]} else {chr="chrmt"} print chr "\t" $2}' | sort
chrI	230218
chrII	813184
chrIII	316620
chrIV	1531933
chrIX	439888
chrmt	85779
chrV	576874
chrVI	270161
chrVII	1090940
chrVIII	562643
chrX	745751
chrXI	666816
chrXII	1078177
chrXIII	924431
chrXIV	784333
chrXV	1091291
chrXVI	948066
```
This command is even more complicated. I am sorry for not being smart enough to figure out a simpler way to do so. I can not explaint it in detail. If you want to know the syntax of `awk`, you can turn to AI. Now we can calculate the density. First we will store the results above into variables for legibility.
```bash
$ gzcat s_cerevisiae_annotation_only.gff.gz | awk -F'\t' '$3=="gene" {print $1}' | sort | uniq -c | awk '{print $2, $1}' OFS='\t' > gene_counts
$ seqkit fx2tab -n -l S288C_reference_sequence_R64-5-1_20240529.fsa.gz | awk -F'\t' '{chr=""; if($0 ~ /chromosome=/) {split($0, a, /chromosome=/); split(a[2], b, /]/); chr="chr"b[1]} else {chr="chrmt"} print chr "\t" $2}' | sort > chr_length
$ join -t $'\t' -1 1 -2 1 gene_counts chr_length > counts_length_pair
$ awk -F'\t' '{printf "%s\t%.1f genes/Mb\n", $0, $2/($3/1000000)}' counts_length_pair
chrI	117	230218	508.2 genes/Mb
chrII	456	813184	560.8 genes/Mb
chrIII	184	316620	581.1 genes/Mb
chrIV	837	1531933	546.4 genes/Mb
chrIX	241	439888	547.9 genes/Mb
chrmt	28	85779	326.4 genes/Mb
chrV	323	576874	559.9 genes/Mb
chrVI	140	270161	518.2 genes/Mb
chrVII	585	1090940	536.2 genes/Mb
chrVIII	323	562643	574.1 genes/Mb
chrX	399	745751	535.0 genes/Mb
chrXI	349	666816	523.4 genes/Mb
chrXII	579	1078177	537.0 genes/Mb
chrXIII	507	924431	548.4 genes/Mb
chrXIV	437	784333	557.2 genes/Mb
chrXV	597	1091291	547.1 genes/Mb
chrXVI	511	948066	539.0 genes/Mb
```
Finally, we have obtained the gene density on each chromosome. We can also sort it.
```bash
$ awk -F'\t' '{printf "%s\t%.1f genes/Mb\n", $0, $2/($3/1000000)}' counts_length_pair | sort -k4,4nr
chrIII	184	316620	581.1 genes/Mb
chrVIII	323	562643	574.1 genes/Mb
chrII	456	813184	560.8 genes/Mb
chrV	323	576874	559.9 genes/Mb
chrXIV	437	784333	557.2 genes/Mb
chrXIII	507	924431	548.4 genes/Mb
chrIX	241	439888	547.9 genes/Mb
chrXV	597	1091291	547.1 genes/Mb
chrIV	837	1531933	546.4 genes/Mb
chrXVI	511	948066	539.0 genes/Mb
chrXII	579	1078177	537.0 genes/Mb
chrVII	585	1090940	536.2 genes/Mb
chrX	399	745751	535.0 genes/Mb
chrXI	349	666816	523.4 genes/Mb
chrVI	140	270161	518.2 genes/Mb
chrI	117	230218	508.2 genes/Mb
chrmt	28	85779	326.4 genes/Mb
```
Note that not the longest chromosome has the highest density. What's the biological significance of high gene density? You can explore it yourself. Maybe it has something to do with evolution. Now let's try to find the longest gene.
```bash
$ gzcat s_cerevisiae_annotation_only.gff.gz | awk -F'\t' '$3=="gene" {len=$5-$4+1; print $1, $3, $4, $5, len, $9}' OFS='\t' | sort -k5,5nr | head -n 3
chrXII	gene	349006	363738	14733	ID=YLR106C;Name=YLR106C;gene=REA1;so_term_name=protein_coding_gene;Alias=REA1,AAA%20family%20ATPase%20midasin,MDN1;Ontology_term=GO:0000027,GO:0005634,GO:0005654,GO:0005739,GO:0005739,GO:0006364,GO:0016887,GO:0016887,GO:0110136,GO:0110136,GO:0110136,GO:2000200,SO:0000236;Note=Huge%20dynein-related%20AAA-type%20ATPase%20%28midasin%29%3B%20forms%20extended%20pre-60S%20particle%20with%20the%20Rix1%20complex%3B%20involved%20with%20interaction%20partners%20Rsa4p%20and%20Ytm1p%2C%20in%20the%20ATP-dependent%20remodeling%20of%20the%20pre-60S%20particle%20at%20successive%20maturation%20steps%20during%20ribosomal%20biogenesis%3B%20involved%20in%20the%20removal%20of%20biogenesis%20factors%20including%20GTPase%20Nog2p%20prior%20to%20nuclear%20export%3B%20contains%20a%20hexameric%20AAA-motor%20head%20domain%20and%20a%20long%20flexible%20tail%20with%20a%20MIDAS%20%28metal%20ion-dependent%20adhesion%20site%29%20domain;display=Huge%20dynein-related%20AAA-type%20ATPase%20%28midasin%29;dbxref=SGD:S000004096;orf_classification=Verified;curie=SGD:S000004096
chrmt	gene	13818	26701	12884	ID=Q0045;Name=Q0045;gene=COX1;so_term_name=protein_coding_gene;Alias=COX1,OXI3,cytochrome%20c%20oxidase%20subunit%201;Ontology_term=GO:0004129,GO:0004129,GO:0005739,GO:0005751,GO:0006123,GO:0009060,SO:0000236;Note=Subunit%20I%20of%20cytochrome%20c%20oxidase%20%28Complex%20IV%29%3B%20Complex%20IV%20is%20the%20terminal%20member%20of%20the%20mitochondrial%20inner%20membrane%20electron%20transport%20chain%3B%20one%20of%20three%20mitochondrially-encoded%20subunits%3B%20number%20of%20introns%20in%20different%20strains%20varies%20from%202%20to%208%2C%20with%20most%20strains%20having%204-6%20introns;display=Subunit%20I%20of%20cytochrome%20c%20oxidase%20%28Complex%20IV%29;dbxref=SGD:S000007260;orf_classification=Verified;curie=SGD:S000007260
chrXI	gene	535647	547925	12279	ID=YKR054C;Name=YKR054C;gene=DYN1;so_term_name=protein_coding_gene;Alias=DYN1,DHC1,PAC6,dynein%20heavy%20chain;Ontology_term=GO:0000070,GO:0000132,GO:0000235,GO:0005524,GO:0005816,GO:0005868,GO:0005881,GO:0005881,GO:0005938,GO:0008569,GO:0008569,GO:0008569,GO:0030473,GO:0030473,GO:0040001,SO:0000236;Note=Cytoplasmic%20heavy%20chain%20dynein%3B%20microtubule%20motor%20protein%3B%20member%20of%20the%20AAA+%20protein%20family%2C%20required%20for%20anaphase%20spindle%20elongation%3B%20involved%20in%20spindle%20assembly%2C%20chromosome%20movement%2C%20and%20spindle%20orientation%20during%20cell%20division%2C%20targeted%20to%20microtubule%20tips%20by%20Pac1p%3B%20motility%20along%20microtubules%20inhibited%20by%20She1p;display=Cytoplasmic%20heavy%20chain%20dynein;dbxref=SGD:S000001762;orf_classification=Verified;curie=SGD:S000001762
```
The fifth field (column) is the length of the gene, ohters are from original GFF. The third and fouth fields are the coordinates. The last field is the attributes. Let's find out what are these genes. Because some special characters in GFF attributes fields are represented using URL encoding, so we need to decode them for legibility. Open your text editor and write some python code.
```python
from urllib.parse import unquote

gff_lines = [
	"chrXII	gene	349006	363738	14733	ID=YLR106C;Name=YLR106C;gene=REA1;so_term_name=protein_coding_gene;Alias=REA1,AAA%20family%20ATPase%20midasin,MDN1;Ontology_term=GO:0000027,GO:0005634,GO:0005654,GO:0005739,GO:0005739,GO:0006364,GO:0016887,GO:0016887,GO:0110136,GO:0110136,GO:0110136,GO:2000200,SO:0000236;Note=Huge%20dynein-related%20AAA-type%20ATPase%20%28midasin%29%3B%20forms%20extended%20pre-60S%20particle%20with%20the%20Rix1%20complex%3B%20involved%20with%20interaction%20partners%20Rsa4p%20and%20Ytm1p%2C%20in%20the%20ATP-dependent%20remodeling%20of%20the%20pre-60S%20particle%20at%20successive%20maturation%20steps%20during%20ribosomal%20biogenesis%3B%20involved%20in%20the%20removal%20of%20biogenesis%20factors%20including%20GTPase%20Nog2p%20prior%20to%20nuclear%20export%3B%20contains%20a%20hexameric%20AAA-motor%20head%20domain%20and%20a%20long%20flexible%20tail%20with%20a%20MIDAS%20%28metal%20ion-dependent%20adhesion%20site%29%20domain;display=Huge%20dynein-related%20AAA-type%20ATPase%20%28midasin%29;dbxref=SGD:S000004096;orf_classification=Verified;curie=SGD:S000004096",
	"chrmt	gene	13818	26701	12884	ID=Q0045;Name=Q0045;gene=COX1;so_term_name=protein_coding_gene;Alias=COX1,OXI3,cytochrome%20c%20oxidase%20subunit%201;Ontology_term=GO:0004129,GO:0004129,GO:0005739,GO:0005751,GO:0006123,GO:0009060,SO:0000236;Note=Subunit%20I%20of%20cytochrome%20c%20oxidase%20%28Complex%20IV%29%3B%20Complex%20IV%20is%20the%20terminal%20member%20of%20the%20mitochondrial%20inner%20membrane%20electron%20transport%20chain%3B%20one%20of%20three%20mitochondrially-encoded%20subunits%3B%20number%20of%20introns%20in%20different%20strains%20varies%20from%202%20to%208%2C%20with%20most%20strains%20having%204-6%20introns;display=Subunit%20I%20of%20cytochrome%20c%20oxidase%20%28Complex%20IV%29;dbxref=SGD:S000007260;orf_classification=Verified;curie=SGD:S000007260",
	"chrXI	gene	535647	547925	12279	ID=YKR054C;Name=YKR054C;gene=DYN1;so_term_name=protein_coding_gene;Alias=DYN1,DHC1,PAC6,dynein%20heavy%20chain;Ontology_term=GO:0000070,GO:0000132,GO:0000235,GO:0005524,GO:0005816,GO:0005868,GO:0005881,GO:0005881,GO:0005938,GO:0008569,GO:0008569,GO:0008569,GO:0030473,GO:0030473,GO:0040001,SO:0000236;Note=Cytoplasmic%20heavy%20chain%20dynein%3B%20microtubule%20motor%20protein%3B%20member%20of%20the%20AAA+%20protein%20family%2C%20required%20for%20anaphase%20spindle%20elongation%3B%20involved%20in%20spindle%20assembly%2C%20chromosome%20movement%2C%20and%20spindle%20orientation%20during%20cell%20division%2C%20targeted%20to%20microtubule%20tips%20by%20Pac1p%3B%20motility%20along%20microtubules%20inhibited%20by%20She1p;display=Cytoplasmic%20heavy%20chain%20dynein;dbxref=SGD:S000001762;orf_classification=Verified;curie=SGD:S000001762"
]

for line in gff_lines:
	parts = line.split("\t")
	parts[5] = unquote(parts[5])
	print("\t".join(parts))
```
Name this file `decode_gff.py`. The items in `gff_lines` are the lines output by the previous command. I simply copy and paste them. Run the script in the terminal.
```bash
$ python3 decode_gff.py
chrXII	gene	349006	363738	14733	ID=YLR106C;Name=YLR106C;gene=REA1;so_term_name=protein_coding_gene;Alias=REA1,AAA family ATPase midasin,MDN1;Ontology_term=GO:0000027,GO:0005634,GO:0005654,GO:0005739,GO:0005739,GO:0006364,GO:0016887,GO:0016887,GO:0110136,GO:0110136,GO:0110136,GO:2000200,SO:0000236;Note=Huge dynein-related AAA-type ATPase (midasin); forms extended pre-60S particle with the Rix1 complex; involved with interaction partners Rsa4p and Ytm1p, in the ATP-dependent remodeling of the pre-60S particle at successive maturation steps during ribosomal biogenesis; involved in the removal of biogenesis factors including GTPase Nog2p prior to nuclear export; contains a hexameric AAA-motor head domain and a long flexible tail with a MIDAS (metal ion-dependent adhesion site) domain;display=Huge dynein-related AAA-type ATPase (midasin);dbxref=SGD:S000004096;orf_classification=Verified;curie=SGD:S000004096
chrmt	gene	13818	26701	12884	ID=Q0045;Name=Q0045;gene=COX1;so_term_name=protein_coding_gene;Alias=COX1,OXI3,cytochrome c oxidase subunit 1;Ontology_term=GO:0004129,GO:0004129,GO:0005739,GO:0005751,GO:0006123,GO:0009060,SO:0000236;Note=Subunit I of cytochrome c oxidase (Complex IV); Complex IV is the terminal member of the mitochondrial inner membrane electron transport chain; one of three mitochondrially-encoded subunits; number of introns in different strains varies from 2 to 8, with most strains having 4-6 introns;display=Subunit I of cytochrome c oxidase (Complex IV);dbxref=SGD:S000007260;orf_classification=Verified;curie=SGD:S000007260
chrXI	gene	535647	547925	12279	ID=YKR054C;Name=YKR054C;gene=DYN1;so_term_name=protein_coding_gene;Alias=DYN1,DHC1,PAC6,dynein heavy chain;Ontology_term=GO:0000070,GO:0000132,GO:0000235,GO:0005524,GO:0005816,GO:0005868,GO:0005881,GO:0005881,GO:0005938,GO:0008569,GO:0008569,GO:0008569,GO:0030473,GO:0030473,GO:0040001,SO:0000236;Note=Cytoplasmic heavy chain dynein; microtubule motor protein; member of the AAA+ protein family, required for anaphase spindle elongation; involved in spindle assembly, chromosome movement, and spindle orientation during cell division, targeted to microtubule tips by Pac1p; motility along microtubules inhibited by She1p;display=Cytoplasmic heavy chain dynein;dbxref=SGD:S000001762;orf_classification=Verified;curie=SGD:S000001762
```
Now you can see that the longest gene (YLR106C/REA1) encodes a ​huge AAA family ATPase​ that plays a ​critical role in ribosomal large subunit assembly. You can explore the functions of the second and third longest genes yourself. Have a good time.

## 6. Gene Association
Now we will explore the `gene_association_R64-5-1_20240529.sgd` file. Take a quick glance first.
```bash
$ head gene_association_R64-5-1_20240529.sgd
!gaf-version: 2.2
!date-generated: 20240529
!generated-by: Saccharomyces Genome Database (SGD)
!URL: https://www.yeastgenome.org/
!Contact Email: sgd-helpdesk@lists.stanford.edu
!Funding: NHGRI at US NIH, grant number U41-HG001315
!
SGD	S000006369	RHO1	involved_in	GO:0090334	PMID:36949198	IDA		P	GTP-binding protein of the rho subfamily of Ras-like small GTPases	YPR165W|Rho family GTPase RHO1	protein	taxon:559292	20230807	SGD		UniProtKB:P06780
SGD	S000006369	RHO1	enables	GO:0008047	PMID:36949198	IDA		F	GTP-binding protein of the rho subfamily of Ras-like small GTPases	YPR165W|Rho family GTPase RHO1	protein	taxon:559292	20230807	SGD	part_of(GO:0090334),has_input(SGD:S000004334)	UniProtKB:P06780
SGD	S000003211	ANK1	is_active_in	GO:0005737	PMID:37531259	IDA		C	Cytoplasmic ankyrin repeat-containing protein	YGL242C	protein	taxon:559292	20230809	SGD		UniProtKB:P53066
```
As is often the first, let's count how many lines do we have.
```bash
$ grep -v '^!' gene_association_R64-5-1_20240529.sgd | wc -l
  201804
```
Wow! 200K, that's a lot. Let's see how many annotations there are for each GO category.
```bash
$ grep -v '^!' gene_association_R64-5-1_20240529.sgd | cut -d $'\t' -f 9 | sort | uniq -c | sort -nr
77832 C
62158 F
61814 P
```
Seems like our scientists know yeast very well, especially the subcellular localization. Do you still remember that we found the longest gene of yeast? We want to find the more detailed funcion of it. Let's try.
```bash
$ awk -F'\t' '$3=="REA1" {print}' gene_association_R64-5-1_20240529.sgd
SGD	S000004096	REA1	enables	GO:0005524	GO_REF:0000043	IEA	UniProtKB-KW:KW-0067	F	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20240408	UniProt		UniProtKB:Q12019
SGD	S000004096	REA1	enables	GO:0005524	GO_REF:0000002	IEA	InterPro:IPR011704	F	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20240408	InterPro		UniProtKB:Q12019
SGD	S000004096	REA1	involved_in	GO:0000027	GO_REF:0000033	IBA	PANTHER:PTN000529705	P	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20240213	GO_Central		UniProtKB:Q12019
SGD	S000004096	REA1	involved_in	GO:0000027	GO_REF:0000033	IBA	SGD:S000004096	P	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20240213	GO_Central		UniProtKB:Q12019
SGD	S000004096	REA1	involved_in	GO:0000027	GO_REF:0000033	IBA	UniProtKB:Q9NU22	P	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20240213	GO_Central		UniProtKB:Q12019
SGD	S000004096	REA1	located_in	GO:0005730	GO_REF:0000044	IEA	UniProtKB-SubCell:SL-0188	C	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20240408	UniProt		UniProtKB:Q12019
SGD	S000004096	REA1	enables	GO:0016887	GO_REF:0000002	IEA	InterPro:IPR003593	F	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20240408	InterPro		UniProtKB:Q12019
SGD	S000004096	REA1	enables	GO:0016887	GO_REF:0000002	IEA	InterPro:IPR011704	F	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20240408	InterPro		UniProtKB:Q12019
SGD	S000004096	REA1	enables	GO:0016887	GO_REF:0000002	IEA	InterPro:IPR012099	F	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20240408	InterPro		UniProtKB:Q12019
SGD	S000004096	REA1	involved_in	GO:0006364	PMID:12837249	IMP		P	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20030805	SGD		UniProtKB:Q12019
SGD	S000004096	REA1	involved_in	GO:0000027	PMID:15528184	IMP		P	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20050310	SGD		UniProtKB:Q12019
SGD	S000004096	REA1	involved_in	GO:0110136	PMID:24240281	IMP		P	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20190129	SGD	has_input(GO:0030687)	UniProtKB:Q12019
SGD	S000004096	REA1	involved_in	GO:0110136	PMID:20542003	IDA		P	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20190129	SGD	has_input(GO:0030687)	UniProtKB:Q12019
SGD	S000004096	REA1	involved_in	GO:0110136	PMID:20542003	IMP		P	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20190129	SGD	has_input(GO:0030687)	UniProtKB:Q12019
SGD	S000004096	REA1	involved_in	GO:2000200	PMID:24240281	IMP		P	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20181205	SGD	has_input(GO:0030687)	UniProtKB:Q12019
SGD	S000004096	REA1	located_in	GO:0005654	PMID:15528184	IDA		C	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20050310	SGD		UniProtKB:Q12019
SGD	S000004096	REA1	located_in	GO:0005739	PMID:14576278	HDA		C	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20181119	SGD		UniProtKB:Q12019
SGD	S000004096	REA1	located_in	GO:0005739	PMID:16823961	HDA		C	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20181119	SGD		UniProtKB:Q12019
SGD	S000004096	REA1	enables	GO:0016887	PMID:30460895	IDA		F	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20190905	SGD		UniProtKB:Q12019
SGD	S000004096	REA1	enables	GO:0016887	PMID:30460895	IMP		F	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20190905	SGD		UniProtKB:Q12019
SGD	S000004096	REA1	involved_in	GO:0000055	GO_REF:0000033	IBA	PANTHER:PTN000529705	P	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20220601	GO_Central		UniProtKB:Q12019
SGD	S000004096	REA1	involved_in	GO:0000055	GO_REF:0000033	IBA	TAIR:locus:2033661	P	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20220601	GO_Central		UniProtKB:Q12019
SGD	S000004096	REA1	is_active_in	GO:0005634	GO_REF:0000033	IBA	PANTHER:PTN000529705	C	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20220309	GO_Central		UniProtKB:Q12019
SGD	S000004096	REA1	is_active_in	GO:0005634	GO_REF:0000033	IBA	SGD:S000004096	C	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20220309	GO_Central		UniProtKB:Q12019
SGD	S000004096	REA1	is_active_in	GO:0005634	GO_REF:0000033	IBA	TAIR:locus:2033661	C	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20220309	GO_Central		UniProtKB:Q12019
SGD	S000004096	REA1	is_active_in	GO:0005634	GO_REF:0000033	IBA	UniProtKB:Q9NU22	C	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20220309	GO_Central		UniProtKB:Q12019
SGD	S000004096	REA1	part_of	GO:0030687	GO_REF:0000033	IBA	PANTHER:PTN000529705	C	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20170228	GO_Central		UniProtKB:Q12019
SGD	S000004096	REA1	part_of	GO:0030687	GO_REF:0000033	IBA	PomBase:SPCC737.08	C	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20170228	GO_Central		UniProtKB:Q12019
SGD	S000004096	REA1	located_in	GO:0005634	PMID:12102729	IDA		C	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20020715	SGD		UniProtKB:Q12019
SGD	S000004096	REA1	located_in	GO:0005634	GO_REF:0000002	IEA	InterPro:IPR012099	C	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20240408	InterPro		UniProtKB:Q12019
SGD	S000004096	REA1	located_in	GO:0005634	GO_REF:0000043	IEA	UniProtKB-KW:KW-0539	C	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20240408	UniProt		UniProtKB:Q12019
SGD	S000004096	REA1	located_in	GO:0005634	GO_REF:0000044	IEA	UniProtKB-SubCell:SL-0191	C	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20240408	UniProt		UniProtKB:Q12019
SGD	S000004096	REA1	located_in	GO:0005654	GO_REF:0000044	IEA	UniProtKB-SubCell:SL-0190	C	Huge dynein-related AAA-type ATPase (midasin)	YLR106C|AAA family ATPase midasin|MDN1	protein	taxon:559292	20240408	UniProt		UniProtKB:Q12019
```
Now we have obtained the lines containing REA1 and GO IDs, and we want to know what are they actually. Extract GO IDs related to REA1 first.
```bash
$ awk -F'\t' '$3=="REA1" {print}' gene_association_R64-5-1_20240529.sgd | cut -d $'\t' -f 5 | sort | uniq > REA1_GO_IDs
```
You can go to [AmiGO](https://amigo.geneontology.org/amigo), click Tools & Resources, click Visualization, copy and paste the IDs we just got into the text box, click visualize, and you can get a fancy hierachical plot. See [REA1_GO_visualize](REA1_GO_visualize.pdf). If you find this PDF very fuzzy, you can try to download it and open again. Thanks to AmiGO.

## 7. 🎉 Well Done 🎉
There are some files I have not covered here. You can play with those files using the techniques you just learned.
```bash
NotFeature_R64-5-1_20240529.fasta.gz
orf_coding_all_R64-5-1_20240529.fasta.gz
orf_trans_all_R64-5-1_20240529.fasta.gz
other_features_genomic_R64-5-1_20240529.fasta.gz
rna_coding_R64-5-1_20240529.fasta.gz
```