# nfcongeneus - nextflow pipeline to create feature or window consensus from whole genome bam files

This pipeline creates multi-fasta files with consensus sequences of multiple individuals. There are two modes. In the feature mode, the consensi for each feature of a given feature type in the reference annotation is created. In the window mode, the consensus for each section of the genome a given window size is created. Both modes work on bam files, for which an index must be available in a bai file.

## Feature mode

The feature mode creates consensus sequences for each feature of a given feature type, for example each gene of the genome. To find the locations of each of those features, the pipeline uses the information given to it by an annotation file. These are the steps of this mode of nfcongeneus.

1. Creating a list with the locations of all feature.
2. Indexing the bam files.
3. Creating bam files for each feature from the feature list.
4. Calling the consensus using the [doFasta 1](http://www.popgen.dk/angsd/index.php/Fasta) algorithm of angsd.
5. Writing the fasta sequences into one file per feature.

**Usage:**

```bash
nextflow run congeneus.nf --bamdir /PATH/TO/BAM/FILES --anno /PATH/TO/ANNOTATION/FILE --feat feat --outdir /PATH/TO/OUTPUT/DIR/ --mindepth INT
```

**Parameters:**

* --bamdir : Path to directory with sub folders for each individual. Each sub folder must hold one bam and one bai file.
* --anno : Path to an annotation file in .gff format.
* --feat : The type of the feature for which consensus sequences should be created (f.e. gene, CDS, etc.). Can be every feature present in the annotation file.
* --mindepth: The minimum depth at which a base should be determined as consensus.
* --outdir : Path to directory were the output of the pipeline will be written into.

**Output:**

Output are multi-fasta files for each feature. Each file contains the sequences of all individuals for which data was available at that specific locus.

## Window mode

In the window mode nfcongenenus creates consensus sequences for all windows of a given size. For this it needs the mapping reference in fasta format and a list of chromosomes, for which window consensus should be created.

**Usage:**

```bash
nextflow run congeneus.nf --bamdir /PATH/TO/BAM/FILES --ref /PATH/TO/REFERENCE/FILE --win_size INT --outdir /PATH/TO/OUTPUT/DIR/ --mindepth INT
```

**Parameters:**

* --bamdir : Path to directory in with sub folders for each individual. Each sub folder must hold one bam and one bai file.
* --ref : Path to the reference genome assembly fasta file.
* --win_size : Size of the windows for which the consensus should be called (f.e. 10000).
* --skip : Size of the part between the windows that is skipped. (f.e. 20000)
* --outdir : Path to directory were the output of the pipeline will be written into.
* --chr :  Path to a file that contains a list chromosome names of from the reference fasta, which should be included in the consensus calling.

Example chr file:

```text
NW_023416346.1
NW_023416390.1
NW_023416330.1
NW_023416351.1
```

**Output:**

Output are multi-fasta files for each feature. Each file contains the sequences of all individuals for which data was available at that specific window.
