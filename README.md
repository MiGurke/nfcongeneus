# nfcongeneus - nextflow pipeline to create feature or window consensus from whole genome bam files

This pipeline creates multi-fasta files with consensus sequences of multiple individuals. There are two modes. In the feature mode, the consensus for each feature of a given feature type in the reference annotation is created. In the window mode, the consensus for each section of the genome given window size is created. Both modes work on bam files, for which an index must be available in a bai file.

## Feature mode

The feature mode creates consensus sequences for each feature of a given feature type, for example each gene of the genome. To find the locations of each of those features, the pipeline uses the information given to it by an annotation file. These are the steps of this mode.

1. Creating a list with the locations of all feature.
2. Indexing the bam files.
3. Creating bam files for each feature from the feature list.
4. Calling the consensus using the [doFasta 1](http://www.popgen.dk/angsd/index.php/Fasta) algorithm of angsd.
5. Writing the fasta sequences into one file per feature.

**Usage:**

```bash
nextflow run congeneus.nf --bamdir /PATH/TO/BAM/FILES/ --anno /PATH/TO/ANNOTATION/FILE.gff --feat feat --outdir /PATH/TO/OUTPUT/DIR/ --mindepth INT -profile mfn
```

**Parameters:**

* --bamdir : Path to directory with sub folders for each individual. Each sub folder must hold one bam and one bai file. Path must end on "/".
* --anno : Path to an annotation file in .gff format.
* --feat : The type of the feature for which consensus sequences should be created (f.e. gene, CDS, etc.). Can be every feature present in the annotation file.
* --mindepth: The minimum depth at which a base should be determined as consensus.
* --outdir : Path to directory were the output of the pipeline will be written into.
* -profile mfn : Include if running on the cluster at the Museum für Naturkunde to use the automatic job submission. If not inlcuded in the command it will run locally. 

**Output:**

Output are multi-fasta files for each feature. Each file contains the sequences of all individuals for which data was available at that specific locus.

## Window mode

In the window mode nfcongenenus creates consensus sequences for all windows of a given size. For this it needs the mapping reference in fasta format and a list of chromosomes, for which window consensus should be created.

**Usage:**

```bash
nextflow run congeneus.nf --bamdir /PATH/TO/BAM/FILES/ --ref /PATH/TO/REFERENCE/FILE --win_size INT --outdir /PATH/TO/OUTPUT/DIR/ --mindepth INT
```

**Parameters:**

* --bamdir : Path to directory in with sub folders for each individual. Each sub folder must hold one bam and one bai file. Path must end on "/".
* --ref : Path to the reference genome assembly fasta file.
* --win_size : Size of the windows for which the consensus should be called (f.e. 10000).
* --skip : Size of the part between the windows that is skipped. (f.e. 20000)
* --outdir : Path to directory were the output of the pipeline will be written into.
* --chr :  Path to a file that contains a list chromosome names of from the reference fasta, which should be included in the consensus calling.
* --mindepth: The minimum depth at which a base should be determined as consensus.
* -profile mfn : Include if running on the cluster at the Museum für Naturkunde to use the automatic job submission. If not inlcuded in the command it will run locally. 

Example chr file:

```text
NW_023416346.1
NW_023416390.1
NW_023416330.1
NW_023416351.1
```

**Output:**

Output are multi-fasta files for each feature. Each file contains the sequences of all individuals for which data was available at that specific window.

## Utils

nfcongeneus inlcudes some python and bash scripts in the folder utils that might be useful in further processing the output of the main pipeline. This section will explain those different scripts. 

### NfracSeqRemover

This script removes sequences from the multi fasta output that have more then a given proportion of missing data and adds the exact proportion of missing data of each sequence to the fasta header.

**Usage:**

```bas
python3 NfracSeqRemover.py -f /PATH/OF/MULTI/FASTA/FOLDER/ -m INT -o /PATH/TO/OUTPUT/DIRECTORY/
```

**Parameters:**

* -f : Path to the folder in which all of the multi fasta consensus file (e.g. the ouput of the main pipeline) are stored. 
* -m : The maximum proportion of missing data (number of N / number of Bases) in one sequence. 
* -o : Path to an output folder were the new multi fasta files with only sequences passing the filter will be written into. 

**Output:**

Filtered multi fasta files. 

### MinIndivRemover

This script filters multi fasta files depending on whether they hold a minimum number of sequences/samples. The exact number of sequences in each file be added to the file names. It should be used after the NfracSeqRemover, which removes sequences. (Before that, all multi fasta files that come out of nfcongeneus should have all sequences/samples given to it as bam files.) 

**Usage:**

```bas
python3 MinIndivRemover.py -f /PATH/OF/MULTI/FASTA/FOLDER/ -m INT -o /PATH/TO/OUTPUT/DIRECTORY/
```

**Parameters:**

* -f : Path to folder with multi fasta files. (Ideally the output directory of the NfracSeqRemover.)
* -m : Minimum number of sequences/samples/individual that must be in on mutli fasta file. 
* -o : Path to an output folder were the new multi fasta files with only sequences passing the filter will be written into. 

**Output:**

Filtered multi fasta files. 

### IntronRemover

When nfcongeneus is used to call the consensus of genes in the feature mode, then this script can be used to remove the introns from the genes and output multi fasta files that contain only the coding regions of the genes. This script must be used before the NfracSeqRemover and MinIndivRemover.

**Usage:**

```bash
python3 IntronRemover.py -a PATH/TO/ANNOTATION/FILE.gff -f /PATH/OF/MULTI/FASTA/FOLDER/ -o /PATH/TO/OUTPUT/DIRECTORY/
```

**Parameters:**

* -a : Path to the annotation file that was already used by the main pipeline to create the gene consensus seqeunces. 
* -f : Path to the folder in which all of the multi fasta consensus file (e.g. the ouput of the main pipeline) are stored. 
* -o : Path to an output folder were the new multi fasta files without introns will be written into. 

**Output:**

Filtered multi fasta files. 

### RaxML_phylogenies

This script automatically and calculates maximum likelihood phylogenies for all multi fasta files in a folder using [RAxML](https://cme.h-its.org/exelixis/web/software/raxml/). It currently does no bootstrapping, but the RAxML command in line 49 can be customized for different analyses. The script also parallelizes the RAxML processes to speed up the whole thing.  

**Usage:**

```bash
bash RaxML_phylogenies.sh -f /PATH/OF/MULTI/FASTA/FOLDER/ -o /PATH/TO/OUTPUT/DIRECTORY/ -n INT -t INT -m INT
```

**Parameters:**

* -f : Path to folder with multi fasta files.
* -o : Path to folder were the RAxML output files will be written into. 
* -n : The number of RAxML jobs that should be started in parallel. 
* -t : The number of threads each RAxML job should use. 
* -m : The number of maximum likelihood searches RAxML should do to find the the best phylogeny. 

**Output**:

RAxML output files for each consensus multi fasta file. The name of the input file will be the file ending of each corresponding RAxML output file.

### ASTRAL_prep

This script prepares the output of the RaxML_phylogenies to be used by [ASTRAL](https://github.com/smirarab/ASTRAL). Practically, it changes the tip names to readable by ASTRAL and writes all phylogenies into one file so that it can be directly used as input for ASTRAL. 

**Usage:**

```bash
python3 ASTRAL_prep.py -f /PATH/TO/INPUT/FOLDER/ -o /PATH/TO/OUTFILE.new
```

**Parameters:**

+ -f : Path to the folder that has all the RAxML output, so the output folder of RaxML_phylogenies. 
+ -o : Path to an output file.

**Output:**

A newick formatted tree file that contains all phylogenies that were in the RAxML output folder. The tip names of each tree are reduced to only the sample names. 

### Nfraxml
Nfraxml is a small nextflow pipeline that does the same job as the RaxML_phylogenies script. I found that for large data sets it is more convienient to use nextflow with multtiple jobs started than one huge job that occupies many ressources for a long time on a cluster. 

**So far this is only tested on the FU curta, not on the mfn cluster** (It is however very likely that it also works there and I inlcuded a config file for it.)

**Usage**

```bash
nextflow run raxml.nf -profile CHAR --fastadir /PATH/TO/DIR/WITH/FASTA/FILES/ --outdir /PATH/TO/OUTPUT/DIR/ --MLsearches INT
```
**Parameters:**

+ -profile : Either local, curta or mfn depending on where the pipeline is supposed to run. 
+ --fastadir : Path to a directory in which the fasta files are located. Each file with the ending '.fasta' will be analyzed. 
+ --outdir : The folder where the raxml output files will stored in. 
+ --MLsearches : And integer telling raxml how many maximum likelihood searches should be carried out for each tree generated.

**Output:**
Output is a folder filled with all raxml output files for each tree generated.
