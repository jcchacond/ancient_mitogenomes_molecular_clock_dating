# Bayesian molecular clock dating using ancient mitogenomes

This repository descibes the steps - including the commands and scripts - to reproduce the methodology for molecular clock dating using ancient mitogenomes published in [Chacón-Duque _et al._, _MBE_ (2025)](https://academic.oup.com/mbe/article/42/4/msaf065/8107989).

### Key methodological and thecnical aspects of our bayesian molecular clock dating approach 

#### Key methodological aspect: "single-sample" dating

In our paper _"we confirm that there is a bias when simultaneously estimating the age of multiple samples at the same time using tip-calibrated dating in_ ```BEAST```_"_ and _"therefore recommend that tip-calibrated dating is performed on one undated sample at a time"_. In other words, if your dataset contains many samples without a finite radiocarbon date (either you got an infinite radiocarbon estimate or you have strong reasons to believe your sample is > 50 000 years old) we strongly recommend you only run one undated sample at a time. In theory you can also date younger samples, but the confidence interval will be so huge that molecular dating can never be considered as a replacement for radiocarbon dating. 

#### Key technical aspect: Using ```BEAST``` efficiently

The software [```BEAST```](https://bmcecolevol.biomedcentral.com/articles/10.1186/1471-2148-7-214) ("Bayesian evolutionary analysis by sampling trees") has been widely used for bayesian molecular clock dating with ancient mitogenomes since the seminal article by [Shapiro _et al_., _MBE_ (2011)](https://academic.oup.com/mbe/article/28/2/879/1212114). Usually, in order to generate the input file (in .XLM format) for this analysis, a Graphic User Interface (GUI) known as ```BEAUti``` needs to be used. This can be particular cumbersome if you want to perform single-sample dating. For this reason, we developed the [```XMLgenerator```](https://github.com/VanssyLi/beastXMLgenerator/tree/main), a ```python``` script to generate the input XML files that significantly reduces the "manual" workload associated with ```BEAUti```. This script was published along with the paper and was written by co-author Wenxi Li under my supervision as a MSc student (for the methodology described here, I attach to this repository the [version 1.0](https://github.com/jcchacond/ancient_mitogenomes_molecular_clock_dating/blob/main/XMLgenerator_v1.py) of the ```XMLgenerator```, also written by Wenxi Li).

> [!NOTE]
> These early versions of the ```XMLgenerator``` were specifically designed and hardcoded for the molecular clock dating of ancient mitogenomes of mammoths. However, we are currently developing a fully automated and customizable pipeline, please stay tuned!

Later on this README file I will show how to use the ```XMLgenerator```. But before I will explain how to prepare the mitogenome sequences and the other input files needed for the script to run.

## Obtaining new mitogenomes from illumina short sequencing reads and merging them with publicly available mitogenomes

### Preparing the sequencing reads

We trimmed adapters, merged forward and reverse reads, and remove merged sequences shorter than 35bp using ```fastp```.

```bash
fastp -i <sampleID>_R1.fastq.gz -I <sampleID>_R1.fastq.gz \
-p -c --merge --merged_out=<sampleID>.trimmed.merged.35bp.fastq.gz -o <sampleID>.trimmed.unmerged.R1.35bp.fastq.gz \
-O <sampleID>.trimmed.unmerged.R2.35bp.fastq.gz -h <sampleID>.trimmed.merged.35bp.html \
-j <sampleID>.trimmed.merged.35bp.json -w <ncores> -l 35
```

If your sample was sequenced across several lanes and/or indexed libraries, you can combine all the merged files using ```zcat```.

### Iterative "_de novo_" assembly of mitogenomes using [```MIA```](https://github.com/mpieva/mapping-iterative-assembler/)

For this purpose we used the reference genome of the closest outgroup to mammoths (Asian elephant) as a guide for the iterative assembly process and all the parameters recommended by ```MIA``` when working with ancient DNA data. We chose this outgroup based on previous recommendations from [van der Valk _et al_. _Nature_ (2021)](https://www.nature.com/articles/s41586-021-03224-9).

```bash
mia -r E.maximus_mito_NC_005129.2.fasta -f <sampleID>.trimmed.merged.35bp.fastq \
-c -C -U -s ancient.submat.txt -i -F -k 14 -m <sampleID>.trimmed.merged.35bp.maln
```

Then, we generated consensus sequences with a sequence agreement of 67% and a minimum per-site depth of coverage of 3x using the [MIA Helper Scripts](https://github.com/aersoares81/mia-helper-scripts).

### Multiple sequence alignment (MSA)

As described in our paper, we used the original alignment published by [van der Valk _et al_. _Nature_ (2021)](https://www.nature.com/articles/s41586-021-03224-9) and we removed a few samples. We ["unaligned"](https://github.com/jcchacond/unalignMSA) this MSA and merged it with all the new consensus mitogenomes.

We aligned all of them together with ```MUSCLE``` v3 using default parameters. We visually inspected the MSA using ```SeaView``` to correct any alignment errors.

## Preparing the input files for the ```XMLgenerator```

### "Lifted" reference mitogenome annotation (GFF3 format)

This is relavant only if you want to use partitions when running your ```BEAST``` analyses.

We "lifted" the reference mitogenome annotation to match the MSA (NC_007596.2.MSA.gff3) using [```Liftoff```](https://github.com/agshumate/Liftoff). In other words, we used the coordinates and features from the reference mitogenome annotation to match these same cooridnates and features to the same reference mitogenome when included in the MSA. We used the woolly mammoth reference mitogenome, to make sure that the lifted coordinates were as accurate as possible. 

```bash
liftoff -g NC_007596.2.gff3 -o NC_007596.2.MSA.gff3 \
-f partition.list.txt target_NC_007596.2.fasta ref_NC_007596.2.fasta
```
The file "NC_007596.2.gff3" (flag -g) contains the original annotation, the file "ref_NC_007596.2.fasta" contains the original reference mitogenome in FASTA format, and the file "target_NC_007596.2.fasta" contains the same reference genome extracted from the MSA (with all the gaps it might include there). Additionally, we provided a plain text file ("partition.list.txt") that contains the names of the partitions to be included in the BEAST analysis, following the same naming conventions of the GFF3 file.

Finally, considering the assembly and alignment limitations on the hypervariable control region, the lower and upper limits of the variable number tandem repeat (VNTR) region (positions 16157 to 16476 in the woolly reference mitogenome; NC_007596.2) were located in the MSA and added to the output annotation file in order to mask it for downstream analyses.

### FASTA files containing the radiocarbon dated samples and one undated sample (FASTA format)

We created subsets of FASTA files each containing one sample with unknown age. This file basically is a multiple sequence alignment containing all the mitogenomes for samples with known radiocarbon ages (the dating reference dataset) and the mitogenome that needs to be dated invidually (sampleID). For this we used the software ```seqkit``` as indicated below:

```bash
# Creating a file containing only the radiocarbon dated samples to be used as reference for dating
seqkit grep -r -p ".*ND" -v aligned.fasta > datedsamples.algined.fasta

# Creating a file that only contains the sample to be dated
seqkit grep -r -p "sampleID" MSA.fasta > sampleID.aligned.fasta

# Concatenaing the files
cat datedsamples.algined.fasta sampleID.aligned.fasta > sampleID.datingset.aligned.fasta
```

### template.xml

This template needs to be created only once using ```BEAUti``` (we are also currently working on scripts to generate the template without using the GUI). I attach the template used for this paper's analyses (created with ```BEAUti```), for reproducibility purposes. Basically, XMLgenerator is going to use this as a draft that contains all the structure needed for the XML file.

### Priors file

A plain text comma-separated (csv) file to set up what tip-date prior distribution will be used for each undated sample. This is an example of how the file looks like for one of the samples:

```csv
taxa,date,prior,mu/lower,sigma/upper,offset
P043_M.pri_NAS_ND,700000,uniform,1000,2000000,0
```

Here we wuse a prior age of 700 thousand years with a range between 1 thousand years and 2 million years, with a uniform (flat) distribution and an offset of 0.

You can make only one file containing all your undated samples, the XMLgenerator will search only for the sample that is included in the alignment 


## How to run  the ```XMLgenerator```

This is a basic example of how this script is used for generating the input XML file for single-sample dating:

```python
python XMLgenerator_v1.py -f sampleID.datingset.aligned.fasta \
 -g NC_007596.2.MSA.gff3 -t template.xml -p priors.csv \
-m 100000000 -l 100000 -o sampleID.uniform.100M.chain1.xml
```

The flags -m and -l indicate the number of iterations (100 million in this case) and burn-in (10%), respectively.

## Now, all ready to run ```BEAST```

This is how I run beast in an HPC environment:

```bash
beast -beagle_cpu -beagle_sse -beagle_double -threads <ncores> sampleID.uniform.100M.chain1.xml
```
