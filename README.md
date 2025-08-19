# Bayesian molecular clock dating using ancient mitogenomes

Methodology for molecular clock dating using ancient mitogenomes published in [Chacón-Duque _et al._, _MBE_ (2025)](https://academic.oup.com/mbe/article/42/4/msaf065/8107989). This repository is complementary to the methodology described in the paper, providing the necessary commands and scripts to reproduce the main analyses.

## How to generate the input files for bayesian molecular clock dating (with ```BEAST```)

The software [```BEAST```](https://bmcecolevol.biomedcentral.com/articles/10.1186/1471-2148-7-214) ("Bayesian evolutionary analysis by sampling trees") has been widely used for bayesian molecular clock dating with ancient mitogenomes since the seminal paper by [Shapiro _et al_., _MBE_ (2011)](https://academic.oup.com/mbe/article/28/2/879/1212114). A central aspect of our methodology is the use of a series of commands and scripts to reduce the "manual" workload required to generate an input file (in .XML format) for running ```BEAST```, which is usually done through a Graphic User Interface (GUI) known as ```BEAUti```. This can be particular cumbersome if you want to perform molecular clock dating in many samples individually ("single-sample" dating), which is the way we recommend to do it based on our results.

For this reason, we developed the [```XMLgenerator```](https://github.com/VanssyLi/beastXMLgenerator/tree/main), a ```python``` script to generate the input XML files that was published along with the paper and was written by co-author Wenxi Li under my supervision as a MSc student (for the methodology described here, I attach to this repository the [version 1.0](https://github.com/jcchacond/ancient_mitogenomes_molecular_clock_dating/blob/main/XMLgenerator_v1.py) of the ```XMLgenerator```, also written by Wenxi Li).

> [!NOTE]
> These early versions of the ```XMLgenerator``` were specifically designed and hardcoded for the molecular clock dating of ancient mitogenomes of mammoths. However, we are currently developing a fully automated and customizable pipeline, please stay tuned!

This is a basic example of how this script is used for generating the input XML file for single-sample dating:

```python
python XMLgenerator_v1.py -f MSA_sampleID.datingset.aligned.fasta \
 -g NC_007596.2.liftoff.gff3 -t template.xml \
-p priors.sampleID.uniform.csv \
-m 100000000 -l 100000 -o sampleID.uniform.chain1.xml
```

But, how can you obtain all these input files?

## MSA_sampleID.datingset.aligned.fasta

This file basically is a multiple sequence alignment containing all the mitogenomes for samples with known radiocarbon ages (the dating reference dataset) and the mitogenome that needs to be dated invidually (sampleID). I attach the alignment for the original samples.

### Preparing the sequencing reads for mitogenome generation

We trimmed adapters, merged forward and reverse reads, and remove merged sequences shorter than 35bp using ```fastp```.

```bash
fastp -i <sampleID>_R1.fastq.gz -I <sampleID>_R1.fastq.gz \
-p -c --merge --merged_out=<sampleID>.trimmed.merged.35bp.fastq.gz -o <sampleID>.trimmed.unmerged.R1.35bp.fastq.gz \
-O <sampleID>.trimmed.unmerged.R2.35bp.fastq.gz -h <sampleID>.trimmed.merged.35bp.html \
-j <sampleID>.trimmed.merged.35bp.json -w <ncores> -l 35
```

If your sample was sequenced across several lanes and/or indexed libraries, you can combine all the merged files using ```zcat```.

### Iterative "_de novo_" assembly of mitogenomes using [```MIA```](https://github.com/mpieva/mapping-iterative-assembler/)

For this purpose we used the reference genome of the closest outgroup to mammoths (Asian elephant) as a guide for the iterative assembly process and all the parameters recommended by MIA for ancient DNA data.

```bash
mia -r E.maximus_mito_NC_005129.2.fasta -f <sampleID>.trimmed.merged.35bp.fastq \
-c -C -U -s ancient.submat.txt -i -F -k 14 -m <sampleID>.trimmed.merged.35bp.maln
```

Then, to generate a consensus sequence filtering with a consensus and a minimum per-site depth of coverage, we used the [MIA Helper Scripts](https://github.com/aersoares81/mia-helper-scripts).

### multiple sequence alignment (MSA)

We had a set of mitogenomes with unknown ages or infinite radiocarbon estimates that we wanted to date. In order to avoid generating MSAs for every sample, we first aligned all of them together with ```MUSCLE``` and then created the subsets containing each individual sample with unknown age.

## Lifted annotation

We ran ```Liftoff``` to create a new annotation based on the MSA. In other words to lift the coordinates from the reference mitogenome annotation (NC_007596.2.gff3), to the coordinates matching this same genome when included in the multiple alignment (NC_007596.2.MSA.gff3). 

```bash
liftoff -g NC_007596.2.gff3 -o NC_007596.2.MSA.gff3 \
-f partition.list.txt target_NC_007596.2.fasta ref_NC_007596.2.fasta
```

## template.xml

This template can be created once using ```BEAUti``` (we are also currently working on scripts to generate the template without using the GUI). I attach the template used for this paper's analyses (created with ```BEAUti```), for reproducibility purposes.

## Priors file

This basically is a plain text file to set up what prior distribution will be used to estimate the tip ages. 
