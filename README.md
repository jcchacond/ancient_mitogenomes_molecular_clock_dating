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
Subsequenty we created subsets containing each one sample with unknown age each (see section XXX below)

## Preparing the input files for the ```XMLgenerator```

### Lifted annotation

We ran ```Liftoff``` to create a new annotation based on the MSA. In other words to lift the coordinates from the reference mitogenome annotation (NC_007596.2.gff3), to the coordinates matching this same genome when included in the multiple alignment (NC_007596.2.MSA.gff3). 

```bash
liftoff -g NC_007596.2.gff3 -o NC_007596.2.MSA.gff3 \
-f partition.list.txt target_NC_007596.2.fasta ref_NC_007596.2.fasta
```

### MSA_sampleID.datingset.aligned.fasta

This file basically is a multiple sequence alignment containing all the mitogenomes for samples with known radiocarbon ages (the dating reference dataset) and the mitogenome that needs to be dated invidually (sampleID). I attach the alignment for the original samples.

### template.xml

This template can be created once using ```BEAUti``` (we are also currently working on scripts to generate the template without using the GUI). I attach the template used for this paper's analyses (created with ```BEAUti```), for reproducibility purposes.

### Priors file

This basically is a plain text file to set up what prior distribution will be used to estimate the tip ages. 



## How to run  the ```XMLgenerator```

This is a basic example of how this script is used for generating the input XML file for single-sample dating:

```python
python XMLgenerator_v1.py -f MSA_sampleID.datingset.aligned.fasta \
 -g NC_007596.2.liftoff.gff3 -t template.xml \
-p priors.sampleID.uniform.csv \
-m 100000000 -l 100000 -o sampleID.uniform.chain1.xml
```


