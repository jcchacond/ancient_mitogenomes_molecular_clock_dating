# Molecular clock dating using ancient mitogenomes

Methodology for molecular clock dating using ancient mitogenomes published in [Chacón-Duque _et al._, _MBE_ (2025)](https://academic.oup.com/mbe/article/42/4/msaf065/8107989). This repository is complementary to the methodology described in the paper, providing the necessary commands and scripts to reproduce the main analyses.

## Generation of XML files (input for running BEAST)

One of the aspects that is central to our methodology is the use of a series of commands and scripts to reduce the "manual" workload usually required to run ```BEAST```, which is usually done through a Graphic User Interface (GUI). This can be particular cumbersome if you want to perform molecular clock dating in every sample individually.

The most important script is the [```XMLgenerator```](https://github.com/VanssyLi/beastXMLgenerator/tree/main), a ```python``` script that was published along with the paper and was written by co-author Wenxi Li under my supervision as a MSc student.

For the methodology described here, the version 1.0 of the ```XMLgenerator``` was used (see script XMLgenerator_v1.py, attached to this repository, also written by Wenxi Li).

> [!NOTE]
> These early versions of the ```XMLgenerator``` were specifically designed for the molecular clock dating of mammoth ancient mitogenomes, however, we are currently developing a fully automated and customizable pipeline, stay tuned!

I first provide a basic example on how this script is run, explaining each of the required input files and how they can be obtained:

```python
python XMLgenerator_v1.py -o sampleID.uniform.chain1.xml \
-t template.xml -f MSA_sampleID.datingset.aligned.fasta \
-g NC_007596.2.liftoff.gff3 -p priors.sampleID.uniform.csv \
-m 100000000 -l 100000
```
- **template.xml** \
This template can be created by using BEAST's GUI (known as BEAUti)

- **Multiple Sequence Alignment**
- **Lifted annotation**
- **Priors file**
