# Molecular clock dating with ancient mitogenomes

Methodology for molecular clock dating using ancient mitogenomes published in [Chacón-Duque _et al._, _MBE_, 2025](https://academic.oup.com/mbe/article/42/4/msaf065/8107989). This repository is complementary to the methodology described in that paper.

One of the novel aspects of this methodology is the use of a series of scripts to reduce the "manual" work that is usually needed to run BEAST.

## Generation of XML files (input for running BEAST)

This is mostly based on the ```python``` script [```XMLgenerator```](https://github.com/VanssyLi/beastXMLgenerator/tree/main) that was published along with the paper and written by co-author Wenxi Li under my supervision as a MSc student. For the methodology described here, the version 1.0 of the ```XMLgenerator``` was used (see script XMLgenerator_v1.py, attached to this repository, also written by Wenxi Li).

> [!NOTE]
> These early versions of the ```XMLgenerator``` were specifically designed for the molecular clock dating of mammoth ancient mitogenomes, however, we are currently developing a fully automated and customizable pipeline, stay tunned!

Example to generate a single file:

```python
python XMLgenerator_v1.py -o sampleID.uniform.chain1.xml \
-t template.xml -f MSA_sampleID.datingset.aligned.fasta \
-g NC_007596.2.liftoff.gff3 -p priors.sampleID.uniform.csv \
-m 10000000 -l 10000
```
