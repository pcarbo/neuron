# neuron

Code and data implementing statistical analyses for *Neuron*
manuscript, "Should single strains suffice?"

### Overview

This repository contains code and data reproducing our analyses of
null allele effects in an F1 panel of 30 inbred strains.

Each F1 panel was obtained by breeding +/- C57BL/6J males with +/+
females from 30 inbred strains. This produced a panel of +/+ and +/-
littermates that are isogenic at genetic loci other than the target
allele. This breeding design was used to independently generate two
cohorts of mice, one for the *Cacna1c* null allele (*n* = 723) and one
for the *Tcf7l2* null allele (*n* = 630).

We collected a variety of physiological and behavioral phenotype data
for each of these panels. Mice in the *Cacna1c* cohort were tested for
anxiety, methamphetamine sensitivity, depression-like behavior, and
acoustic startle response. Mice in the *Tcf7l2* cohort were tested for
several behavioral traits: anxiety, fear conditioning, and
sensorimotor gating, as well as several metabolic traits: body weight,
baseline blood glucose and fasting blood glucose. In total, we
obtained data for 15 phenotypes, 12 of which are behavioral. These
data are stored in CSV files in the [data](data) folder.

The R scripts in the [code](code) folder implement analyses to assess
and summarize the phenotypic effects of two null alleles expressed
from crosses of different inbred strains.

### Citing this repository

If the data or code in this repository is useful for your research
project, please cite our paper:

L. J. Sittig, P. Carbonetto, K. A. Engel, K. S. Krauss, C.
Barrios-Camacho & A. A. Palmer. Should single strains suffice?
Submitted to *Neuron*.

###Code and data licenses

All files in the [code](code) folder are free software: you can
redistribute them under the terms of the
[GNU General Public License](http://www.gnu.org/licenses/gpl.html). All
the files in this folder are part of the source code. This source code
is distributed in the hope that it will be useful, but **without any
warranty**; without even the implied warranty of **merchantability or
fitness for a particular purpose**. See file [LICENSE](code/LICENSE)
for the full text of the license.

All files in the [data](data) folder are released to the public domain
under
[Creative Commons Zero](http://creativecommons.org/publicdomain/zero)
(CC0) license. To the extent possible under law, the authors have
waived all copyright and related or neighboring rights to these data.

### Contents

1. [README.md](README.md): this file.

2. [data/pheno.tcf7l2.csv](pheno.tcf7l2.csv): physiological and
behavioural phenotype data collected in *Tcf7l2* cohort.

3. [data/pheno.cacna1c.csv](pheno.cacna1c.csv): physiological and
behavioural phenotype data collected in *Cacna1c* cohort.

### Physiological and behavioral trait data

*Description of data files goes here.*

### R scripts implementing data analyses

*Description of R scripts goes here.*

### Contact info

For questions and feedback, please contact:

Laura Sittig<br>
Department of Psychiatry<br>
University of California, San Diego<br>
9500 Gilman Drive<br>
La Jolla, California, USA<br>
lsittig@ucsd.edu

### Contributors

Laura J. Sittig<br>
Peter Carbonetto<br>
Kyle A. Engel<br>
Kathleen S. Krauss<br>
Camila Barrios-Camacho<br>
Abraham A. Palmer

