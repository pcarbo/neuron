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
[Creative Commons Zero](http://creativecommons.org/publicdomain/zero/1.0)
(CC0) license. To the extent possible under law, the authors have
waived all copyright and related or neighboring rights to these data.

### Contents

1. [README.md](README.md): this file.

2. [pheno.tcf7l2.csv](data/pheno.tcf7l2.csv): physiological and
behavioural phenotype data collected in *Tcf7l2* cohort, stored in CSV
text format. Read below for details.

3. [pheno.cacna1c.csv](data/pheno.cacna1c.csv): physiological and
behavioural phenotype data collected in *Cacna1c* cohort, stored in
CSV text format. Read below for details

4. [read.data.R](code/read.data.R): R source code defining functions
to read and process data stored in text files.

5. [transformation.functions.R](code/transformation.functions.R): R
source code defining some functions used to apply various
transformations to the data.

6. [data.manip.R](code/data.manip.R): R source code defining functions
to process and manipulate the phenotype data.

7. [misc.R](code/misc.R): R source code defining additional functions
used in the statistical analyses.

8. [defaults.tcf7l2.R](code/defaults.tcf7l2.R): R source code
specifying the default analysis settings for each phenotype analyzed
in the *Tcf7l2* cohort.

9. [defaults.cacna1c.R](code/defaults.cacna1c.R): R source code
specifying the default analysis settings for each phenotype analyzed
in the *Cacna1c* cohort.

10. [run.all.anova.R](code/run.all.anova.R): This R script generates
the ANOVA results for all phenotypes analyzed in the *Tcf7l2* and
*Cacna1c* cohorts.

11. [run.anova.analysis.R](code/run.anova.analysis.R): An R Script
used by [run.all.anova.R](code/run.all.anova.R) to implement 3-way
ANOVA for a single phenotype.

### Physiological and behavioral trait data

File [pheno.tcf7l2.csv](data/pheno.tcf7l2.csv) contains physiological
and behavioural phenotype data collected for 630 F1 mice from the
*Tcf7l2* cohort. File [pheno.cacna1c.csv](data/pheno.cacna1c.csv)
contains physiological and behavioral phenotype data for 723 F1 mice
from the *Cacna1c* cohort. These data are stored in comma-delimited
("csv") format, with one line per sample. Missing entries are marked
as "NA", following the convention used in R. Use R function read.pheno
in file [read.data.R](code/read.data.R) to read these data into an R
data frame.

This table includes the following columns:

+ id: Unique number assigned to each mouse.

+ strain: The mother of this mouse is from this inbred mouse stock. In
all cases, the father is a C57BL/6J (B6) mouse. This is an
abbreviation of standard inbred strain id.

+ sex: gender of mouse (M = male, F = female)

+ TCF7L2: whether the *Tcf7l2* gene functions normally (WT), or only
partially functions (HET).

+ CACNA1C: whether the *Cacna1c* gene functions normally (WT), or only
partially functions (HET).

+ bw2, d50bw: Body weight (in g) measured on day 50 of testing.

+ bw3: Body weight (in g) measured on day 100 of testing.

+ fastglucose: Glucose level measured after fasting for 16 hours.

+ baseglucose: Baseline glucose level.

+ centerduration: Duration in center of field during open field
testing.

+ d1centerdur: *Description of this phenotype goes here.*

+ pctdurlight: *Description of this phenotype goes here.*

+ totalactivity: Total activity measured on day 2 during the 30-min
open field test.

+ d1totalactivity: *Description of this phenotype goes here.*

+ d2totalactivity: *Description of this phenotype goes here.*

+ d3.d2totalactivity: *Description of this phenotype goes here.*

+ FCtimeofday: Time of day in which fear conditioning tests were
conducted (either "AM" or "PM").

+ d1tone: Average proportion of freezing on day 1 of the conditioned
fear tests during the two 30-second intervals (180-210 seconds and
240-270 seconds) after exposure to the conditioned stimulus.

+ d2context: Average proportion of freezing over the 30-180 second
interval in response to the test chamber on day 2 of the conditioned
fear tests.

+ d3tone: Average proportion of freezing on day 3 of the conditioned
fear tests during the two 30-second intervals (180-210 seconds and
240-270 seconds) after exposure to the conditioned stimulus.

+ PPIbox: Apparatus used for PPI behavioral testing.

+ PPI6: Average of the inhibition intensity, taken as the ratio of the
prepulse response during 6-dB prepulse trials over the pulse-alone
startle amplitude.

+ PPIstartle, startle: Average startle response during the pulse-alone
trials (with 120-dB pulses).

+ immobdur: Trait measured during the forced swim tests. Defined as
the amount of time in which the mouse is stationary in the water
("immobility").

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

