[![DOI](https://zenodo.org/badge/41293576.svg)](https://zenodo.org/badge/latestdoi/41293576)
# Supporting data and code for: Host plant species and insecticides shape the evolution of genetic and clonal diversity in a major aphid crop pest.
*This repository contains the R code used for the data analyses and production of the figures of the related article*

![alt text](https://j2ejmg.db.files.1drv.com/y4mfs0HpAp-0lm3RXzqAl_6ox6ANJQa-eeY3mIva0J6-lCC_iOKhirczqHbvFa1CbVb0zPHC62CYNYdRDSlUcYTQsepfEoC7Rmwm5mL_yKFWTqgLlbRiQ8RWuDxwEzTYUQqne5s6Sj7aI_ky82MSBhwN4rsbfdgoEmAVv7WUUCsUatxVesPePWVoVl-Sv0hMsYnAh5W2h4q5jLprGqbSMofWQ?width=1584&height=588&cropmode=none)


## Context
Here a small text on why and how we did this study 


## Datasets
In this section, you will find the list of the datasets used in this study. The data files can be found in the "data" folder. For the data tables, the name of the different variables are listed and explained as well. There are 5 data sets used in this study.  

+ **AgrAph5.dat:** the first data set contains the data for all the indivudals analyzed. Each line correspond to one individuals and the following information for each individuals can be found in this table: 
  + *indiv_ID*:
  + *data_batch*:
  + *country	year*:
  + *sampling_date*:
  + *patch_ID*:
  + *host*:
  + *host_corrected*:
  + *longitude*:
  + *latitude*:
  + *MP_27, MP_39, MP_44, MP_5,	MP_7,	MP_23, MP_45,	MP_28,	MP_9,	MP_13,	MP_2,	MP_38,	MP_4,	MP_46*:
  + *MLG*:
  + *MLG_ID*:
  + *missing_data*:
  + *KDR*:
  + *sKDR*:
  + *MACE*:
  + *R81T*:
  + *pot_PB*:
  + *repeated*:
  + *one_MLG*:
  + *several_hosts*:
  + *one_MLG_host*:
  + *several_years*:
  + *one_MLG_year*:
  + *K3_Q1,	K3_Q2, K3_Q3*:
  + *Clust_K3*:
  + *K4_Q1,	K4_Q2,	K4_Q3,	K4_Q4*:
  + *Clust_K4*:

+ **AgrAphout2.str:** the data set summarizing the STRUCTURE runs performed to analyze individual based genetic clusterisation. This data set allows to perform the Delta-K analysis and plot

+ **100str:** a folder containing 15 files that allow the plotting of 100 STRUCTURE runs for each K of interest (from K=2 to K=6)

+ **onebyhost.txt:** Microsatellite results for Myzus grouped by sampling host. For multicopies MLG, only one copy of each MLG by host was kept. The data file is formated as a Genepop file. 

+ **oneclustrap4.txt:** Microsatellite results for Myzus grouped by genetic clusters for K=4. For multicopies MLG, only one copy of each MLG by genetic cluster was kept. The data file is formated as a Genepop file.  


## R scripts
In this section, you will find the list of the different scripts used in the article with a brief description of their purpose.

+ **Agra_load.R:** the script to load the different data sets, functions and packages that are useful for the data analyses and representation in the R environment
+ **Agra_deltaK_fun.R:** functions to perform the delta-K analysis
+ **Agra_div_fun.R:** a script that includes several function to compute a set of different genetic diversity indices
+ **Agra_strplot_fun.R:** a function to plot beautiful STRUCTURE-like plot with several parameters to control the output
+ **Agra_str_100plot.R:** the code to plot the multiple runs of STRUCTURE for K ranging from 2 to 6. 
+ **Agra_str_finalplot.R:** the code to produce Figure 1
+ **Agra_dapc_ana.R:** a script for DAPC analysis of the microsatellite data set. 
+ **Agra_DAStree.R:** a script to compute a dissimilarity matrix as well as to build a neighbour joining tree
+ **Agra_genepop_ana.R:** a script that use some of the functions of the genepop package to compute Fst, Linkage disequilibrium between loci and related tests of significance
+ **Agra_netw.R:** a script to compute dissimilarity matrix between MLG and to build a network of MLG


## Citation
You will soon be able (hopefully) to cite the related study as follow: 
+ Lise Roy, Benoit Barrès, Cécile Capderrey, Frédérique Mahéo, Annie Micoud, Maurice Hullé, Jean-Christophe Simon
[Host plant species and insecticides shape the evolution of genetic and clonal diversity in a major aphid crop pest. *In preparation*.]()

If you want to use (some of) the code found on this page or if you want to cite this repository:
+ Benoit Barrès. bbarres/AgrAphid: [Supporting data and code for: Host plant species and insecticides shape the evolution of genetic and clonal diversity in a major aphid crop pest. Zenodo; 2021.](https://zenodo.org/badge/latestdoi/41293576)
