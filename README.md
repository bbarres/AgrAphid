[![DOI](https://zenodo.org/badge/41293576.svg)](https://zenodo.org/badge/latestdoi/41293576)
# Supporting data and code for: Host plant species and insecticides shape the evolution of genetic and clonal diversity in a major aphid crop pest.
*This repository contains the R code used for the data analyses and production of the figures of the related article*

![alt text](https://j2ejmg.db.files.1drv.com/y4mfs0HpAp-0lm3RXzqAl_6ox6ANJQa-eeY3mIva0J6-lCC_iOKhirczqHbvFa1CbVb0zPHC62CYNYdRDSlUcYTQsepfEoC7Rmwm5mL_yKFWTqgLlbRiQ8RWuDxwEzTYUQqne5s6Sj7aI_ky82MSBhwN4rsbfdgoEmAVv7WUUCsUatxVesPePWVoVl-Sv0hMsYnAh5W2h4q5jLprGqbSMofWQ?width=1584&height=588&cropmode=none)


## Context
In this study, we wanted to explore the genetic diversity and structure of *Myzus persicae* populations sampled from different host. We also studied the evolution through time of the different genetic group as well as the variation of target site resistance against 3 different classes of insecticide (pyrethroids, carbamates and neonicotinoids). 


## Datasets
In this section, you will find the list of the data sets used in this study. The data files can be found in the "data" folder. For the data tables, the name of the different variables are listed and explained as well. There are 5 data sets used in this study.  

+ **AgrAph5.dat:** the first data set contains the data for all the indivudals analyzed. Each line correspond to one individuals and the following information for each individuals can be found in this table: 
  + *indiv_ID*: individual's ID, this is a unique string of character
  + *data_batch*: the global dataset can be divided in several sub-projects (AgrAphid, RAtransect, Rpp and Zepeda)
  + *country*: country in which the sample has been collected
  + *year*: year when the sample was collected
  + *sampling_date*: sampling date for the AgrAphid aerial trap
  + *patch_ID*:the ID of the field/orchard/region/… from which the sample has been collected
  + *host*: host from which the individual was collected
  + *host_corrected*: for the table in which only one MLG was kept, an additional column for the host is added. It is the same host if the MLG has been
found only on one host, but it is turned to "several_host" if the MLG was found on different hosts
  + *longitude*: longitude coordinates in the WGS84 geodesic system
  + *latitude*: latitude coordinates in the WGS84 geodesic system
  + *MP_27,	MP_39,	MP_44,	MP_5,	MP_7,	MP_23,	MP_45,	MP_28,	MP_9,	MP_13,	MP_2,	MP_38,	MP_4,	MP_46*: Microsatellite data (coded as the length of the alleles with 3 digits for each allele)
  + *MLG*: MultiLocus Genotypes
  + *MLG_ID*: MLG's ID, a unique string of character for each different MLG
  + *missing_data*: number of null allele for microsatellites data
  + *KDR*: Genotype at the KDR position (coded as the sequence of the codon of interest)
  + *sKDR*: Genotype at the sKDR position (coded as the sequence of the codon of interest)
  + *MACE*: Genotype at the MACE position (coded as the sequence of the codon of interest)
  + *R81T*: Genotype at the R81T position (coded as the sequence of the codon of interest)
  + *RG3loci*: 
  + *pot_PB*: indicates if we suspect a problem with this individual
  + *repeated*: Is the MLG found several times in the global dataset (1=yes, 0=no)
  + *one_MLG*: Individuals to keep in order to have only one MLG in the dataset
  + *several_hosts*: Is the MLG found on several hosts? (1=yes, 0=no)
  + *one_MLG_host*: Individuals to keep in order to have only one MLG for each host (ie there could be several times the same MLG, but sampled on different host)
  + *several_years*: Is the MLG found on several years? (1=yes, 0=no)
  + *one_MLG_year*: Individuals to keep in order to have only one MLG for each year (ie there could be several times the same MLG, but sampled on different year)
  + *K3_Q1,	K3_Q2,	K3_Q3*: Q-matrix, obtained by averaging the runs belonging to the major solution over 100 runs with K=3
  + *Clust_K3*: which cluster does the individual belongs too when K=3. A threshold of 0.7 was applied (Clust_primary = red cluster, Clust_second = green cluster, Clust_Wild = blue cluster, NA = not assign to a cluster)
  + *K4_Q1,	K4_Q2,	K4_Q3,	K4_Q4*: Q-matrix, obtained by averaging the runs belonging to the major solution over 100 runs with K=4
  + *Clust_K4*: which cluster does the individual belongs too when K=4. A threshold of 0.7 was applied (Clust_primary = red cluster, Clust_Oil = green cluster, Clust_tobacoil = yellow cluster, Clust_Wild = blue cluster, NA = not assign to a cluster)
  + *cor_MP_27,	cor_MP_39,	cor_MP_44,	cor_MP_5,	cor_MP_7,	cor_MP_23,	cor_MP_45,	cor_MP_28,	cor_MP_9,	cor_MP_13,	cor_MP_2,	cor_MP_38,	cor_MP_4,	cor_MP_46*: Microsatellite corrected data (coded as the length of the alleles with 3 digits for each allele). The "correction" consist of homogenizing the allele for individuals identified as belonging to the same MLG 

+ **AgrAphout2.str:** the data set summarizing the STRUCTURE runs performed to analyze individual based genetic clusterisation. This data set allows to perform the Delta-K analysis and plot

+ **100str:** a folder containing 15 files that allow the plotting of 100 STRUCTURE runs for each K of interest (from K=2 to K=6)

+ **onebyhost.txt:** Microsatellite results for *Myzus persicae* grouped by sampling host. For multicopies MLG, only one copy of each MLG by host was kept. The data file is formatted as a Genepop file. 

+ **oneclustrap4.txt:** Microsatellite results for *Myzus persicae* grouped by genetic clusters for K=4. For multicopies MLG, only one copy of each MLG by genetic cluster was kept. The data file is formatted as a Genepop file.  


## R scripts
In this section, you will find the list of the different scripts used in the article with a brief description of their purpose.

+ **Agra_load.R:** the script to load the different data sets, functions and packages that are necessary for the data analyses and representation in the R environment. 
+ **Agra_deltaK_fun.R:** functions to perform the delta-K analysis. 
+ **Agra_div_fun.R:** a script that includes several function to compute a set of different genetic diversity indices. 
+ **Agra_strplot_fun.R:** a function to plot beautiful STRUCTURE-like plot with several parameters to control the output. 
+ **Agra_str_100plot.R:** the code to plot the multiple runs of STRUCTURE for K ranging from 2 to 6. 
+ **Agra_str_finalplot.R:** the code to produce Figure 1. 
+ **Agra_dapc_ana.R:** a script for DAPC analysis of the microsatellite data set. 
+ **Agra_DAStree.R:** a script to compute a dissimilarity matrix as well as to build a neighbour joining tree. 
+ **Agra_genepop_ana.R:** a script that use some of the functions of the GENEPOP package to compute Fst, Linkage disequilibrium between loci and related tests of significance. 
+ **Agra_netw.R:** a script to compute dissimilarity matrix between MLG and to build a network of MLG. 


## Citation
You will soon be able (hopefully) to cite the related study as follow: 
+ Lise Roy, Benoit Barrès, Cécile Capderrey, Frédérique Mahéo, Annie Micoud, Maurice Hullé, Jean-Christophe Simon
[Host plant species and insecticides shape the evolution of genetic and clonal diversity in a major aphid crop pest. *submitted*.]()

If you want to use (some of) the code found on this page or if you want to cite this repository:
+ Benoit Barrès. bbarres/AgrAphid: [Supporting data and code for: Host plant species and insecticides shape the evolution of genetic and clonal diversity in a major aphid crop pest. Zenodo; 2021.](https://zenodo.org/badge/latestdoi/41293576)
