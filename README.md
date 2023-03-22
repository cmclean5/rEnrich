# rEnrich

##      Package Description

This package was built to wrap underlying C/C++ code used in the down-stream analysis results found in several proteomics papers [1-7]. The package contains a number of statistical tests to test the enrichment of a network, or clustered network, given a set of node annotated data.

### Network Enrichment given two annotation sets 

Initally the package was constructed to calculate the probability of gene intersection between two annotation sets, at the network level, using hypergeometric distribution [8]:

```math
P\left(X=\mu_{AB}; \mu_{AB},A[a],B[b],N \right) = \frac{ \binom{A[a]}{\mu_{AB}} \binom{N-A[a]}{B[b]-\mu_{AB}} } { \binom{N}{B[b]} }
```
Where $N$ is taken as the network size, $A[a]$ and $B[b]$ the number of annotations of types $a$ and $b$ in annotation sets $A$ and $B$ respectively, and $\mu_{AB}$ the number of genes (i.e. network nodes) overlapping between the two annotations sets.

#### One-sided Enrichment

```math
\text{p.value$^{enr}_{1}$($\mu_{AB}$)} =
 \displaystyle\sum^{\mu_{AB}}_{i=0}
  \begin{cases}
    P\left( X=i \right)       & P\left(X=i\right) \leq P\left(X=\mu_{AB}\right)\\
    0                         & P\left(X=i\right) > P\left(X=\mu_{AB}\right)
  \end{cases}
```

#### One-Sided Depletion

```math
\text{p.value$^{dep}_{1}$($\mu_{AB}$)} =
 \displaystyle\sum^{\mu_{AB}}_{i=0}
  \begin{cases}
    P\left( X=i \right)       & P\left(X=i\right) \geq P\left(X=\mu_{AB}\right)\\
    0                         & P\left(X=i\right) < P\left(X=\mu_{AB}\right)
  \end{cases}
```

#### Two-sided Enrichment

```math
   \text{p.value$^{enr}_{2}$($\mu_{AB}$)} = 2 \times \text{p.value$^{enr}_{1}$($\mu_{AB}$)} 
```

#### Two-sided Depletion

```math
   \text{p.value$^{dep}_{2}$($\mu_{AB}$)} = 2 \times \text{p.value$^{dep}_{1}$($\mu_{AB}$)} 
```

#### P-values

The code will calculate both one- and two-side p-values for enrichment and depletion. By default the two-sided p-values for enrichment and depletion are returned to the user: `useTwoSided=1` and `useOneSided=0`.         

#### Hypergeometric mean

```math
   \text{mean$_{AB}$} = \frac{A[a] \times B[b]}{N} 
```
<!--Where $A[a]$ is the number of annotation types $a$ in annotation set $A$, $B[b]$ is the number of annotation types $b$ in annotation set $B$.-->

#### Odds Ratio

```math
   \text{OR$_{AB}$} = \frac{ (\mu_{AB} \times (N-A[a] + \mu_{AB} - B[b]) }{ (B[b] - \mu_{AB}) \times (A[a] - \mu{AB}) } 
```

Where the 95% Confidence Intervals (CI) are calculated as [9]:

```math
\text{95\% CI} = \log(\text{OR$_{AB}$}) \pm 1.96 \times \left(\frac{1}{\mu_{AB}} + \frac{1}{(B[b]-\mu_{AB})} + \frac{1}{(A[a]-\mu_{AB})} + \frac{1}{(N-A[a]-\mu_{AB}-B[b])} \right)^{1/2}
```

### False Discovery Rate

Each p-value is corrected for multiple hypothesis testing by selecting one of the following methods: Benjamini and Hochberg FDR (BH) [10], Benjamini and Liu (BL) [11] or Benjamini and Yekutieli (BY) [12]. The default used is (BH): `FDRmeth="BY"`

### Clustered Network Enrichment given one annotation set

The hypergeometric distribution was also used to calculate the significance of enrichment a clustered network given an annotation type:

```math
P\left(X=\mu_{a}; \mu_{a},A,n ,N \right) = \frac{ \binom{A}{\mu_{a}} \binom{N-A}{n-\mu_{a}} } { \binom{N}{n} }
```
Where $N$ is the total number of genes in the network; $n$ the number of nodes in the community; $A[a]$ the total number of annotation types $a$ in annotation set $A$ in the network, and $\mu_{a}$ the number of annotated nodes per community.

#### One-sided Enrichment

```math
\text{p.value$^{enr}_{1}$($\mu_{a}$)} =
 \displaystyle\sum^{\mu_{a}}_{i=0}
  \begin{cases}
    P\left( X=i \right)       & P\left(X=i\right) \leq P\left(X=\mu_{a}\right)\\
    0                         & P\left(X=i\right) > P\left(X=\mu_{a}\right)
  \end{cases}
```
#### One-Sided Depletion

```math
\text{p.value$^{dep}_{1}$($\mu_{a}$)} =
 \displaystyle\sum^{\mu_{a}}_{i=0}
  \begin{cases}
    P\left( X=i \right)       & P\left(X=i\right) \geq P\left(X=\mu_{a}\right)\\
    0                         & P\left(X=i\right) < P\left(X=\mu_{a}\right)
  \end{cases}
```

#### Two-sided Enrichment

```math
   \text{p.value$^{enr}_{2}$($\mu_{a}$)} = 2 \times \text{p.value$^{enr}_{1}$($\mu_{a}$)} 
```

#### Two-sided Depletion

```math
   \text{p.value$^{dep}_{2}$($\mu_{a}$)} = 2 \times \text{p.value$^{dep}_{1}$($\mu_{a}$)} 
```

#### P-values

The code will calculate both one- and two-side p-values for enrichment and depletion. By default the two-sided p-values for enrichment and depletion are returned to the user: `useTwoSided=1` and `useOneSided=0`.    

#### Hypergeometric mean

```math
   \text{mean$_{nA}$} = \frac{ n \times A[a]}{N} 
```
<!--Where $A[a]$ is the number of annotation types $a$ in annotation set $A$, $n$ the number of nodes in the community.-->

#### Odds Ratio

```math
   \text{OR$_{nA}$} = \frac{ (\mu_{a} \times (N-A[a] + \mu_{a} - n) }{ (n - \mu_{a}) \times (A[a] - \mu_{a}) } 
```

Where the 95% Confidence Intervals are calculated as [9]:

```math
\text{95\% CI} = \log(\text{OR$_{nA}$}) \pm 1.96 \times \left(\frac{1}{\mu_{a}} + \frac{1}{(n-\mu_{a})} + \frac{1}{(A[a]-\mu_{a})} + \frac{1}{(N-A[a]-\mu_{a}-n)} \right)^{1/2}
```

### False Discovery Rate

Each p-value is corrected for multiple hypothesis testing by selecting one of the following methods: Benjamini and Hochberg FDR (BH) [10], Benjamini and Liu (BL) [11] or Benjamini and Yekutieli (BY) [12]. The default used is (BH): `FDRmeth="BY"`


### Permuted Clustered Network Enrichment given one annotation set

It's common to calculate the p-values for enrichment of a clustered network given one annotation set. For this reason we provide code to perform a permutation study on such calculated p-values. For example to test those p-values $\leq$ 10-2 for their strength of significance (sig). Our permution study works by recording the percentage of permutated p-values found from every community/annotation combination, lower than or equal to the observed p-value, when $np$ random permutations ($np$ default to 1000 random iteration) of the annotation labels are made. A pseudocount of 1 was added to avoid permutated p-values of zero. P-values found with a strength of significance < 1% are considered statistically significant. Those p-values values are also tested against the more stringent Bonferroni correction at the 0.05 (\*), 0.01 (\*\*) and 0.001 (\*\*\*) significance levels. The default options for the permutation study code is: `runPerm=FALSE`, `setNOP=1000` and `pseudoCount=1.0`.

### Clustered Network Enrichment given two annotation sets

We also tested the significance of the overlap between two annotation sets within a community relative to the annotation set sizes at the network level:

```math
P\left(X=\mu_{ab}; \mu_{ab},n_a, n_b, n, A[a],B[b],N \right) =
\frac{ \binom{A[a] \cup B[b]}{\mu_{ab}} \binom{N-A[a] \cup B[b]}{n-\mu_{ab}} \binom{A[a]}{n_a} \binom{N-A[a]}{n-n_a} \binom{B[b]}{n_b} \binom{N-B[b]}{n-n_b} } {3 \binom{N}{n} }
```
Where $n_a$ and $n_b$ are the number of annotations of types $a$ and $b$, and $\mu_{ab}$ the number of nodes overlapping between the two annotation sets in a community of size $n$. This is similar in spirit to calculating the probability of the intersection distance between two distributions given in eqn (13) pg 8 in [13]. Where we have set $v1 = v2$, and where we have focused on the population overlap relative to the the size of the community, and overlap found in it.

#### One-sided Enrichment

```math
\text{p.value$^{enr}_{1}$($\mu_{ab}$)} =
 \displaystyle\sum^{\mu_{ab}}_{i=0}
  \begin{cases}
    P\left( X=i \right)       & P\left(X=i\right) \leq P\left(X=\mu_{ab}\right)\\
    0                         & P\left(X=i\right) > P\left(X=\mu_{ab}\right)
  \end{cases}
```
#### One-Sided Depletion

```math
\text{p.value$^{dep}_{1}$($\mu_{ab}$)} =
 \displaystyle\sum^{\mu_{ab}}_{i=0}
  \begin{cases}
    P\left( X=i \right)       & P\left(X=i\right) \geq P\left(X=\mu_{ab}\right)\\
    0                         & P\left(X=i\right) < P\left(X=\mu_{ab}\right)
  \end{cases}
```

#### Two-sided Enrichment

```math
   \text{p.value$^{enr}_{2}$($\mu_{ab}$)} = 2 \times \text{p.value$^{enr}_{1}$($\mu_{ab}$)} 
```

#### Two-sided Depletion

```math
   \text{p.value$^{dep}_{2}$($\mu_{ab}$)} = 2 \times \text{p.value$^{dep}_{1}$($\mu_{ab}$)} 
```
#### P-values

The code will calculate both one- and two-side p-values for enrichment and depletion. By default the two-sided p-values for enrichment and depletion are returned to the user: `useTwoSided=1` and `useOneSided=0`. 

#### Hypergeometric mean

```math
   \text{mean$_{ab}$} = \frac{ n \times A[a] \cup B[b]}{N} 
```
<!--Where $A[a]$ is the number of annotation types $a$ in annotation set $A$, $n$ the number of nodes in the community.-->

#### Odds Ratio

```math
   \text{OR$_{ab}$} = \frac{ (\mu_{a} \times (N- A[a] \cup B[b] + \mu_{ab} - n) }{ (n - \mu_{ab}) \times (A[a] \cup B[b] - \mu_{ab}) } 
```

Where the 95% Confidence Intervals are calculated as [9]:

```math
\text{95\% CI} = \log(\text{OR$_{ab}$}) \pm 1.96 \times \left(\frac{1}{\mu_{ab}} + \frac{1}{(n-\mu_{ab})} + \frac{1}{(A[a] \cup B[b]-\mu_{ab})} + \frac{1}{(N-A[a] \cup B[b]-\mu_{ab}-n)} \right)^{1/2}
```

### False Discovery Rate

Each p-value is corrected for multiple hypothesis testing by selecting one of the following methods: Benjamini and Hochberg FDR (BH) [10], Benjamini and Liu (BL) [11] or Benjamini and Yekutieli (BY) [12]. The default used is (BH): `useFDR="BH"`


#### Relative Distance

As part of the clustered network enrichment given two annotation sets analysis, we provided to use to calculate the probability of the intersection distance between two distributions given in eqn (13) pg 8 [13] `relDist=TRUE`.

#### $\chi^2$ Distribution

We first construct the $2\times2$ contingency table (CT):

|                            |                                    |              |                |            |                |         |
| -------------------------- | ---------------------------------- | ------------ | -------------- | ---------- | -------------- | ------- |
| $\mu_{ab}$                 | $n - \mu_{ab}$                     | $n_a$        | $n-n_a$        | $n_b$      | $n-n_b$        | $3n$    |
| $A[a] \cup B[b]-\mu_{ab}$  | $N - A[a] \cup B[b] + \mu_{ab}- n$ | $A[a] - n_a$ | $N-n-n_a-A[a]$ | $B[b]-n_b$ | $N-n-n_b-B[b]$ | $3(N-n)$|
| $A[a] \cup B[b]$           | $N-A[a] \cup B[b]$                 | $A$          | $N-A[a]$       | $B[b]$     | $N-B[b]$       | $3N$    |

The $\chi^2$-(chi-squared) test statistic is then:

```math
\chi^2 = \displaystyle\sum^{Nr}_{i=0} \displaystyle\sum^{Nc}_{j=0} 
\frac{(O_{ij} - E_{ij})^2}{E_{ij}}
```
Where $O_{ij}$ is the observed entry in the contingency table, and $E_{ij}$ the expected value given by $\frac{\sum_{i} O_{i.}}{\sum_{ij} O_{ij}}$

<!--$ \frac{ \sum_{i} O_{i.} \times \sum_{j} O_{.j} } { \sum_{ij} O_{ij} }$ --> 

### Notation

```math
\begin{cases}
N  & \text{Number of nodes in the network}\\
A  & \text{Number of annotation types in set $A$}\\
A[a]  & \text{number of annotation types $a$ in annotation set $A$}\\
B  & \text{Number of annotation types $B$}\\
B[b]  & \text{number of annotation types $b$ in annotation set $B$}\\
\mu_{AB} & \text{Overlap of annotation types $a$ in set $A$ \& annotations types $b$ in set $B$}\\
C  & \text{Number of annotation types C}\\
M  & \text{Number of communities}\\
n  & \text{Number of nodes in a community}\\
n_a  & \text{Number of annotation types $a$ in a community}\\
n_b  & \text{Number of annotation types $b$ in a community}\\
n_c  & \text{Number of annotation types $c$ in a community}\\
\mu_{ab} & \text{Overlap of annotation types $a$ in set $A$ \& annotations types $b$ in set $B$ in a community}\\
np & \text{number of random permutations}
\end{cases}
```

### References

[1] Sorokina, O., Mclean, C., Croning, M.D.R. et al. A unified resource and configurable model of the synapse proteome and its role in disease. Sci Rep 11, 9967 (2021).

[2] Kanellopoulos, A. et al: Aralar Sequesters GABA into Hyperactive Mitochondria, Causing Social Behaviour Deficits, Cell, 180, 1-20 (2020).

[3] Chapelle, J. et al: Dissecting the shared and context-dependent pathways mediated by the p140Cap adaptor protein in cancer and in neurons, Front. Cell Dev. Biol., 15 (2019).

[4] R. Roy, et al: Regional Diversity in the Postsynaptic Proteome of the Mouse Brain, Proteomes, 6, 31, (2018).

[5] A. Alfieri, et al: Synaptic Interactome Mining Reveals p140Cap as a New Hub for PSD Proteins Involved in Psychiatric and Neurological Disorders, Front. Mol. Neurosci., 10, (2017).

[6] E. Fernandez, et al: Arc requires PSD95 for assembly into postsynaptic complexes involved with neural dysfunction and intelligence, Cell Reports, 21, 679-691, (2017).

[7] C. Mclean, X. He, I.T Simpson, D.J Armstrong: Improved Functional Enrichment Analysis of Biological Networks using Scalable Modularity Based Clustering, (2016), J Proteomics Bioinformatics, 9:9-18, doi:10.4172/jpb.1000383.

[8] Pocklington A, Cumiskey D, Armstrong D, Grant S: The proteomes of neurotransmitter receptor complexes from modular networks
with distributed functionality underlying plasticity and behaviour, MSB, 2, (2006).

[9] Szumilas, M. Explaining Odds Ratios, J Can Acad Child Adolesc Psychiatry. 2010 Aug; 19(3): 227–229.

[10] Benjamini, Y., and Hochberg, Y. Controlling the false discovery rate:  a practical and powerful approach to multiple testing.
Journal of the Royal Statistical Society Series B 57 (1995), 289–300.

[11] Benjamini, Y., and Liu, W. A step-down multiple hypotheses testing procedure that controls the false discovery rate under independence. Journal of Statistical Planning and Inference 82 (1999), 163–170.

[12] Benjamini, Y., and Yekutieli, D. (2001). The control of the false discovery rate in multiple testing under dependency. Annals of Statistics, 29, 1165-1188.

[13] Alex T. Kalinka, The probablility of drawing intersections: extending the hypergeometric distribution, arXiv:1305.0717v5, (2014).

[15] M. Galassi et al, GNU Scientific Library Reference Manual (3rd Ed.), ISBN 0954612078.

### TO INSTALL AND BUILD

(1) This package makes use of the GNU Scientific Library (GSL); upon request, we can provide an implementation which is independent of libraries, making use of the subroutines from Numerical Recipes We use the default location of the gsl directory at: /usr/local/include/gsl. If the directory is not installed on the standard search path of your compiler, you will also need to provide its location in the Makefile at: GSLCFLAGS = -I/location/to/your/gsl

## Compile & GSL environment variables setup

1) Install GSL:

1.1) Using mac can run: > brew install gsl

2) Set environment variables (replacing your GSL path in GSL_HOME):

2.1) in .bashrc:    
      
      export GSL_HOME="/opt/homebrew/Cellar/gsl/2.7.1"
      
      export GSL_CFLAGS="${GSL_HOME}/include"
      
      export GSL_LIBS="${GSL_HOME}/lib"
      
      export GSL_CONFIG="${GSL_HOME}/bin/gsl-config"
      
      export PATH="${GSL_HOME}/bin:${PATH}"
      
      export LD_LIBRARY_PATH="${LD_LIBRARY_PATH}:${GSL_HOME}/lib"

3) Run autoconfig.ac:

3.1) > autoconf

4) Build and install:

4.1) > R CMD build rEnrich

4.2) > R CMD INSTALL rEnrich_1.0.tar.gz

5) Run example:

5.1) > cd rEnrich/example/

5.2) > R

    
### TO RUN

The package makes use of the clustering results (of the respective graphs located in the Graphs directory) found the directory 'Clustering', and the annotation files stored in directory 'Annotation'. The directory 'parameterFiles' 

To run the clustering package at the command line type the following:

EXAMPLE 1:

To run cluster enrichment of the Spectral method on the Presynaptic network for disease annotation:  
> ./run -opt 1 -Comfile ../Clustering/PPI_Presynaptic_Published/Spectral_communities_cytoscape.csv -Annofile ../Annotations/flatfile_human_gene2HDO.csv

EXAMPLE 2:

To run cluster enrichment of the Spectral method on the Presynaptic network for the overlap of disease and synaptic functional annotation: 
> ./run -opt 2 -Comfile  ../Clustering/PPI_Presynaptic_Published/Spectral_communities_cytoscape.csv -Annofile ../Annotations/flatfile_human_gene2HDO.csv -Annofile ../Annotations/SCH_flatfile.csv

EXAMPLE 3:

To run network level enrichment for disease, synaptic function and cell type annotation sets:  	
> ./run -opt 4 -Comfile  ../Clustering/PPI_Presynaptic_Published/Spectral_communities_cytoscape.csv -Annofile ../Annotations/flatfile_human_gene2HDO.csv -Annofile ../Annotations/SCH_flatfile.csv -Annofile ../Annotations/celltypes_PMID27991900_L2.csv


##      GNU General Public Licenses v3 

This program is free software: you can redistribute it and/or modify it under the terms of the GNU General Public License as published by the Free Software Foundation, either version 3 of the License, or (at your option) any later version.

This program is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License (GNU_GPL_v3)  along with this program.  If not, see
<http://www.gnu.org/licenses/>.


##      Funding Acknowledgement


This open source software code was developed in part or in whole in the Human Brain Project, funded from the European Union’s Horizon 2020 Framework Programme for Research and Innovation under the Specific Grant Agreement No. 720270 (Human Brain Project SGA1).
