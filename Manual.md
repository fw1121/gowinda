

# Introduction #

Gowinda is a multi-threaded Java application that allows an unbiased analysis of gene set enrichment for Genome Wide Association Studies. Classical analysis of gene set (e.g.: Gene Ontology) enrichment assumes that all genes are sampled independently from each other with the same probability. These assumptions are violated in Genome Wide Association (GWA) studies since (i) longer genes typically have more SNPs resulting in a higher probability of being sampled and (ii) overlapping genes are sampled in clusters. Gowinda has been specifically designed to test for enrichment of gene sets in GWA studies. We show that Gene Ontology (GO) tests on GWA data could result in a substantial number of false positive GO terms. Permutation tests implemented in Gowinda eliminate these biases, but maintain sufficient power to detect enrichment of GO terms.


For a validation of Gowinda please see: http://code.google.com/p/gowinda/wiki/Validation
In the validation we show that Gowinda yields highly reliable results (by comparision with GoMiner) efficiently corrects for the gene length bias while still identifying significantly overrepresented GO categories. We also demonstrate that the gene length bias has a tremendous influence on the GO analysis potentially causing an substantial amount of false positive GO categories.

# Requirements #
Gowinda requires the following software:
  * Java 6 or higher
Furthermore the following input files are required
  * a file containing the annotation of the genome in .gtf
  * a gene set file, containing for every gene set (e.g.: GO category) a list of the associated gene IDs
  * a file containing the total set of SNPs
  * a file containing the candidate SNPs

## Annotation (.gtf) ##
The annotation must be in the `.gtf` format: http://mblab.wustl.edu/GTF22.html
Gowinda however only requires the attribute `gene_id`, the attribute `transcript_id` is optional. Following an example of a minimal annotation file:

```
2L      FlyBase exon    8193    9484    .       +       .       gene_id "CG11023";
2L      FlyBase exon    8193    8589    .       +       .       gene_id "CG11023";
2L      FlyBase exon    8668    9484    .       +       .       gene_id "CG11023";
2L      FlyBase exon    9839    11344   .       -       .       gene_id "CG2671";
2L      FlyBase exon    11410   11518   .       -       .       gene_id "CG2671";
```

Furthermore the feature column (column 3) must either contain `exon` or `CDS`, all other features will be ignored. Column 2 (source), column 6 (score), and column 8 (offset)  have no influence on the results.

It is not necessary that the entries are unique, Gowinda is internally reducing multiple copies of the same entry to a single one.

Gowinda does not support the `.gff` format as it less transparent as the `.gtf` format (e.g.: in order to obtain the gene ID for a given exon it is necessary to traverse the hierarchy exon -> mRNa -> gene) and has several sources of inconsistencies (exons without parent (mRNA), exons with several parents).
We however provide a script for converting a `.gff` file into a `.gtf` file:

http://gowinda.googlecode.com/files/Gff2Gtf.py

The user may thus assure that the annotation is correctly converted before providing it to Gowinda.

## Gene set ##
The gene set file for Gene Ontology terms may be directly obtained from the download section of FuncAssociate2 (http://llama.mshri.on.ca/funcassociate/download_go_associations) or indirectly from GoMiner (see tutorial: http://code.google.com/p/gowinda/wiki/Tutorial)

Following an example of a gene set file for Gene Ontology
```
GO:0000002      mitochondrial genome maintenance        CG11077 CG33650 CG4337 CG5924
GO:0000003      reproduction    CG10112 CG10128 CG1262 CG13873 CG14034 CG15117 CG15616 CG1656 CG17011 CG17097 CG17673 CG17799 CG17843 CG1803 CG2665 CG2668 CG2852 CG30450 CG30473 CG31680 CG31
704 CG31872 CG31883 CG31941 CG32203 CG32498 CG32667 CG33943 CG34033 CG34034 CG34102 CG3662 CG3801 CG42461 CG42462 CG42466 CG42468 CG42469 CG42472 CG42474 CG42475 CG42477 CG42478 CG42479 CG42
480 CG42482 CG42483 CG42485 CG42564 CG42602 CG42603 CG42604 CG42605 CG42606 CG42607 CG42608 CG42609 CG4546 CG4706 CG4847 CG4986 CG6555 CG6690 CG6917 CG7157 CG8137 CG8462 CG8622 CG8626 CG8982
 CG9024 CG9029 CG9074 CG9111 CG9334 CG9997
GO:0000009      alpha-1,6-mannosyltransferase activity  CG8412
GO:0000010      trans-hexaprenyltranstransferase activity       CG10585 CG31005
GO:0000012      single strand break repair      CG4208 CG5316
GO:0000014      single-stranded DNA specific endodeoxyribonuclease activity     CG10215 CG10670 CG10890 CG2990
GO:0000015      phosphopyruvate hydratase complex       CG17654
GO:0000017      alpha-glucoside transport       CG30035
```
The file consists of three tab-delimited columns
  * column 1: the GO category
  * column 2: the description of the GO category. Note that spaces are allowed but no tabs
  * column 3: a space separated list of gene ids for the given GO category

**NOTE:** the `gene_id`s of the annotation need to be identical to the `gene_id`s of the GO association file. The case is not considered as Gowinda internally converts all `gene_id`s to lower-case.
However a problem may arise as many `gene_id`s have synonyms. One strategy to deal with this problem is to use the GO association file from GoMiner (see tutorial http://code.google.com/p/gowinda/wiki/Tutorial)

## Total SNP file ##

This file must contain all the SNPs used for the GWAS, in a simple tab-delimited file format. Following you can see the simplest example
```
2L	117081
2L	117082
2L	144234
2L	252591
2L	283388
2L	318365
2L	320282
2L	378118
2L	378119
2L	476447
```
  * column 1: the chromosome
  * column 2: the position

Gowinda ignores all additional columns after column 2, thus it is for example also possible to provide a '.mpileup' file (http://samtools.sourceforge.net/) as shown in the followinge example
```
2R	2299	N	4	TTT^FT	AAAA	4	TTT^FT	AAAA
2R	2300	N	5	AAAA^FA	AAAAA	5	AAAA^FA	AAAAA
2R	2301	N	6	TTTTT^FT	AAAAAA	6	TTTTT^FT	AAAAAA
2R	2302	N	7	TTTTTT^FT	AAAAAAA	7	TTTTTT^FT	AAAAAAA
2R	2303	N	8	TTTTTTT^FT	AAAAAAAA	8	TTTTTTT^FT	AAAAAAAA
```

Gowinda also ignores entries starting with a '#' thus a '.vcf' file (http://vcftools.sourceforge.net/specs.html) may also be provided as shown in the following example:
```
##fileformat=VCFv4.0
##fileDate=20090805
##...
##...
#CHROM POS     ID        REF ALT    QUAL FILTER INFO                              FORMAT      NA00001        NA00002        NA00003
20     14370   rs6054257 G      A       29   PASS   NS=3;DP=14;AF=0.5;DB;H2           GT:GQ:DP:HQ 0|0:48:1:51,51 1|0:48:8:51,51 1/1:43:5:.,.
20     17330   .         T      A       3    q10    NS=3;DP=11;AF=0.017               GT:GQ:DP:HQ 0|0:49:3:58,50 0|1:3:5:65,3   0/0:41:3
20     1110696 rs6040355 A      G,T     67   PASS   NS=2;DP=10;AF=0.333,0.667;AA=T;DB GT:GQ:DP:HQ 1|2:21:6:23,27 2|1:2:0:18,2   2/2:35:4
```


## Candidate SNP file ##

the same applies as for the total SNP file, except that the candidate SNPs must be a subset of the total SNPs.

# Gowinda #

For an example of how to use Gowinda with a sample dataset please see the tutorial: http://code.google.com/p/gowinda/wiki/Tutorial

## Parameters ##

Gowinda has the following input parameters

  * `--annotation-file`: a file containing the annotation for the species of interest. Only the `.gtf` format is accepted (see above). Mandatory parameter
  * `--gene-set-file`: a file containing for every gene set (e.g.: Gene Ontology term) the associated genes; For the file format see above; Mandatory parameter
  * `--snp-file`: a file containing the total set of SNPs that were used for the GWAS. For the file format see above; Mandatory parameter
  * `--candidate-snp-file`: a file containing the candidate SNPs that show some association with the trait of interest. For the file format see above; Mandatory parameter
  * `--output-file`: where to store the output. Mandatory parameter
  * `--mode`: As a major feature Gowinda offers two main analysis modes either `snp` or `gene`; For description see below; Optional parameter; default=`gene`
  * `--gene-definition`: As another major feature Gowinda allows to adjust the SNP to gene mapping, i.e.: to decide which genes are associated with a given SNP. For example the user may decide that only SNPs being located in an exon are associated with the corresponding gene. See detailed description below; possible arguments: `exon`, `cds`, `utr`,`gene`, `upstreamDDDD`, `downstreamDDDD`, `updownstreamDDDD`; Mandatory parameter
  * `--simulations`: the number of simulations that should be performed. For more information on the number of simulation see below
  * `--min-genes`: filter for GO categories having at least `--min-genes` number of genes; This parameter is for example useful to remove small GO categories having only one associated gene; Optional parameter; default=1
  * `--min-significance`: only report GO categories having  after FDR correction `p-value <= --min-significance`; Optional parameter; default=1.0
  * `--detailed-log`: switch to the detailed log mode. The IDs of genes present in the GO association file but not present in the annotation will be displayed. Also the progress of the simulations will be shown in steps of 10.000; Optional parameter
  * `--threads`: the simulations of Gowinda utilize multi-threading. Adjust the number of threads to use. Optional parameter; default=1
  * `--help`: show the help for Gowinda; Optional parameter

### SNP to gene mapping ###

SNP to gene mapping refers to the assignment of genes to SNPs. For example it is possible to associate only genes with a SNP, when an exon is overlapping with the SNP. As another example - when regulatory regions should be considered - a SNP may be considered as associated with a gene when it's position is less than 500 bp upstream of the gene. Gowinda only loads the `.gtf` features `exon` and `cds`, all other features are ignored. The following features are thus internally computed from the features `exon` and `cds`:
  * `exon`: SNPs within exons are associated with genes
  * `cds`: SNPs within CDS are associated with genes
  * `utr`: SNPs within 5'-UTR and 3'-UTR are associated with genes. Caluclated as `exon - cds`
  * `gene`: SNPs within exons or introns are associated with genes. Internally the distance from the start position of the first exon to the end position of the last exon is computed
  * `upstreamDDDD`: in addition to exons and introns also the DDDD bases upstream the start of a gene are considered for mapping a SNP to a gene. DDDD must be replaced with an arbitrary number. This method requires the strand information.
  * `downstreamDDDD`: in addition to exons and introns also the DDDD bases downstream of the end of the gene are considered for mapping a SNP to a gene. DDDD must be replaced with an arbitrary number. This method requires the strand information.
  * `updownstreamDDDD`: in addition to  exons and introns also the DDDD bases upstream the start of a gene and the DDDD bases downstream of the end of a gene are considered for mapping a SNP to a gene. DDDD must be replaced with an arbitrary number. This method does NOT require strand information.


**NOTE:** We found examples (in the annotation of _D. melanogaster_) where different exons from the same gene are located on opposite strands, in this case Gowinda uses the majority vote, i.e.: the strand being supported by the most exons


### simulations ###
The number of simulations has a direct influence on the minimum achievable p-value; For example with 1.000.000 simulations the minimum achievable p-value is 1.0e-6; Increase the number of simulations if you need a higher accuracy.

**Note** increasing the number of simulations will not result in higher numbers of significant GO categories! It will just increase the accuracy.

We recommend to start with 100.000 simulations, this should take about 4 minutes (with 8 CPUs).
When a higher accuracy is required the number of simulations may be increased to 1 million (~30 minutes) or 10 million (~5h).


## Algorithm ##
Gowinda does not reproduce the exact pattern of linkage disequilibrium (LD) between SNPs but offers two complementary test strategies making two extreme assumptions about LD.
The first strategy (`--mode snp`) assumes linkage equilibrim between SNPs and the second strategy (`--mode gene`) assumes complete linkage disequilibrium between SNPs within a gene, thus basically treating every gene as a large haplotype.

### Analysis mode SNP ###
This strategy is based on the extreme assumption that all SNPs are independent, i.e.: in linkage equilibrium. In this mode of operation, Gowinda randomly samples the same number of SNPs as candidate SNPs. Subsequently the genes corresponding to the sampled SNPs are identified. Finally the significance of enrichment is estimated from the empirical null distribution and FDR is calculated (see below). The number of SNPs is thus constant (=number of candidate SNPs) in each simulation. However the number of genes associated with the randomly sampled SNPs may vary between the simulations. Furthermore, when using this strategy, genes will be counted multiple times according to the number of candidate/random SNPs. This is based on the rationale that - when SNPs are in linkage equilibrium - every candidate SNP constitutes an independent observation. However any remaining linkage between variants will cause a bias in this strategy, that is the significance of enrichment will be overestimated.


### Analysis mode gene ###
This strategy is based on the assumption that all SNPs within a gene are completely linked. In this mode of operation, Gowinda first computes the number of genes that are corresponding to the candidate SNPs (=candidate genes). Subsequently, Gowinda randomly samples SNPs until the corresponding number of genes equals the number of candidate genes. Finally the significance of overrepresentation is estimated from the empirical null distribution and FDR is calculated (see below). In this strategy the number of randomly sampled genes is constant in each simulation whereas the number of randomly sampled SNPs may vary. This strategy basically treats genes as large haplotype blocks. However, if SNPs are not in complete linkage this analysis strategy may result in underestimating the significance of enrichment.

**Note**: This mode assumes complete linkage of SNPs within genes and does not account for LD between genes.

**Note**: It is however possible to extent the size of the haplotype block of genes by using the option `--gene-definition updownstream`  (see above).

We recommend this mode per default.


### Calculate the p-value of gene set enrichment ###

After the simulations have been completed, Gowinda derives an empirical null distribution of gene abundance for every gene set. An example of such an empirical null distribution for a single gene set and 1 million simulations may look like in the following example:

![http://gowinda.googlecode.com/files/gowinda_edist.png](http://gowinda.googlecode.com/files/gowinda_edist.png)

Let the number of simulations be `S`, the gene set for which the p-value should be calculated be `g`, the number of genes found for the given gene set within a single simulation be c<sup>g</sup><sub>s</sub> and the number of candidate genes be c<sup>g</sup><sub>cand</sub> than the p-value of enrichment (p<sup>g</sup>) for the given gene set can be calculated according to the following equation:

![http://gowinda.googlecode.com/files/equation_pvalue.png](http://gowinda.googlecode.com/files/equation_pvalue.png)

In the example shown above c<sup>g</sup><sub>cand</sub> would be 18 genes, `S` would be 1 million simulations and the sum of c<sup>g</sup><sub>s</sub> larger or equal than 18 (c<sup>g</sup><sub>cand</sub>) about 47.000. Thus the p-value p<sup>g</sup> for `intracellular protein transport` would be 47.000/1.000.000 = 0.047

This p-value is calculated for every gene set.


### FDR correction ###

We used the empirical FDR correction described in 'Elements of Statistical Learning' (2009): Trevor Hastie, Robert Tibshirani and Jerome Friedman, 2nd edition, pp687-690, (http://www-stat.stanford.edu/~tibs/ElemStatLearn/)

In detail, Gowinda implements the following algorithm for calculating the FDR:
Let `G` be the number of gene sets and `S` the number of simulations. Than an observed (R<sub>obs</sub>) and an expected (R<sub>exp</sub>) count of gene sets having p-values smaller or equal than a threshold (`P` = the p-value which should be FDR corrected) can be computed.  The FDR can subsequently simply be calculated by dividing R<sub>exp</sub> by R<sub>obs</sub>.


![http://gowinda.googlecode.com/files/fdr_equations.png](http://gowinda.googlecode.com/files/fdr_equations.png)


## Output ##

Following an example output:
```
GO:0045155      0.050   2       0.0004800000    0.1648850000    2       2       2       electron transporter, transferring electrons from CoQH2cytochrome c reductase complex and cytochrome c oxidase complex activity        cg13263,cg17903
GO:0006119      0.050   2       0.0004800000    0.1648850000    2       2       8       oxidative phosphorylation       cg13263,cg17903
GO:0009112      0.066   2       0.0009500000    0.2321100000    2       2       16      nucleobase metabolic process    cg7811,cg7171
GO:0006725      0.308   3       0.0027700000    0.3813620000    3       11      73      cellular aromatic compound metabolic process    cg7811,cg7171,cg10501
GO:0005759      0.348   3       0.0031100000    0.3813620000    3       9       55      mitochondrial matrix    cg13263,cg7235,cg17903
```

  * column 1: the GO term
  * column 2: on the average this number of genes are found per simulation for the given GO category. In `--mode gene` every gene is only counted once whereas in `--mode snp` a single gene may be counted several times dependent on the SNP
  * column 3: using the candidate SNPs this number of genes was found for the given GO category. In `--mode gene` every gene is only counted once whereas in `--mode snp` a single gene may be counted several times dependent on the SNP
  * column 4: p-value (uncorrected for multiple testing)
  * column 5: FDR (p-value after adjustment for multiple testing)
  * column 6: the number of genes (uniq) found for the given GO category
  * column 7: the number of genes that could at most be found for the given GO category, i.e.: genes of the given GO category that have an corresponding entry in the annotation file and contain at least one SNP
  * column 8: total number of genes for the given GO category in the GO association file
  * column 9: description of the given GO term
  * column 10: comma separated list of the gene\_ids found for the given GO category


# Scripts #

## Gff2Gtf.py ##

The script `Gff2Gtf.py` may be used to convert a `.gff` file into a `.gtf` file. The script scans the `.gff` file for the features exon, cds, and mRNA all other features are ignored. These features require the fields 'ID=' and 'Parent='. The ID of the corresponding gene will be  infered for every exon/cds from the hierarchy exon->mRNA->gene.
Exons not having a parent mRNA will be ignored. Furthermore if exons have several parent mRNAs only the first one will be used.

Following are the parameters:
  * --input: the path to the gtf file
  * --help: Display help for the script

Redirect the output (`.gtf`) to any file of your choice (see tutorial)

## Gominer2FuncAssociate.py ##

This script converts a `.gce` file from GoMiner to a GO association file compatible with Gowinda and FuncAssociate. See the tutorial for how to obtain a `.gce` file and how to use this script.

Following a small example from a `.gce` file:
```
GO:0000003_reproduction	FBGN0033033	137	137	1.000000	0.000000	4633	4633.0	1.000000
GO:0000003_reproduction	FBGN0250842	137	137	1.000000	0.000000	4633	4633.0	1.000000
GO:0000003_reproduction	FBGN0260934	137	137	1.000000	0.000000	4633	4633.0	1.000000
GO:0000003_reproduction	FBGN0027538	137	137	1.000000	0.000000	4633	4633.0	1.000000
GO:0000003_reproduction	FBGN0015584	137	137	1.000000	0.000000	4633	4633.0	1.000000
GO:0000003_reproduction	FBGN0023172	137	137	1.000000	0.000000	4633	4633.0	1.000000
GO:0000003_reproduction	FBGN0016053	137	137	1.000000	0.000000	4633	4633.0	1.000000
GO:0000003_reproduction	FBGN0011236	137	137	1.000000	0.000000	4633	4633.0	1.000000
GO:0000003_reproduction	FBGN0005655	137	137	1.000000	0.000000	4633	4633.0	1.000000
GO:0000003_reproduction	FBGN0034430	137	137	1.000000	0.000000	4633	4633.0	1.000000
```

Following an example of the resulting GO assocation file
```
GO:0030381      chorion-containing_eggshell_pattern_formation   FBGN0003731 FBGN0015754 FBGN0004400 FBGN0003733
GO:0045887      positive_regulation_of_synaptic_growth_at_neuromuscular_junction        FBGN0003317
GO:0005795      Golgi_stack     FBGN0033500 FBGN0027558 FBGN0034697 FBGN0034025 FBGN0033339
GO:0034660      ncRNA_metabolic_process FBGN0034401 FBGN0028426 FBGN0034177 FBGN0003062 FBGN0034735 FBGN0026722 FBGN0033741 FBGN0002069 FBGN0035064 FBGN0262560 FBGN0033454 
FBGN0053505 FBGN0033375 FBGN0020766 FBGN0033900 FBGN0027091 FBGN0027079 FBGN0028744 FBGN0259937 FBGN0011824
GO:0042826      histone_deacetylase_binding     FBGN0033748
GO:0043139      5'-3'_DNA_helicase_activity     FBGN0261850
```

The parameters of this script are:
  * --input: The `.gce` file which should be converted
  * --help: Display help for the script

# Use a custom pathway database #

Although Gowinda will probably mostly be used to test for overrepresentation of genes in GO categories, it is also possible to use different databases. In this case it is only necessary to provide a custom 'gene set' file like in the following example:
```
category_A    description of category A    gene1 gene2 gene3 gene4 gene5
category_B    description of category B    gene1 gene6 gene7 gene8
category_C    description of category C    gene10 gene11 gene12
```

  * column1: name of the custom category
  * column2: description of the custom categories
  * column3: list of genes associated with the custom category
**Note** that columns are separated by a tab whereas entries within columns need to be separated by a space!

# Links #
Several tools have been developed for the analysis of pathway enrichement in GWA. For an excellent review see Wang et. al. (2010) "Analysing biological pathways in genome-wide association studies" (http://www.nature.com/nrg/journal/v11/n12/abs/nrg2884.html)

  * ALIGATOR http://www.ncbi.nlm.nih.gov/pubmed/19539887
  * i-GSEA4GWAS http://www.ncbi.nlm.nih.gov/pubmed/20435672
  * GenGen http://www.ncbi.nlm.nih.gov/pubmed/17966091
  * GeSBAP http://www.ncbi.nlm.nih.gov/pubmed/19502494/
  * GRASS http://www.ncbi.nlm.nih.gov/pubmed/20560206
  * GSA-SNP http://www.ncbi.nlm.nih.gov/pubmed/20501604
  * GSEA-SNP http://www.ncbi.nlm.nih.gov/pubmed/18854360
  * PLINK set-test http://www.ncbi.nlm.nih.gov/pubmed/17701901
  * SNP ratio test http://www.ncbi.nlm.nih.gov/pubmed/19620097