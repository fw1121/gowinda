

# Introduction #
The following section provides a quick and fast introduction into the analysis of gene set enrichment with Gowinda. For the sake of this tutorial we provide a small artificial data set from _Drosophila melanogaster_ chromosome 2R. In these sample data we identified 348.701 SNPs and a GWAS-like study showed that 279 of those are associated with the trait of interest, lets say body weight. Furthermore the dataset contains an annotation of _D. melanogaster_ from FlyBase and the GO associations from FuncAssociate2. In the following section we demonstrate how to test, whether our set of candidate SNPs shows a significant enrichment for any GO terms.


# Data #

## Tutorial Gowinda ##

The following sample data set contains SNPs, annotation and the GO associations for _D.melanogaster_ chromosome arm 2R.
First, download the sample data for the Walkthrough: http://gowinda.googlecode.com/files/data_gowinda_tutorial.zip

After unzipping the archive you should have the following files:
  * annotation.gtf: Annotation of the genome in the .gtf file format
  * annotation.gff: Annotation of the genome in the .gff file format
  * goassociations\_cg.txt: A GO annotation file obtained from the FuncAssociate2 web page
  * total\_snps.txt: the set of all SNPs used in the GWAS-like  study
  * cand\_snps.txt: the candidate SNPs. Some GWAS-like study showed that these SNPs are associated with a trait of interest, for example wing length.


# Tutorial Gowinda #

**NOTE: Please follow the examples in the given order**

## Example 1: Basic analysis ##

First, let's test whether the candidate SNPs show any enrichment for any GO category when counting every gene only once (`--mode gene`). Furthermore we want a SNP to be associated with a gene if the SNP is overlapping with exons  or introns (`--gene-definition gene`). To speed up the tutorial let's use 8 threads for the simulations and perform 100.000 simulations.

```
java -Xmx4g -jar <path-to-gowinda>/Gowinda.jar --snp-file total_snps.txt --candidate-snp-file cand_snps.txt --gene-set-file goassociations_cg.txt --annotation-file annotation.gtf --simulations 100000 --min-significance 1 --gene-definition gene --threads 8 --output-file results_gene_gene.txt --mode gene --min-genes 1
```

After approximately ~1-2 minutes the results should be  available:
```
GO:0045155      0.050   2       0.0004800000    0.1648850000    2       2       2       electron transporter, transferring electrons from CoQH2-cytochrome c reductase compl
ex and cytochrome c oxidase complex activity        cg13263,cg17903
GO:0006119      0.050   2       0.0004800000    0.1648850000    2       2       8       oxidative phosphorylation       cg13263,cg17903
GO:0009112      0.066   2       0.0009500000    0.2321100000    2       2       16      nucleobase metabolic process    cg7811,cg7171
GO:0006725      0.308   3       0.0027700000    0.3813620000    3       11      73      cellular aromatic compound metabolic process    cg7811,cg7171,cg10501
GO:0005759      0.348   3       0.0031100000    0.3813620000    3       9       55      mitochondrial matrix    cg13263,cg7235,cg17903
```
Please see the manual for a detailed description of the output: http://code.google.com/p/gowinda/wiki/Manual


However, in this example, with an FDR <0.05, no GO categories show a significant enrichment.

## Example 2: Including regulatory regions ##

In the previous analysis only SNPs overlapping with exons or introns were associated with a gene. In order to account for regulatory regions, let's also include 2000 bp upstream and 2000 bp downstream of a gene (`--gene-definition updownstream2000`), for the next analysis.

```
java -Xmx4g -jar <path-to-gowinda>/Gowinda.jar --snp-file total_snps.txt --candidate-snp-file cand_snps.txt --gene-set-file goassociations_cg.txt --annotation-file annotation.gtf --simulations 100000 --min-significance 1 --gene-definition updownstream2000 --threads 8 --output-file results_gene_2000ud.txt --mode gene --min-genes 1
```

Following an example of the results.  Again, no significantly enriched GO categories were found with an FDR <0.05.
```
GO:0005516      1.019   4       0.0027100000    0.8419150000    4       4       29      calmodulin binding      cg5125,cg8475,cg14472,cg6713
GO:0016831      0.301   3       0.0027100000    0.8419150000    3       4       25      carboxy-lyase activity  cg5029,cg7811,cg10501
GO:0022618      0.389   3       0.0047600000    0.8419150000    3       6       33      ribonucleoprotein complex assembly      cg7138,cg4602,cg6137
GO:0016830      0.370   3       0.0048500000    0.8419150000    3       6       34      carbon-carbon lyase activity    cg5029,cg7811,cg10501
GO:0007349      3.040   8       0.0048600000    0.8419150000    8       24      91      cellularization cg9506,cg13263,cg14472,cg10728,cg8049,cg5020,cg17161,cg6137
```

## Example 3: high resolution GO term enrichment ##
As we are getting quite desperate and we know that linkage disequilibrium is decaying quite rapidly in our organism we want to perform a higher resolution analysis in which genes may be counted several times, dependent on the number of candidate SNPs (`--mode snp`). For example a gene containing 5 candidate SNPs is counted 5 times (with the simulations as well as with the candidate SNPs).

```
java -Xmx4g -jar <path-to-gowinda>/Gowinda.jar --snp-file total_snps.txt --candidate-snp-file cand_snps.txt --gene-set-file goassociations_cg.txt --annotation-file annotation.gtf --simulations 100000 --min-significance 1 --gene-definition updownstream2000 --threads 8 --output-file results_snp_2000ud.txt --mode snp --min-genes 1
```

The results should look like in the following example:
```
GO:0006497      0.606   8       0.0000100000    0.0014018182    1       7       45      protein lipidation      cg18810
GO:0043543      0.797   8       0.0000100000    0.0014018182    1       8       36      protein acylation       cg18810
GO:0006119      0.167   7       0.0000100000    0.0014018182    2       2       8       oxidative phosphorylation       cg13263,cg17903
GO:0005794      1.865   12      0.0000100000    0.0014018182    2       12      65      Golgi apparatus cg18810,cg5020
GO:0045155      0.167   7       0.0000100000    0.0014018182    2       2       2       electron transporter, transferring electrons from CoQH2cytochrome c reductase complex and cytochrome c oxidase complex activity        cg13263,cg17903
GO:0005212      0.118   4       0.0000100000    0.0014018182    1       1       1       structural constituent of eye lens      cg16963
GO:0019707      0.464   8       0.0000100000    0.0014018182    1       5       21      protein-cysteine S-acyltransferase activity     cg18810
GO:0019706      0.464   8       0.0000100000    0.0014018182    1       5       21      protein-cysteine S-palmitoleyltransferase activity      cg18810
GO:0016409      0.640   8       0.0000100000    0.0014018182    1       6       25      palmitoyltransferase activity   cg18810
GO:0007349      4.461   16      0.0000100000    0.0014018182    8       24      91      cellularization cg9506,cg13263,cg14472,cg10728,cg8049,cg5020,cg17161,cg6137
GO:0018345      0.464   8       0.0000100000    0.0014018182    1       5       22      protein palmitoylation  cg18810
```
Finally, we have found a significantly enriched GO categories with an FDR <0.05 . For example, for the first GO term (GO:0006497) we find on the average only one candidate SNP, but in our analysis we found 8. Note that we have a significant enrichment for this category although candidate SNPs were only found for a single gene. In this mode we assume that selection is acting on all 8 SNPs independently. However this could also be due to selection acting on a single SNP and the others hitch-hicking along. Therefore results of this analysis mode must be interpreted with care.


# Convert a `.gff` file into a `.gtf` #

In this section we demonstrate how to convert a `.gff` file into a `.gtf` file using the script `Gff2Gtf.py`:

```
python ~/dev/Gowinda/src/scripts/Gff2Gtf.py --input annotation.gff > converted.gtf
```

The output should be similar to the following:

```
2R      FlyBase exon    18442   18629   .       +       .       gene_id "FBgn0262115";
2R      FlyBase exon    18487   18629   .       +       .       gene_id "FBgn0262115";
2R      FlyBase exon    18498   18773   .       +       .       gene_id "FBgn0262115";
2R      FlyBase exon    18681   19484   .       +       .       gene_id "FBgn0262115";
2R      FlyBase exon    18681   18773   .       +       .       gene_id "FBgn0262115";
```

This converted annotation has the `FBgn` gene IDs whereas the GO association file has the `CG` gene IDs. It is thus necessary to obtain a different GO association file. Most conveniently such a GO association file can be directly obtained from FuncAssociate2 (http://llama.mshri.on.ca/funcassociate/download_go_associations). However if many genes in the association file have an unmatched entry in the annotation, or vice versa, we recommend to obtain the GO association file from GoMiner (see following tutorial):

# Obtain a GO association file from GoMiner #
Usually we recommend to obtain a GO association file from FuncAssociate2. However in some cases it may be advantageous to obtain this file from GoMiner instead. We recommend to obtain an GO association file from GoMiner when you have many unmatched gene IDs, between the annotation and the GO association file (the log of Gowinda contains information about the number of matching gene IDs). GoMiner uses a database of synonymous gene IDs and this can be utilized to produce a GO association file suitable for Gowinda. However, GoMiner does not directly offer such a download, therefore we have to perform a 'fake' query for GO term enrichment and subsequently extract the GO association file from the results.

First we need a unique set of the gene IDs (using the `converted.gtf` file from the previous tutorial)
```
cat converted.gtf|awk '{print $10}'|sed 's/";//'|sed 's/"//'|sort |uniq > total_genes.txt
```

Than we start a query at High Throughput GoMiner (http://discover.nci.nih.gov/gominer/relationship.jsp):
  * Step 1: provide the total\_genes.txt file shown above
  * Step 2: again provide the total\_genes.txt file shown above
  * Step 3,4,5: choose your species and evidence level
  * Step 6: make sure that `Cross Reference` and `Synonym` are checked
  * Step 7,8: leave the devaults
  * Step 9: set smalles category to 1
  * Step 10: set largest category to 10.000
  * Step 11: Choose all/gene ontology
  * Step 13: Provide email address and submit query

Your parameters should be similar to the following:
```
totalfile=total_genes.txt
changedfile=total_genes.txt
datasource=FB
organism=7227
evidencecode=all
crossref=true
synonym=true
thresholdtype=BOTH
timeseriesthreshold=0.05
randomization=100
categorymin=1
cim=1
categorymax=10000
rootcategory=all
tf=0
```


Finally convert the resulting `.gce` file into the file format compatible with Gowinda:
```
python <path-to-conversiontool>/GoMiner2FuncAssociate.py --input work1338826265/total_genes.txt1338826265.dir/total_genes.txt.dir/total_genes.txt.change.gce > association_gominer.txt
```

# Dealing with synonymous gene IDs #
Many genes have a tremendous amount of aliases i.e.: alternative gene IDs. As a requirement for using Gowinda the gene IDs of the annotation have to match the gene IDs from the gene set file.
Due to an extensive amount of synonymous gene IDs this may not always be the case.
Therefore we provide a script which allowes to synchronized the gene IDs of the annotation and the  gene set file given a list of synonymous gene IDs. The script basically updates the gene IDs of the annotation ('.gtf') with the corresponding gene ID from the gene set. A list of synonymous gene IDs is required. Such a list of synonymous gene IDs may for example be obtained from a .gff file, see next example:

## Obtain a list of synonymous gene IDs from a 'gff' file ##
```
python ExtractGeneAliasesFromGff.py --input annotation.gff > synonymous.genids.txt
```

Which should result in a file similar to the following example
```
FBgn0259968	Sfp60F	CG42478	Seminal fluid protein 60F	Seminal fluid protein 60F
FBgn0011694	PebII	PEB	PEBmeII	CG2665	PEBII	PEB2	Protein ejaculatory bulb II
FBgn0004181	Peb	FBgn0028493	FBgn0044740	CG2668	PEB-me	PEB	P38	peb	cDNA:GH06048	anon-WO0140519.140	ejaculatory bulb protein	Protein ejaculatory bulb
FBgn0035094	CG9380
FBgn0001325	Kr	FBgn0001251	FBgn0015732	Kuppel	kruppel	kr	CG3340	Krpple	If	Irregular facets	Kruppel	Dm-Kr	krueppel	Lruppel	KR	KRUPPEL	Krueppel	Kruppel
FBgn0050429	CG30429
FBgn0053680	CG33680
FBgn0050428	CG30428	FBgn0064772	BcDNA:SD21514
```

Every line contains all synonymous gene IDs for a single gene

## Synchronize the gene IDs of the annotation with the gene sets ##

The following example used the file 'converted.gtf' from a previous tutorial. Converted.gtf contains the gene IDs as FBGN. After updating they will be replaced with the corresponding 'CG' gene IDs.
```
python SynchronizeGtfGeneIDs.py --gtf converted.gtf --gene-set goassociations_cg.txt --gene-aliases synonymous.genids.txt > converted_cg.gtf
```

Example of the '.gtf' before updating:
```
2R      FlyBase exon    209942  210220  .       -       .       gene_id "FBgn0260798";
2R      FlyBase exon    210281  210443  .       -       .       gene_id "FBgn0260798";
2R      FlyBase exon    223882  224088  .       -       .       gene_id "FBgn0260798";
```

Example of the '.gtf' after updating:
```
2R      FlyBase exon    209942  210220  .       -       .       gene_id "CG40129";
2R      FlyBase exon    210281  210443  .       -       .       gene_id "CG40129";
2R      FlyBase exon    223882  224088  .       -       .       gene_id "CG40129";
```

Note that the IDs of genes not having a synonymous gene ID will not be updated!