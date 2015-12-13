

# Introduction #

We validated Gowinda with the following approaches:
  * first we compared the results of Gowinda with GoMiner using an unbiased dataset (no gene length bias and no overlapping genes). We found that the p-values reported by Gowinda and GoMiner are highly correlated (Spearman's rank correlation; rho=0.99).
  * second we assessed the impact of the gene length bias and estimated whether Gowinda efficiently deals with the bias. We found that with a biased dataset GoMiner reports an significant enrichment for ~300 GO categories, whereas Gowinda correctly reports none.
  * and third we tested whether Gowinda correctly identifies a significant enrichment for 5 preselected GO terms and found that Gowinda indeed correctly identified these preselected GO categories (FDR<0.01).

## Files ##
The data and scripts of following file allow to repeat the whole validation of Gowinda:
> http://gowinda.googlecode.com/files/validation_files.zip

This archive contains:
  * the Python scripts used in the validation
  * the annotation of _D. melanogaster_ (v5.43)
  * the GO associations obtained from GoMiner for 'CG' gene\_ids

### Results for validating Gowinda ###
As our approach for validating Gowinda contains several random drawing steps we also provide the results we obtained for these random procedures, in order to allow an exact reproduction of our results.
  * biased dataset: http://gowinda.googlecode.com/files/step1_valGoMiner.zip
  * unbiased dataset: http://gowinda.googlecode.com/files/step2_lengthBias.zip
  * known recovery: http://gowinda.googlecode.com/files/step3_knownRecovery.zip




# Validating Gowinda by comparing the results to GoMiner #
As GoMiner does not correct for the gene length bias nor consider overlapping genes it was necessary to use a dataset not showing these two biases. We therefore (i) filtered for not-overlapping genes (ii) filtered for genes that have an associated GO term and (iii) introduced exactly 5 SNPs in each of the thus obtained genes.
Subsequently we randomly picked 1000 SNPs and compared the p-values obtained with Gowinda to the p-values obtained with GoMiner (without FDR correction).

## Obtain a set of genes that are not overlapping and have an associated GO term ##

First obtain a set of non-overlapping genes (gene IDs)
```
python ../scripts/get_nonoverlapping_geneids.py --gtf Flybase.gtf > nonoverlapping_genes.txt
```
Than filter this set for genes that have an associated GO term
```
python ../scripts/get_geneids_havinggocategory.py --go association_gominer.txt --genelist nonoverlapping_genes.txt > geneids_withgoterm_nooverlap.txt
```
## Introduce five SNPs into each of the non-overlapping genes ##
```
python ../scripts/create_snps_for_genes.py --genelist geneids_withgoterm_nooverlap.txt --gtf Flybase.gtf > snps_5pgene.txt
```
## Choose 1000 random SNPs ##
```
cat snps_5pgene.txt|perl -ne 'print if rand()<0.03'|head -1000 > rand_snps_1k.txt
```
## Perform GO term analysis with Gowinda ##
```
/System/Library/Frameworks/JavaVM.framework/Versions/1.6.0/Commands/java -Xmx4g -jar /Users/robertkofler/dev/PopGenTools/gowinda/Gowinda.jar --snp-file snps_5pgene.txt  --candidate-snp-file rand_snps_1k.txt --gene-set-file association_gominer.txt --annotation-file Flybase.gtf --simulations 1000000 --min-significance 1 --gene-definition exon --threads 8 --output-file gowinda_res.txt --mode gene --min-genes 5
```
## Perform GO term analysis with GoMiner ##
For the analysis with GoMiner we need two lists of genes (gene IDs). The first list with must contain the total set of genes and the second list must contain the subset of genes that should be tested for enrichment.
```
cat snps_5pgene.txt| awk '{print $3}'|sort |uniq > gominer_total.txt
cat rand_snps_1k.txt| awk '{print $3}'|sort |uniq > gominer_totest.txt
```

Subsequently use High-Troughput GoMiner: http://discover.nci.nih.gov/gominer/GoCommandWebInterface.jsp

Following the settings that we used with GoMiner
```
totalfile=gominer_total.txt
changedfile=gominer_totest.txt
datasource=FB
organism=7227
evidencecode=all
crossref=true
synonym=true
thresholdtype=BOTH
timeseriesthreshold=0.05
randomization=100
categorymin=5
cim=1
categorymax=10000
rootcategory=all
tf=0
email=xxx@xxx
```

## Merge the results of Gowinda and GoMiner ##
From the GoMiner results you need to extract the file: gominer\_totest.txt.change

```
cat gowinda_res.txt| awk '{print $1,$4}'|sort -k1,1 > merged/tm_gowinda.txt
cat GoMiner/gominer_totest.txt.change |awk 'NF>5{print $1,10^$5}'|perl -pe 's/_[^ ]+//' |sort -k1,1> merged/tm_gominer.txt
join -1 1 -2 1 tm_gominer.txt tm_gowinda.txt > merged_min_win.txt
```
The resulting file contains for every GO category (column 1) the p-value obtained with GoMiner (column 2) and Gowinda (column 3) and may easily be analyzed with R.

## Results for unbiased analysis ##
With an FDR of <0.05 neither GoMiner nor Gowinda found a significant overrepresentation for any GO categories. The p-values (not FDR corrected) for the different GO categories obtained by Gowinda and GoMiner are highly correlated (Spearman's rank correlation; rho=0.9999005; p-value < 2.2e-16), see graph below:

![http://gowinda.googlecode.com/files/cor_gowinda_gominer.png](http://gowinda.googlecode.com/files/cor_gowinda_gominer.png)


Our results (including our random SNPs) for the comparision of Gowinda with GoMiner can be found here: http://gowinda.googlecode.com/files/step1_valGoMiner.zip


# Validating whether Gowinda eliminates the gene length bias #
In GWAS studies longer genes typically have more SNPs and thus also have a higher probability of containing a candidate SNP, which could lead to an overrepresentation of GO categories having longer genes. To test the extent of the bias of GWA datasets and whether Gowinda efficiently corrects for this bias, we introduced one SNP every 100bp into all genes (having a GO category). As overlapping genes may be an additional source of bias in GWA datesets,  overlapping genes were **not** discarded in this analysis. Finally we randomly picked 1000 SNPs and compared the results obtained with GoMiner to the results obtained with Gowinda.

## Extract the genes having an associated GO term ##
```
python ../scripts/extract_geneids_fromgoaassociation.py --go association_gominer.txt > total_genes.txt
```

## Introduce one SNP every 100 bp ##
```
python ../scripts/create_length_biased_genes.py --gtf Flybase.gtf --genelist total_genes.txt > snps_lengthbiased.txt
```

## Randomly draw 1000 SNPs ##
```
cat snps_lengthbiased.txt|perl -ne 'print if rand()<0.004'|head -1000> rand_1k_snps.txt
```
## Perform GO term analysis with Gowinda ##
```
/System/Library/Frameworks/JavaVM.framework/Versions/1.6.0/Commands/java -Xmx4g -jar /Users/robertkofler/dev/PopGenTools/gowinda/Gowinda.jar --snp-file snps_lengthbiased.txt  --candidate-snp-file rand_1k_snps.txt --gene-set-file association_gominer.txt --annotation-file Flybase.gtf --simulations 1000000 --min-significance 1 --gene-definition exon --threads 8 --output-file gowinda_res.txt --mode gene --min-genes 5
```

## Perform GO term analysis with GoMiner ##
First obtain the gene\_ids from the randomly drawn SNPs
```
cat rand_1k_snps.txt|awk '{print $3}'|sort |uniq > subset_genes.txt
```

Than perform the analysis for GO term enrichment with GoMiner: http://discover.nci.nih.gov/gominer/htgm.jsp

We used the following parameters
```
totalfile=total_genes.txt
changedfile=subset_genes.txt
datasource=FB
organism=7227
evidencecode=all
crossref=true
synonym=true
thresholdtype=BOTH
timeseriesthreshold=0.05
randomization=100
categorymin=5
cim=1
categorymax=10000
rootcategory=all
tf=0
email=rokofler@gmail.com
```

## Merge the results of Gowinda and GoMiner ##
From the GoMiner results you need to extract the file subset\_genes.txt.change
```
cat gowinda_res.txt| awk '{print $1,$4}'|sort -k1,1 > merged/tm_gowinda.txt 
cat GoMiner/subset_genes.txt.change |awk 'NF>5{print $1,10^$5}'|perl -pe 's/_[^ ]+//' |sort -k1,1> merged/tm_gominer.txt
join -1 1 -2 1 tm_gominer.txt tm_gowinda.txt > merged_min_win.txt
```

## Results for length biased analysis ##
Although we randomly drew SNPs and would thus not expect an overrepresentation of any GO term, GoMiner reports a significant overrepresentation for 341 GO terms  (FDR <0.05). This demonstrates that the biases of GWA datasets substantially affect GO analysis. By contrast, Gowinda correctly reports that none of the GO terms show any significant overrepresentation (FDR <0.05). The correlation between the p-values obtained with Gowinda and GoMiner (Spearman's rank correlation; rho=0.5610065; p-value < 2.2e-16) is much worse as compared to the unbiased data set. The following graph shows the correlation of the p-values obtained with Gowinda and GoMiner using the length biased dataset.

[![http://gowinda.googlecode.com/files/corBiased_Gominer_Gowinda.png](http://gowinda.googlecode.com/files/corBiased_Gominer_Gowinda.png)]

Our results of this analysis can be found here: http://gowinda.googlecode.com/files/step2_lengthBias.zip

# Validating whether Gowinda identifies preselected GO categories #
To validate that Gowinda correctly identifies overrepresentation of GO terms, we proceeded in two steps. First we randomly picked five GO categories with a small number of genes (between 5 and 10). We choose small GO categories in order to get conservative estimates of the reliability of Gowinda. We than introduced one candidate SNP for ever gene associated with these preselected GO categories. We added randomly drawn candidate SNPs from a length biased SNP dataset (see above) until a total of 1000 candidate SNPs was obtained. Finally we performed an analysis for gene set enrichment using Gowinda and checked whether the five randomly preselected GO categories correctly showed a significant overrepresentation.

## Identify genes having an associated GO term ##
```
python ../scripts/extract_geneids_fromgoaassociation.py --go association_gominer.txt > total_genes.txt
```
## Introduce one SNP every 100bp for genes with an associated GO term ##
```
python ../scripts/create_length_biased_genes.py --gtf Flybase.gtf --genelist total_genes.txt > snps_lengthbiased.txt
```

## Randomly pick five small GO categories ##

We randomly picked five small GO categories having 5 - 10 genes.
```
python ../scripts/pick_random_gocats.py --input association_gominer.txt > random_gocategories.txt
```

In our example we picked the following categories, where the first column is the number of genes in this category:
```
cat random_gocategories.txt|awk '{print $2}' |sort |uniq -c
   5 GO:0007289
   6 GO:0009074
   5 GO:0030173
   5 GO:0040020
   8 GO:0046364
```

## Create 1000 candidate SNPs ##

We than introduce candidate SNPs into the genes associated with the preselected GO categories and add randomly drawn candidate SNPs from the biased dataset until a total of 1000 candidate SNPs was obtained


```
python ../scripts/pick_random_snps_excepttargetlist.py --snp snps_lengthbiased.txt --genelist random_gocategories.txt > random_targetsnps.txt
```

## Perform GO term analysis with Gowinda ##

```
/System/Library/Frameworks/JavaVM.framework/Versions/1.6.0/Commands/java -Xmx4g -jar /Users/robertkofler/dev/PopGenTools/gowinda/Gowinda.jar --snp-file snps_lengthbiased.txt --candidate-snp-file random_targetsnps.txt --gene-set-file association_gominer.txt --annotation-file Flybase.gtf --simulations 1000000 --min-significance 1 --gene-definition exon --threads 8 --output-file gowinda.res --mode gene --min-genes 1
```
## Results for recovery of preselected GO categories ##
With an FDR <0.01 we identified 19 significantly enriched GO categories (see below). All five preselected GO categories were successfully identified (marked with an arrow). However we also found an additional 14 not targeted GO categories, which can be explained by the nesting of GO categories.
For example the not preselected GO category GO:0046165 contains the genes cg32220,cg17725,cg8890,cg10688,cg3495,cg10924,cg8251,cg1516 which are all associated with the preselected GO category: GO:0046364
```
GO:0046165      1.278   9       0.0000010000    0.0002654000    9       14      14      alcohol_biosynthetic_process    cg32220,cg17725,cg10118,cg8890,cg10688,cg3495,cg10924,cg8251,cg1516
GO:0019319      0.714   7       0.0000010000    0.0002654000    7       7       7       hexose_biosynthetic_process     cg17725,cg8890,cg10688,cg3495,cg10924,cg8251,cg1516
=>GO:0009074      0.339   6       0.0000010000    0.0002654000    6       6       6       aromatic_amino_acid_family_catabolic_process    cg7399,cg1555,cg4779,cg9363,cg9362,cg2155
=>GO:0046364      0.740   8       0.0000010000    0.0002654000    8       8       8       monosaccharide_biosynthetic_process     cg32220,cg17725,cg8890,cg10688,cg3495,cg10924,cg8251,cg1516
=>GO:0007289      0.483   5       0.0000010000    0.0002654000    5       5       5       spermatid_nucleus_differentiation       cg3354,cg6998,cg5648,cg12284,cg8827
GO:0009225      0.152   4       0.0000020000    0.0003176667    4       4       4       nucleotide-sugar_metabolic_process      cg32220,cg8890,cg10688,cg3495
GO:0009226      0.152   4       0.0000020000    0.0003176667    4       4       4       nucleotide-sugar_biosynthetic_process   cg32220,cg8890,cg10688,cg3495
GO:0019439      0.464   6       0.0000020000    0.0003176667    6       7       7       aromatic_compound_catabolic_process     cg7399,cg1555,cg4779,cg9363,cg9362,cg2155
GO:0051445      0.966   7       0.0000020000    0.0003176667    7       9       9       regulation_of_meiotic_cell_cycle        cg18543,cg6513,cg9900,cg13584,cg4727,cg4336,cg7719
=>GO:0030173      0.517   5       0.0000040000    0.0005274545    5       5       5       integral_to_Golgi_membrane      cg10580,cg4871,cg10772,cg2448,cg12366
GO:0031228      0.517   5       0.0000040000    0.0005274545    5       5       5       intrinsic_to_Golgi_membrane     cg10580,cg4871,cg10772,cg2448,cg12366
GO:0006558      0.227   4       0.0000070000    0.0008063077    4       4       4       L-phenylalanine_metabolic_process       cg7399,cg4779,cg9363,cg9362
GO:0006559      0.227   4       0.0000070000    0.0008063077    4       4       4       L-phenylalanine_catabolic_process       cg7399,cg4779,cg9363,cg9362
=>GO:0040020      0.622   5       0.0000110000    0.0011997857    5       5       5       regulation_of_meiosis   cg9900,cg13584,cg4727,cg4336,cg7719
GO:0009072      1.411   8       0.0000230000    0.0024266667    8       18      18      aromatic_amino_acid_family_metabolic_process    cg7399,cg10118,cg1555,cg4779,cg17870,cg9363,cg9362,cg2155
GO:0034637      2.181   9       0.0000870000    0.0094014375    9       21      21      cellular_carbohydrate_biosynthetic_process      cg32220,cg17725,cg8890,cg9485,cg10688,cg3495,cg10924,cg8251,cg1516
GO:0006572      0.142   3       0.0000970000    0.0096698421    3       3       3       tyrosine_catabolic_process      cg4779,cg9363,cg9362
GO:0046395      3.068   11      0.0001050000    0.0096698421    11      38      38      carboxylic_acid_catabolic_process       cg7399,cg4586,cg1555,cg4779,cg3626,cg6638,cg9709,cg9527,cg9363,cg9362,cg2155
GO:0016054      3.068   11      0.0001050000    0.0096698421    11      38      38      organic_acid_catabolic_process  cg7399,cg4586,cg1555,cg4779,cg3626,cg6638,cg9709,cg9527,cg9363,cg9362,cg2155
```

Our results for this analysis can be found here: http://gowinda.googlecode.com/files/step3_knownRecovery.zip