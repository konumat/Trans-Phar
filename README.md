# Trans-Phar
Trans-Phar (integration of **Trans**criptome-wide association study and **Phar**macological database)

This software achieves *in silico* screening of chemical compounds, which have inverse effects in expression profiles compared with genetically regulated gene expression of common diseases, from large-scale pharmacological database (Connectivity Map [CMap] L1000 library).

## Overview
![Graphical_abstract](https://user-images.githubusercontent.com/69625255/90997292-21371a00-e5fc-11ea-9c92-8980b999419d.png)

## Publication/Citation

*Human Molecular Genetics* 2021 Konuma T, Ogawa K, Okada Y. "Integration of genetically regulated gene expression and pharmacological library provides therapeutic drug candidates."
https://doi.org/10.1093/hmg/ddab049


## Requirements
- R
- dichromat (R package) (installing by install.packages("dichromat"))
- python 3.X
- scipy
- numpy
- pandas
- math
- cycler
- kiwisolver
- FOCUS (Fine-mapping Of CaUsal gene Sets) as TWAS software


## Installation (Trans-Phar)
In order to get started with **Trans-Phar**, you can just clone this repo as follows;
```bash
git lfs clone https://github.com/konumat/Trans-Phar.git
cd ./Trans-Phar

#unzip QCed Cmap L1000 data
cd ./Cmap_QCeddata
for filename in $( ls *.gz ); do
echo ${filename}
gunzip ${filename}
done

cd ../
```

## Installation (FOCUS)
You have to install [FOCUS (Fine-mapping Of CaUsal gene Sets) soft ware] (https://github.com/bogdanlab/focus) as follows.
For detailed explanations, please visit [the original repository and installing tutorial] (https://github.com/bogdanlab/focus) and [wiki] (https://github.com/bogdanlab/focus/wiki).

When installing FOCUS, please make focus folder under the Trans-Phar folder.

```bash
git clone https://github.com/bogdanlab/focus.git
cd ./focus
python setup.py install
cd ../
```

## Usage
### Step 1: Prepare your input 
All you need is a text file with GWAS summary statistics. (A file extension is .sumstats)

| Column | Column name | Descriptions |
|:-----------:|:-----------:|:------------|
|1|CHR|Chromosome|
|2|SNP|rsID|
|3|BP|BP position|
|4|A1|Effect allele|
|5|A2|Other allele|
|6|MAF|Minor allele frequency (*optional*)|
|7|N|#Samples|
|8|BETA|Beta (effect allele)|
|9|P|P-value|

Please have a look at an example input at `./tutorial_input/Schizo.sumstats`.


### Step 2: Put your input data to predetermined folder (named as Input_GWASsummary)


```bash
mkdir ./Input_GWASsummary
mkdir ./Input_GWASsummary_done
mkdir ./Output

#if you use tutorial GWAS summary data;
gunzip ./tutorial_input/Schizo.sumstats.gz
cp ./tutorial_input/Schizo.sumstats ./Input_GWASsummary
```

### Step 3: Trans-Phar from GWAS summary to chemical compounds in all-in-one script

1) If you input ICD-10 code (for example, F20 for Schizophrenia as below), you will also get gold-standard drug (approved drugs for ICD-10 F20 in ChEMBL and TTD [Therapeutic Target Database] in an output Q-Q plot data.
ICD-10 codes which are not listed in ChEMBL and TTD are not applicable.
The example command is as follows;
```bash
cd ./script
./Trans-Phar.sh F20
```

2) If you need not get gold-standard Q-Q plot, you only enter the example command as follows;
```bash
cd ./script
./Trans-Phar.sh
```


## Output

1) The example TWAS result outputs are as follows (if you use tutorial GWAS data);

```bash
#TWAS results according to each 29 GTEx (v7) tissue and combined files from all 29 tissues at Output/Schizo/TWASresults.

#For Example
cd ../Output/Schizo/TWASresults/ALLTISSUE
less GTEx_Adipose_Subcutaneous.chr_all.focus_shaped.tsv #TWAS result file (shaped), file format is described in https://github.com/bogdanlab/focus/wiki/Fine-mapping-TWAS-associations
#TWAS result png files are also in Output/Schizo/TWASresults/ALLTISSUE
```

2) The example Spearman result outputs are as follows (if you use tutorial GWAS data);

```bash
#Output p-values for Negative Spearmans's correlation tests according to total 308,872 pairs of TWAS tissue - CMap cell - Compunds
#Data of TWAS tissue - CMap cell - Compunds whose P-value < 0.0001 are in Output/Schizo/Spearmanresults/spearman_eachpair_results and Output/Schizo/Spearmanresults/spearman_eachpair_coplots
#For Example
cd ../../Spearmanresults/spearman_totalresults
less ALLpairs_spearmanresults.txt
#Q-Q plot for distribution of these P-value is also in Output/Schizo/Spearmanresults/spearman_totalresults
```

## Acknowledgements
* The original [FOCUS](https://github.com/bogdanlab/focus) was written by Nicholas Mancuso et al.

## Licence
This software is freely available for academic users. Usage for commercial purposes is not allowed.
Please refer to the [LICENCE](https://github.com/konumat/Trans-Phar/blob/master/LICENSE.md/LICENSE.md) page.

