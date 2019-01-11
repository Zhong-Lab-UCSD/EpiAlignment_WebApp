## Using the python template program to access the EpiAlignment web service

We provided EpiAlignment_PyClient.py, a python template program for users to access EpiAlignment from their own server, submit multiple jobs and parse the output. The user may also modify the code and incorporate it into their own programs.

### Python platform and dependencies
The current verion of EpiAlignment_PyClient.py only supports Python 2. The HTTP library [Requests] (http://docs.python-requests.org/en/master/) is required.

### Usage
The basic usage of the template program is:

```
python EpiAlignment_PyClient.py sampleSheet.txt
```
where sampleSheet.txt is a tab-delimited text file with input information.

A sample input file named sampleSheet.txt can be found in the zip file. This sample uses the enhancer mode and will take 1 to 2 minutes to finish.

The template program also provides two auxiliary functions that allow users to view preset data in EpiAlignment:

(1) View all paired ENCODE / public epigenomic datasets in EpiAlignment

```
python EpiAlignment_PyClient.py --public_data
```
A tab-delimited file named "EpiAlign_publicData.txt" will be generated. The file contains details of preset epigenomic datasets in EpiAlignment. 

(2) Search for gene clusters

```
python EpiAlignment_PyClient.py --find_gene_cluster geneid
```
where gene id can be either a gene symbol/partial gene symbol, or an Ensembl id. The command will write gene clusters with fully-matched and partially-matched names to the standard output (stdout). You may try the following commands to see the results:

```python
# find gene clusters with partially-matched names.
python EpiAlignment_PyClient.py --find_gene_cluster GNG
# find gene clusters with fully-matched names.
python EpiAlignment_PyClient.py --find_gene_cluster GNG7
``` 

### Submitting jobs to EpiAlignment

#### Input
The template program parses a sample sheet to get input data of each job. The sample sheet starts with a header, specifying the cotent that should be put in each column. There are 24 fields in the sample sheet. Typically, the user only need to fill in part of the first 11 fields, whereas the others are parameters with default values.

A file named sampleSheet_headerOnly.txt with the header and default parameter values is provided in the same zip file with EpiAlignment\_PyClient.py.

The 24 fields are:
>alignMode, searchRegionMode, genomeAssembly\_query, genomeAssembly\_target,
>encodeData\_query, encodeData\_target, speciesPeak\_query, speciesPeak\_target,
>speciesInput\_1, speciesInput\_2,clusters
>promoterUp, promoterDown, enhancerUp, enhancerDown,
>epiweight, paras, paramu,parak,
>piA,piC,piG,piT,pi1

The program will parse the information following the same order.

##### alignMode and searchRegionMode
These two fields specified the alignment mode and submode to be used. Valid values and their correspondences with the modes on the website is shown below:

| alignMode value| searchRegionMode value| alignment mode on the website  | submode on the website |
| --------- |:-----------:| :-----:|:-----:|
| promoter  | genecluster | promoter |Search a gene cluster|
| promoter  | genomeregion| promoter |Define target regions with a BED file / a gene list|
| enhancer  | homoregion  | enhancer |Use homologous regions in this species|
| enhancer  | genomeregion| enhancer | Define target regions with a BED file|

##### genomeAssembly
The genomeAssembly\_query and genomeAssembly\_target fields can be "hg38" or "mm10". 

##### encodeData and speciesPeak
These fields specify the input peak files. If you would like to use a preset ENCODE/public data, please put the ENCODE/GEO id in the encodeData\_query and encodeData\_target fields and leave the speciesPeak\_query, speciesPeak\_target blank. For example: 

| encodeData\_query | encodeData\_target| 
| :---------: |:-----------:| 
| GSM1673960  | GSM1674003 |

You may also provide your own input files. To do this, put the file names in speciesPeak\_query, speciesPeak\_target and leave the encodeData\_query and encodeData\_target fields empty. Please note that if your files are not in the same folder as the program, you need to put the absolute of our files here.

##### speciesInput\_1 and speciesInput\_2
These fields specify the input query and target files. Please put the file names in these two fields. For promoter mode, the two files can be either BED6 files or gene name lists. For enhancer mode, only BED6 files are acceptable. 

Please note that if you are using the "promoter-genecluster" or "enhancer-homoregion" mode, only the speciesInput\_1 is required. You don't have to provide the second input file in this case.

##### clusters
This field only needs to be filled when alignMode is "promoter" and searchRegionMode is "genecluster". A valid gene cluster name in the format of "Cluster_XXX" needs to be entered.

##### fields 16 to 24: parameters
These fields specify the parameters to be used for the alignments. Default values are available in the downloaded sample sheet. You may keep them as they are unless you'd like to adjust them.  

#### Output
One output file will be generated for one job and put in the same folder as the program. The output file is tab-delimited, containing alignment scores and region coordinates as the results provided on the website.

 