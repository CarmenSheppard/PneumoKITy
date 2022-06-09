# README for PneumoKITy

## Table of content

* [Introduction](https://github.com/CarmenSheppard/PneumoKITy#Introduction)
* [Dependencies](https://github.com/CarmenSheppard/PneumoKITy#dependencies)
* [Running PneumoKITy](https://github.com/CarmenSheppard/PneumoKITy#running-PneumoKITy)
* [User customisable options](https://github.com/CarmenSheppard/PneumoKITy#User-customisable-options)
* [Example command lines](https://github.com/CarmenSheppard/PneumoKITy#Example-command-lines)
* [How PneumoKITy works](https://github.com/CarmenSheppard/PneumoKITy#How-PneumoKITy-works)
* [Quality checks](https://github.com/CarmenSheppard/PneumoKITy#Quality-checks)
* [PneumoKITy Output](https://github.com/CarmenSheppard/PneumoKITy#Output-files)
* [Interpretation of results](https://github.com/CarmenSheppard/PneumoKITy#Interpretation-of-results)
* [Licence agreement](https://github.com/CarmenSheppard/PneumoKITy#licence-agreement)

## Introduction

PneumoKITy (**Pneumo**coccal **K**mer **I**ntegrated **Ty**ping) is a 
lite serotyping tool, developed for speed and sensitivity for mixed serotype calling and using experience of development and use of PneumoCaT to inform it's development. It is written for **Python 3.7+**. PneumoKITy uses the excellent tool [MASH](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x) for 
kmer based analysis. We decided to create a lite version (PneumoKITy) as this his tool has significant value 
for both sensitive detection of mixed serotypes in fastQ data and is also significantly faster for serotype screening than existing methods which may mean
it has utility  either as a quick assessment method which has input type flexibility for running PURE cultures, and when using fastq data it can be run on
expected mixed serotype samples (eg. Carriage study culture -enriched samples or even deep sequenced primary swabs) in order to attempt to assess serotype composition in 
these more complex situations.


PneumoKITy, like the original PneumoCaT tool assigns capsular types to
*S.pneumoniae* genomic data  using a using a two step approach, But PneumoKITy, is limited only to second stage determinations that can be assessed using presence or absence or gene allele variants.
However due to the extra specificity of the MASH method, PneumoKITy fully  serotype about 58% of the serotypes defined 
by the SSI Diagnostica typing sera, and can get to serogroup, or subgroup level for the remaining types. In addition 
PneumoKITy and provides some useful information regarding  some described genetic subtypes and genetic types. 
Serotypes that require determination using alternative variants such as SNPs cannot be distinquished using PneumoKITy.
PneumoKITy has the advantage that it can be used on assembly files OR illumina fastq read files and it is incredibly fast.

PneumoKITy has a useful RAG (red - amber- green) status flag which alerts users to the status of the results, with green 
flagged samples passing within expected metrics, amber flagged metrics indicating that a result was obtained but that the results
should be checked as some outcomes were not within range, or RED indicating a failure to serotype. 


**PneumoKITy has two run options:**

`pure`  suitable for analysis of expected pure cultures. For this run method both fastq and assembly inputs can be used for analysis. This method which will increase specificity and 
reduce over reporting of mixed serotypes which may be present due to high levels of DNA from pure culture isolates. Although
it will report genogroups of any mixed samples, PneumoKITy will assign these as an AMBER RAG status due to them being unexpected in a 
pure culture. 


`mix` PneumoKITy is particularly sensitive to mixed serotypes when run with fastq data - which is useful for determination of types in multiple carriage.
When run in mix mode PneumoKITy reports stage 2 variants where possible within mixtures of multiple 
serotype carriage and reports the results slightly differently.  Mixtures are not automatically assigned an AMBER rag status. MIX analysismethod can ONLY be used with fastq inputs. 


   
PneumoKITy uses very different methods to previous PneumoCaT versions and requires a new running environment.

 The **C**apsular **T**ype **V**ariant database (CTVdb) used in PneumoKITy is now a real database, running in SQLite3,
 allowing  for much easier updating of information and better scope for the database to grow in the future as more variants and serotype 
 determinants are added. The use of this format will also allow us to store  extra related information about the serotypes, such as subtype information 
 (subtyping etc is planned in a future update). 

For PneumoKITy specifically any snp or gene_function variants have been removed from the CTVdb as it is not possible to perform
these determinations using kmer screening alone, so the database is very small. 


The database is populated from the included Excel template in database_tools. 

At present the database import script can only import new data and has not yet been programmed with any functions to update existing records (future update)

## Dependencies and getting up and running 

PneumoKITy is written for Python 3.7+ and is **NOT** compatible with earlier versions of Python. In particular the
method of running the mash subprocess requires python 3.7+ (not 3.6)


PneumoKITy requires the following packages installed before running:
* Mash version 2.3 (or 2.0+) [https://github.com/marbl/Mash](https://github.com/marbl/Mash)
* numpy [http://www.scipy.org/scipylib/download.html](http://www.scipy.org/scipylib/download.html)
* pandas [https://pandas.pydata.org/](https://pandas.pydata.org/)
* SQLite3 [https://www.sqlite.org/index.html](https://www.sqlite.org/index)
* SQLalchemy [https://www.sqlalchemy.org/](https://www.sqlalchemy.org/)

Due to the dependencies PneumoKITy can only be run on Linux based operating systems,
however the software can be run on Windows 10 using the Windows Subsystem for Linux 
[WSL](https://docs.microsoft.com/en-us/windows/wsl/). 

An easy way to install the dependencies is to use a Python 3 conda or venv environment.  

Install numpy, pandas and SQLalchemy in the environment (SQLite3 is likely to be bundled anyway).

Mash 2.3 can be installed from conda using the bioconda channel - please ensure that the channel priorities in your
conda environment are set up in the order as follows (later added channels have greater priority).
If installed from default conda channels the version of mash is NOT compatible with
PneumoKITy. 

`conda config --add channels defaults `   
`conda config --add channels bioconda`  
`conda config --add channels conda-forge`

Then run `conda install mash` to install mash. Check the mash version with `mash --version`
PneumoKITy has been tested with versions 2.0 to 2.3.

 

Once this is working you should be able to run PneumoKITy as detailed below.

## Running PneumoKITy
#### Mandatory commandline inputs

The first stage is to select run type as detailed above, either `pure` for runs expecting pure culture input (eg. reference
samples from colony picks).  or `mix` for expected mixed samples eg. from carriage studies.

For PURE culture analysis PneumoKITy accepts 3 input options, two for read input and one for assembly
 input. It mandatory to give at least one of these options. 
 
**Option 1 -i**: a folder containing two fastq files (argument: `-i 
folder_path`)
 
**Option 2 -f**: Specified paths to forward and reverse read files with input 
paths for two fastq files (argument: `-f read1_path read2_path`)

**Option 3 -a**: A specified assembly file (argument: `-a path_to_assembly`)

For MIX analysis PneumoKITy will only accept fastq files and therefore only Option 1 and Option 2 are available.


### Customisable settings (with default behaviour)
 
**-o** (output directory): An output directory can be specified using `-o 
output_dir`. PneumoKITy will place a subfolder `pneumo_capsular_typing` 
within any specified output directory.  If no output dir is specified 
PneumoKITy defaults to using the input directory for either read1 if using 
fastq input or the assembly in this case ensure the input directory is writable or an
error will occur.
 
 
**-t** (threads) Number of threads to use for subprocesses (default = 1) eg: `-t 8`

**-m** (mash): path to mash software (If necessary). By default PneumoKITy will use the 
command  `mash` which will only work if the mash software is recognised as included in 
 `PATH` variable (eg if installed from conda and runnign in conda environment).
Otherwise the path to the mash software file location
  **must** be provided eg: `-m path_to_mash/mash`.

**-p** (minpercent): Alternative filter cut-off value for kmer percentage, ie 
percentage of kmer hits to reference. eg:  `-p 80` (default = 90) NB: The software will automatically 
incrementally drop the kmer percentage if no hits are initially obtained - and report serotype with an amber 
RAG status if determined with the dropped back cut off (auto drop min= 70%). 

**-n** (minmulti): minimum, median-multiplicity value cut off (relevant for fastq 
read input only). eg: `-n 10` (default = 4) Used to minimise the kmers present due to 
sequencing errors only.  We have set the default to 4 to allow sensitivity for mixed serotype detection, while
retaining some specficity to ignore novel reads caused by sequencing errors. However if you are running expected pure
cultures we recommend setting this value to 10 to avoid over sensitivity of detection of reads caused by contamination for eg.
during pure culture DNA extraction. Recommended minimum 4. 

**-s** (sampleid): Specify a sample ID for output files. If not specified 
PneumoKITy will default to the assembly file name or the fastq file name 
(split on character as specified in -S or on "." default). eg: `-s sample-name`

**-S** (split): Specify a split character to split filename on for sample ID for output files. If not specified 
PneumoKITy will default to using ".". Eg if using `-S _`, filename `sampleid_1.fastq.gz` becomes sampleID `sampleid`
(split on first "_").

**-c** (collate): Specify a folder for PneumoKITy to collate results from the run into a file called "Collated_result_data.csv" (folder *MUST* already exist).
This is useful when running multiple PneumoKITy jobs for a particular project, for example via a queue submission system or Bash loop command. The basic result data will be appended to this file until either the flag is not specified, a different folder is specified or the resulting file is moved or renamed. In rare instances multiple processing MAY result in this file not being writable, and a result beng missed from the collation. The original data files from the run will be saved in their output location.

**-d** (database): path to capsular type variant database (ctvdb). By default this 
will be the ctvdb that is contained within the PneumoKITy program. 


An alternative ctvdb folder can be specified (ADVANCED USERS ONLY). If you wish to use an 
 alternative ctvdb we recommend that this is stored in a separate folder location
 and called using the -d option. To avoid confusion when troubleshooting,
  please do **NOT**  overwrite the included ctvdb and run without specifying the option -d.
 eg: `-d path_to_alternativectvdb

## Example command lines

1.  Pure Culture analysis, input folder containing read 1 and read 2, 8 threads, path to mash, custom sample name.

`python pneumokity.py pure -i my_input_folder -t 8 -m path_to_mash/mash -s my_name`

2. Mix analysis, input read files, custom output_dir, collate file folder location

`python pneumokity.py mix -f path_to_fastq/fastq1 path_to_fastq/fastq2 -o my_output_dir
 -c path_to_collate_folder `

3. Mix analysis ,i nput fastq file paths, custom output directory,  custom kmer percentage cut off

`python pneumokity.py mix -f path_to_fastq/fastq1 path_to_fastq/fastq2 -o my_output_dir
 -k 75`

4. Pure analysis input assembly, custom initial kmer percentage cut off, specified sampleid for filename
`python pneumokity.py pure -a path_to_assembly/assembly  -m path_to_mash/mash -k 75 -s mysampleid`

5. Pure analysis input assembly, collate file folder location

`python pneumokity.py pure -a path_to_assembly/assembly  -c path_to_collate_folder`

6. mix analysis input fastq folder, path to mash ,specified split character for filename

`python pneumokity.py  mix -i my_input_folder -m path_to_mash/mash -S _`


## How PneumoKITy works  

PneumoKITy  uses a kmer based approach to 
screen the read files or input assembly file against the capsular operon 
references (references.fasta) or gene presence absence references is stored as a mash sketch file (references.msh).

In the first stage, the input sequences are screened against a file of the
 capsular operon references (references.msh). The resulting data is analysed for hit
  percentage (kmers in sample / kmers in reference x 100) and then filtered according to
  the cut off values specified.

-p kmer percentage - Percentage of hit kmers against total reference kmers for a given
serotype reference sequence, calculated by PneumoKITy.

-n median multiplicity -the median number of multiples of a given kmer in the dataset,
given in mash screen output. Low kmer multiplicity could be caused by 
sequencing errors or mixed samples - only applicable for input read files if
 assembly files are input this value is automatically set to 1.  By default for 
 fastq input data  this value is set to 10 for expected PURE culture runs and 4 for expected MIXED culture runs.


This results in several stage 1 outcome categories and also come with a RAG 
status (see Quality and Error checking):

- `type` - A serotype that can be determined in stage 1 only. Resulted as final result 
when run in normal mode 
- `subtype` - a serotype that determined in stage 1 only but also has subtypes which 
can be determined to add further information (*future update*)
- `variants`  - a serotype that cannot determined in stage 1 only - needs further 
determination against capsular type variant database.
PneumoKITy will attempt to subtype where variants exist in the limited database.
- `mix` - Mixed serotypes are present
- `mix_variants`  - mixed serotype analysis selected, and mix serotypes found which include some that cannot determined in stage 1 only - needs further 
determination against capsular type variant database of subtypes determined by PneumoKITy. 
- PneumoKITy will attempt to subtype where variants exist in the limited database.
- `acapsular` - the kmer percentages against all capsular operon references were 20% and 
suggests that the sample does NOT have a capsular operon present and may be an acapsular 
organism. Check phenotype and species ID.
- `No hits` - No hits were determined but the kmer percentage to the reference operons
 was higher than that that might suggest an acapsular organism. Could be due to serotype
variant or due to poor sequencing quality. 

If the `variants` category is returned in stage 1 and the serotype has presence absence or allele variants
 then PneumoKITy proceeds to stage 2.

**Stage 2 - Capsular Type Variants** -  only two categories of variants are available in stage 2, these are
gene presence/absence and allele variants. Both of these methods are implemented using MASH in
a similar procedure to described above. Only genogroups 15F_15A and 19A_19AF(variant) are able to fully
determine in stage 2. However some genogroups contain types that can be separated using the kmer screening determinations
eg 12F from 12A_12B_44_46 and 38 from 25F_25A. So depending on the serotype with in a genogroup you
may get either a type specific or genogroup specific result.

Additional reference files have been added to the CTVdb for serotypes 24F and 33F and referecen sequence for 18F has been replaced 
due to the extra specificity of the kmer screening method causing many test samples to give low top hits and AMBER calls, due to divergence between the old reference strains
used for the original PneumoCaT reference.fasta and isolates now circulating (in the UK). Once more data has been collected, potentially 
it may be possible to call other serotypes from groups to serotype level due to the greater specificity of the kmer method and increased
match with the new references (eg 24F from B), however additional testing is needed to validate this.

All other serotypes requiring variant analysis will return a genogroup as the methods for determining SNPS and gene function/non-function
variants are not available. A partial gene profile will be present for those where some of the determinations are allele or presence/absence based.


## Quality checks

PneumoKITy is written for use in an accredited laboratory and to aid this,
stores various metrics which are reported in the final 
results report file created at the end of the run.  This includes the 
versions of software used (PneumoKITy itself, Mash etc), paths to input 
files,  information about the reference database used and the cut offs 
specified. 

Outputs from PneumoKITy are automatically assigned a RAG status, please take care to note the status assigned when interpreting the results of the analysis.

![#4CAF50](https://via.placeholder.com/15/4CAF50/000000?text=+) GREEN: Analysis passed within expected cut-off

![#F39C12](https://via.placeholder.com/15/F39C12/000000?text=+) AMBER: Result obtained but caution advised, check top hit percentages possible variant or low sequence quality.

![#f03c15](https://via.placeholder.com/15/f03c15/000000?text=+) RED: Analysis failed.

Amber result status is assigned when no-hits occur in stage 1 analysis, the 
kmer percentage cut off is automatically dropped by 10% from the initial cut 
off (so usually reset to 81% unless the user has input a custom kmer percent
 cut off) and analysis is run again. This was implemented to help avoid 
 those annoying situations when a sample would miss just below the cut off as sometimes happened with PneumoCaT1 while also warning the user that something may be wrong with the result. 
If an AMBER result is obtained it could be due to either poor sequence quality, or a variant of the sequence which does not match very well with the reference sequences available in the CTVdb - please check the results. 

Amber result status also occurs in stage 2 for unrecognised variant profiles, or if a mixed result is obtained in an 
expected pure culture run. 

   
RED rag status alerts the user to failure of the serotyping. This could be 
due to an unexpected pattern of results, mixed serotypes or no-hits in the 
analysis.

**IMPORTANT:** Once software is released publicly - each time an update is added to the included ctvdb we increment the 
version of the overall PneumoKITy software (eg. from version 2.0 to 2.0.1). 
However if an alternative version of the ctvdb is used it is up to the user to record 
which version is used for their analysis.

## Output Files

PneumoKIty produces several output files. 

`SAMPLEID_serotyping_results.txt` *(All serotypes)*

A text file containing human-readable formatted information about the run metrics and sample results

`SAMPLEID_quality_system_data.csv` *(All serotypes)*

A csv file containing information about the run metrics (folder and file locations, software versions, cut offs etc)

`SAMPLEID_result_data.csv` *(All serotypes)*

csv file containing final result data from run

`SAMPLEID_stage1_screen.csv` *(All serotypes)*

csv file containing all of the stage1 screen MASH run hits.


`SAMPLEID_GENE_screen.csv` *(Serotypes resulting at stage 2 only)*

csv files containing MASH screen run hits information for the relevant variant genes.

`SAMPLEID_mixed_serotypes.csv` *(Mix run and when Mixed serotypes are detected only)*

csv file containing details of mixed serotypes found in the sample.


## Interpretation of results

**`predicted serotype`** -  This is the predicted PHENOTYPICAL type of the organism if characterised using the commercially available [SSI Diagnostica serotyping sera](https://www.ssidiagnostica.com/antisera/pneumococcus-antisera/). However, the organism could have a specific underlying genetic type as in the case of 23B1 and 19A/F for example. 

Some previously described important genetic subtypes **are** represented in the ctvdb and can be determined by looking at the stage 1 hits. Eg. in the case of the variant 19A/F isolates that have the [genetic background of a 19A but produce a 19F capsular polysaccharide](https://pubmed.ncbi.nlm.nih.gov/19439547/) , the predicted serotype result would be 19F, but the stage 1 result is recorded as 19AF and the result is determined from 19A in stage 2 via wzy analysis, whereas a standard 19F isolate would be determine in stage 1 analysis alone by hit to the 19F reference. For 23B1 the top hit in stage 1 will be the 23B1 reference. At present there is no specific "genetic subtype" result field implemented (this is planned in a future version), so the user needs to look at the stage 1 hits to determine these.

New genetic types are being discovered all the time and it is not possible for us to keep up with them, however we hope that the outputs obtained from PneumoKITy will help the user to determine if their isolate may be a novel genetic serotype for which a reference is not available in the ctvdb, and can then be further investigated.

if the sample has failed to hit a serotype a description of the result is output in this field.

**`top hits`** - this is a list of the top 5 hits in the stage 1 analysis, with their MASH hit kmer percentage score in order highest to lowest. For serogroups 6 and 19, subtype references from [Elberse et al](https://journals.plos.org/plosone/article?id=10.1371/journal.pone.0025018_) are included in the stage 1 references (eg 19F-III, 19A-IV).  Subtypes are not yet interpreted in the program but the closest hit subtype can be determined by reference to the top hit.

**`max percent`** - the kmer percentage of the maximum hit in the MASH screen stage 1 analysis.

**`folder`** - the folder location within the \ctvdb folder used for the stage 2 analysis/

**`stage 1 result`** - the outcome of stage 1 analysis.

**`mix_mm`** - the estimated % abundance of each serotype in a mixture (calculated from the median multiplicity values >cut off)  
- Only applicable if input is fastq and mixed serotypes are found

**`stage 2 varids`** - the variant ID (keys) in the ctvdb sql database of the variants used for determination in stage2.

**`stage 2 hits`** - the variant genes and results of stage 2 analysis, eg for those variants determined using MASH screen (gene presence/absence and allele) the result will be a hit % of kmers from sample vs kmers in the gene. 

**`stage 2 result`** - outcome of any stage 2 analysis, eg hit variant determined (eg detected, not detected)

**`rag status`** - the overall quality status (traffic light system) of the run as described above.

### Outputs specific for mixed culture analysis

When mix detection is selected on run, PneumoKity will attempt to subtype any serotypes within the mix that are subtypeable using the limited collection in the PneumoKITy CTVdb. PneumoKITy will alos output a csv file contining details of the mixtures for each sample, and add a table into the report.txt file. 


### Licence Agreement 

This software is covered by GNU General Public License, version 3 (GPL-3.0).




