# README for PneumoCaT2 

## Table of content

* Introduction
* [Dependencies](https://github.com/CarmenSheppard/pneumocat2#dependencies)
* [Running PneumoCaT](https://github.com/CarmenSheppard/pneumocat2#running-pneumocat2)
* [User customisable options](https://github.com/CarmenSheppard/pneumocat2#User-customisable-options)
* [Example command lines](https://github.com/CarmenSheppard/pneumocat2#Example-command-lines)
* [How PneumoCaT2 works](https://github.com/CarmenSheppard/pneumocat2#How-PneumoCaT2-works)
* [Quality checks](https://github.com/CarmenSheppard/pneumocat2#Quality-checks)
* [PneumoCaT Output](https://github.com/CarmenSheppard/pneumocat2#Output-files)
* CTV database
* Examples
* Troubleshooting
* Contact Information
* Licence Agreement


PneumoCaT2 (**Pneumo**coccal **Ca**psular **T**yping version 2) is a 
complete from the ground up, redevelopment of the original [PneumoCaT](https://github.com/phe-bioinformatics/PneumoCaT) capsular typing tool, written for **Python 3.6+**, using different methods. Stage 1 uses the excellent tool [MASH](https://genomebiology.biomedcentral.com/articles/10.1186/s13059-016-0997-x) for 
kmer based analysis.  

PneumoCaT2, like the original PneumoCaT tool assigns capsular types to
*S.pneumoniae* genomic data  using a using a two step 
approach. PneumoCaT2 has the advantage that it can be used on assembly files
 as well as illumina fastq read files. 
  
This version incorporates developments that we hope are useful, drawing on 
experience of running PneumoCaT routinely for over 2 years on thousands of 
isolates, and from feedback from external users. Our priority however is (as always)
our own clinical reference service. PneumoCaT has always been written primarily 
for our own use in the PHE pipeline but we have now tried to incorporate changes to make running and interpreting 
 outputs more user-friendly for external users. 
 
PneumoCaT2 aims to improve the accuracy, sensitivity and efficiency of 
serotype calling and generally improve PneumoCaT. This version uses very different
 methods to previous versions and requires a new running environment.

 The **C**apsular **T**ype **V**ariant database (CTVdb) used in 
PneumoCaT 2 is now a real database, running in SQLite3. We felt that, though
 the database is currently small, the use of a true database structure
 allows for much easier updating of information and gives better 
 scope for the database to grow in the future as more variants and serotype 
 determinants are added. The use of this format will also allow us to store 
 extra related information about the serotypes, such as subtype information 
 (subtyping is planned in a future update). 



## Dependencies


PneumoCaT2 is written for Python 3.6+ and is **NOT** compatible with earlier versions of Python.

PneumoCaT2 requires the following packages installed before running:
* Mash version 2.2 (or 2.1) [https://github.com/marbl/Mash](https://github.com/marbl/Mash)
* numpy [http://www.scipy.org/scipylib/download.html](http://www.scipy.org/scipylib/download.html)
* pandas [https://pandas.pydata.org/](https://pandas.pydata.org/)
* SQLite3 [https://www.sqlite.org/index.html](https://www.sqlite.org/index)
* SQLalchemy [https://www.sqlalchemy.org/](https://www.sqlalchemy.org/)

Due to the dependencies PneumoCaT2 can only be run on Linux based operating systems. Please note if using conda enviroments the version of mash installed from Conda (1.X) is NOT Compatible with PneumoCaT2. Please use 2.2 (or 2.1 - 2.2 recommended)

## Running PneumoCaT2
#### Mandatory commandline inputs

PneumoCaT 2 accepts 3 input options, two for read input and one for assembly
 input. It mandatory to give at least one of these options. 
 
**Option 1 -i**: a folder containing two fastq files (argument: `-i 
folder_path`)
 
**Option 2 -f**: Specified paths to forward and reverse read files with input 
paths for two fastq files (argument: `-f read1_path read2_path`)

**Option 3 -a**: A specified assembly file (argument: `-a path_to_assembly`)

 Please note: for serotypes that result in stage 1 only, when using assembly input files,
 the analysis is EXTREMELY fast, do not assume it has failed if the job completes
 seemingly instantly ;-)

### Customisable settings (with default behaviour)
 
**-o** (output directory): An output directory can be specified using `-o 
output_dir`. PneumoCaT2 will place a subfolder `pneumo_capsular_typing` 
within any specified output directory.  If no output dir is specified 
PneumoCaT2 defaults to using the input directory for either read1 if using 
fastq input or the assembly in this case ensure the input directory is writable or an
error will occur.
 
**-d** (database): path to capsular type variant database (ctvdb). By default this 
will be the ctvdb that is contained within the PneumoCaT2 program. 


An alternative ctvdb folder can be specified (ADVANCED USERS ONLY). If you wish to use an 
 alternative ctvdb we recommend that this is stored in a separate folder location
 and called using the -d option. To avoid confusion when troubleshooting,
  please do **NOT**  overwrite the included ctvdb and run without specifying the option -d.
 eg: `-d path_to_alternative_CTV.db` 
  
 The alternative ctvdb folder location **MUST** contain a `references.msh` file 
 created using [mash](https://mash.readthedocs.io/en/latest/) sketch 
 **with sketch parameters -k 31 and (default) -s 1000** run on all the references that you wish to include,
  and a `CTV.db` which is an sqlite database file containing information about all serotypes for inclusion 
  in the serotype determination. The database structure (table and field names) should
   match that defined in the sqlalchemydeclarative.py script. Matching foldernames 
   for each serogroup contained with the CTV.db must be added within the
    ctvdb root folder.
     Further information on creation of custom db will be added in the future.
 
**-t** (threads) Number of threads to use for subprocesses (default = 1) eg: `-t 8`

**-m** (mash): path to mash software. By default PneumoCaT2 will use the 
command  `mash` which will only work if the mash software is included in the
 `PATH` variable. Otherwise the path to the mash software file location
  **must** be provided eg: `-m path_to_mash\mash`.

**-p** (minpercent): Alternative filter cut-off value for kmer percentage, ie 
percentage of kmer hits to reference. eg:  `-p 80` (default = 90) NB: The software will automatically 
incrementally drop the kmer percentage if no hits are initially obtained - and report serotype with an amber 
RAG status if determined with the dropped back cut off (auto drop min= 70%). 

**-n** (minmulti): minimum, median-multiplicity value cut off (relevant for fastq 
read input only). eg: `-n 2` (default = 4) Used to minimise the kmers present due to 
sequencing errors only. 

**-s** (sampleid): Specify a sample ID for output files. If not specified 
PneumoCaT2 will default to the assembly file name or the fastq file name 
(split on first "."). eg: `-s sample-name`

**-c** (collate): Specify a folder for PneumoCaT2 to collate results from the run into a file called "Collated_result_data.csv"
This is useful when running multiple PneumoCaT2 jobs for a particular project, for example via a queue submission system or Bash loop command. The basic result data will be appended to this file until either the flag is not specified, a different folder is specified or the resulting file is moved or renamed.

## Example command lines

1. Input folder containing read 1 and read 2, 8 threads, path to mash, custom sample name.

`python pneumocat2.py -i my_input_folder -t 8 -m path_to_mash/mash -s my_name`

2. Input read files, path to mash, custom output_dir, collate file folder location
`python pneumocat2.py -f path_to_fastq/fastq1 path_to_fastq/fastq2 -o my_output_dir
 -c path_to_collate_folder `

3. Input fastq file paths, custom output directory, path to mash, custom kmer percentage cut off

`python pneumocat2.py -f path_to_fastq/fastq1 path_to_fastq/fastq2 -o my_output_dir
 -m path_to_mash/mash -k 75`

4. Input assembly, path to mash, custom initial kmer percentage cut off

`python pneumocat2.py -a path_to_assembly/assembly  -m path_to_mash/mash -k 75`

5. Input assembly, path to mash ,collate file folder location

`python pneumocat2.py -a path_to_assembly/assembly  -m path_to_mash/mash -c path_to_collate_folder`


## How PneumoCaT2 works  

Stage 1 of PneumoCaT2 (and some calls in stage 2) use a kmer based approach to 
screen the read files or input assembly file against the capsular operon 
references (references.fasta) which is stored as a mash sketch file (references.msh).

In the first stage, the input sequences are screened against a file of the
 capsular operon references (references.msh). The resulting data is analysed for hit
  percentage (kmers in sample / kmers in reference x 100) and then filtered according to
  the cut off values specified.

-p kmer percentage - Percentage of hit kmers against total reference kmers for a given
serotype reference sequence, calculated by PneumoCaT2.

-n median multiplicity -the median number of multiples of a given kmer in the dataset,
given in mash screen output. Low kmer multiplicity could be caused by 
sequencing errors or mixed samples - only applicable for input read files if
 assembly files are input this value is automatically set to 1.


This results in several stage 1 outcome categories and also come with a RAG 
status (see Quality and Error checking):

- `type` - A serotype that can be determined in stage 1 only. Resulted as final result 
when run in normal mode 
- `subtype` - a serotype that determined in stage 1 only but also has subtypes which 
can be determined to add further information (*future update*)
- `variants`  - a serotype that cannot determined in stage 1 only - needs further 
determination against capsular type variant database.
- `mix` - Mixed serotypes are present
- `acapsular` - the kmer percentages against all capsular operon references were 20% and 
suggests that the sample does NOT have a capsular operon present and may be an acapsular 
organism. Check phenotype and species ID.
- `No hits` - No hits were determined but the kmer percentage to the reference operons
 was higher than that that might suggest an acapsular organism. Could be due to serotype
variant or due to poor sequencing quality. 

If the `variants` category is returned in stage 1 then PneumoCaT2 proceeds to stage 2.

**Stage 2** - currently only two categories of variants are available in stage 2, these are
gene presence/absence and allele variants. Both of these methods are implemented using MASH in
a similar procedure to described above. Only serotype 15A and 19A / variant 19F(A) are able to
determine in stage 2. All other serotypes requiring variant analysiswill return an error message
regarding unrecognised variant profile as the methods for determining SNPS and gene function/non-function
variants are not yet implemented.


## Quality checks

PneumoCaT2 is written for use in an accredited laboratory and to aid this,
stores various metrics which are reported in the final 
results report file created at the end of the run.  This includes the 
versions of software used (PneumoCaT itself, Mash etc), paths to input 
files,  information about the reference database used and the cut offs 
specified. 

Outputs from PneumoCaT2 are automatically assigned a RAG status:

* GREEN: Analysis passed within expected cut-off
* AMBER: Result obtained but caution advised, check top hit percentages.
* RED: Analysis failed

Amber result status is assigned when no-hits occur in stage 1 analysis, the 
kmer percentage cut of is automatically dropped by 10% from the initial cut 
off (so usually reset to 81% unless the user has input a custom kmer percent
 cut off) and analysis is run again. This was implemented to help avoid 
 those annoying situations when a sample would miss just below the cut off.

Amber result status also occurs in stage 2 for unrecognised variant profiles.

The analysis is given the RAG status AMBER to flag to the user that 
  they should check the results. This type of result may occur in situations
   when the local pneumococcal population has differing clones than that to 
   which the CTVdb was created, and may help avoid many fails in this 
   situation while still alerting users to the issue. Local procedures can then be 
   implemented to interpret these results.
   
RED rag status alerts the user to failure of the serotyping. This could be 
due to an unexpected pattern of results, mixed serotypes or no-hits in the 
analysis.

**IMPORTANT:** Each time an update is added to the included ctvdb we increment the 
version of the overall PneumoCaT software (eg. from version 2.0 to 2.0.1). 
However if an alternative version of the ctvdb is used it is up to the user to record 
which version is used for their analysis.

## Output Files

PneumoCaT2 produces several output files. 

`SAMPLEID_serotyping_results.txt` *(All serotypes)*

A text file contaning human-readable formatted information about the run metrics and sample results

`SAMPLEID_quality_system_data.csv` *(All serotypes)*

A csv file containing information about the run metrics (folder and file locations, software versions, cut offs etc)

`SAMPLEID_result_data.csv` *(All serotypes)*

csv file containing final result data from run

`SAMPLEID_stage1_screen.csv` *(All serotypes)*

csv file containing all of the stage1 screen MASH run hits.


`SAMPLEID_GENE_screen.csv` *(Serotypes resulting at stage 2 only)*

csv files containing MASH screen run hits information for the relevant variant genes.
