# Creating a new pipeline for long-read RNA-seq using StringTie v2


This repository contains data and Nextflow workflow for the tests of the Google Summer of Code .


#### Creating a new pipeline for long-read RNA-seq using StringTie v2
The majority of contemporary genomics and transcriptomics research is carried out using short-read technology such as the output of Illumina sequencers. However, newer long-read technologies such as PacBio and Oxford Nanopore (ONT) are becoming more prevalent due to the advantages they offer over short reads. The main advantage of long reads is that they span a much larger portion of the genome or transcriptome, making it easier to detect events such as structural variants or isoforms. As more researchers begin to take advantage of long reads, GenPipes needs to evolve and support long-read technology in its main pipelines.

The objective of this project is to create a new RNA-seq pipeline (rnaseq_longreads), that supports long-read inputs. The new pipeline would be created based on the current versions of the RNA-seq pipeline, with the addition of the minimap2 aligner as well as StringTie v2 and Ballgown for generating transcript counts.

If time allows, in addition to implementing a long-read RNA-seq pipeline within the GenPipes framework, the student will help update the RNA-seq pipeline to also use StringTie2 as well as test both pipelines to ensure they are working properly.

####Selection tests
Install the following software in your computing environment:

a. StringTie v2 (run test suite to ensure it is working).
b. The Ballgown R-Bioconductor library (If required, install R and any additional dependencies).
c. Tablemaker.

Download the following nanopore mouse RNA-seq datasets from this recent study:

a. Brain C1: BAM [file](http://www.genoscope.cns.fr/externe/ONT_mouse_RNA/data/transcriptome/RNA_nanopore.brain.C1R1_mapping_E94_minimap2_primary_no_read_less_than_80QC.bam)
b. Liver C1: BAM [file](http://www.genoscope.cns.fr/externe/ONT_mouse_RNA/data/transcriptome/RNA_nanopore.liver.C1R1_mapping_E94_minimap2_primary_no_read_less_than_80QC.bam)

Using StringTie v2 long-read mode, generate a GTF file for each of these samples.

Run Tablemaker on the two samples above (use the BAM and GTF files used above).

Import the output of Tablemaker into R using Ballgown.

Save your work in an external repository and send us a link to it. The repository must contain:

- A bash script with your StringTie v2 and Tablemaker commands.
- An R script with your Ballgown commands.
- The first 1000 lines of both GTF files generated using StringTie2.
- The first 1000 lines of the FPKM matrix produced by Ballgown, saved as a csv.
- Hint: use the texpr function in Ballgown to obtain this matrix.



### Pipeline
The pipeline is built using [Nextflow](https://www.nextflow.io), a workflow tool to run tasks across 
multiple compute infrastructures in a very portable manner. It comes with a docker container making installation trivial and results highly reproducible.
This pipeline has 4 different processes:
- stringtie
- tablemaker
- ballgown
- postProduceOutputs

### Pipeline Quick Start
Make sure you have either docker/singularity installed or the required dependencies listed 
in the last section.

Install the Nextflow runtime by running the following command:

    $ curl -fsSL get.nextflow.io | bash


When done, you can launch the pipeline execution by entering the command shown below:

    $ nextflow run edgano/C3G_longRead2020
    

By default the pipeline is executed against the provided example dataset. 
Check the *Pipeline parameters*  section below to see how enter your data on the program 
command line.     
  

### Containers

Stringtie v2 and Tablemaker above are available in a [Docker](http://www.docker.com) image on DockerHub [here](https://hub.docker.com/r/edgano/c3g) and the Ballgown is in a native container from Dockerhub [here](https://hub.docker.com/r/kapeel/ballgown-r-package) and the images are tested to be compatible with the [Singularity](http://singularity.lbl.gov/).

The container also contains test data consisting of bam files in the directory `/data`.

To launch the container interactively with Docker run:

`docker run edgano/c3g`

To launch the container interactivly with Singularity run:

`singularity shell docker://edgano/c3g`


### Pipeline parameters

#### `--bam` 
   
* Specifies the location of the input *bam* file(s).
* Multiple files can be specified using the usual wildcards (*, ?), in this case make sure to surround the parameter string
  value by single quote characters (see the example below)

Example: 

    $ nextflow run edgano/C3G_longRead2020 --bam '/home/files/*.bam'

This will handle each bam file as a seperate sample.

#### `--outdir`

* Location of the results.
* Default locations is `results` directory.

## Command lines
The command lines run inside the pipeline by each process are:

 *stringTie* :
```
stringtie -L -o ${pair_id}.out.gtf ${bam_file}
```

*tablemaker* :
```
tablemaker -W -G ${gtf_file} -o ${pair_id}.tablemaker ${bam_file}
```

*ballgown*:
```
bg <- ballgown(dataDir = "!{baseDir}/results/table", samplePattern="!{pair_id}.tablemaker", meas='all')
transcript_fpkm = texpr(bg, 'FPKM')
write.csv(transcript_fpkm,file="!{pair_id}_fpkm.csv")
```

And the *postProduceOutputs* its just to have the output files in the format required by CCG.
```
    head -1000 ${gtf_file} > ${pair_id}_trim.gtf
    head -1000 ${fkm_file} > ${pair_id}_trim.csv
```




