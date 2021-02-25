# Nextflow pipleine for preprocessing of ribosomal profiling sequencing

[[_TOC_]]


The pipeline was developed by 
[Independent Data Lab](https://www.independentdatalab.com/) 
in collaboration with 
[Immagina Biotechnology, The Ribosome Company ](https://www.immaginabiotech.com/).

The pipeline was adopted to work on, and tested on samples prepared and 
sequenced by Immagina.

In order to run on a local computer, following software should be installed:

*  nextflow
*  fastqc
*  cutadapt
*  bowtie2
*  samtools
*  picard
*  multiqc
*  R and [riboWaltz R package](https://github.com/LabTranslationalArchitectomics/riboWaltz)

If any utilized files are located in AWS S3 bucket, or to run the pipeline on AWS,
aws cli tool needs to be installed and configured with access key, secret key, and region

To do that, install AWS CLI, and run
```
aws configure
```

 and fill in the fields as following:

```
AWS Access Key ID :  <YOUR KEY>
AWS Secret Access Key : <YOUR SECRET KEY>
Default region name [eu-central-1]: eu-central-1
Default output format [json]: json
```

## Configure the run to specific samples

The key parameters like the sample fastq files location, command line arguments for 
all the tools, reference files locations are listed in `nextflow.config`, with 
symbolic example pf parameter values. The best way to run it is to copy 
`nextflow.config` into `my.config`


```
cp nextflow.config my.config
```

and modify values in `my.config`

There one can also change the cutadapt adapter sequences and other parameters.

## Running the pipeline

To run the pipeline locally, navigate to directory where the pipeline is located 
(where the `main.nf` is) and run:

```
nextflow run main.nf -c my.config -profile local
```

To run on AWS:

```
nextflow run main.nf -c my.config -profile awsbatch
```

## Troubleshooting

### Failed alignment or indexing

If you get an error during alignment or indexing, check that the file locations 
are set-up correctly, remove bowtie index directory and run the pipeline again without
`-resume` option. 

It is possible that during previous run indexing completed only partially, but the directory 
with the index was already created. This can happen when the pipeline was interrupted 
due to lost internet connection, or lack of available memory or killed by the user.
The pipeline then identifies the presence of index directory, and assumes that 
indexing can be skipped. Removing index directory will force the pipeline to index again. 

# License

 This software is distributed under the BSD-style License. See LICENSE for details.

