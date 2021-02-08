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
all the tools, reference files locations are specified in 

```
nextflow.config
```

There one can also change the cutadapt adapter sequences and other parameters.

## Running the pipeline

To run the pipeline locally, navigate to directory where the pipeline is located 
(where the `main.nf` is) and run:

```
nextflow run main.nf -profile local
```

To run on AWS:

```
nextflow run main.nf -profile awsbatch
```

# License

 This software is distributed under the BSD-style License. See LICENSE for details.

