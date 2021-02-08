#!/usr/bin/env nextflow
/*
################################################################################ 
##   Program: Nextflow riboSeq preprocessing pipeline  
##                                                                              
##   Copyright (c) Independent Data Lab UG, All rights reserved.                
##   The use of this source code is governed by a BSD-style license which can be
##   found in the LICENSE file in the root directory of this repository.        
##
##   This software is distributed without any warranty; without even the implied
##   warranty of merchantability or fitness for a particular purpose.  See the
##   above copyright notice for more information.
################################################################################
*/

// Channel for reading raw fastq files
Channel
    .fromFilePairs(params.raw_fastq, size: 1)
    .ifEmpty { exit 1, "Cannot find any fastq files matching: ${params.raw_fastq}" }
    .into { raw_fastq_to_fastqc; raw_fastq_to_cutadapt }

// for frequently used parameters, define short name
mode = params.publish_dir_mode
reference = defineReference()
run_bowtie2_indexing = !reference.bowtie2_ref_dir.exists()

align_wait_4_index_1 = ( run_bowtie2_indexing 
    ? Channel.empty() 
    : Channel.from(1) )

// Printing logs
log.info idlHeader()
def summary = [:]
summary['Reference'] = params.reference_name
summary['Reads'] = params.raw_fastq
summary['Output directory'] = params.outdir
summary['run_bowtie2_indexing'] = run_bowtie2_indexing
log.info summary.collect { k,v -> "${k.padRight(18)}: $v" }.join("\n")

// Starting the pipeline


process fastqc{
  
  tag "$sampleID"
  publishDir "$params.outdir", mode:"$mode"  
  
  input:
    set sampleID, file(raw_fastq) from raw_fastq_to_fastqc

  output:
    file "qc/fastq/*_fastqc.{zip, html}" into fastqc_to_multiqc  

  script:
  """
    mkdir -p qc/fastq
    fastqc ${params.fastqc_cmd_args} $raw_fastq -o qc/fastq
  """
}


process bowtie2_index_ref{

  publishDir "$params.references.bowtie2_ref_dir", mode:"$mode"  

  input:
      file rRNA_fasta from reference.rRNA_fasta
      file tRNA_fasta from reference.tRNA_fasta
      file cDNA_fasta from reference.cDNA_fasta

  output:
      file ("**")
      val 1 into align_wait_4_index_2
  
  when: run_bowtie2_indexing

  script:
  """
    bowtie2-build ${params.bowtie2_index_cmd_args} $rRNA_fasta rRNA
    bowtie2-build ${params.bowtie2_index_cmd_args} $tRNA_fasta tRNA
    bowtie2-build ${params.bowtie2_index_cmd_args} $cDNA_fasta cDNA
  """
}  

process bowtie2_index_ncRNA {
  publishDir "$params.references.bowtie2_ref_dir", mode:"$mode"  

  input:
    file ncRNA_fasta from reference.ncRNA_fasta 

  output:
    file ("**")
  
  when: run_bowtie2_indexing_1 && ncRNA_fasta.name != 'NO_FILE'

  script:
  """
    bowtie2-build ${params.bowtie2_index_cmd_args} $ncRNA_fasta ncRNA
  """
}


process cutadapt3 {
    
  tag "$sampleID"

  input:
    set sampleID, file(raw_fastq) from raw_fastq_to_cutadapt

  output:
    set (
      sampleID,
      file("${sampleID}_trimmed3.fastq.gz")
    ) into trimmed_3_to_5
    file "${sampleID}.trim3.log" into cutadapt3_to_multiqc  

  script:
  """
    cutadapt ${params.cutadapt3_cmd_args} \
      -o ${sampleID}_trimmed3.fastq.gz \
      $raw_fastq > ${sampleID}.trim3.log
  """
}

process cutadapt5 {
    
  tag "$sampleID"
  publishDir "$params.outdir/fastq/trimmed", mode:"$mode"

  input:
    set sampleID, file(trimmed3_fastq) from trimmed_3_to_5

  output:
    set (
      sampleID,
      file("${sampleID}_trimmed.fastq.gz")
    ) into (trimmed_fastq_to_bowtie2, trimmed_fastq_to_qc)
    file "${sampleID}.trim.log" into cutadapt5_to_multiqc  

  script:
  """
    cutadapt ${params.cutadapt5_cmd_args} \
      -o ${sampleID}_trimmed.fastq.gz \
      $trimmed3_fastq > ${sampleID}.trim.log
  """
}


process trimmed_fastqc{
  
  tag "$sampleID"
  publishDir "$params.outdir", mode:"$mode"  
  
  input:
    set sampleID, file(trimmed_fastq) from trimmed_fastq_to_qc

  output:
    file "qc/trimmed_fastq/*_fastqc.{zip, html}" into trimmed_fastqc_to_multiqc  

  script:
  """
    mkdir -p qc/trimmed_fastq
    fastqc ${params.fastqc_cmd_args} $trimmed_fastq -o qc/trimmed_fastq
  """
}

process bowtie2_align {
  tag "$sampleID"
  publishDir "$params.outdir/qc", mode:"$mode", pattern:"bowtie2_logs"  
  
  input:
    set sampleID, file(trimmed_fastq) from trimmed_fastq_to_bowtie2
    file bowtie2_ref from file(reference.bowtie2_ref_dir)
    val flag from align_wait_4_index_1.mix(align_wait_4_index_2)

  output:
    set (
      sampleID,
      file ("${sampleID}.sam")
    ) into sam_to_sort
    file("bowtie2_logs/*") into bowtie2_logs_to_multiqc

  script:
  """
    mkdir bowtie2_logs
    gunzip -c ${trimmed_fastq} > unzipped.fq
    bowtie2 ${params.bowtie2_contaminants_cmd_args} \
        -q unzipped.fq \
        --un no_rRNA.fq \
        -x ${bowtie2_ref}/rRNA > /dev/null 2> bowtie2_logs/${sampleID}_1_rRNA.log
    bowtie2 ${params.bowtie2_contaminants_cmd_args} \
        -q no_rRNA.fq \
        --un=no_rRNA_tRNA.fq \
        -x ${bowtie2_ref}/tRNA > /dev/null 2> bowtie2_logs/${sampleID}_2_tRNA.log
    if [[ $reference.ncRNA_fasta.name != "NO_FILE" ]]; then
        echo "Align to ncRNA"
        bowtie2 ${params.bowtie2_contaminants_cmd_args} \
            -q no_rRNA_tRNA.fq \
            --un=no_rRNA_tRNA_ncRNA.fq \
            -x ${bowtie2_ref}/ncRNA > /dev/null 2> bowtie2_logs/${sampleID}_3_ncRNA.log
        bowtie2 ${params.bowtie2_cdna_cmd_args} \
            -q no_rRNA_tRNA_ncRNA.fq \
            --un=${sampleID}_unaligned.fq \
            -x ${bowtie2_ref}/cDNA --no-unal  > ${sampleID}.sam 2> bowtie2_logs/${sampleID}_4_cDNA.log
    fi
    if [[ $reference.ncRNA_fasta.name == "NO_FILE"  ]]; then
        echo "No ncRNA specified, align to cDNA"
        bowtie2 ${params.bowtie2_cdna_cmd_args} \
            -q no_rRNA_tRNA.fq \
            --un=${sampleID}_unaligned.fq \
            -x ${bowtie2_ref}/cDNA --no-unal  > ${sampleID}.sam 2> bowtie2_logs/${sampleID}_3_cDNA.log
    fi
  """
}

process samtools_sam_sort{
  
  tag "$sampleID"
  
  input:
    set sampleID, file(sam) from sam_to_sort

  output:
    set (
      sampleID,
      file("${sampleID}.srtd.bam"),
    ) into bam_sorted_to_markdups

  script:
  """
    samtools view -F 16 -b -o ${sampleID}.bam $sam     
    samtools sort ${sampleID}.bam -o ${sampleID}.srtd.bam -m 2GiB -@ 8
  """
}

process picard_MarkDuplicates{
  
  tag "$sampleID"
  publishDir "$params.outdir", mode:"$mode" 
  
  input:
    set sampleID, file(bam) from bam_sorted_to_markdups

  output:
    set (
      sampleID,
      file("bam/${sampleID}/${sampleID}.srtd.markdup.bam"),
    ) into (final_bam_to_index, final_bam_qc, final_bam_to_ribowaltz)
    file("qc/markDuplicates/${sampleID}.markDuplicates") into picard_markdups_to_multiqc

  script:
  """
    mkdir -p qc/markDuplicates
    mkdir -p bam/${sampleID}/
    java -jar \${PICARD_JAR} MarkDuplicates ${picard_markdups_cmd_args} \
                -I ${bam} \
                -O bam/${sampleID}/${sampleID}.srtd.markdup.bam \
                -M qc/markDuplicates/${sampleID}.markDuplicates
  """
}

process samtools_bam_index{
  
  tag "$sampleID"
  publishDir "$params.outdir/bam/${sampleID}", mode:"$mode" 
  
  input:
    set sampleID, file(bam) from final_bam_to_index

  output:
    set (
      sampleID,
      file("${bam}.bai")
    ) 

  script:
  """
    samtools index ${bam} ${bam}.bai
  """
}


process multiqc {
  publishDir "$params.outdir/qc", mode:"$mode"
  
  input:
    file ('fastqc/*') from fastqc_to_multiqc.collect().ifEmpty([])
    file ('cutadapt3/*') from cutadapt3_to_multiqc.collect().ifEmpty([])
    file ('cutadapt5/*') from cutadapt5_to_multiqc.collect().ifEmpty([])
    file ('trimmed_fastqc/*') from trimmed_fastqc_to_multiqc.collect().ifEmpty([])
    file ('bowtie2/*') from bowtie2_logs_to_multiqc.collect().ifEmpty([])
    file ('markDuplicates/*') from picard_markdups_to_multiqc.collect().ifEmpty([])

  output:
    file ("**")

  script:  
  """   
    multiqc .                                                                                                   
  """
}

process ribowaltz {

  tag "$sampleID"
  publishDir "$params.outdir/riboWaltz", mode:"$mode"

  input:
    set (
      sampleID, 
      file(bam)
    ) from final_bam_to_ribowaltz 
    file cdna_fasta from file(reference.cDNA_fasta)
    file gtf from file(reference.gtf)
    file rmd from file(file("${baseDir}/rmd_reports/riboWaltz_report.rmd"))

  output:
    file("**")
  
  script:  
  """  
    cp -L ${rmd} ${sampleID}_riboWaltz_report.rmd
    R -e "rmarkdown::render('${sampleID}_riboWaltz_report.rmd', output_dir = '${sampleID}', params = list('sample_name' = '${sampleID}', 'bam_dir' = '.', 'gtf_file' = '${gtf}', 'cdna_fasta_file' = '${cdna_fasta}'))" 
  """

}



/*
________________________________________________________________________________

                            F U N C T I O N S
________________________________________________________________________________

*/

// some references contain a separate file for ncRNA (like Ensembl), while others
// don't (like gencode)
def checkParamReturnFileReferences(item) {
    if( params.references."${item}")
      return file(params.references."${item}")
    else
      return(file('NO_FILE'))
}

def defineReference() {
  return [
    'rRNA_fasta' : checkParamReturnFileReferences("rRNA_fasta"),
    'tRNA_fasta' : checkParamReturnFileReferences("tRNA_fasta"),
    'ncRNA_fasta' : checkParamReturnFileReferences("ncRNA_fasta"),
    'cDNA_fasta' : checkParamReturnFileReferences("cDNA_fasta"),
    'bowtie2_ref_dir' : checkParamReturnFileReferences("bowtie2_ref_dir"),
    'gtf'      : checkParamReturnFileReferences("gtf"),
  ]
}


def idlHeader() {
    // Log colors ANSI codes
    c_reset = "\033[0m";
    c_dim = "\033[2m";
    c_black = "\033[0;30m";
    c_green = "\033[0;32m";
    c_yellow = "\033[0;33m";
    c_blue = "\033[0;34m";
    c_purple = "\033[0;35m";
    c_cyan = "\033[0;36m";
    c_white = "\033[0;37m";

    return """    -${c_dim}--------------------------------------------------${c_reset}-
                                           
    ${c_blue}  ${workflow.manifest.name} ${c_reset}

    ${c_purple}  Independent Data Lab ${c_reset}
                                             
    ${c_purple}  Author: ${workflow.manifest.author}${c_reset}
    ${c_purple}  Home Page: ${workflow.manifest.homePage}${c_reset}
    
    -${c_dim}--------------------------------------------------${c_reset}-
    """.stripIndent()
}

