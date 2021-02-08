################################################################################
#   Program:   Nextflow riboSeq preprocessing pipeline
#
#   Copyright (c) Independent Data Lab UG, All rights reserved.
#   The use of this source code is governed by a BSD-style license which can be
#   found in the LICENSE file in the root directory of this repository.
#
#   This software is distributed without any warranty; without even the implied
#   warranty of merchantability or fitness for a particular purpose.  See the
#   above copyright notice for more information.
################################################################################

import argparse
import os
import subprocess
from Bio import SeqIO

def extract_rRNA_tRNA(reference_fasta_file):
    # first handle all contaminants that are present in reference transcriptome
    ref = SeqIO.parse(open(reference_fasta_file, "r"), "fasta")
    gene_biotype = ["Mt_rRNA","Mt_tRNA","rRNA"]

    contaminants = []
    for record in ref:
        if any(x in record.description for x in gene_biotype):
            contaminants.append(record)
  
    output_file_name = os.path.splitext(reference_fasta_file)[0] + ".rRNA_tRNA.fa"
    output_file=open(output_file_name,"w")

    SeqIO.write(contaminants, output_file, "fasta")
    output_file.close()
    
    return output_file_name


def main():
    argparser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter)

    filepath_group = argparser.add_argument_group(title="File paths")
         
    filepath_group.add_argument(
        "--reference_fasta_file",
        type=str,
        help="reference fasta file containing non-coding RNA",
        required=True,
    )
    filepath_group.add_argument(
        "--tRNA_fasta",
        type=str,
        help="tRNA fasta file from http://gtrnadb.ucsc.edu/",
        required=True,
    )

    args, unknown = argparser.parse_known_args()

    print("Extracting Mt_tRNA, Mt_rRNA and rRNA from the reference...")
    contaminants_from_ref_file = extract_rRNA_tRNA(args.reference_fasta_file)

    print("Concatenating contaminants from reference and tRNA")
    fin = open(args.tRNA_fasta, "r")
    tRNAs = fin.read()
    fin.close()
    fout = open(contaminants_from_ref_file, "a")
    fout.write(tRNAs)
    fout.close()

    print("Final file with all contaminants: %s" % contaminants_from_ref_file)

    print("Completed !!")



if __name__ == '__main__':
    main()
