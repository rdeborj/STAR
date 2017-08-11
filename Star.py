"""
A Python3 wrapper for the STAR RNASeq aligner.
"""

import os
import sys
import yaml


def align_rnaseq(config, star="STAR"):
    """
    Produce the command for aligning FASTQ reads to a reference using the
    STAR aligner.

    Usage:
        align_rnaseq(config, star)

    Input:
          * config: dictionary containing STAR options to use in alignment (required)
          * star: full path to the STAR aligner (default: STAR)

    Output:
        Returns a dictionary containing:
          * cmd: command to execute for alignment
          * output: name of output file alignment is written to
    """
    # flatten the config dictionary into an options string that can be
    # combined with the STAR program to create an executable command line
    # string
    options = ' '.join([' '.join([''.join(['--', key]), str(value)]) for key, value in config.items()])
    cmd = " ".join([star, options])
    return cmd


def get_default_config(file):
    """
    Get the default configuration for the STAR aligner as a dictionary.  Note
    that the values can be overwritten for custom configurations.

    Usage:
        get_default_config(file)

    Input:
        * file: default configuration file in YAML format containing all STAR
                options (required)

    Output:
        Returns a dictionary containing standard STAR options
    """
    with open(file, 'r') as stream:
        try:
            config = yaml.load(stream)
        except yaml.YAMLError as exc:
            print(exc)
    return config


def get_fusions():
    """
    STAR --runMode  alignReads --readFilesIn ./FCRL4_neg_Day_0_CGGCTATG-TAAGATTA_R1.fastq.gz  ./FCRL4_neg_Day_0_CGGCTATG-TAAGATTA_R2.fastq.gz --outFileNamePrefix STAR/FCRL4_neg_Day_0_CGGCTATG-TAAGATTA --genomeDir /cluster/projects/carlsgroup/workinprogress/abdel/20170418_INSPIRE_RNA/STARIndex_hg38_gencode_genomeSAsparseD2/  --runThreadN 4 --chimSegmentMin 20 --readFilesCommand zcat --twopassMode Basic --outSAMprimaryFlag AllBestScore --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outSAMunmapped Within --genomeSAsparseD 2 --limitBAMsortRAM  35000000000
    STAR-Fusion \
        --chimeric_junction  STAR/FCRL4_neg_Day_0_CGGCTATG-TAAGATTAChimeric.out.junction \
        --chimeric_out_sam STAR/FCRL4_neg_Day_0_CGGCTATG-TAAGATTAChimeric.out.sam \
        --ref_GTF /cluster/projects/carlsgroup/workinprogress/abdel/20170418_INSPIRE_RNA//gencode.v26.annotation.gtf \
        --out_prefix STAR_fusion/FCRL4_neg_Day_0_CGGCTATG-TAAGATTA        
    """


def get_star_aligner_output_files(directory):
    """
    Get the output

    Usage:
        get_star_output_files(directory)

    Input:
        * directory: path to the directory containing STAR output
    
    Output:
        Returns a dictionary of output files from the STAR aligner

    NOTES:
        From a standard pipeline run, here's the output from 

    """
    files = os.listdir(directory)

    # what are the default files