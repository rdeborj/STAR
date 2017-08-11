"""
A Python3 wrapper for the STAR fusion program.
"""

import sys
import os
import yaml


def run_fusion(star="STAR-Fusion"):
    """
    NAME
        run_fusion -- Run the STAR fusion program.
    SYNOPSIS
        run_fusion()
    DESCRIPTION
    STAR --runMode  alignReads --readFilesIn ./FCRL4_neg_Day_0_CGGCTATG-TAAGATTA_R1.fastq.gz  ./FCRL4_neg_Day_0_CGGCTATG-TAAGATTA_R2.fastq.gz --outFileNamePrefix STAR/FCRL4_neg_Day_0_CGGCTATG-TAAGATTA --genomeDir /cluster/projects/carlsgroup/workinprogress/abdel/20170418_INSPIRE_RNA/STARIndex_hg38_gencode_genomeSAsparseD2/  --runThreadN 4 --chimSegmentMin 20 --readFilesCommand zcat --twopassMode Basic --outSAMprimaryFlag AllBestScore --outFilterIntronMotifs RemoveNoncanonical --outSAMtype BAM SortedByCoordinate --quantMode TranscriptomeSAM GeneCounts --outSAMunmapped Within --genomeSAsparseD 2 --limitBAMsortRAM  35000000000
    STAR-Fusion --chimeric_junction  STAR/FCRL4_neg_Day_0_CGGCTATG-TAAGATTAChimeric.out.junction --chimeric_out_sam STAR/FCRL4_neg_Day_0_CGGCTATG-TAAGATTAChimeric.out.sam --ref_GTF /cluster/projects/carlsgroup/workinprogress/abdel/20170418_INSPIRE_RNA//gencode.v26.annotation.gtf --out_prefix STAR_fusion/FCRL4_neg_Day_0_CGGCTATG-TAAGATTA        
    """

    

    cmd = " ".join([star, options])

    return cmd
