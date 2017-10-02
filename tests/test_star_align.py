import unittest
from star import align_rnaseq

class TestAlign(unittest.TestCase):
    def test_align_command(self):
        self.maxDiff = None
        command = align_rnaseq(
            read1='./FCRL4_neg_Day_0_CGGCTATG-TAAGATTA_R1.fastq.gz',
            read2='./FCRL4_neg_Day_0_CGGCTATG-TAAGATTA_R2.fastq.gz',
            prefix='STAR/FCRL4_neg_Day_0_CGGCTATG-TAAGATTA',
            genome_dir='/cluster/projects/carlsgroup/workinprogress/abdel/20170418_INSPIRE_RNA/STARIndex_hg38_gencode_genomeSAsparseD2')
        print(command)
        expected_command = " ".join([
            'STAR',
            '--runMode alignReads',
            '--readFilesIn ./FCRL4_neg_Day_0_CGGCTATG-TAAGATTA_R1.fastq.gz ./FCRL4_neg_Day_0_CGGCTATG-TAAGATTA_R2.fastq.gz',
            '--outFileNamePrefix STAR/FCRL4_neg_Day_0_CGGCTATG-TAAGATTA',
            '--genomeDir /cluster/projects/carlsgroup/workinprogress/abdel/20170418_INSPIRE_RNA/STARIndex_hg38_gencode_genomeSAsparseD2',
            '--runThreadN 4',
            '--chimSegmentMin 20',
            '--readFilesCommand zcat',
            '--twopassMode Basic',
            '--outSAMprimaryFlag AllBestScore',
            '--outFilterIntronMotifs RemoveNoncanonical',
            '--outSAMtype BAM SortedByCoordinate',
            '--quantMode TranscriptomeSAM GeneCounts',
            '--outSAMunmapped Within',
            '--genomeSAsparseD 2',
            '--limitBAMsortRAM 35000000000'])
        print(expected_command)
        self.assertEqual(command, expected_command)
