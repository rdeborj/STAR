import unittest
from star import get_fusions

class TestFusion(unittest.TestCase):
    def test_fusion_command(self):
        fusions = get_fusions(
            junction='STAR/LTS-035_T_RNA_GCCAAT_L007Chimeric.out.junction',
            prefix='STAR_fusion/LTS-035_T_RNA_GCCAAT_L007',
            out_sam='STAR/LTS-035_T_RNA_GCCAAT_L007Chimeric.out.sam',
            gtf='/cluster/projects/carlsgroup/workinprogress/abdel/20170418_INSPIRE_RNA/gencode.v26.annotation.gtf')
        expected_command = " ".join([
            'STAR-Fusion',
            '--chimeric_junction STAR/LTS-035_T_RNA_GCCAAT_L007Chimeric.out.junction',
            '--chimeric_out_sam STAR/LTS-035_T_RNA_GCCAAT_L007Chimeric.out.sam',
            '--ref_GTF /cluster/projects/carlsgroup/workinprogress/abdel/20170418_INSPIRE_RNA/gencode.v26.annotation.gtf',
            '--out_prefix STAR_fusion/LTS-035_T_RNA_GCCAAT_L007'])
        self.assertEqual(fusions['command'], expected_command)
