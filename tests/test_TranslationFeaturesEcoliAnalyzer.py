"""
Created on Nov 1, 2012

@author: jcg
@author: Shyam Saladi (saladi@caltech.edu)
"""

from SequenceAnalyzer import SequenceAnalyzer
from Features import CAI, Structure, RNADuplex
from Functions import validateCDS
from Data import project_dir


class TranslationFeaturesEcoliAnalyzer(SequenceAnalyzer):

    """
    Initializes class that analyzes sequence features
    """

    def __init__(self, input_file, input_type):
        SequenceAnalyzer.__init__(self, input_file, input_type)

    def configureSolution(self, solution):
        solution.valid = validateCDS(solution.sequence[49:])

        if solution.valid:
            # CAI
            cai_obj = CAI.CAI(
                solution=solution,
                label="cds",
                args={'cai_range': (49, len(solution.sequence))})

            # Look for RBS
            dup_obj1 = RNADuplex.RNADuplexRibosome(
                solution1=solution,
                label="sd16s",
                args={'rnaMolecule1region': (25, 48)})
            dup_mfe = RNADuplex.RNADuplexMFE(dup_obj1)
            dup_obj1.add_subfeature(dup_mfe)

            #MFE [-30,30]
            st1_obj = Structure.Structure(
                solution=solution,
                label="utr", \
                args={'structure_range': (49 - 30, 49 + 30)})
            st_mfe = Structure.StructureMFE(st1_obj)
            st1_obj.add_subfeature(st_mfe)

            solution.add_feature(cai_obj)
            solution.add_feature(dup_obj1)
            solution.add_feature(st1_obj)

    def outputStart(self):
        print("gene_name,sd_hyb_energy,mfe_structure,cai")

    def output(self, solution):
        if solution.valid:
            print(solution.solid,
                  solution.scores['sd16sRNADuplexMFE'],
                  solution.scores['utrStructureMFE'],
                  solution.scores['cdsCAI'],
                  sep=",")

def test_TranslationFeaturesEcoliAnalyzer():
    seqAnalyzerTest = TranslationFeaturesEcoliAnalyzer(
        "genomes/partial_ecoli_genome.csv", "CSV")
    seqAnalyzerTest.run()

if __name__ == '__main__':
    test_TranslationFeaturesEcoliAnalyzer()
