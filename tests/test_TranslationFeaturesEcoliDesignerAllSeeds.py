'''
Created on Dec 22, 2012

@author: jcg
'''

import sys

import Functions
from SequenceDesigner import SequenceDesigner
from Features.Structure import Structure, StructureMFE
from Features import CAI, RNADuplex
from DesignOfExperiments.Design import RandomSampling, Optimization, CustomDesign, FullFactorial
import Data

class TranslationFeaturesEcoliDesigner(SequenceDesigner):

    def __init__(self, name, seed, design, dbfile, createDB=True):
        SequenceDesigner.__init__(self, name, seed, design, dbfile, createDB)

    def configureSolution(self, solution):
        '''
        Solution configuration
        '''

        if solution.sequence is None:
            return 0

        # Populate solution with desired features
        solution.mutable_region = list(
            range(0, len(solution.sequence)))  # whole region
        solution.cds_region = (49, len(solution.sequence))
        solution.keep_aa = True

        cai_obj = CAI.CAI(
            solution=solution, label="cds", args={
                'cai_range': (
                    49, len(
                        solution.sequence)), 'mutable_region': list(
                    range(
                        49, len(
                            solution.sequence)))})

        # Look for RBS
        dup_obj1 = RNADuplex.RNADuplexRibosome(
            solution1=solution, label="sd16s", args={
                'rnaMolecule1region': (
                    25, 48), 'mutable_region': list(
                    range(
                        25, 48))})
        dup_mfe = RNADuplex.RNADuplexMFE(dup_obj1)
        dup_obj1.add_subfeature(dup_mfe)

        #MFE [-30,30]
        st1_obj = Structure(
            solution=solution,
            label="utr",
            args={
                'structure_range': (
                    49 - 30,
                    49 + 30),
                'mutable_region': list(
                    range(
                        49 - 30,
                        49 + 30))})
        st_mfe = StructureMFE(st1_obj)
        st1_obj.add_subfeature(st_mfe)

        solution.add_feature(cai_obj)
        solution.add_feature(dup_obj1)
        solution.add_feature(st1_obj)

    def validateSolution(self, solution):
        '''
        Solution validation tests
        '''
        if solution.sequence is None or (
                '?' in list(solution.levels.values())):
            sys.stderr.write(
                "SolutionValidator: Level unknown - " + str(solution.levels) + "\n")
            solution.valid = False
            return 0

        # check if solution is valid
        valid = True

        designed_region = solution.sequence

        # No internal Promoters
        (score, position, spacer) = Functions.look_for_promoters(designed_region)
        if score >= 15.3990166:  # 0.95 percentile for Promoter PWM scores
            valid = False
            sys.stderr.write(
                "SolutionValidator: High Promoter score: " +
                str(score) +
                "\n")

        # No internal Terminator
        score = Functions.look_for_terminators(designed_region)
        if score >= 90:  # 90% confidence from transtermHP
            valid = False
            sys.stderr.write("SolutionValidator: High Terminator score\n")

        # No restriction enzymes
        if 'ggtctc' in designed_region or 'gagacc' in designed_region:
            sys.stderr.write("SolutionValidator: Restriction enzyme found\n")
            valid = False

        solution.valid = valid

        return valid

def test_TranslationFeaturesEcoliDesignerAllSeeds():
    # Seed sequence from which mutants will be derived
    seeds = [
    'aacattcttgttaagattatgtgatctttagcgcgggaggaaaatattgatgaaacagcctgcgcccgtttatcagagaattgcgggtcatcaatggcgacatatctggctttctggcgatatacacggttgtcttgagcagttgcgccgcaaattatggcattgtcgttttgatccgtggcgagatttacttatctcagtgggagacgttatcgatcgtgggccgcaaagtttacgttgtctgcagttactggaacaacattgggtttgtgcggtaagaggcaatcatgaacagatggcgatggatgcgctggcatcccagcagatgtctttgtggttgatgaatggcggcgactggtttattgcgctggcagataatcaacagaaacaagcgaaaacggcgctggaaaaatgtcagcatttgccctttattcttgaagtacacagtcgcaccggcaagcatgttattgctcatgccgattatccagatgatgtttatgaatggcaaaaggacgttgatttgcatcaggtcttgtggagccgctcgcgattaggtgaacgccaaaaagggcagggaattacaggtgctgatcatttctggtttggtcatacaccgttgcgacatcgcgtggatattggcaacctgcattatattgataccggtgctgtctttgggggcgaactgactcttgtgcaattgcaataa',
    'ctgaaatatgaattttaacttttagtcattttataaagaggacattttcatgaatcgtattgaacattatcatgactggttacgtgacgcccacgcaatggaaaagcaagccgaatctatgcttgaatccatggccagccgtatagataattatcctgaactacgcgctcgtattgaacaacatcttagtgaaaccaaaaaccagattgttcaactggaaactattcttgatcgtaatgacatttcacgttcagtcattaaagattccatgagtaaaatggctgcgcttgggcagtcaatcggtggtatattcccttctgatgaaatagtcaaaggctctattagcggatatgtcttcgagcaatttgaaatcgcctgttacacctcactattagcagcagcaaaaaatgccggtgatacagcttcaattccaaccatcgaagcgattttaaatgaggaaaagcaaatggccgactggctgattcagaatattccgcaaacaactgagaaatttttaattcgctctgaaactgatggcgtagaagcgaagaaataa',
    'cgctggtgatgggcttaagtatcctgctgctgaaaaaacaggagggatgatgcgccatttacgcaatatttttaatctgggtatcaaagagttgcgcagtctgctcggtgataaagcgatgctgacgctgattgtcttctcgtttacggtgtcggtgtattcgtcagcgaccgttacgccaggatcgttgaacctcgcgccgatcgccattgccgatatggatcaatcgcagttatcgaaccggatcgttaacagcttctatcgtccgtggtttttgccaccggagatgatcaccgccgatgagatggatgccggactggacgccggacgctataccttcgcgataaatattccgcctaattttcagcgtgatgtcctcgccggacgccagccggatattcaggtgaacgtcgatgccacgcgcatgagccaggcatttaccggcaatgggtatatccagaatattatcaacggtgaagtgaacagctttgtcgcgcgctaccgtgataacagcgaaccgttggtatcgctggaaacccggatgcgctttaacccgaacctcgatcccgcgtggtttggcggggtgatggcgatcatcaacaacattaccatgctggcgattgtattgaccggatcggcgctgatccgcgagcgtgaacacggcacggtggaacacttactggtgatgccgataacgccgtttgagatcatgatggcgaagatctggtcgatggggctggtggtgctggtggtatcgggattatcgctggtgctgatggtgaaaggtgtactgggcgtaccgattgaaggctcgatcccgctgtttatgctgggcgtggcgctcagtctgtttgccaccacgtcaatcggcatttttatggggacgatagcgcgttcaatgccgcaactggggctgctggtgattctggtgctgctgccgctgcaaatgctttccggtggttccacgccgcgcgaaagtatgccgcagatggtgcaggacattatgctgaccatgccgacgacacactttgttagcctcgcgcaggccatcctctaccggggtgccggattcgaaatcgtctggccgcagtttctgacgctgatggcaattggcggcgcatttttcaccattgcgctgctgcgattcaggaagacgattgggacaatggcgtaa',
    'ttacaacgataaaaggctgtactttttctttagctcatggattaacacaatgaaattaatcactgcaccatgcagagcattacttgctctgccgttttgctacgccttttctgcggcaggagaagaagcacgtccggcagaacatgacgacacaaaaacacccgcaattacctcgacatcttctccttcatttcgtttttacggcgaattaggggttggtggatatatggatttagagggtgagaataaacataaatacagcgacggtacctatattgaaggtggcctggagatgaagtacggctcctggttcggcctgatttacggcgaaggctggaccgtgcaggccgaccacgacggcaatgcctgggtgccagaccatagctggggtggtttcgagggcggaattaaccgtttctatggcggttatcgtaccaatgatggcaccgaaatcatgctcagtctgcgtcaggattcctcgctggatgacctgcaatggtggggcgatttcacccccgatctgggctacgtcattcccaatacccgcgacattatgactgcgctgaaggtacagaacttaagcggcaactttcgttatagcgtcaccgcgactcctgccggacatcatgatgaaagcaaagcctggctacattttggcaaatacgatcgctatgacgacaaatacacctatccggcaatgatgaacggttacatccagtatgaccttgccgaaggcatcacctggatgaacggtctggaaatcaccgacggcacaggacagctctatctcacgggcctgctaactcctaactttgccgctcgcgcctggcaccataccggacgcgccgacgggctggacgtaccgggaagtgaaagtgggatgatggtgagcgccatgtatgaagcgttaaagggcgtttatctctccaccgcttacacctacgccaaacatcgccctgaccacgctgacgatgaaaccacctctttcatgcagtttggtatctggtacgaatacggcggcggacgtttcgccacggcttttgatagccgcttctacatgaaaaatgcctctcacgatcccagcgaccaaatcttcctgatgcaatatttctactggtaa',
    'tgggtcgtcagcaaaagaaagataacgctgactcacggggagaataaccgtgaacagacgtaattttattaaagcagcctcctgcggggcattgctgacgggcgcgctgccgtctgtcagtcatgcggctgctgaaaaccgcccgccaattccgggatcgctggggatgttgtacgactcgaccttgtgcgtaggctgccaggcttgcgtcaccaagtgtcaggatatcaatttccctgaacgtaacccgcaaggggaacagacctggtcgaacaacgacaaactgtcgccgtataccaataacatcattcaggtgtggaccagcggcacaggggtcaacaaagaccaggaggagaacggctacgcgtacattaagaaacagtgtatgcactgcgtcgatccgaactgtgtctctgtgtgcccggtctctgcactgaaaaaagatccgaaaaccggcattgtccattacgacaaagatgtgtgcaccggctgccgttactgcatggtcgcctgtccgtacaacgtgccgaagtacgactacaacaacccgtttggtgcgctgcataagtgcgagctgtgcaaccagaaaggtgtggaacgtctcgataaaggcggtctacctggctgcgtagaagtgtgcccggcgggcgcggtgattttcggtacgcgtgaagagctgatggcggaggcgaaaaaacgtctggcgctgaagcctggcagcgaataccactatccgcgtcagacgctgaaatctggcgacacttacctgcatacggtgccgaaatattatccgcatctgtacggcgagaaagagggcggcggtactcaggttctggtactgacgggtgtgccttatgaaaatctcgacctgccgaaactggacgatctttctaccggtgcgcgttccgaaaatattcaacacaccctgtataaaggcatgatgctaccactggctgtgctggcgggcttaaccgtgctggttcgtcgcaacaccaaaaacgaccatcacgacggaggagacgatcatgagtcatga',
    'ccagtagcactggctgctggggtgcgttttattcataaagcaaggctgtatgagcgagaaattaaagatagtctatcgcccattacaagaattgtcaccgtatgcgcacaacgccaggacgcacagtactgagcaggtggcacaactggtagaaagtattaagcaattcggctggactaatccggtgctgattgacgaaaagggcgaaattattgcgggtcacggtcgtgttatggcggctgaaatgctcaaaatggattctgttccggtcattgttctgtctggcctgacggatgagcagaagcagcgataa',
    'gtataaacttgacgctttcaaaataaaaagaaaatcgaagcattcacacatgaataaaaaattaatgtatatattcgcaatttttatagttgcagcaattacctgtattagccaacccaagaaaacgacgttgcgtgataaagccatggtgaattatgcctttgattatttaagctcaccgggcagtcttccattcaccacggcagccacggagctttccgcgattcatggtcactcaacgtcgcaatatcgccttggagaattttatcttcatggtagcgacggtaaaccactggattatacacaggcgagatactggtatgagcaatcagcggaacaggaaaatccacgcgcgcaaagtaaactggggtggatctacctcaaaggtctgggggtcaaacccgacacccgtaaagcaattctctggtataaggaagcagctgaacaagggtatgctcatgctcaatatactttaggtttgatctacagaaatggctcaggtattaatgttaaccattatgaatctcaaaaatggttaaaactgaccgccaaacaacattacaaaaatgcggaaagattacttgccgggcttcccgcacattaa',
    'tctcctttgttattactgtcgtgctttcacttctcgcaggagtcctcgtatggtaagcaacgcctccgcattaggacgcaatggcgtacatgatttcatcctcgttcgcgctaccgctatcgtcctgacgctctacatcatttatatggtcggttttttcgctaccagtggcgagctgacatatgaagtctggatcggtttcttcgcctctgcgttcaccaaagtgttcaccctgctggcgctgttttctatcttgatccatgcctggatcggcatgtggcaggtgttgaccgactacgttaaaccgctggctttgcgcctgatgctgcaactggtgattgtcgttgcactggtggtttacgtgatttatggattcgttgtggtgtggggtgtgtga',
    'ggttgttgcagaatatgcaaggatgttgtttttcgttaacggagctgccatgaatctgcctgtaaaaatccgccgtgactggcactactatgcgttcgccattggccttatattcattcttaatggcgtggtggggttactgggatttgaagcaaaaggttggcagacctatgccgtcggtctggtgacgtgggtgattagtttctggctggcggggttgattattcgtcgtcgcgatgaagaaactgaaaacgcccaataa',
    'caattagcaagacatctttttagaacacgctgaataaattgaggttgctatgtctattgtggtgaaaaataacattcattgggttggtcaacgtgactgggaagtgcgtgattttcacggcacggaatataaaacgctgcgcggcagcagctacaatagctacctcatccgcgaagaaaaaaacgtgctgatcgacaccgtcgaccataaattcagccgcgaatttgtgcagaacctgcgtaatgaaatcgatctggcggatatcgattacatcgtgattaaccatgcagaagaggaccacgctggggcgctgaccgaactgatggcacaaattcccgatacgccgatctactgtacagccaacgctatcgactcgataaatggtcatcaccatcatccggagtggaattttaatgtggtgaaaactggcgacacgctggatatcggcaacggcaaacagctcatttttgtcgaaacaccaatgctgcactggccggacagcatgatgacttacctgacaggcgacgcggtgctgttcagtaacgatgctttcggtcaacactactgcgacgagcatctgttcaacgatgaagtggatcagacggagcttttcgagcagtgccagcgttactacgccaatatcctgacgccgttcagccgcctggtaacaccgaaaattaccgagatcctgggctttaacttaccagtcgatatgatagccacttcccacggcgtggtatggcgcgataacccgacgcaaattgtcgagctgtacctgaaatgggcggctgattatcaggaagacagaatcaccattttctacgacaccatgtcgaataacacccgcatgatggctgacgctatcgcccaggggattgcggaaaccgacccacgcgtggcggtgaaaattttcaacgtcgcccgaagcgataaaaacgaaatcctgactaatgtcttccgctcaaaaggcgtgctggtcggcacttcgacgatgaataacgtgatgatgccgaaaatcgccgggctggtggaggagatgactggtttacgcttccgtaacaaacgcgccagtgctttcggctctcacggctggagcggcggtgcggtggatcgtctttccacgcgcctgcaggatgcgggtttcgaaatgtcgcttagcctgaaagcgaaatggcgaccagaccaggacgctctgaagttatgccgtgaacacggtcgcgaaatcgcccgtcagtgggcgctcgcgccgctgccgcagagcacggtgaatacggtagttaaagaagaaacctctgccaccacgacggctgacctcggcccacggatgcagtgcagcgtctgccagtggatttacgatccggcaaaaggcgagccaatgcaggacgttgcgccaggaacgccgtggagtgaagtcccggataacttcctctgcccggaatgctccctcggcaaagacgtctttgaagaactggcatcggaggcaaaatga',
    'gaatcaggctgttaatcataaataagaccacgggccacggaggctatcaatgttgagtatttttaaaccagcgccacacaaagcgcgcttacctgccgcggagatcgatccgacttatcgtcgattgcgctggcaaattttcctggggatattctttggctatgcggcttactatttggttcgtaagaactttgcgcttgctatgccttatctggttgagcagggattctcacgcggtgatttaggttttgccctttcggggatctcgattgcttatggattttcgaaattcatcatgggttcggtatcggatcgctcgaatccgcgcgttttcctgcccgcaggtttgattctggcggcggcagtgatgttgtttatgggctttgtgccatgggcgacgtcgagcattgcggtgatgtttgtactgttgttcctctgcggttggttccaggggatggggtggccgccgtgtggtcgtactatggtgcactggtggtcgcagaaagaacgtggcggcattgtgtcagtgtggaactgtgcgcacaacgtcggtggtggtattccgccgctgctgttcctgctggggatggcctggttcaatgactggcatgcggcgctctatatgcctgctttctgcgccattctggtggcattattcgcctttgcgatgatgcgcgataccccgcaatcctgtggcttgccgccgatcgaagagtacaaaaatgattatccggacgactataacgaaaaagcggaacaggagctgacggcgaagcaaatcttcatgcagtacgtactgccgaacaaactgctgtggtatatcgccatcgccaacgtgttcgtttatctgctgcgttacggcatcctcgactggtcaccgacttatctgaaagaggttaagcatttcgcgctagataaatcctcctgggcctacttcctttatgaatatgcaggtattccgggcactctgctgtgcggctggatgtcggataaagtcttccgtggcaaccgtggggcaaccggcgttttctttatgacactggtgaccatcgcgactatcgtttactggatgaacccggcaggtaacccaaccgtcgatatgatttgtatgattgttatcggcttcctgatctacggtcctgtgatgctgatcggtctgcatgcgctggaactggcaccgaaaaaagcggcaggtacggcagcgggctttaccgggctgtttggttacctgggcggttcggtggcggcgagcgcgattgttggctacaccgtggacttcttcggctgggatggcggctttatggtaatgattggcggcagcattctggcggttatcttgttgattgttgtgatgattggcgaaaaacgtcgccatgaacaattactgcaagaacgcaacggaggctaa',
    'acaattcaagaatagccgcaaaatgttgtcattacaacaaggcggctatatgacgctcgcgcagtttgccatgattttctggcacgacctggcagcaccgatcctggcgggaattattaccgcagcgattgtcagctggtggcgtaaccggaagtaa',
    'atttataaagattaagtaaacacgcaaacacaacaataacggagccgtgatggctggaaacacaattggacaactctttcgcgtaaccaccttcggcgaatcgcacgggctggcgctcggctgcatcgtcgatggtgttccgccaggcattccgctgacggaagcggacctgcaacatgacctcgaccgtcgtcgccctgggacatcgcgctataccacccagcgccgcgagccggatcaggtcaaaattctctccggtgtttttgaaggcgttactaccggcaccagcattggcttgttgatcgaaaacactgaccagcgctctcaggattacagtgcgattaaggacgttttccgtccaggccatgccgattacacctacgaacaaaaatacggtctgcgcgattatcgcggcggtggacgttcttccgcccgcgaaaccgccatgcgcgtggcggcaggagctattgccaaaaaatatctcgccgagaaatttggtattgaaatccgtggctgcctgacccagatgggcgacattccgctggatatcaaagactggtcgcaggtcgagcaaaatccgtttttttgcccggaccccgacaaaatcgacgcgttagacgagttgatgcgtgcgctgaaaaaagagggcgactccatcggcgctaaagtcaccgttgttgccagtggcgttcctgccggacttggcgagccggtctttgaccgcctggatgctgacatcgcccatgcgctgatgagcatcaacgcggtgaaaggcgtggaaattggcgacggctttgacgtggtggcgctgcgcggcagccagaaccgcgatgaaatcaccaaagacggtttccagagcaaccatgcgggcggcattctcggcggtatcagcagcgggcagcaaatcattgcccatatggcgctgaaaccgacctccagcattaccgtgccgggtcgtaccattaaccgctttggcgaagaagttgagatgatcaccaaaggccgtcacgatccctgtgtcgggatccgcgcagtgccgatcgcagaagcgatgctggcgatcgttttaatggatcacctgttacggcaacgggcgcaaaatgccgatgtgaagactgatattccacgctggtaa',
    'taaaaacagcgcggtgtattgtgacgtttttatatctaccgtgaatgttatgaacactatcgtatttgtggaagatgatgcggaagtcggttcactgattgccgcgtacctggcaaaacatgatatgcaggttaccgtagagccgcgcggcgaccaggccgaagaaaccattttgcgagaaaatccggatttggtgttactcgacatcatgctaccaggcaaggacggcatgaccatttgtcgtgatttacgcgcaaagtggtctggaccgattgttcttctaacctctctcgatagcgatatgaaccacatcctggcactggaaatgggtgcctgcgactatattctcaaaacgacgccccctgctgttttgctagcgcgtttacgtttgcatttgcgtcagaatgagcaagccacactgaccaaaggtcttcaggaaacgtctctgactccctacaaagccctgcatttcggcacgttgaccatcgatcccatcaaccgcgtagtcaccctggctaacactgaaatctcgctctcgacagctgatttcgaattattgtgggaattagctacccatgccgggcaaatcatggaccgcgatgcattgctgaaaaatttacgcggcgtcagttatgacggactggatcgtagcgtggacgtggctatttcgcggttaagaaaaaaactgctcgataacgccgcagaaccttatcgcattaaaactgtgcgtaacaaaggctatctttttgcgcctcatgcatgggaataa',
    'cgtcgcactcgatgcttagcaagcgataaacacattgtaaggataacttatgaacaagactcaactgattgatgtaattgcagagaaagcagaactgtccaaaacccaggctaaagctgctctggagtccactctggctgcaattactgagtctctgaaagaaggcgatgctgtacaactggttggtttcggtaccttcaaagtgaaccaccgcgctgagcgtactggccgcaacccgcagaccggtaaagaaatcaaaattgccgcagctaacgtaccggcatttgtttctggcaaggcactgaaagacgcagttaagtaa',
    'gccgaataatcgtcgttggcgaattttacgactctgacaggaggtggcaatgctggttgccgcaggacagtttgctgttacatctgtgtgggaaaagaacgctgagatttgtgcctcgttgatggcgcaggcggcggaaaacgacgcatcgctgtttgccctgccggaagcattgctggcgcgcgatgatcatgatgcagatctatcggttaaatcagcacagctgctggaaggcgaattcctcggactttacggcgagaaagtaaacgtaacatgatgacgacaattctgacgattcatgttccttcaacgccggggcgcgcatggaatatgctggtggcacttcaggcaggaaacatcgtcgcccgttatgccaaactgcatctctatgatgcatttgccattcaggaatcacgccgtgttgatgctggtaatgaaatcgctccgttactggaggtggaagggatgaaggtcggtctgatgacctgttatgacttacgctttccagagctggcgctggcacaggcattacagggagctgaaatcctggtacttcctgccgcctgggttcgcgggccgctcaaagagcatcactggtcaacgttgcttgccgctcgtgcgctggataccacctgttatatggtggcggcgggggagtgcgggaacaaaaatatcggtcaaagccggattatagatccctttggcgtcaccattgcggcagcgtcagaaatgcctgcactcattatggcggaagtgacgcccgaacgtgtgcgtcaggtgcgcgcgcaactgcccgtcttaaacaaccgtcgctttgcgccgccgcaattattatga',
    'aacgactgcccatgtcgatttagaaatagttttttgaaaggaaagcagcatgaaaattaaaactctggcaatcgttgttctgtcggctctgtccctcagttctacagcggctctggccgctgccacgacggttaatggtgggaccgttcactttaaaggggaagttgttaacgccgcttgcgcagttgatgcaggctctgttgatcaaaccgttcagttaggacaggttcgtaccgcatcgctggcacaggaaggagcaaccagttctgctgtcggttttaacattcagctgaatgattgcgataccaatgttgcatctaaagccgctgttgcctttttaggtacggcgattgatgcgggtcataccaacgttctggctctgcagagttcagctgcgggtagcgcaacaaacgttggtgtgcagatcctggacagaacgggtgctgcgctgacgctggatggtgcgacatttagttcagaaacaaccctgaataacggaaccaataccattccgttccaggcgcgttattttgcaaccggggccgcaaccccgggtgctgctaatgcggatgcgaccttcaaggttcagtatcaataa',
    'ttttgcagggatgttgtcgtccctgaaaaagcaaaaatggagaaaaggaatgagtgaatcattacatctgacccgcaatggatcaattctggaaattacccttgatcgtccaaaagcgaatgctattgatgcaaaaaccagctttgaaatgggcgaagtatttctaaatttccgtgacgatccgcaattacgtgtcgccattattaccggtgccggagagaagttcttttccgcgggctgggatttaaaagcggcagcagaaggcgaagcaccggatgctgactttggtccgggtggttttgcgggattaaccgaaattttcaatctcgacaaaccggttatcgcagctgtgaacggctatgcctttggcggcggctttgaactggcgctggcggcagattttattgtttgtgccgataacgccagcttcgccctgccggaagccaaactgggcatcgttcctgacagcggcggtgtgctgcgtctgccgaagatcctgccgcctgccatcgtcaatgaaatggtgatgaccggcagacgaatgggcgcagaagaggcgctgcgttgggggatagtcaaccgcgtggttagccaggcggaactgatggataacgcccgcgaactggctcagcagctggttaacagcgccccgctggcgattgcggcgctgaaagagatctaccgcaccaccagcgaaatgccggtagaagaagcgtatcgctatattcgcagcggcgtgttgaaacactatccatcggttctgcattcggaagatgccattgaagggccgctggcgtttgccgagaagcgcgatccggtgtggaaaggacgttaa',
    'caaccagtaaactacgcgccagttatgtacacactcaggacaaaaaaacgtgacgattaaattgattgtcggcctggcgaaccccggtgctgaatacgccgcaacgcgacataatgctggtgcctggttcgttgacttactggcagagcgtttgcgcgctccgctgcgcgaagaggctaaattctttggttatacttcgcgagtcactcttggaggcgaagatgtccgcctgttagtcccgactacatttatgaatctcagcggcaaagccgttgcggcgatggccagttttttccgcattaatccggacgaaattctggtggcccacgacgaactggatctgcctcctggcgtcgccaaatttaaattgggcggtggccatggtggtcacaatggactgaaagacatcatcagtaaattgggtaataaccctaactttcaccgtttacgcatcggaatcggtcatccgggcgataaaaataaagttgtcggttttgtgttaggcaaaccgcctgttagtgaacagaagttaattgatgaagccattgacgaagcggcgcgttgtactgaaatgtggtttacagatggcttgaccaaagcaacgaaccgattgcacgcctttaaagcgcaataa',
    'atcacaaatgttttttgattgtgaagttttgcacggacggggaagatgaatgaaaaagattgcatttggctgtgatcatgtcggtttcattttaaaacatgaaatagtggcacatttagttgagcgtggcgttgaagtgattgataaaggaacctggtcgtcagagcgtactgattatccacattacgccagtcaagtcgcactggctgttgctggcggagaggttgatggcgggattttgatttgtggtactggcgtcggtatttcgatagcggcgaacaagtttgccggaattcgcgcggtcgtctgtagcgaaccttattccgcgcaactttcgcggcagcataacgacaccaacgtgctggcttttggttcacgagtggttggcctcgaactggcaaaaatgattgtggatgcgtggctgggcgcacagtacgaaggcggtcgtcatcaacaacgcgtggaggcgattacggcaatagagcagcggagaaattga',
    'gccccggcacaggctgcccaggccgttgcgactctataaggacacgataatgacgatttttgataattatgaagtgtggtttgtcattggcagccagcatctgtatggcccggaaaccctgcgtcaggtcacccaacatgccgagcacgtcgttaatgcgctgaatacggaagcgaaactgccctgcaaactggtgttgaaaccgctgggcaccacgccggatgaaatcaccgctatttgccgcgacgcgaattacgacgatcgttgcgctggtctggtggtgtggctgcacaccttctccccggccaaaatgtggatcaacggcctgaccatgctcaacaaaccgttgctgcaattccacacccagttcaacgcggcgctgccgtgggacagtatcgatatggactttatgaacctgaaccagactgcacatggcggtcgcgagttcggcttcattggcgcgcgtatgcgtcagcaacatgccgtggttaccggtcactggcaggataaacaagcccatgagcgtatcggctcctggatgcgtcaggcggtctctaaacaggatacccgtcatctgaaagtctgccgatttggcgataacatgcgtgaagtggcggtcaccgatggcgataaagttgccgcacagatcaagttcggtttctccgtcaatacctgggcggttggcgatctggtgcaggtggtgaactccatcagcgacggcgatgttaacgcgctggtcgatgagtacgaaagctgctacaccatgacgcctgccacacaaatccacggcaaaaaacgacagaacgtgctggaagcggcgcgtattgagctggggatgaagcgtttcctggaacaaggtggcttccacgcgttcaccaccacctttgaagatttgcacggtctgaaacagcttcctggtctggccgtacagcgtctgatgcagcagggttacggctttgcgggcgaaggcgactggaaaactgccgccctgcttcgcatcatgaaggtgatgtcaaccggtctgcagggcggcacctcctttatggaggactacacctatcacttcgagaaaggtaatgacctggtgctcggctcccatatgctggaagtctgcccgtcgatcgccgcagaagagaaaccgatcctcgacgttcagcatctcggtattggtggtaaggacgatcctgcccgcctgatcttcaatacccaaaccggcccagcgattgtcgccagcttgattgatctcggcgatcgttaccgtctactggttaactgcatcgacacggtgaaaacaccgcactccctgccgaaactgccggtggcgaatgcgctgtggaaagcgcaaccggatctgccaactgcttccgaagcgtggatcctcgctggtggcgcgcaccataccgtcttcagccatgcactgaacctcaacgatatgcgccaattcgccgagatgcacgacattgaaatcacggtgattgataacgacacacgcctgccagcgtttaaagacgcgctgcgctggaacgaagtgtattacgggtttcgtcgctaa',
    'catgatattcatcaggaaaacgccatctatttgatggtgaggagactgcgtgactgacgttttactctgtgttggcaatagcatgatgggcgatgatggcgcaggtccgctgctggcggaaaagtgcgccgccgcgccgaaaggtaactgggtggtgattgacggcggtagcgcaccggaaaacgacatcgtcgctatccgtgaactgcgcccgacacgactgctgattgtcgacgccacggatatggggctaaaccccggcgagatccgcatcatcgacccggatgatatcgccgagatgtttatgatgactacccataacatgccgttgaattaccttatcgaccagttgaaagaagatattggcgaagtgattttcctcggcattcagccggatatcgtcggcttttactacccgatgacccagccgattaaagatgcggtagaaaccgtttatcaacgactggaaggctgggaaggaaatggcggcttcgcgcagttagcggtggaagaagagtag',
    'ccgacgatgattacggcctcaggcgacaggcacaaatcggagagaaactatgtttgaaccaatggaacttaccaatgacgcggtgattaaagtcatcggcgtcggcggcggcggcggtaatgctgttgaacacatggtgcgcgagcgcattgaaggtgttgaattcttcgcggtaaataccgatgcacaagcgctgcgtaaaacagcggttggacagacgattcaaatcggtagcggtatcaccaaaggactgggcgctggcgctaatccagaagttggccgcaatgcggctgatgaggatcgcgatgcattgcgtgcggcgctggaaggtgcagacatggtctttattgctgcgggtatgggtggtggtaccggtacaggtgcagcaccagtcgtcgctgaagtggcaaaagatttgggtatcctgaccgttgctgtcgtcactaagcctttcaactttgaaggcaagaagcgtatggcattcgcggagcaggggatcactgaactgtccaagcatgtggactctctgatcactatcccgaacgacaaactgctgaaagttctgggccgcggtatctccctgctggatgcgtttggcgcagcgaacgatgtactgaaaggcgctgtgcaaggtatcgctgaactgattactcgtccgggtttgatgaacgtggactttgcagacgtacgcaccgtaatgtctgagatgggctacgcaatgatgggttctggcgtggcgagcggtgaagaccgtgcggaagaagctgctgaaatggctatctcttctccgctgctggaagatatcgacctgtctggcgcgcgcggcgtgctggttaacatcacggcgggcttcgacctgcgtctggatgagttcgaaacggtaggtaacaccatccgtgcatttgcttccgacaacgcgactgtggttatcggtacttctcttgacccggatatgaatgacgagctgcgcgtaaccgttgttgcgacaggtatcggcatggacaaacgtcctgaaatcactctggtgaccaataagcaggttcagcagccagtgatggatcgctaccagcagcatgggatggctccgctgacccaggagcagaagccggttgctaaagtcgtgaatgacaatgcgccgcaaactgcgaaagagccggattatctggatatcccagcattcctgcgtaagcaagctgattaa',
    'cgggggtgggggtataatgaccattctgttattgcatagagtagttaacatgaagcggagtagaacggaagtggggcgctggcgcatgcagcgtcaggctagccgacgtaaatcgcgttggcttgaggggcaatcgcgccgaaatatgcgtatccacagcatcaggaagtgcattctaaacaaacagcgtaactcgttattgtttgcgatctacaatatctaa',
    'gatgagttatgtagactggccgccattaattttgaggcacacgtactacatggctgaattcgaaaccacttttgcagatctgggcctgaaggctcctatccttgaagcccttaacgatctgggttacgaaaaaccatctccaattcaggcagagtgtattccacatctgctgaatggccgcgacgttctgggtatggcccagacggggagcggaaaaactgcagcattctctttacctctgttgcagaatcttgatcctgagctgaaagcaccacagattctggtgctggcaccgacccgcgaactggcggtacaggttgctgaagcaatgacggatttctctaaacacatgcgcggcgtaaatgtggttgctctgtacggcggccagcgttatgacgtgcaattacgcgccctgcgtcaggggccgcagatcgttgtcggtactccgggccgtctgctggaccacctgaaacgtggcactctggacctctctaaactgagcggtctggttctggatgaagctgacgaaatgctgcgcatgggcttcatcgaagacgttgaaaccattatggcgcagatcccggaaggtcatcagaccgctctgttctctgcaaccatgccggaagcgattcgtcgcattacccgccgctttatgaaagagccgcaggaagtgcgcattcagtccagcgtgactacccgtcctgacatcagccagagctactggactgtctggggtatgcgcaaaaacgaagcactggtacgtttcctggaagcggaagattttgatgcggcgattatcttcgttcgtaccaaaaacgcgactctggaagtggctgaagctcttgagcgtaacggctacaacagcgccgcgctgaacggtgacatgaaccaggcgctgcgtgaacagacactggaacgcctgaaagatggtcgtctggacatcctgattgcgaccgacgttgcagcccgtggcctggacgttgagcgtatcagcctggtagttaactacgatatcccgatggattctgagtcttacgttcaccgtatcggtcgtaccggtcgtgcgggtcgtgctggccgcgcgctgctgttcgttgagaaccgcgagcgtcgtctgctgcgcaacattgaacgtactatgaagctgactattccggaagtagaactgccgaacgcagaactgctaggcaaacgccgtctggaaaaattcgccgctaaagtacagcagcagctggaaagcagcgatctggatcaataccgcgcactgctgagcaaaattcagccgactgctgaaggtgaagagctggatctcgaaactctggctgcggcactgctgaaaatggcacagggtgaacgtactctgatcgtaccgccagatgcgccgatgcgtccgaaacgtgaattccgtgaccgtgatgaccgtggtccgcgcgatcgtaacgaccgtggcccgcgtggtgaccgtgaagatcgtccgcgtcgtgaacgtcgtgatgttggcgatatgcagctgtaccgcattgaagtgggccgcgatgatggtgttgaagttcgtcatatcgttggtgcgattgctaacgaaggcgacatcagcagccgttacattggtaacatcaagctgtttgcttctcactccaccatcgaactgccgaaaggtatgccgggtgaagtgctgcaacactttacgcgcactcgcattctcaacaagccgatgaacatgcagttactgggcgatgcacagccgcatactggcggtgagcgtcgtggcggtggtcgtggtttcggtggcgaacgtcgtgaaggcggtcgtaacttcagcggtgaacgccgtgaaggtggccgtggtgatggtcgtcgttttagcggcgaacgtcgtgaaggccgcgctccgcgtcgtgatgattctaccggtcgtcgtcgtttcggtggtgatgcgtaa',
    'ccggggctatgcttatagcgataatcatactgatgagagagggaaggtcatggatcaggcgctactggacgggggttatcgctgttataccggcgaaaagatcgatgtctatttcaacactgcgatatgtcagcattctggcaattgcgtacgtggcaacggcaagttatttaatctcaaacgaaagccgtggatcatgccggatgaagtcgacgtcgccactgtggttaaagtgattgatacgtgcccgagcggcgcgctgaaataccgtcataaataa',
    'ctggcgagggtttccagatatcatgagttctgattacgcaggagaactcatgatctggataatgctcgccacgctggcggtagtgtttgtggttggttttcgggtgctgacatccggggccagaaaagcgattcgccgtctcagcgatcggctgaacatcgatgtcgtacccgtggagtcgatggtcgatcaaatgggaaagtcagctggtgacgaatttttacgttatttgcatcgtccggatgagtcgcacctgcaaaacgccgcgcaggtgttgctcatctggcaaattgtcattgtcgatggtagcgaacagaacctgctgcaatggcatcggattttacaaaaagctcgcctcgccgcgccgattaccgacgctcaggtcaggctggcgctaggttttctgcgcgaaaccgaacctgaaatgcaggatattaatgcttttcagatgcgctataacgcgttctttcagcctgccgagggcgttcactggttgcattga',
    'ctgcattatttctggcgtcgaatagctattccttaagcaggagcttgtcatggaattcttaatggacccctcaatttgggcggggctactcacgcttgttgttctcgaaattgtgctgggtatcgataacctggtcttcatcgccattcttgctgacaaactgccgccaaaacaacgcgataaagcgcgtttgctggggttatcactggcgctgattatgcgtctggggctgctgtcgctgatttcatggatggtcacgctgaccaaaccgctatttaccgtcatggatttctccttctccggacgcgacctgattatgttgttcggggggatattcttgctgttcaaagcaacaaccgaactgcatgaacggctggaaaaccgcgatcatgattccggccacggtaaaggctacgccagtttctgggtggtcgtcacacagattgtcatccttgacgccgtcttctcgttggatgcggtaattactgcagtagggatggttaaccatctgccggtgatgatggcggcggtagtgattgcgatggcggttatgttgctggcatccaaaccgctgacgcgattcgttaaccagcaccccacggtggtggtgctctgtctgagcttcctgttaatgattggtctgagtctggtggcagaaggtttcggtttccacattccgaaaggttacctgtatgccgcgattggcttctcgatcatcatcgaagtgtttaaccagattgcgcgtcgcaactttattcgccaccagtcgactttgccgctgcgagcgcgtactgccgatgccatcctgcgtttgatgggcgggaaacgtcaggccaatgttcagcacgatgccgataacccgatgccgatgccgatcccggaaggtgcatttgccgaagaagaacgttacatgattaacggcgtactgacgctggcgtcgcgttctctgcgcgggatcatgacgccgcgcggtgaaataagctgggttgacgctaatctcggggtcgatgaaatccgcgagcaactgctctcttcaccgcacagtctgttcccggtatgtcgcggtgaactggatgaaatcatcggtattgtacgtgctaaagaactgctggtggcgctggaagagggcgttgatgtggcggcgattgcttcggcgtctccggcgattatcgtcccggaaaccctcgatccgatcaacctgttgggcgtgctgcgtcgtgctcgcgggagctttgttatcgtgaccaacgagtttggtgtggtacaaggtctggtcacgccgctggatgtgctggaagccattgcgggtgaattcccggacgctgacgaaacgccggaaatcattactgatggtgacggctggctggtaaaaggcggtacagatttgcatgccttgcagcaggcgcttgatgttgagcaccttgccgatgacgatgatatcgcgacggtcgcgggcctcgtgatctcggcaaatggtcacattccccgtgtgggcgatgtgattgatgtagggccactgcatatcaccatcattgaagccaatgattatcgtgttgatctggttcgcattgttaaagagcaaccggcgcacgatgaagatgagtaa',
    'ctaacgcatgctagtttaatgacataaggtaggtgaaacggagattggagtgaaaaagtttcgatgggtcgttctggttgtcgtggtgttggcttgcttgctgctttgggcgcaggtattcaacatgatgtgcgatcaggatgtacaatttttcagcggaatttgtgccattaaccagtttatcccgtggtga',
    'ataaaagttatctcccttctcgttcatcgttccatatttgagaaacagtatgtcttccagagttttgaccccggacgtcgttggtattgacgccctggtacacgatcaccaaaccgttctggcaaaagctgaaggcggtgtggttgccgtatttgctaacaatgccccggcgttttatgccgtcacgcctgcacgcctggctgaactgctggcgctggaagaaaagctggcgcgtccgggaagcgatgtcgctctggacgatcaactctatcaggaaccgcaagccgctcccgttgctgtacccatggggaaattcgccatgtatccggactggcaacccgatgccgattttatccgcctggcggcgctatggggcgtggcgctaagagagccggtgaccaccgaagaactggcctcattcattgcctactggcaggcggaaggtaaagtctttcaccatgtgcagtggcaacaaaaactggcgcgcagcctgcaaatcggtcgtgccagcaacggcggactgccgaaacgagatgtgaatacggtcagcgaacctgacagccaaattccaccaggattcagagggtaa'
]
    # Design Methodology and thresholds
    design_param = {"sd16sRNADuplexMFE": {'type': 'REAL',
                                          'thresholds': {'1': (-12.7, -7.3),
                                                         '2': (-7.3, -5.8),
                                                         '3': (-5.8, -5.2),
                                                         '4': (-5.2, -3.3),
                                                         '5': (-3.3, 2.0)}},
                    "utrStructureMFE": {'type': 'REAL',
                                        'thresholds': {'1': (-29.2, -12.2),
                                                       '2': (-12.2, -9.95),
                                                       '3': (-9.95, -8.4),
                                                       '4': (-8.4, -6.73),
                                                       '5': (-6.73, 0.65)}},
                    "cdsCAI": {'type': 'REAL',
                               'thresholds': {'1': (0.13, 0.29),
                                              '2': (0.29, 0.33),
                                              '3': (0.33, 0.37),
                                              '4': (0.37, 0.42),
                                              '5': (0.42, 0.86)}}
                    }

    for i, seed in enumerate(seeds):
        print("seed=", i)

        design = FullFactorial(
            ["sd16sRNADuplexMFE", "utrStructureMFE", "cdsCAI"], design_param)
        #design = RandomSampling(["sd16sRNADuplexMFE","utrStructureMFE","cdsCAI"],design_param, 3000)
        #tfec_designer = TranslationFeaturesEcoliDesigner("tfec", seed, design, "/Users/jcg/Documents/workspace/D-Tailor/testFiles/outputFiles/tfec_ff_rnd_"+str(i), createDB=True)
        tfec_designer = TranslationFeaturesEcoliDesigner(
            "tfec",
            seed,
            design,
            Data.project_dir + "/testFiles/outputFiles/tfec_ff_neutral_s" + str(i),
            createDB=True)
        tfec_designer.run(selection="directional")

if __name__ == '__main__':
    test_TranslationFeaturesEcoliDesignerAllSeeds()
