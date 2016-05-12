"""
Created on Dec 22, 2012

@author: jcg
@author: Shyam Saladi (saladi@caltech.edu)

"""

import os
import sys

import Bio.SeqIO

import Data
import Functions
import SequenceDesigner
import Features
import Features.CAI
import Features.RNADuplex
import Features.Structure
import DesignOfExperiments.Design

class TranslationFeaturesEcoliDesigner(SequenceDesigner.SequenceDesigner):

    def __init__(self, name, seed, design, dbfile, createDB=True):
        SequenceDesigner.SequenceDesigner.__init__(
            self, name, seed, design, dbfile, createDB)

    def configureSolution(self, solution):
        """
        Solution configuration
        """

        if solution.sequence is None:
            return 0

        # Populate solution with desired features
        solution.mutable_region = list(
            range(0, len(solution.sequence)))  # whole region
        solution.cds_region = (49, len(solution.sequence))
        solution.keep_aa = True

        cai_obj = Features.CAI.CAI(
            solution=solution,
            label="cds",
            args={'cai_range': (49, len(solution.sequence)),
                  'mutable_region': list(range(49, len(solution.sequence)))
                  })

        # Look for RBS
        dup_obj1 = Features.RNADuplex.RNADuplexRibosome(
            solution1=solution,
            label="sd16s",
            args={'rnaMolecule1region': (25, 48),
                  'mutable_region': list(range(25, 48))
                  })
        dup_mfe = Features.RNADuplex.RNADuplexMFE(dup_obj1)
        dup_obj1.add_subfeature(dup_mfe)

        #MFE [-30,30]
        st1_obj = Features.Structure.Structure(
            solution=solution,
            label="utr",
            args={'structure_range': (49 - 30, 49 + 30),
                  'mutable_region': list(range(49 - 30, 49 + 30))
                  })
        st_mfe = Features.Structure.StructureMFE(st1_obj)
        st1_obj.add_subfeature(st_mfe)

        solution.add_feature(cai_obj)
        solution.add_feature(dup_obj1)
        solution.add_feature(st1_obj)

    def validateSolution(self, solution):
        """
        Solution validation tests
        """
        if solution.sequence is None or ('?' in list(solution.levels.values())):
            print("SolutionValidator: Level unknown - ", solution.levels,
                  file=sys.stderr)
            solution.valid = False
            return 0

        # check if solution is valid
        valid = True

        designed_region = solution.sequence

        # No internal Promoters
        (score, position, spacer) = Functions.look_for_promoters(designed_region)
        if score >= 15.3990166:  # 0.95 percentile for Promoter PWM scores
            valid = False
            print("SolutionValidator: High Promoter score: ", score,
                  file=sys.stderr)

        # No internal Terminator
        score = Functions.look_for_terminators(designed_region)
        if score >= 90:  # 90% confidence from transtermHP
            valid = False
            print("SolutionValidator: High Terminator score", file=sys.stderr)

        # No restriction enzymes
        if 'ggtctc' in designed_region or 'gagacc' in designed_region:
            print("SolutionValidator: Restriction enzyme found",
                  file=sys.stderr)
            valid = False

        solution.valid = valid

        return valid

def test_TranslationFeaturesEcoliDesignerAllSeeds():
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

    for i, record in enumerate(Bio.SeqIO.parse("tests/seeds.fna", "fasta")):
        print("seed=", i)

        design = DesignOfExperiments.Design.FullFactorial(
            ["sd16sRNADuplexMFE", "utrStructureMFE", "cdsCAI"], design_param)
        # design = DesignOfExperiments.Design.RandomSampling(
        #     ["sd16sRNADuplexMFE","utrStructureMFE","cdsCAI"], design_param, 3000)

        os.makedirs("output/tfec_ff_neutral_s%d" % i, exist_ok=True)

        tfec_designer = TranslationFeaturesEcoliDesigner(
            "tfec",
            str(record.seq),
            design,
            "output/tfec_ff_neutral_s%d" % i,
            createDB=True)
        tfec_designer.run(selection="directional")

if __name__ == '__main__':
    test_TranslationFeaturesEcoliDesignerAllSeeds()
