'''
Created on Nov 16, 2011

@author: jcg
'''

from Features.Feature import Feature
import Functions
from uuid import uuid4


class PWMScore(Feature):
    """
    PWMScore Feature
        solution - solution where PWM score should be computed
        label - some label to append to the name
        pwm_range - start and end position to calculate PWM score - a tuple in the form (start, end)
        mutable_region - a list with all bases that can be mutated
        cds_region - a pair with begin and end of CDSs - example: (0,100)
        keep_aa - boolean option indicating if in the design mode amino acids should be kept
    """

    def __init__(self, solution=None, label="", args={'pwm': None,
                                                      'pwm_range': (0, 9),
                                                      'mutable_region': None,
                                                      'cds_regions': None,
                                                      'keep_aa': True}):
        # General properties of feature
        Feature.__init__(self, solution=solution, label=label)
        # Specifics of this Feature
        self.pwm = args['pwm']
        self.pwm_range = args['pwm_range']
        self.sequence = solution.sequence[
            self.pwm_range[0]:self.pwm_range[1] + 1]
        self.mutable_region = [
            i -
            self.pwm_range[0] for i in set(
                range(
                    self.pwm_range[0],
                    self.pwm_range[1] +
                    1)) & set(
                args['mutable_region'])] if 'mutable_region' in args else solution.mutable_region
        self.cds_region = args[
            'cds_region'] if 'cds_region' in args else solution.cds_region
        self.keep_aa = args[
            'keep_aa'] if 'keep_aa' in args else solution.keep_aa
        self.set_scores()
        self.set_level()

    def set_scores(self, scoring_function=Functions.analyze_pwm_score):
        self.scores = Functions.appendLabelToDict(
            scoring_function(self.sequence, self.pwm), self.label)

    def mutate(self, operator=Functions.SimplePWMScoreOperator):
        if not self.targetInstructions:
            return None

        new_seq = list(self.solution.sequence)
        mutated_seq = operator(
            self.sequence,
            self.pwm,
            self.targetInstructions['direction'],
            self.mutable_region,
            keep_aa=self.keep_aa)
        if mutated_seq is None:
            return None
        else:
            new_seq[
                self.ntusage_range[0]:self.ntusage_range[1] +
                1] = list(mutated_seq)
        new_seq = "".join(new_seq)

        return Solution.Solution(sol_id=str(uuid4().int), sequence=new_seq, cds_regions=self.cds_regions,
                                 mutable_region=self.mutable_region, parent=self.solution, design=self.solution.designMethod)

import Solution
