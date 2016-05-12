"""
Created on Nov 1, 2011

@author: jcg
@author: Shyam Saladi (saladi@caltech.edu)

"""

import sys
import uuid
import random

import Functions

class Solution:
    """
    A Solution encapsulates a sequence and their inherent attributes:
        sol_id - ID for Solution
        seqence - sequence for Solution
        cds_region - a tuple indicating the location of (Begin,End) of CDS sequence (this will be necessary in the design mode if one want to contrain mutations).
        mutable_region - a list with all positions that can be mutated
        parent - Solution from which the current Solution was derived

    """

    def __init__(self, sol_id=0, sequence="", cds_region=None,
                 keep_aa=False, mutable_region=None, parent=None, design=None):

        if sequence is None:
            sys.stderr.write("Tried to create a solution with sequence NULL\n")
            self.sequence = None
            return None

        # check if solution is in DB
        self.mutable_region = mutable_region
        self.cds_region = cds_region
        self.keep_aa = keep_aa
        self.solid = sol_id
        self.parent = parent
        self.sequence = sequence.lower()
        self.scores = {}
        self.levels = {}
        self.features = {}
        self.designMethod = design
        self.valid = True

    def add_feature(self, feature):
        featureLabel = feature.label + feature.__class__.__name__
        if featureLabel not in self.features:
            self.features[featureLabel] = feature
            # update scores
            self.scores.update(feature.scores)
            # update levels
            if feature.level is not None:
                self.levels[featureLabel + "Level"] = feature.level
            for subfeature in list(feature.subfeatures.values()):
                self.add_feature(subfeature)
        else:
            sys.stderr.write("Feature label already exists!")

        return

    def checkSolution(self, desiredSolution):

        if desiredSolution is None:
            return False

        same = True
        for feature in list(self.designMethod.features.keys()):
            key = feature + "Level"
            same = same & (desiredSolution[key] == 0 or desiredSolution[
                           key] == self.levels[key])

        return same

    def mutate(self, desiredSolution=None, arg_random=False):

        if desiredSolution is None or \
           arg_random or \
           self.designMethod.listDesigns == [] or \
           self.features == {}:
            return self.randomMutation()
        else:
            # get features with targets
            mutable = []
            for feature in list(self.features.values()):
                if feature.defineTarget(desiredSolution):
                    mutable.append(feature)

            if mutable == []:
                return None

            return random.choice(mutable).mutate()
            # return random.choice(mutable).randomMutation()
            # return self.randomMutation()

    def randomMutation(self, pos=None, n_mut=[1, 2]):
        new_seq = Functions.randomMutationOperator(self.sequence,
                                                   self.keep_aa,
                                                   self.mutable_region,
                                                   self.cds_region,
                                                   pos,
                                                   n_mut)
        return Solution(sol_id=str(uuid.uuid4().int),
                        sequence=new_seq,
                        cds_region=self.cds_region,
                        keep_aa=self.keep_aa,
                        mutable_region=self.mutable_region,
                        parent=self,
                        design=self.designMethod)
