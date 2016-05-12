"""
Created on Nov 1, 2012

@author: jcg
@author: Shyam Saladi (saladi@caltech.edu)

"""

import sys
import csv

import Bio.SeqIO

import Solution

class SequenceAnalyzer(object):
    """
    Initializes class that analyzes sequence features
    """

    def __init__(self, input_file, input_type="fasta", sep=","):
        if input_type.upper() == "CSV":
            self.readCSV(input_file, sep=sep)
        else:
            self.readSeqIO(input_file, type=input_type)

    def readCSV(self, input_file, sep):
        list_seq = []

        reader = csv.DictReader(open(input_file), delimiter=sep, quotechar='"')

        for l in reader:
            list_seq.append(l)

        self.list_of_input_sequences = list_seq
        return

    def readSeqIO(self, input_file, input_type):
        """
        Use Bio.SeqIO to try parsing records
        """

        with open(input_file, 'r+') as fh:
            record_list = list(SeqIO.parse(fh, input_type))

        for i, record in record_list:
            record_list[i] = record.seq

        self.list_of_input_sequences = record_list
        return

    def configureSolution(self, solution):
        pass

    def outputStart(self):
        pass

    def output(self, solution):
        pass

    def run(self):

        self.outputStart()

        for sequence in self.list_of_input_sequences:
            sol_id = sequence['name']
            seq = sequence['sequence']

            solution = Solution.Solution(sol_id=sol_id, sequence=seq)
            self.configureSolution(solution)

            self.output(solution)
