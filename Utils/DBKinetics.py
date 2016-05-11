'''
Created on Apr 4, 2012

@author: jcg
'''

from os import path
from sqlite3 import connect, Row
import sys


def DBKinetics(db_file, resolution=50):
    """
        Prints a time series of desired solutions found as a function of generated solutions generated
    """
    # Create connection to DB
    con = connect(db_file)
    con.isolation_level = None
    con.row_factory = Row
    cur = con.cursor()

    cur.execute("select count(1) from generated_solution")
    total_sol = cur.fetchone()[0]

    # print "Generated Solutions \t Desired Solutions Found"

    #it = range(0,total_sol,total_sol/resolution)
    it = list(range(0, total_sol, 10))
    it.append(total_sol)

    # get statistics
    for i in it:
        cur.execute(
            "select count(DISTINCT code) as c from (select des_solution_id as code from generated_solution LIMIT " +
            str(i) +
            ")")
        result = cur.fetchone()

        # print i,"\t",result['c']
        print(result['c'])

    con.close()


if __name__ == '__main__':
    pass
    if len(sys.argv) == 2:
        db_file = sys.argv[1]
    else:
        db_file = "../testFiles/outputFiles/tfec_ff_1.sqlite"

    for i in range(1, 31):
        db_file = "../testFiles/outputFiles/tfec_ff_rnd_samp_s" + \
            str(i) + ".sqlite"
        print("seed" + str(i) + ",end")
        DBKinetics(db_file)

    # One-time
    # DBKinetics(db_file)
