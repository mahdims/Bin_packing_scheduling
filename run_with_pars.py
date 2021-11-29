import argparse
import logging
import sys
from runner import Parameters
from Input import Input
import pickle as Pick
import time
from GA_AL import GA


def main(data, parameters):
    # just a test
    start = time.time()
    Solver = GA(data, parameters)
    Best_Sol = Solver.run()
    Runtime = time.time() - start

    # save the fo values in DATFILE
    print(Best_Sol.total_cost)
    # with open("./tunning/" + DATFILE, 'w') as f:
    #    f.write(Best_Sol.total_cost)


if __name__ == "__main__":

    ap = argparse.ArgumentParser(description='Feature Selection using GA with DecisionTreeClassifier')

    ap.add_argument('--nPop', dest='nPop', type=int, required=True, help='Population size')
    ap.add_argument('--Maxit', dest="Maxit", type=int, required=True)
    ap.add_argument('--Noimp', dest='Noimp', type= float, required=True)
    ap.add_argument('--CRate', dest='CRate', type=float, required=True, help='Crossover probability')
    ap.add_argument('--MRate', dest='MRate', type=float, required=True, help='Mutation probability')
    ap.add_argument('--Mu', dest='Mu', type=float, required=True)
    ap.add_argument('--IM', dest='IM', type=float, required=True)
    ap.add_argument('--Sp', dest='Sp', type=float, required=True)
    ap.add_argument('--BR', dest='BR', type=int, required=True)
    ap.add_argument('--i', dest='InputFile', nargs='+', type=str, required=True)
    # 1 arg file name to save and load fo value
    # ap.add_argument('--datfile', dest='datfile', type=str, required=True,
    #                help='File where it will be save the score (result)')
    args = ap.parse_known_args()

    rest_of_input = " ".join(args[1])

    args = args[0]
    # args = ap.parse_args()
    logging.debug(args)

    args.InputFile = args.InputFile[0] + ' ' + rest_of_input

    with open(args.InputFile, 'rb') as input_obj:
        Data = Pick.load(input_obj, encoding="latin1")

    # call main function passing args
    parameters = Parameters(Data)
    parameters.nPop = 20 + args.nPop * int(Data.N)
    parameters.Maxit = 100 + int(Data.N) * args.Maxit
    parameters.Max_noimprove = int(args.Noimp * parameters.Maxit)
    parameters.CRate = args.CRate
    parameters.MRate = args.MRate
    parameters.mu = args.Mu
    parameters.IM = args.IM
    parameters.sp = args.Sp
    parameters.BR_rep =args.BR

    parameters.MaxTime = 170
    parameters.RGA_flag = 1

    main(Data, parameters)
