from os import path
import argparse
import time
import pickle as Pick
from GA import GA
from utils.Chromo import inbalance_measure
from utils.Input import Input
import warnings
warnings.filterwarnings("error")


BASEDIR = path.dirname(path.realpath(__file__))


class Parameters:
    def __init__(self, N):
        self.RGA_flag = 1  # 1 -> Will run RGA | 0 -> Will run UGA
        self.MaxTime = 3600 # Termination criteria (Total run time limit)
        self.nPop = 30 + int(N) # Number of the solution in Genetic alg population
        self.Maxit = 100 + 3 * int(N) # Termination criteria (maximum number of generation)
        self.Max_noimprove = int(self.Maxit * 0.2) # Termination criteria (maximum number of generation without improvemnt)
        self.CRate = 0.7 # Percentage of new solution created by crossover
        self.MRate = 0.3 # Percentage of new solutions created by mutation
        self.mu = 0.3  # Mutation Rate (percentage of genes changed in each mutation)
        self.IM = 0.2  # Minimum printing quantity balance
        self.sp = 1.5  # parameter in parents selection
        self.BR_rep = 1  # Bin reduction applied after % iterations


def read_object(FileName, folder):

    if folder == "Input":
        path = BASEDIR + f"/Data/{folder}/{FileName}"
    else:
        path = BASEDIR + f"/Data/{folder}/{FileName}_ModelSol"

    with open(path, 'rb') as input:
        obj = Pick.load(input, encoding="latin1")
    return obj


def solution_display(sol, runtime):
    print("\033[94m------------ Best Solution ------------\033[0m")
    print("Total cost of printing all bins: %s" %sol.total_cost)
    print("Total lateness and earliness: %s" %sol.early_lateness)
    print("Number of Bins: %s" %len(sol.Bins))
    print("Algorithm Run Time: %s" %runtime)
    print("Bin to days %s" %[(i,len(d.bin2print)) for i,d in enumerate(sol.days)])
    for i, b in enumerate(sol.Bins):
        b.Draw()
        print("Printing Quantity: %s" %b.quantity)
        # print("Items quantity-due: ", [(it.name, it.q, it.d) for it in b.items] )
        print("IM: %s" %inbalance_measure(b.items))
        print("Is a reversal Bin= %d" % b.Revolta)


def item_name_correction(data):

    for it in data.items.values():
        it.name = str(it.ID)
        it.e = data.Earliness
        it.l = data.Lateness


if __name__ == "__main__":

    parser = argparse.ArgumentParser(description='')
    # Required positional argument

    parser.add_argument('which_GA', type=str,
                        help="Specify the GA version you want to use: RGA or UGA")
    parser.add_argument('inputFile', type=str,
                        help='Name of the input file located in Data/Input directory')

    # Optional argument
    parser.add_argument('--rep', type=int,
                        help='How many time you want to run GA on the input file? (Default: 1)')

    parser.add_argument('--draw', action='store_true',
                        help='Flag to draw the final solution')
    args = parser.parse_args()

    FileName = args.inputFile
    print(f"We are solving {FileName}")

    if args.which_GA not in ['RGA', 'UGA']:
        exit("The GA version is not correct enter RGA or UGA")

    if args.rep:
        No_reps = args.rep
    else:
        No_reps = 1

    Data = read_object(FileName, "Input")
    item_name_correction(Data)
    Pars = Parameters(Data.N)
    Pars.RGA_flag = args.which_GA == 'RGA'
    if Pars.RGA_flag:
        GA_type = 'RGA'
    else:
        GA_type = 'UGA'

    best_cost = 0
    Avg_cost = 0
    Avg_time = 0
    for runs in range(No_reps):
        start = time.time()
        Solver = GA(Data, Pars)
        Best_Sol = Solver.run()
        Runtime = time.time()-start

        if runs == 0 or Best_Sol.total_cost < best_cost:
            best_cost = Best_Sol.total_cost
            BB_Sol = Best_Sol
        Avg_cost += Best_Sol.total_cost / No_reps
        Avg_time += Runtime / No_reps

    if args.draw:
        solution_display(BB_Sol, Avg_time)

    print("\033[94m---------Summery of Results--------\033[0m")
    print(f"We solved {FileName} with {GA_type} for {No_reps} times")
    print(f"Best Cost: {round(best_cost, 2)}")
    print(f"Avg. Cost {str(round(Avg_cost, 2))}")
    print(f"Avg. Time {str(round(Avg_time, 2))}")



