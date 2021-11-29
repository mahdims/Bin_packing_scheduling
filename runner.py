import os
import pandas as pd
import time
import pickle as Pick
from GA_AL import GA
from Chromo import inbalance_measure


class Parameters:
    def __init__(self, data):
        self.RGA_flag = 1
        self.MaxTime = 3600
        self.nPop = 20 + 3 * int(data.N)
        self.Maxit = 100 + int(data.N)
        self.Max_noimprove = int(self.Maxit * 0.2)
        self.CRate = 0.7
        self.MRate = 0.3
        self.mu = 0.25  # Mutation Rate
        self.IM = 0.3  # Minimum printing quantity balance
        self.sp = 1.8  # parameter in parents selection
        self.BR_rep = 1  # Bin reduction applied after % iterations


def read_object(FileName, folder):
    WD = os.getcwd()
    if folder == "Input":
        path = WD + f"/Model/{folder}/{FileName}"
    else:
        path = WD + f"/Model/{folder}/{FileName}_ModelSol"

    with open(path, 'rb') as input:
        obj = Pick.load(input, encoding="latin1")
    return obj


def solution_display(sol, runtime):
    print("############### Results ################")

    # print("Lower Bound: %s" %lower_bound(Data))
    print("Total cost of printing all bins: %s" %sol.total_cost)
    print("Total lateness and earliness: %s" %sol.early_lateness)
    print("Number of Bins: %s" %len(sol.Bins))
    print("Algorithm Run Time: %s" %runtime)
    for i, b in enumerate(sol.Bins):
        b.Draw()
        print("Printing Quqntity: %s" %b.quantity)
        print("Items quantity: ", [it.q for it in b.items] )
        print("IM: %s" %inbalance_measure(b.items) )
        print("Revolting Bin= %d" % b.Revolta)


def item_name_correction(data):

    for it in data.items.values():
        it.name = str(it.ID)


if __name__ == "__main__":
    results = []
    # N = 12
    T = 3
    No_reps = 1
    for N in [8, 11, 14]:  # range(15, 17):
        for name in ['Tstr']:
            for rep in [0]:  # range(1,2):
                FileName = "Data_%d_%d_%d_%s" %(N, T, rep, name)
                Data = read_object(FileName, "Input")
                item_name_correction(Data)
                Pars = Parameters(Data)
                Pars.RGA_flag = 1

                best_cost = 0
                Avg_cost = 0
                Avg_time = 0
                for runs in range(No_reps):
                    Listofsolutions = []
                    pop = []
                    start = time.time()
                    Solver = GA(Data, Pars)
                    Best_Sol = Solver.run()
                    Runtime = time.time()-start

                    if runs == 0 or Best_Sol.total_cost < best_cost:
                        best_cost = Best_Sol.total_cost
                        BB_Sol = Best_Sol
                    Avg_cost += Best_Sol.total_cost / No_reps
                    Avg_time += Runtime / No_reps

                print('%d_%d %s %d %s %s %s' % (
                    N, T, name, rep, str(round(best_cost, 1)), str(round(Avg_cost, 1)), str(round(Avg_time, 1))))
                results.append([N, T, name, rep, round(best_cost, 1), round(Avg_cost, 1), round(Avg_time, 1)])

                # display the results
                # solution_display(BB_Sol, Runtime)

    results = pd.DataFrame(results, columns=['N', 'T', 'Mode', 'index', 'B.obj', 'A.obj', 'A.time'])
    if Pars.RGA_flag:
        print("RGA")
        results.to_csv("RGA_New.csv")
    else:
        print("UGA")
        results.to_csv("UGA_New.csv")