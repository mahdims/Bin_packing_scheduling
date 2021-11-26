import os
import pandas as pd
import time
import pickle as Pick
from GA_AL import GA


class Parameters:
    def __init__(self, data):
        self.RGA_flag = 0
        self.MaxTime = 3600
        self.nPop = 20 + 3 * int(data.N)
        self.Maxit = 100 + int(data.N)
        self.Max_noimprove = int(self.Maxit * 0.2)
        self.CRate = 0.7
        self.MRate = 0.3
        self.mu = 0.25  # Mutation Rate
        self.IM = 0.3
        self.sp = 1.8  # parameter in parents selection


def read_object(FileName, folder):
    WD = os.getcwd()
    if folder == "Input":
        path = WD + f"/Model/{folder}/{FileName}"
    else:
        path = WD + f"/Model/{folder}/{FileName}_ModelSol"

    with open(path, 'rb') as input:
        obj = Pick.load(input, encoding="latin1")
    return obj



if __name__ == "__main__":


    results = []
    N = 12
    T = 3
    No_reps = 1
    for N in range(6, 17):
        for name in ['Tlos', 'Tstr']:
            for rep in [0]: # range(1,2):
                FileName = "Data_%d_%d_%d_%s" %(N,T,rep,name)
                Data = read_object(FileName,"Input")
                Pars = Parameters(Data)

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
                    Avg_cost += Best_Sol.total_cost / No_reps
                    Avg_time += Runtime / No_reps

                    print('%d_%d %s %d %s %s %s' % (
                    N, T, name, rep, str(round(best_cost, 1)), str(round(Avg_cost, 1)), str(round(Avg_time, 1))))
                    results.append([N, T, name, rep, round(best_cost, 1), round(Avg_cost, 1), round(Avg_time, 1)])

                # display the results
                    #print("############### Results ################")
                    #
                    #print("Lower Bound: %s" %lower_bound(Data))
                    #print("Total cost of printing all bins: %s" %Best_Sol.total_cost)
                    #print("Total lateness and earliness: %s" %Best_Sol.early_lateness)
                    #print("Number of Bins: %s" %len(Best_Sol.Bins))
                    #print ("Algorithm Run Time: %s" %Runtime)
                    # for i,b in enumerate(pop[0].Bins):
                        # Draw_the_Bin(Data,b)
                        # print("Printing Quqntity: %s" %b.quantity)
                        # print("Items quantity: ", [it.q for it in b.items] )
                        # print("IM: %s" %inbalance_measure(b.items) )
                        # print("Revolting Bin= %d" % b.Revolta)

    results = pd.DataFrame(results, columns=['N', 'T', 'Mode', 'index', 'B.obj', 'A.obj', 'A.time'])
    results.to_csv("RGA_New.csv")
