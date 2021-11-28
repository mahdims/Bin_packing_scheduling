from Input import Input
from Chromo import Chromo, inbalance_measure,ItemSet
from lowerbound import lower_bound
import random
import os
import pandas as pd
import time
import math
import numpy as np
import pickle as Pick


class GA:
    def __init__(self, Data, Pars):
        self.Data = Data
        self.Pars = Pars
        self.RGA_flag = Pars.RGA_flag
        self.pop = []
        self.Listofsolutions = []
        self.iterationNO = 0

        self.Maxit = Pars.Maxit
        self.Max_noimprove = Pars.Max_noimprove
        self.MaxRunTime = Pars.MaxTime
        self.nPop = Pars.nPop
        self.Varsize = Data.N
        self.IM = Pars.IM
        self.mu = Pars.mu
        self.mutation_rate = Pars.MRate
        self.crossover_rate = Pars.CRate
        self.sp = Pars.sp  # parameter in parents selection

    def initialpop(self):
        # First calculate the minimum bin number needed to allocate all itmes
        Bmin = self.Data.MinBinNo + math.ceil(2.5*self.Data.MinBinNo)

        Value = [0] * self.Varsize
        items = ItemSet(self.Data, list(self.Data.items.values()) )
        items.Qsort()
        tempBinIt = ItemSet(self.Data, [items[0]])
        BN = 0
        while len(items.list) != 0:
            it = items.list[0] # select the unassigned item with the largnest quantity
            # if adding this item to current bin is not going to change the IM
            if inbalance_measure(tempBinIt.list+[it]) <= self.IM:
                del items.list[0]
                tempBinIt.add([it])
            else:
                # set value of all items in current bin
                for i in tempBinIt.list:
                    Value[i.ID] = BN
                BN += 1  # create a new bin
                tempBinIt = ItemSet(self.Data, [items[0]])
                del items.list[0]

        (Num_Bin, Value) = self.Calc_Bin_No(Value)
        # randomly decide on revolting option
        Revolting = list(np.random.randint(2, size=(1, Num_Bin))[0])
        # check if there is two-sided option
        Revolting = [a*b for a, b in zip(Revolting, self.Bin_can_revolt(Value))]
        sol = Chromo(self.Pars, self.iterationNO, Value, Revolting)
        self.pop.append(sol.Fitness_Calc(self.Data, self.pop))
        while len(self.pop) < self.nPop:

            # Generating Value genes
            Value=[]
            for _ in range(self.Varsize):
                Value.append(random.randint(0, Bmin))

            # Update the "Value" with ordered bin numbers
            (Num_Bin,Value) = self.Calc_Bin_No(Value)

            ### Genearting revolting genes
            Revolting = list(np.random.randint(2, size=(1, Num_Bin))[0])
            Revolting = [a*b for a, b in zip(Revolting, self.Bin_can_revolt(Value))]

            if self.Is_solution_new(Value, Revolting):
                sol = Chromo(self.Pars,self.iterationNO, Value, Revolting)
                self.pop.append(sol.Fitness_Calc(self.Data, self.pop))

    def mutation(self, sol):
        ## Mutation for value ##
        gen2change = int(math.ceil(self.mu*self.Varsize))
        Value = sol.value
        (Bin_no, Value) = self.Calc_Bin_No(Value)
        rnd_value = random.random()

        if rnd_value <= 0.5:
            # Swap two items between their bins
            for _ in range(gen2change):
                rep = np.random.choice(self.Varsize, 2, replace=False)
                rep.sort()
                Value = Value[:rep[0]]+[Value[rep[1]]] + Value[rep[0]+1:rep[1]] + [Value[rep[0]]]+Value[rep[1]+1:]

        elif rnd_value <= 0.5:
            ## change one item bin ##
            for _ in range(gen2change):
                item2change = np.random.choice(self.Varsize,replace=False)
                newbin = np.random.randint(Bin_no)
                Value[item2change] = newbin
        else:

            Value = sol.Random_bin_no_change(self.Data)

       ## Mutation for revolting ###
        Revolting = sol.Revolting
        (Bin_no, Value) = self.Calc_Bin_No(Value)

        if Bin_no <= len(Revolting):
            Revolting = Revolting[:Bin_no]
        else:
            Revolting = Revolting + list(np.random.randint(2, size=(1, Bin_no-len(Revolting)))[0])
        gen2change = int(math.ceil(self.mu * Bin_no))
        rep = np.random.choice(Bin_no, gen2change, replace=False)
        for i in rep:
            Revolting[i] = 1-Revolting[i]
        Revolting = [a*b for a, b in zip(Revolting, self.Bin_can_revolt(Value))]

        return Chromo(self.Pars, self.iterationNO,Value, Revolting)

    def crossover(self, DadSol, MomSol):

        Dad_value = DadSol.value
        Mom_value = MomSol.value
        child_Val = [[], []]

        if random.random() <= 0.5:
            """ one-point crossover """
            x = np.random.randint(0, self.Varsize)
            child_Val[0] = Dad_value[:x]+Mom_value[x:]
            child_Val[1] = Mom_value[:x]+Dad_value[x:]
        else:
            """ two-point crossover """
            (x, y) = np.random.choice(self.Varsize, 2, False)
            if x > y:
                x, y = y, x
            child_Val[0] = Mom_value[:x]+Dad_value[x:y]+Mom_value[y:]
            child_Val[1] = Dad_value[:x]+Mom_value[x:y]+Dad_value[y:]

        """ Revolting part crossover  """
        Dad_Revolting = DadSol.Revolting
        Mom_Revolting = MomSol.Revolting
        Dad_Bin_no = len(Dad_Revolting)
        Mom_Bin_no = len(Mom_Revolting)
        # Varsize = min(Dad_Bin_no, Mom_Bin_no)
        Child_Bin_no = [[], []]
        child_Rev = [[], []]
        child = [[], []]
        # calculate children bin number
        for i in [0, 1]:
            (Child_Bin_no[i], child_Val[i]) = self.Calc_Bin_No(child_Val[i])
            for _ in range(Child_Bin_no[i]):
                if random.random() <= 0.5:
                    x = np.random.randint(0, Dad_Bin_no)
                    child_Rev[i].append(Dad_Revolting[x])
                else:
                    x = np.random.randint(0, Mom_Bin_no)
                    child_Rev[i].append(Mom_Revolting[x])

            child_Rev[i] = [a*b for a, b in zip(child_Rev[i], self.Bin_can_revolt(child_Val[i]))]

            child[i] = Chromo(self.Pars, self.iterationNO, child_Val[i], child_Rev[i])

        return child

    def Bin_can_revolt(self, Value):
        Num_Bin = max(Value)
        Revolta = np.zeros((Num_Bin+1))
        for b in range(Num_Bin+1):
            Bin_Items = np.where(np.array(Value) == b)[0]
            for it in Bin_Items:
                if self.Data.items[it].two_side == 1:
                   Revolta[b] = 1
                   break

        return Revolta

    def Is_solution_new(self, Value, Revolting):
        # Identify if it ios new and add it to solution pool if it is.
        indicator = 0
        a = Value+Revolting
        if a not in self.Listofsolutions:
            self.Listofsolutions.append(a)
            indicator = 1
        return indicator

    def Calc_Bin_No(self, Value):
        Item_Bin = []
        i = 0
        Num_Bin = max(Value)
        Value=np.array(Value)
        for b in range(Num_Bin+1):
            if b in Value:
                Item_Bin.append(np.array(list(self.Data.items.values()))[np.where(Value == b)] )
                Value[[it.ID for it in Item_Bin[-1]]] = i
                i+=1
        Value=list(Value)
        Num_Bin=len(Item_Bin)

        return Num_Bin, Value

    def evolve(self):
        ############################# parents Selection ################################
        # number of the parents
        parents_length = int(self.nPop/2)
        pv = []
        # calculate the parents selection probability
        for r, individual in enumerate(self.pop):
            rank = float(self.nPop-r-1)
            pv.append(round((2-self.sp)/self.nPop+2*rank*(self.sp-1)/(self.nPop*(self.nPop-1)),5))
        pv = np.array(pv)
        pv /= pv.sum()
        # selecting the parents
        parents = roulette_wheel_pop(self.pop, pv, parents_length)
        #################### mutate some individuals###########################
        Mutation_number = math.ceil(self.mutation_rate*self.nPop)
        counter = 1
        Mutants = []
        Mut_inner_counter = 0
        while counter <= Mutation_number and Mut_inner_counter <= 20*Mutation_number:
            # Select the individual
            individual = parents[np.random.randint(len(parents))]
            individual = self.mutation(individual)

            if self.Is_solution_new(individual.value, individual.Revolting):
                Mutants.append(individual)
                counter += 1
                Mut_inner_counter = 0
            else:
                Mut_inner_counter += 1

        ########################### crossover ##############################
        Crossover_number = int(self.crossover_rate * self.nPop)
        crosscounter = 0
        children = []
        while len(children) <= Crossover_number and crosscounter <= 2*Crossover_number:
            crosscounter += 1

            (male, female) = np.random.choice(parents_length,2,False)

            (child1, child2) = self.crossover(parents[male], parents[female])

            if self.Is_solution_new(child1.value, child1.Revolting):
                children.append(child1)

            if self.Is_solution_new(child2.value, child2.Revolting):
                children.append(child2)

        children = [x.Fitness_Calc(self.Data, self.pop) for x in children]
        Mutants = [x.Fitness_Calc(self.Data, self.pop) for x in Mutants]

        # create the pool
        pool = self.pop[:parents_length]
        # Since th fitness depend on the iteration number we have to re calculate the fitness for old solutions
        for sol in pool:
            sol.Fitness_measure(self.Data.N, self.iterationNO, self.pop)

        pool.extend(children)
        pool.extend(Mutants)
        # evaluate the pool
        pool = sorted(pool, key=lambda pool: pool.Fitness_Value, reverse=False)
        # truncate the pool and create the new generation
        self.pop=pool[0:self.nPop]

        return

    def run(self):

        start = time.time()

        self.iterationNO = 1
        self.initialpop()  # generate initial solutions
        self.pop = sorted(self.pop, key=lambda x: x.Fitness_Value, reverse=False)
        current_bestsol = self.pop[0]
        noimprove = 0
        while self.iterationNO <= self.Maxit and noimprove < self.Max_noimprove and time.time() <= start+self.MaxRunTime:

            self.evolve()

            last_bestsol = current_bestsol

            current_bestsol = self.pop[0]

            if current_bestsol.Fitness_Value == last_bestsol.Fitness_Value:
                noimprove += 1
            else:
                noimprove = 0
            #print("Iteration %s" % iterationNO + "--#searched solutions= %s" % len(Listofsolutions))
            #print("Best objective: %s" %(current_bestsol.Fitness_Value) )

            #AvgNOBins=np.average([len(a.Bins) for a in pop]    )
            #print("Average number of bins = %s" % AvgNOBins)
            #AvgEL=np.average([a.early_lateness for a in pop]    )
            #print("Average lateness earliness measure = %s" % AvgEL)

            self.iterationNO += 1
        return current_bestsol


def roulette_wheel_pop(pop, p, number):
    chosen=np.random.choice(len(pop),number,False,p)
    chosen = [pop[a] for a in chosen]
    return chosen




