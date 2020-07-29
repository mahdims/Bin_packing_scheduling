
"""
Created on Wed Sep 26 12:43:52 2018

@author: mahdi
"""

from Draw_the_Bin import Draw_the_Bin
from Input import Input
from Chromo import Chromo
import random
import time
import math
import numpy as np
  


def initialpop(nPop,Data):
   global pop
   # First calculate the minimum bin number needed to allocate all itmes 
   Bmin=Data.MinBinNo+ math.ceil(2.5*Data.MinBinNo)
   while len(pop) < nPop:
        
        Value=[]
        for _ in range(Varsize):
            Value.append(random.randint(0,Bmin))
        
        if Value not in Listofsolutions:
            sol = Chromo(Value)
            Listofsolutions.append(Value)
            pop.append(sol.Fintess_Calc(Data,1,pop))


  
def mutation(sol,mu):
    gen2change=int(math.ceil(mu*Varsize))
    Value=sol.value
    if np.random.rand()>=0.25:
        for  _ in xrange(gen2change):
            rep=np.random.choice(Varsize,2,replace=False)
            rep.sort()
            
            Value= Value[:rep[0]]+[Value[rep[1]]] +Value[rep[0]+1:rep[1]]  +[Value[rep[0]]]+Value[rep[1]+1:]
        
    else:
        Value=sol.Random_bin_no_change(Data)

    return Chromo(Value)
    
    
def crossover(Dad_value, Mom_value):


    child_Val_1=[]
    child_Val_2=[]
    if random.random()<=0.5:
        """ one-point crossover """
        x=np.random.randint(0,Varsize) 
        child_Val_1=Dad_value[:x]+Mom_value[x:]
        child_Val_2=Mom_value[:x]+Dad_value[x:]
    else:
        """ two-point crossover """   
        (x,y)=np.random.choice(Varsize,2,False)        
        if x > y: x,y = y,x
        child_Val_1 = Mom_value[:x]+Dad_value[x:y]+Mom_value[y:]
        child_Val_2 = Dad_value[:x]+Mom_value[x:y]+Dad_value[y:]
   
        
        
    child1 = Chromo(child_Val_1)
    child2 = Chromo(child_Val_2)
    
    return (child1, child2)
    

 
def evolve(nPop,mutation_rate,crossover_rate,mu):

    global pop
    
    sp=1.8 # parameter in parents selection
    parents_length = int(nPop/1.5) # number of the parents
    pv=[]
    ############################# parents ################################
    # calculate the parents selection probablity
    for r,individual in enumerate(pop):
        rank=float(nPop-r-1)
        pv.append(round((2-sp)/nPop+2*rank*(sp-1)/(nPop*(nPop-1)),5))
    pv = np.array(pv)    
    pv /= pv.sum()
    # selecting the parents 
    parents= roulette_wheel_pop(pop, pv, parents_length)
    #################### mutate some individuals###########################
    Mutation_number=math.ceil(mutation_rate*nPop)
    counter=1
    Mutants=[]
    Mut_inner_counter=0
    while counter<=Mutation_number and Mut_inner_counter<=20*Mutation_number:
        #Select the individual
        individual=parents[np.random.randint(len(parents))]
        individual=mutation(individual,mu)

        if individual.value not in Listofsolutions:
            Listofsolutions.append(individual.value)
            Mutants.append(individual)
            counter+=1
            Mut_inner_counter=0
        else:
            Mut_inner_counter+=1

    ########################### crossover ##############################
    Crossover_number = int(crossover_rate*nPop)
    crosscounter=0
    children = []
    while len(children) <= Crossover_number and crosscounter<=2*Crossover_number:
        crosscounter+=1
        
        (male,female)=np.random.choice(parents_length,2,False)
        
        (child1, child2)=crossover( parents[male].value[:] , parents[female].value[:] )
  
        if child1.value not in Listofsolutions : 
            Listofsolutions.append(child1.value)
            children.append(child1)
            
        if child2.value not in Listofsolutions :
            Listofsolutions.append(child2.value)
            children.append(child2)
            
           
   ###################################################################         
    children=[x.Fintess_Calc(Data,0,pop) for x in children]
    Mutants=[x.Fintess_Calc(Data,0,pop) for x in Mutants]
    
    # create the pool
    pool=pop[:parents_length]
    pool.extend(children)
    pool.extend(Mutants)
    # evaluate the pool
    pool=sorted(pool,key=lambda pool:pool.Fitness_Value,reverse=False)
    # truncate the pool and create the new generation
    pop=pool[0:nPop]
    
    return
  


def roulette_wheel_pop(pop, p, number):    
    chosen=np.random.choice(len(pop),number,False,p)
    chosen = [pop[a] for a in chosen]
    return chosen 



def GA(Data):
    global pop
    MaxRunTime=300
    nPop=3*int(Varsize )  #Population Size 
    MaxIt=100  # Maximum Number of Iterations
    Max_noimprove=35 # Maximum number of iterations without improvement before termination
    crossover_rate=0.7 # 
    mutation_rate=0.3 # Mutation Percentage
    mu=0.2  # Mutation Rate 
    
    start=time.time()
    
    iterationNO=1
    initialpop(nPop,Data) # generate intial solution 
    pop=sorted(pop,key=lambda x:x.Fitness_Value,reverse=False) # sort population base on fitness value
    
    current_bestsol=pop[0]
    noimprove=0
    while iterationNO<=MaxIt and noimprove<Max_noimprove and time.time()<=start+MaxRunTime:
        evolve(nPop,mutation_rate,crossover_rate,mu)
        last_bestsol=current_bestsol
        current_bestsol=pop[0]
        if current_bestsol.Fitness_Value==last_bestsol.Fitness_Value:
            noimprove+=1
        else: noimprove=0
        print("iteration%s" % iterationNO + "--#searched solutions= %s" % len(Listofsolutions))
        print("Best objective: %s" %(current_bestsol.Fitness_Value) )
        
        #AvgNOBins=np.average([len(a.Bins) for a in pop]    )   
        #print("Average number of bins = %s" % AvgNOBins)
        #AvgEL=np.average([a.early_lateness for a in pop]    )  
        #print("Average lateness earliness measure = %s" % AvgEL)
        
        
        
        iterationNO+=1
    
    return current_bestsol


global pop
global Listofsolutions
Listofsolutions=[]
pop=[]  

### read Execl Data file ####
### First number in "Input(29,6)"  is the number of item to read from Data2 excel file
### Second number is the number of planning period days
Data=Input(13,3)
Varsize=Data.N

# Run the Genetic algorithm 
start=time.time()

Best_Sol=GA(Data)

Runtime=time.time()-start

# display the results
print("######################## Results #############################")
print("Total cost of printing all bins: %s" %Best_Sol.total_cost)
print("Total lateness and earliness: %s" %Best_Sol.early_lateness)
print("Number of Bins: %s" %len(Best_Sol.Bins))
print ("Algorithm Run Time: %s" %Runtime)
for i,b in enumerate(pop[0].Bins):
    Draw_the_Bin(Data,b)
    print("Printing Quqntity: %s" %b.quantity)



