from Draw_the_Bin import Draw_the_Bin
from Input import Input
from Chromo import Chromo, inbalance_measure,ItemSet
from lowerbound import lower_bound
import random
import os
import time
import math
import numpy as np
import pickle as Pick

def initialpop(nPop,Data):
   global pop
   # First calculate the minimum bin number needed to allocate all itmes 
   Bmin=Data.MinBinNo + math.ceil(2.5*Data.MinBinNo)
   
  
   ### initial heurestics ##
   Value=[0]*Data.N
   items = ItemSet(Data, list(Data.items.values()) )
   items.Qsort()
   tempBinIt=ItemSet(Data,[ items[0]])
   BN=0
   while  len(items.list)!=0 : 
       it=items.list[0] # select the unassigned item with the largnest quantity 
       if inbalance_measure(tempBinIt.list+[it]) <= 0.3: #if adding this item to current bin is not going to change the IM
          del items.list[0]
          tempBinIt.add([it])
       else:
          for i in tempBinIt.list: # set value of all items in current bin 
              Value[i.ID]=BN
          BN += 1 # create a new bin
          tempBinIt=ItemSet(Data, [items[0]])
          del items.list[0]
    
   (Num_Bin,Value)=Calc_Bin_No(Value,Data) 
   Revolting=list(np.random.randint(2, size=(1, Num_Bin))[0]) # randomly decide on revolting option
   Revolting=[a*b for a,b in zip(Revolting,Bin_can_revolt(Data,Value) ) ] # check if there is two-sided option
   sol = Chromo(Value,Revolting)
   pop.append(sol.Fitness_Calc(Data,pop,iterationNO,Maxit))
   
   
   
   while len(pop) < nPop:
        
        # Generating Value genes
        Value=[]
        for _ in range(Varsize):
            Value.append(random.randint(0,Bmin))
        
        # Update the "Value" with ordered bin numbers         
        (Num_Bin,Value)=Calc_Bin_No(Value,Data)
        
        ### Genearting revolting genes
        Revolting=list(np.random.randint(2, size=(1, Num_Bin))[0])
        ### CHANGE ###
        Revolting=[a*b for a,b in zip(Revolting,Bin_can_revolt(Data,Value) ) ] 
      
        if Is_solution_new(Value,Revolting):
            sol = Chromo(Value,Revolting)
            pop.append(sol.Fitness_Calc(Data,pop,iterationNO,Maxit))
            
def mutation(sol,mu):
    ## Mutation for value ##
    gen2change=int(math.ceil(mu*Varsize))
    Value=sol.value
    (Bin_no,Value)=Calc_Bin_No( Value , Data )
    rnd_value = random.random()
    if rnd_value <= 0.3 :
        # Swap two items between their bins
        for  _ in range(gen2change):
            rep=np.random.choice(Varsize,2,replace=False)
            rep.sort()
            Value= Value[:rep[0]]+[Value[rep[1]]] +Value[rep[0]+1:rep[1]]  +[Value[rep[0]]]+Value[rep[1]+1:]
            
    elif rnd_value <= 0.7 :
        ## change one item bin ##
        for  _ in range(gen2change):
            item2change = np.random.choice(Varsize,replace=False)
            newbin = np.random.randint( Bin_no )
            Value[ item2change ] = newbin
    else:
        
        Value=sol.Random_bin_no_change(Data)
        
   ## Mutation for revolting ###
    Revolting=sol.Revolting 
    (Bin_no,Value)=Calc_Bin_No( Value , Data )
    
    if Bin_no<=len(Revolting):
        Revolting=Revolting[:Bin_no]
    else:
        Revolting=Revolting+ list(np.random.randint(2, size=(1,Bin_no-len(Revolting) ))[0])
    gen2change=int(math.ceil(mu*Bin_no))
    rep=np.random.choice(Bin_no,gen2change,replace=False)
    for i in rep:
        Revolting[i]=1-Revolting[i]
    ### CHANGE ###    
    Revolting=[a*b for a,b in zip(Revolting,Bin_can_revolt(Data,Value) ) ] 
    
    return Chromo(Value,Revolting)
    
    
def crossover(DadSol, MomSol):
    Varsize=Data.N
    Dad_value=DadSol.value
    Mom_value=MomSol.value
    child_Val=[[],[]]
    
    if random.random()<=0.5:
        """ one-point crossover """
        x=np.random.randint(0,Varsize) 
        child_Val[0]=Dad_value[:x]+Mom_value[x:]
        child_Val[1]=Mom_value[:x]+Dad_value[x:]
    else:
        """ two-point crossover """   
        (x,y)=np.random.choice(Varsize,2,False)        
        if x > y: x,y = y,x
        child_Val[0] = Mom_value[:x]+Dad_value[x:y]+Mom_value[y:]
        child_Val[1]= Dad_value[:x]+Mom_value[x:y]+Dad_value[y:]
   
     
    
    ### crossover for revolting ###
    Dad_Revolting=DadSol.Revolting
    Mom_Revolting=MomSol.Revolting
    Dad_Bin_no=len(Dad_Revolting)
    Mom_Bin_no=len(Mom_Revolting)
    Varsize=min(Dad_Bin_no,Mom_Bin_no)
    Child_Bin_no=[[],[]]
    child_Rev=[[],[]]
    child=[[],[]]
    # calculate childs bin number
    
    """ Revolting part crossover  """
    for i in [0,1]:   
        (Child_Bin_no[i],child_Val[i])=Calc_Bin_No( child_Val[i] , Data )
        for _ in range(Child_Bin_no[i]):
            if random.random()<=0.5:
                x=np.random.randint(0,Dad_Bin_no)
                child_Rev[i].append(Dad_Revolting[x])
            else:
                x=np.random.randint(0,Mom_Bin_no)
                child_Rev[i].append(Mom_Revolting[x])

        ### CHANGE ###
        child_Rev[i]=[a*b for a,b in zip(child_Rev[i],Bin_can_revolt(Data,child_Val[i]) )]
        
        child[i] = Chromo(child_Val[i],child_Rev[i])

    return child
    
def Bin_can_revolt(Data,Value):
    ### CHANGE ###
    Num_Bin = max(Value)
    Revolta=np.zeros((Num_Bin+1))
    for b in range(Num_Bin+1):
        Bin_Items=np.where(np.array(Value)==b)[0]
        for it in Bin_Items:
            if Data.items[it].two_side==1:
               Revolta[b]=1
               break
           
    return Revolta

def Is_solution_new(Value,Revolting):
    indicator=0
    a=  Value+Revolting
    if a not in Listofsolutions:
        Listofsolutions.append(a)
        indicator=1
    return indicator
    
def Calc_Bin_No(Value,Data):
    Item_Bin=[]
    i=0            
    Num_Bin = max(Value)
    Value=np.array(Value)
    for b in range(Num_Bin+1):
        if b in Value:
            Item_Bin.append( np.array(list(Data.items.values()))[np.where(Value == b)] )
            Value[[it.ID for it in Item_Bin[-1]]]=i
            i+=1   
    Value=list(Value)
    Num_Bin=len(Item_Bin)
    return (Num_Bin,Value)
 
def evolve(nPop,mutation_rate,crossover_rate,mu):

    global pop
    
    sp=1.8 # parameter in parents selection
    parents_length = int(nPop/2) # number of the parents
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

        if Is_solution_new(individual.value,individual.Revolting):
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
        
        (child1, child2)=crossover( parents[male] , parents[female] )
  
        if Is_solution_new(child1.value,child1.Revolting):
            children.append(child1)
            
        if Is_solution_new(child2.value,child2.Revolting) :
            children.append(child2)
            
           
   ###################################################################         
    children=[x.Fitness_Calc(Data,pop,iterationNO,Maxit) for x in children]
    Mutants=[x.Fitness_Calc(Data,pop,iterationNO,Maxit) for x in Mutants]
    
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
    global Maxit
    global iterationNO
    
    MaxRunTime=3600
    nPop=20+3*int(Varsize )  #Population Size 
    Maxit=100 + int(Varsize )  ### CHANGE ### # Maximum Number of Iterations
    Max_noimprove=int(Maxit*0.2) ### CHANGE ###  Maximum number of iterations without improvement before termination
    crossover_rate=0.7 # 
    mutation_rate=0.3 # Mutation Percentage
    mu=0.25  # Mutation Rate 
    
    start=time.time()
    
    iterationNO=1
    initialpop(nPop,Data) # generate intial solution 
    pop=sorted(pop,key=lambda x:x.Fitness_Value,reverse=False) # sort population base on fitness value
    
    current_bestsol=pop[0]
    noimprove=0
    while iterationNO<=Maxit and noimprove<Max_noimprove and time.time()<=start+MaxRunTime:
        evolve(nPop,mutation_rate,crossover_rate,mu)
        last_bestsol=current_bestsol
        current_bestsol=pop[0]
        if current_bestsol.Fitness_Value==last_bestsol.Fitness_Value:
            noimprove+=1
        else: noimprove=0
        #print("iteration%s" % iterationNO + "--#searched solutions= %s" % len(Listofsolutions))
        #print("Best objective: %s" %(current_bestsol.Fitness_Value) )
        
        #AvgNOBins=np.average([len(a.Bins) for a in pop]    )   
        #print("Average number of bins = %s" % AvgNOBins)
        #AvgEL=np.average([a.early_lateness for a in pop]    )  
        #print("Average lateness earliness measure = %s" % AvgEL)
        
        
        
        iterationNO+=1
    
    return current_bestsol




def read_object(FileName,folder):
    WD = os.getcwd()
    if folder == "Input":
        path = WD + f"/Model/{folder}/{FileName}"
    else:
        path = WD + "\Model/%s/%s_ModelSol" %(folder,FileName)
    
    with open(path , 'rb') as input:
        obj= Pick.load(input, encoding="latin1")
    return  obj 


if __name__ == "__main__":

    global pop
    global Listofsolutions


    ### read Execl Data file ####
    ### First number in "Input(29,6)"  is the number of item to read from Data2 excel file
    ### Second number is the number of planning period days
    #Data=Input( 10  , 2 )
    N = 12
    T=5
    for _ in range(1):
        for name in ['ST']:
            print(name)
            for rep in range(4,6):

                FileName="Data_%d_%d_%d_%s" %(N,T,rep,name)
                Data= read_object(FileName,"Input")

                Varsize=Data.N
                Listofsolutions=[]
                pop=[]
                # Run the Genetic algorithm
                start=time.time()

                Best_Sol=GA(Data)

                Runtime=time.time()-start
                print('%d_%d %s %d %s  %s'%( N,T,name,rep , str(round(Best_Sol.total_cost,1)), str(round(Runtime,1))) )
                #print('%s  %s'%( str(round(Best_Sol.total_cost,1)), str(round(Runtime,1))) )
                # display the results
                #print("############### Results ################")
                #
                #print("Lower Bound: %s" %lower_bound(Data))
                #print("Total cost of printing all bins: %s" %Best_Sol.total_cost)
                #print("Total lateness and earliness: %s" %Best_Sol.early_lateness)
                #print("Number of Bins: %s" %len(Best_Sol.Bins))
                #print ("Algorithm Run Time: %s" %Runtime)

                for i,b in enumerate(pop[0].Bins):
                    Draw_the_Bin(Data,b)
                    print("Printing Quqntity: %s" %b.quantity)
                    #print("Items quantity: ", [it.q for it in b.items] )
                    print("IM: %s" %inbalance_measure(b.items) )
                    #print("Revolting Bin= %d" % b.Revolta)



