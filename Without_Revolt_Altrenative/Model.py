from gurobipy import *
from itertools import combinations 
import pandas as pd
import numpy as np
from Input import Input
from copy import copy
from Chromo import stack
from Chromo import level
from Chromo import Bin
from Draw_the_Bin import Draw_the_Bin



def Solve():

    N=Data.N
    T=Data.T  
    d=[Data.items[i].d for i in range(N)]
    Qmax=max([Data.items[a].q for a in range(N)])
    NN=tuplelist( [(i,j) for i in range(N) for j in range(N) if j>=i ] )
    
    
    
    items=copy(Data.items.values())
    for it in items:
        if it.two_side==0:
            it.q=it.q/2
    
    alphaNN=[]
    
    for i,j in NN:
        if (Data.items[i].w==Data.items[j].w and Data.items[i].h+Data.items[j].h<=Data.H) or i==j:
            alphaNN.append((i,j))
            
    alphaNN=tuplelist( alphaNN )
    
    
    LKJI=[]
    bindex=[]
    for (j,i) in alphaNN:
        for k in range(N):
            for l in range(k+1):
                LKJI.append((l,k,j,i))
                if (i,l) not in bindex: bindex.append((i,l))
    bindex2=[]
    for i,l in bindex:
        for t in range(T):
            bindex2.append((i,l,t))
    
    deltaindex=[]
    alphadelta=[]
    for l,i in combinations(range(N),2):
        for j in range(l,i):
            deltaindex.append((l,i,j))
            if (j,i) in alphaNN: alphadelta.append( (l,i,j) )
    
    MIP=Model("packing")
    alpha=MIP.addVars(alphaNN,name="alpha",vtype=GRB.BINARY)
    beta=MIP.addVars(N,N,name="beta",vtype=GRB.BINARY)
    gamma=MIP.addVars(NN,name="gamma",vtype=GRB.BINARY)
    delta=MIP.addVars(deltaindex,name="delta",vtype=GRB.BINARY)    
    
    
    a=MIP.addVars(LKJI,name="a",vtype=GRB.BINARY)
    b=MIP.addVars(bindex2,name="b",vtype=GRB.BINARY)    
    Q=MIP.addVars(N,name="Q",vtype=GRB.INTEGER)
    #D=MIP.addVars(N,name="D",vtype=GRB.INTEGER)
    z=MIP.addVars(N,T,name="z",vtype=GRB.INTEGER)
    x=MIP.addVars(N,T,name="x",vtype=GRB.BINARY)
    tp=MIP.addVars(N,name="tp",vtype=GRB.INTEGER)    
    tn=MIP.addVars(N,name="tn",vtype=GRB.INTEGER)
    
    MIP.update()

    MIP.addConstrs(quicksum(alpha.select('*',i) )==1 for i in range(N))
    MIP.addConstrs(quicksum(alpha.select(j,'*')  ) <= (N-j-1)*alpha[j,j] for j in range(N-1) )
    
    MIP.addConstrs(quicksum(beta[k,j] for k in range(N) )==alpha[j,j] for j in range(N) )
    
    MIP.addConstrs(quicksum(Data.items[i].h*alpha[j,i] for _,i in alphaNN.select(j,'*') )
        <= quicksum(Data.items[i].h*alpha[k,i] for _,i in alphaNN.select(k,'*')) + (Data.H+1)*(1-beta[k,j]) -0.003
        for k in range(1,N) for j in range(k)   )
        
        
    MIP.addConstrs(quicksum(Data.items[i].h*alpha[j,i] for _,i in alphaNN.select(j,'*') )
        <= quicksum(Data.items[i].h*alpha[k,i] for _,i in alphaNN.select(k,'*')) + Data.H*(1-beta[k,j])
        for k in range(N-1) for j in range(k,N)   ) 
    
    MIP.addConstrs(quicksum(Data.items[j].w*beta[k,j] for j in range(N))<=Data.W*beta[k,k]/2 for k in range(N))
    
    MIP.addConstrs(quicksum( gamma.select('*',k) )==beta[k,k] for k in range(N))
    
    MIP.addConstrs(quicksum(Data.items[i].h*gamma[l,i] for i in range(l,N))
                    +quicksum(Data.items[i].h*delta[l,i,j] for i in range(l+1,N) for j in range(l,i))
                    <= Data.H*gamma[l,l] for l in range(N-1))
    
    
    MIP.addConstrs(alpha[j,i]+gamma[l,j]-1<=delta[l,i,j] for l,i,j in alphadelta)
    MIP.addConstrs((alpha[j,i]+gamma[l,j])/2 >= delta[l,i,j] for l,i,j in alphadelta)
    
    MIP.addConstrs(quicksum(gamma[l,k] for k in range(l+1,N) ) <= (N-l)*gamma[l,l] for l in range(N-1) )
    
    
    
    #new constraints
    MIP.addConstrs(quicksum(x[l,t]for t in range(T))==gamma[l,l] for l in range(N))

    #Item-Bin relationship     
    MIP.addConstrs(alpha[j,i]+beta[k,j]+gamma[l,k] <= 2 + a[l,k,j,i] for l,k,j,i in LKJI )
    MIP.addConstrs(alpha[j,i]+beta[k,j]+gamma[l,k] >= 3*a[l,k,j,i] for l,k,j,i in LKJI )
    #bin quantity
    MIP.addConstrs(quicksum(a.select(l,'*','*',i))*items[i].q <= Q[l] for i in range(N) for l in range(N) )
    #item processing date
    MIP.addConstrs(2*b[i,l,t]<=quicksum(a.select(l,'*','*',i))+x[l,t] for i,l,t in bindex2 )    
    MIP.addConstrs(1+b[i,l,t]>=quicksum(a.select(l,'*','*',i))+x[l,t] for i,l,t in bindex2 )
    #Production capacity linearization  
    MIP.addConstrs(z[l,t]<=Q[l]  for t in range(T) for l in range(N) )
    MIP.addConstrs(z[l,t]<= Qmax * x[l,t]   for t in range(T) for l in range(N) )
    MIP.addConstrs(z[l,t]>=Q[l]-Qmax*(1-x[l,t])  for t in range(T) for l in range(N) )
    # production capacity
    MIP.addConstrs(quicksum(z.select('*',t)) <= Data.proCap for t in range(T) )
    #tardiness calculation
    MIP.addConstrs( tp[i] >= quicksum(t*b[i,l,t] for l in range(N) for t in range(T)) -d[i] for i in range(N) )
    #earliness calculation
    MIP.addConstrs( tn[i] >= d[i] - quicksum(t*b[i,l,t] for l in range(N) for t in range(T)) for i in range(N) )
    #tardiness limit
    MIP.addConstrs(tp[i] <= 1 for i in range(N) )
    #earliness limit
    MIP.addConstrs(tn[i] <=2 for i in range (N) )
    
    
    MIP.setObjective(quicksum(Q[l]*0.1325808/2 for l in range(N)) +60*quicksum(gamma[l,l] for l in range(N)) , GRB.MINIMIZE )

   
    MIP.update()

    MIP.optimize()
    MIP.write('out.lp')
    
    print(MIP.Runtime)
    if MIP.status==2:
        MIP.write('soll.sol')
        print(MIP.objVal)
        #Dv=MIP.getAttr('x',D)
        Xv=MIP.getAttr('x',x)
        alphav=MIP.getAttr('x',alpha)
        betav=MIP.getAttr('x',beta)
        gammav=MIP.getAttr('x',gamma)
        Qv=MIP.getAttr('x',Q)

        
        alphavv=[(i,j) for (i,j) in alphav.keys() if alphav[(i,j)]==1]
        
        print("Item to stack assignment")
        print(alphavv)
        betavv=[(i,j) for (i,j) in betav.keys() if betav[(i,j)]==1]
        print("Stack to stripe assignment")
        print(betavv)
        gammavv=[(i,j) for (i,j) in gammav.keys() if gammav[(i,j)]==1]
        print("Stripe to bin assignment")
        print(gammavv)
        Xvv=[(i,j) for (i,j) in Xv.keys() if Xv[(i,j)]==1]
        print("Bin to day assignment")
        print(Xvv)
        Qvv=[(i,Q[i]) for i in Qv.keys() if Qv[i]!=0]
        print(Qvv)
    return (alphavv,betavv,gammavv)            




def Draw_Bins(Data,alphavv,betavv,gammavv):
    NN=Data.N
    
    Bin_list=[0 for _ in range(NN)]
    stack_list=[0 for _ in range(NN)]
    level_list=[0 for _ in range(NN)]
    
    alpha=np.array(alphavv)
    beta=np.array(betavv)
    gamma= np.array(gammavv)
    
    for st in range(NN):
        a,b=np.where(alpha==st)
        stock=[alpha[a[i]][1] for i in range(len(a)) if b[i]==0]
        if len(stock)!=0:
            st_h=sum([Data.items[i].h for i in stock])
            stack_list[st]=stack(st_h,Data.items[stock[0]])
            for it in stock[1:]:
                stack_list[st].add( Data.items[it] )
        
    for le in range(NN):
        a,b=np.where(beta==le)
        stock_in_le_inx=[beta[a[i]][1] for i in range(len(a)) if b[i]==0]
        
        if len(stock_in_le_inx)!=0:
            stock_in_le=np.array(stack_list)[stock_in_le_inx]
            level_h=max([st.level_h for st in stock_in_le])
            level_list[le]=level(level_h,Data.W/2)
            for st in stock_in_le: 
                for it in st.items:
                    indi=level_list[le].add(it)    
                    if indi==0:
                        print("Whaaat!!")
    
    for bi in range(NN):
        a,b=np.where(gamma==bi)
        level_in_bin=[gamma[a[i]][1] for i in range(len(a)) if b[i]==0]
        if len(level_in_bin)!=0:
           Bin_list[bi]=Bin(level_list[level_in_bin[0]],Data,1) 
           for le in level_in_bin[1:]:
               Bin_list[bi].add(level_list[le])
    
    for i,bi in enumerate(Bin_list):
        if bi!=0:
            Draw_the_Bin(Data,bi)
            print("Printing Quqntity: %s" %bi.quantity)
            print("Revolting Bin= %d" % bi.Revolta)  
    return     
    
    
    
    

Data=Input(10,3)
alphavv,betavv,gammavv=Solve()
Draw_Bins(Data,alphavv,betavv,gammavv)


