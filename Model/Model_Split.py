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
    Qmin=min([Data.items[a].q for a in range(N)])
    NN=tuplelist( [(i,j) for i in range(N) for j in range(N) if j>=i ] )
    M = [int((Data.H*Data.W)/Data.items[i].area) for i in range(N)  ]
    Mmax=int(max(M))
    Mmin=int(min(M))
    
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
    
    
    # breaking points in X^2 function
    rplow=1/2  
    rphigh=(Mmax+Qmax)/2
    BP_rp= np.linspace(0, rphigh, 150 ,dtype=int)
    f_rp=[ it**2 for it in BP_rp ]    
    rnlow=(1-Qmax)/2
    rnhigh=(Mmax-0)/2
    BP_rn=np.linspace(rnlow, rnhigh, 150 ,dtype=int)
    #BP_rn= np.array( [rnlow , -1500, -1000, -500, -100, 0, 10, rnhigh] )
    f_rn=[ it**2 for it in BP_rn ]  
    
    MIP=Model("packing")   
    alpha = MIP.addVars(alphaNN,name="alpha",vtype=GRB.BINARY)
    o = MIP.addVars(alpha,name="o",vtype=GRB.INTEGER) # number of item i in stock j
    y = MIP.addVars(N,name="y",lb=1, ub=M, vtype=GRB.INTEGER) # number of item i splists
    
    beta = MIP.addVars(N,N,name="beta",vtype=GRB.BINARY)
    gamma = MIP.addVars(NN,name="gamma",vtype=GRB.BINARY)
    delta = MIP.addVars(deltaindex,name="delta",vtype=GRB.INTEGER)    
    
    
    f=MIP.addVars(gamma,name="f")
    a=MIP.addVars(LKJI,name="a",vtype=GRB.BINARY)
    b=MIP.addVars(bindex2,name="b",vtype=GRB.BINARY)    
    Q=MIP.addVars(N,name="Q",lb=0)
    
    rp=MIP.addVars(N,N,name="rp")
    rn=MIP.addVars(N,N,name="rn",lb=-GRB.INFINITY)
    lambdaP=MIP.addVars(N,N,len(BP_rp),name="lambdaP",lb=0,vtype=GRB.CONTINUOUS)
    lambdaN=MIP.addVars(N,N,len(BP_rn),name="lambdaN",lb=0,vtype=GRB.CONTINUOUS)



    z=MIP.addVars(N,T,name="z",vtype=GRB.INTEGER)
    x=MIP.addVars(N,T,name="x",vtype=GRB.BINARY)
    tp=MIP.addVars(N,name="tp",vtype=GRB.INTEGER)    
    tn=MIP.addVars(N,name="tn",vtype=GRB.INTEGER)
    
    MIP.update()
    #c1
    MIP.addConstrs(quicksum( o.select('*',i) ) == y[i] for i in range(N))
    #c2
    MIP.addConstrs( o[j,i] <= M[i]*alpha[j,i] for j,i in alphaNN )
    #c3
    MIP.addConstrs(quicksum(alpha.select(j,'*') ) <= (N-j-1)*alpha[j,j] for j in range(N-1) )
    #c4
    MIP.addConstrs(quicksum(beta[k,j] for k in range(N) ) == alpha[j,j] for j in range(N) )
    #c5
    MIP.addConstrs(quicksum(Data.items[i].h*o[j,i] for _,i in alphaNN.select(j,'*') )
        <= quicksum(Data.items[i].h*o[k,i] for _,i in alphaNN.select(k,'*')) + (Data.H+1)*(1-beta[k,j]) -0.003
        for k in range(1,N) for j in range(k)   )
    
    #c6        
    MIP.addConstrs(quicksum(Data.items[i].h*o[j,i] for _,i in alphaNN.select(j,'*') )
        <= quicksum(Data.items[i].h*o[k,i] for _,i in alphaNN.select(k,'*')) + Data.H*(1-beta[k,j])
        for k in range(N-1) for j in range(k,N)   ) 
    #c7
    MIP.addConstrs(quicksum(Data.items[j].w*beta[k,j] for j in range(N)) <= Data.W*beta[k,k]/2 for k in range(N))
    #c8
    MIP.addConstrs(quicksum( gamma.select('*',k) ) == beta[k,k] for k in range(N))
    
    #MIP.addConstrs(f[l,i]<=o[i,i] for l,i in f.keys())
    #MIP.addConstrs(f[l,i] <= M[i]*gamma[l,i] for l,i in gamma.keys())
    #MIP.addConstrs(o[i,i]-M[i]*(1-gamma[l,i]) <= f[l,i] for l,i in gamma.keys() )
    # f[l,i]= delta[l,i,i]
    
    
    #c9-c11
    MIP.addConstrs(delta[l,i,j] <= o[j,i] for l,i,j in alphadelta)
    MIP.addConstrs(delta[l,i,j] <= M[i]*gamma[l,j] for l,i,j in alphadelta)
    MIP.addConstrs(o[j,i]-M[i]*(1-gamma[l,j]) <= delta[l,i,j] for l,i,j in alphadelta)    
    
    #c12
    MIP.addConstrs(quicksum(Data.items[i].h*delta[l,i,i] for i in range(l,N))
                    +quicksum(Data.items[i].h*delta[l,i,j] for i in range(l+1,N) for j in range(l,i))
                    <= Data.H*gamma[l,l] for l in range(N-1))
    #c13
    MIP.addConstrs(quicksum(gamma[l,k] for k in range(l+1,N) ) <= (N-l)*gamma[l,l] for l in range(N-1) )
    
    #new constraints
    #c14
    MIP.addConstrs(quicksum(x[l,t]for t in range(T))==gamma[l,l] for l in range(N))

    #Item-Bin relationship     # c15
    MIP.addConstrs(alpha[j,i]+beta[k,j]+gamma[l,k] <= 2 + a[l,k,j,i] for l,k,j,i in LKJI )
    MIP.addConstrs(alpha[j,i]+beta[k,j]+gamma[l,k] >= 3*a[l,k,j,i] for l,k,j,i in LKJI )
    # All splits should be in one bin     
    #MIP.addConstrs(quicksum(a.select('*','*','*',i)) ==1 for i in range(N) ) 
    #bin quantity
    #for i in range(N):
    #    for l in range(N):
    #        MIP.addQConstr(quicksum(a.select(l,'*','*',i))*items[i].q <= y[i]*Q[l])    
    #C16
    MIP.addConstrs(quicksum(a.select(l,'*','*',i))*items[i].q <= 
                    quicksum( lambdaP[i,l,k]*f_rp[k] for k in range(len(BP_rp)) ) 
                    - quicksum( lambdaN[i,l,k]*f_rn[k] for k in range(len(BP_rn)) )
                    for i in range(N) for l in range(N) )
    
    # linearization of y*Q    
    MIP.addConstrs(rp[i,l] == 0.5* (y[i] + Q[l]) for i in range(N) for l in range(N) )
    MIP.addConstrs(rn[i,l] == 0.5* (y[i] - Q[l]) for i in range(N) for l in range(N) )
    
    MIP.addConstrs(rp[i,l] == quicksum(lambdaP[i,l,k]*BP_rp[k] for k in range( len(BP_rp) )) for i in range(N) for l in range(N) )
    MIP.addConstrs(quicksum(lambdaP.select(i,l,'*') ) == 1 for i in range(N) for l in range(N))   
    for i in range(N):
        for l in range(N):
            MIP.addSOS(GRB.SOS_TYPE2, lambdaP.select(i,l,'*'))    
    
    
    MIP.addConstrs(rn[i,l] == quicksum(lambdaN[i,l,k]*BP_rn[k] for k in range(len(BP_rn)) ) for i in range(N) for l in range(N) )
    MIP.addConstrs(quicksum(lambdaN.select(i,l,'*')) == 1 for i in range(N) for l in range(N))   
    for i in range(N):
        for l in range(N):
            MIP.addSOS(GRB.SOS_TYPE2, lambdaN.select(i,l,'*')) 
    
    #item processing date
    
    #c17
    MIP.addConstrs(2*b[i,l,t]<=quicksum(a.select(l,'*','*',i))+x[l,t] for i,l,t in bindex2 )    
    MIP.addConstrs(1+b[i,l,t]>=quicksum(a.select(l,'*','*',i))+x[l,t] for i,l,t in bindex2 )
    #Production capacity linearization  
    MIP.addConstrs(z[l,t]<=Q[l]  for t in range(T) for l in range(N) )
    MIP.addConstrs(z[l,t]<= Qmax * x[l,t]   for t in range(T) for l in range(N) )
    MIP.addConstrs(z[l,t]>=Q[l]-Qmax*(1-x[l,t])  for t in range(T) for l in range(N) )
    # production capacity
    #c18
    MIP.addConstrs(quicksum(z.select('*',t)) <= Data.proCap for t in range(T) )
    #tardiness calculation
    #c19
    MIP.addConstrs( tp[i] >= quicksum(t*b[i,l,t] for l in range(N) for t in range(T)) -d[i] for i in range(N) )
    #earliness calculation
    #c20
    MIP.addConstrs( tn[i] >= d[i] - quicksum(t*b[i,l,t] for l in range(N) for t in range(T)) for i in range(N) )
    #tardiness limit
    #c21
    MIP.addConstrs(tp[i] <= 1 for i in range(N) )
    #earliness limit
    #c22
    MIP.addConstrs(tn[i] <=2 for i in range (N) )
    
    
    MIP.setObjective(quicksum(Q[l]*0.1325808/2 for l in range(N)) +60*quicksum(gamma[l,l] for l in range(N)) , GRB.MINIMIZE )
    
   
    MIP.update()

    MIP.optimize()
    MIP.write('out.lp')
    
    print(MIP.Runtime)
    if MIP.status==2:
        MIP.write('soll.sol')
        print(MIP.objVal)
        Xv=MIP.getAttr('x',x)
        Yv=MIP.getAttr('x',y)
        Ov=MIP.getAttr('x',o)
        
        alphav=MIP.getAttr('x',alpha)
        betav=MIP.getAttr('x',beta)
        gammav=MIP.getAttr('x',gamma)
        deltav=MIP.getAttr('x',delta)
        
        av=MIP.getAttr('x',a)
        Qv=MIP.getAttr('x',Q)
        lambdaNv=MIP.getAttr('x',lambdaN)
        lambdaPv=MIP.getAttr('x',lambdaP)
        
        alphavv=[(i,j) for (i,j) in alphav.keys() if alphav[(i,j)]>0]
        alphavv2=[(i,o[i].X) for i in alphav.keys() if alphav[i]>0]
        print("Item to stack assignment")
        print(alphavv2)
        betavv=[(i,j) for (i,j) in betav.keys() if betav[(i,j)]==1]
        print("Stack to stripe assignment")
        print(betavv)
        gammavv=[(i,j) for (i,j) in gammav.keys() if gammav[(i,j)]==1]
        print("Stripe to bin assignment")
        print(gammavv)
        Xvv=[(i,j) for (i,j) in Xv.keys() if Xv[(i,j)]==1]
        print("Bin to day assignment")
        print(Xvv)
        Qvv=[(i,Q[i].X) for i in Qv.keys() if Qv[i]!=0]
        print(Qvv)
        print ('Real rp: ', rp[0,0].X)
        print('Real rp^2: ' , rp[0,0].X**2)
        print('Approximate rp : ' , sum( [lambdaP[0,0,k].X*BP_rp[k] for k in range( len(BP_rp) )] ) )
        print('Approximate rp^2 : ' ,sum( [lambdaP[0,0,k].X*f_rp[k] for k in range(len(BP_rp))] ))
        
        
        print ('Real rn: ', rn[0,0].X)        
        print('Real rn^2: ' ,rn[0,0].X**2)
        print('Approximate rn : ' , sum( [lambdaN[0,0,k].X*BP_rn[k] for k in range( len(BP_rn) )] ) )
        print(sum( [lambdaN[0,0,k].X*f_rn[k] for k in range(len(BP_rn))] ))
        
                
        #deltav=[(i,j,k) for (i,j,k) in deltav.keys() if deltav[(i,j,k)]==1]
        #print(deltav)
        #avv=[(i,j,k,r) for (i,j,k,r) in av.keys() if av[(i,j,k,r)]==1]
        #print (avv)
    return (alphavv,betavv,gammavv,Yv)            




def Draw_Bins(Data,alphavv,betavv,gammavv,Yv):
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
    
    
    
    

Data=Input(5,1)
alphavv,betavv,gammavv, Yv=Solve()
Draw_Bins(Data,alphavv,betavv,gammavv,Yv)


