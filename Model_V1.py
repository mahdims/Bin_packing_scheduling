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
import pickle as Pick
import time
from output import Output
import os


def save_object(obj, filename):
    WD = os.getcwd()
    with open(WD + f"/Model/Output/{filename}_ModelSol", 'wb') as out:  # Overwrites any existing file.
        Pick.dump(obj, out, Pick.HIGHEST_PROTOCOL)
    out.close()


def read_object(filename, folder):
    WD = os.getcwd()
    if folder == "Input":
        path = WD + f"/Model/%s/%s" %(folder, filename)
    else:
        path = WD + f"/Model/%s/%s_ModelSol" %(folder, filename)
    
    with open(path, 'rb') as input:
        obj = Pick.load(input, encoding="latin1")
    input.close()
    return obj


def data_cb(model, where):
    if where == GRB.Callback.MIP:
        cur_obj = model.cbGet(GRB.Callback.MIP_OBJBST)
        cur_bd = model.cbGet(GRB.Callback.MIP_OBJBND)

        # Did objective value or best bound change?
        if model._obj != cur_obj or model._bd != cur_bd:
            model._obj = cur_obj
            model._bd = cur_bd
            model._data.append([time.time() - model._start, cur_obj, cur_bd])

# Build model m here


def Solve():

    N=Data.N
    T=Data.T  
    d=[Data.items[i].d for i in range(N)]
    Qmax=max([Data.items[a].q for a in range(N)])
    NN= tuplelist( [(i,j) for i in range(N) for j in range(N) if j>=i ] )
    
    items=copy(list(Data.items.values()))
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
    delta=MIP.addVars(deltaindex,name="delta",vtype=GRB.BINARY) #     
    
    
    a=MIP.addVars(LKJI,name="a",vtype=GRB.BINARY) # 
    b=MIP.addVars(bindex2,name="b",vtype=GRB.BINARY)    
    Q=MIP.addVars(N,name="Q",vtype=GRB.INTEGER)
    P=MIP.addVars(N,name="P",vtype=GRB.INTEGER) #Printing quantity after the revolting option.
    r=MIP.addVars(N,name="r",vtype=GRB.BINARY) # if bin l in revolting or not?
    w=MIP.addVars(gamma,name="w",vtype=GRB.BINARY)
    s=MIP.addVars(N,name="s",vtype=GRB.BINARY) # if two sided items are exist in the bin
    sr=MIP.addVars(N,name="sr",vtype=GRB.BINARY)
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
        
        
    # revolting an width of a revolting bin 
    MIP.addConstrs(gamma[l,k]+r[l] <= 1+w[l,k] for l,k in gamma.keys())
    MIP.addConstrs(gamma[l,k]+r[l] >= 2*w[l,k] for l,k in gamma.keys())
    MIP.addConstrs(quicksum(Data.items[j].w*beta[k,j] for j in range(N))<=Data.W*beta[k,k]-Data.W*0.5*quicksum(w.select('*',k)) for k in range(N))
    
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
    MIP.addConstrs( Q[l]*0.5-Qmax*(1-r[l]) <=P[l] for l in range(N) )   
    MIP.addConstrs( Q[l]-Qmax*r[l] <=P[l] for l in range(N) )  
    # revolting varaiable and two sided items relation  
    MIP.addConstrs(r[l] <= quicksum( quicksum(a.select(l,'*','*',i))*items[i].two_side for i in range(N) ) for l in range(N) )
    #item processing date
    MIP.addConstrs(2*b[i,l,t]<=quicksum(a.select(l,'*','*',i))+x[l,t] for i,l,t in bindex2 )    
    MIP.addConstrs(1+b[i,l,t]>=quicksum(a.select(l,'*','*',i))+x[l,t] for i,l,t in bindex2 )
    #Production capacity linearization  
    MIP.addConstrs(z[l,t]<=P[l]  for t in range(T) for l in range(N) )
    MIP.addConstrs(z[l,t]<= Qmax * x[l,t]   for t in range(T) for l in range(N) )
    MIP.addConstrs(z[l,t]>=P[l]-Qmax*(1-x[l,t])  for t in range(T) for l in range(N) )
    # production capacity
    MIP.addConstrs(quicksum(z.select('*',t)) <= Data.proCap for t in range(T) )
    #tardiness calculation
    MIP.addConstrs( tp[i] >= quicksum(t*b[i,l,t] for l in range(N) for t in range(T)) -d[i] for i in range(N) )
    #earliness calculation
    MIP.addConstrs( tn[i] >= d[i] - quicksum(t*b[i,l,t] for l in range(N) for t in range(T)) for i in range(N) )
    #tardiness limit
    MIP.addConstrs(tp[i] <= Data.Lateness for i in range(N) )
    #earliness limit
    MIP.addConstrs(tn[i] <=Data.Earliness for i in range (N) )
    
    MIP.addConstrs( s[l]*N>= quicksum( quicksum(a.select(l,'*','*',i))*items[i].two_side for i in range(N) ) for l in range(N) )
    MIP.addConstrs(s[l]+r[l ]<= 1+sr[l] for l in range(N))
    MIP.addConstrs(s[l]+r[l] >= 2*sr[l] for l in range(N))
    MIP.setObjective(quicksum(P[l]*0.1325808 for l in range(N)) +60*quicksum( gamma[l,l] +(s[l]-sr[l]) for l in range(N) ) , GRB.MINIMIZE )

    MIP.params.outputflag = 0
    MIP.params.TimeLimit = 1200
    MIP.update()

    MIP._obj = None
    MIP._bd = None
    MIP._data = [0]
    MIP._start = time.time()

    MIP.optimize(callback=data_cb)
    MIP.write('out.lp')

    #print(MIP.Runtime)
    if MIP.status==2 or  MIP.status == 9:
        #MIP.write('soll.sol')
        #Dv=MIP.getAttr('x',D)
        try :
            Xv=MIP.getAttr('x',x)
            alphav=MIP.getAttr('x',alpha)
            betav=MIP.getAttr('x',beta)
            gammav=MIP.getAttr('x',gamma)
            rv=MIP.getAttr('x',r)
            Pv=MIP.getAttr('x',P)
            tnv=MIP.getAttr('x',tn)
            tpv=MIP.getAttr('x',tp)
            alphavv=[(i,j) for (i,j) in alphav.keys() if alphav[(i,j)]==1]

            betavv=[(i,j) for (i,j) in betav.keys() if betav[(i,j)]==1]
            gammavv=[(i,j) for (i,j) in gammav.keys() if gammav[(i,j)]==1]
            Xvv=[(i,j) for (i,j) in Xv.keys() if Xv[(i,j)]==1]
            Pvv=[(i,P[i].X) for i in Pv.keys() if Pv[i]!=0]
            tpvv=[(i,tp[i].X) for i in tpv.keys() if tpv[i]!=0]
            tnvv=[(i,tn[i].X) for i in tnv.keys() if tnv[i]!=0]
            rvv=[(i,r[i].X) for i in rv.keys() if rv[i]!=0]

            return MIP,alphavv,betavv,gammavv,Xvv,Pvv,tpvv,tnvv,rvv
        except:
            pass

    return MIP, [], [], [], [], [], [], [], []


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


def check_if_there(FName):
    WD = os.getcwd()
    Output_path = WD + "\Model\Output"
    FName = FName + "_ModelSol"
    flag = os.path.isfile(Output_path + FName)
    return flag
    
N = 12
T=3
results= []
Ns =[13, 15]#, 7, 10, 13, 15]
for T in [ 1]:
    for N in Ns:  # range(15, 17):
        for TW in ['WL', 'WS']:
            for PC in ["PL", "PS"]:
                for rep in range(5):
                    FileName = 'Data_%d_%d_%s_%s_%d' %(T, N, TW, PC, rep)
                    Data = read_object(FileName, "Input")
                    FileName += "_NoSplit"
                    if check_if_there(FileName):
                        print("I solved %s before" % FileName)
                        continue
                    start = time.time()
                    MIP, alpha, beta, gamma, x, P, tp, tn, Revolt = Solve()
                    # Draw_Bins(Data,alpha,beta,gamma)
                    runtime = time.time()-start
                    try:
                        objval = MIP.objVal
                        lowerbound = MIP.ObjBound
                        lowerbound2 = []

                        print("%s %s %s %s" % (FileName, lowerbound, objval, runtime))
                        out=Output([], [], alpha, beta, gamma, x, P, Revolt, objval, lowerbound, lowerbound2, runtime, tn, tp)
                        save_object(out, FileName)
                    except:
                        print("Last UB or LB Update: %s" %MIP._data[-1][0])
                        objval = -1
                        lowerbound = -1
                    results.append([N, T, TW, PC, rep, round(objval, 1), round(lowerbound,1),
                                    round(runtime, 1), round(MIP._data[-1][0],1)])
            #print("Lower Bound: %s" %lower_bound(Data))
            #print("Total cost of printing all bins: %s" %Best_Sol.total_cost)
            #print("Total lateness and earliness: %s" %Best_Sol.early_lateness)
            #print("Number of Bins: %s" %len(Best_Sol.Bins))
            #print ("Algorithm Run Time: %s" %Runtime)

results = pd.DataFrame(results, columns=['N', 'T', 'TW', 'PC', 'index', 'UB', 'LB', 'Time', 'L-time'])
results.to_csv("Model_NoSplit_13_15-T1.csv")
