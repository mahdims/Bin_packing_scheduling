
import pandas as pd
import numpy as np
from Input import Input
from gams import *
import cPickle as Pick

def read_object(FileName,folder):
    if folder=="Input":
        path = "G:\My Drive\Side Projects\Cutting stock (APA)\Code\Model\%s\%s" %(folder,FileName)
    else:
        path = "G:\My Drive\Side Projects\Cutting stock (APA)\Code\Model\%s\%s_ModelSol" %(folder,FileName)
    
    with open(path , 'rb') as input:
        obj= Pick.load(input)
    return  obj 
    
def save_object(obj, filename):
    with open(filename, 'wb') as output:  # Overwrites any existing file.
        Pick.dump(obj, output, Pick.HIGHEST_PROTOCOL)

def Export_Data_to_Gdx(FileName,Data):
    N=Data.N
    T=Data.T  
    d=[Data.items[i].d for i in range(N)]
    Qmax=max([Data.items[a].q for a in range(N)])
    Qmin=min([Data.items[a].q for a in range(N)])
    NN= [(i,j) for i in range(N) for j in range(N) if j>=i ] 
    M = [int((Data.H*Data.W)/Data.items[i].area) for i in range(N)  ]
    Mmax=int(max(M))
    Mmin=int(min(M))
    alphaNN=[]
    for i,j in NN:
      if (Data.items[i].w==Data.items[j].w and Data.items[i].h+Data.items[j].h<=Data.H) or i==j:
         alphaNN.append( [ str(i) , str(j) ] )
         
      
    
    ws = GamsWorkspace( 	working_directory = "G:\My Drive\Side Projects\Cutting stock (APA)\Code\Model\Input")
    db = ws.add_database()
    
    i = db.add_set("i", 1, "Items")
    for p in range(N):
        i.add_record( str(p) )
        
    t = db.add_set("t", 1, "Time")
    for p in range(T):
        t.add_record( str(p) )
    
    due= db.add_parameter_dc("due", [i], "Due date of item i")
    for j,du in enumerate(d):
        due.add_record(str(j)).value = int(du)
    
    MaxS= db.add_parameter_dc("M", [i], "Max split of item i")
    for j,m in enumerate(M):
        MaxS.add_record(str(j)).value = int(m)  
    
    itemH=db.add_parameter_dc("h", [i], "hight of item i")
    itemW=db.add_parameter_dc("w", [i], "wide of item i")
    for j in range(N):
        itemH.add_record(str(j)).value=Data.items[j].h
        itemW.add_record(str(j)).value=Data.items[j].w
    
    quantity=db.add_parameter_dc("q", [i], "quantity of item i")  
    for j in range(N):
        quantity.add_record(str(j)).value=float(Data.items[j].q)
    
    Two_side= db.add_parameter_dc("S", [i], "If item i is two sided or not")
    for  j in range(N):
        Two_side.add_record(str(j)).value = int(Data.items[j].two_side)
        
    NN = db.add_parameter("N", 0, "number of items")
    NN.add_record().value = N
    
    Procap= db.add_parameter("Cap", 0, "Production capacity")
    Procap.add_record().value = Data.proCap
    
    BinH = db.add_parameter("binH", 0, "Bin H")
    BinH.add_record().value = Data.H
    
    BinW = db.add_parameter("binW", 0, "Bin W")
    BinW.add_record().value = Data.W
    
    Pcost = db.add_parameter("P_cost", 0, "paper cost")
    Pcost.add_record().value = Data.papercost
    
    BinCost = db.add_parameter("B_cost", 0, "one bin cost")
    BinCost.add_record().value = Data.BinCost
    
    Earliness_allos = db.add_parameter("Early",0,"Maximum erliness allowed")
    Earliness_allos.add_record().value = Earliness
    Lateness_allos = db.add_parameter("Late",0,"Maximum latness allowed")
    Lateness_allos.add_record().value = Lateness
    
    
    db.export("%s.gdx" %FileName)


#Read the big data at first just once!
df = pd.read_excel('Big_Data.xlsx', sheetname='Sheet1')  
n=5
rep=0
for n in [12]:
    for rep in [1]:
    #n=7
        Pro_Cap_indi=1 # 1 means tight production capacity and 0 means ample capacity
        Lateness = 2
        Earliness = 1
        T=3
        N_of_items=n
        
        Data=Input(N_of_items,T)
        Data.randomly_select(df,Pro_Cap_indi,Lateness,Earliness)
        
        FileName="Data_%d_%d_%d_LO" %(N_of_items,T,rep)
        #Data  =read_object(FileName,"Input")
        #Data.Latness=Lateness
        #FileName="Data_%d_%d_%d_ST" %(N_of_items,T,rep)
        save_object(Data,FileName)
        Export_Data_to_Gdx(FileName,Data)