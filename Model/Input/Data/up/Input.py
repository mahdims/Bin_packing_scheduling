import pandas as pd
import numpy as np
import math
import copy
class Item:
    df=[]
    def __init__(self,j,i,df):
        self.name="%d" %(i)
        self.ID=j
        self.h=df['H'][i]
        self.w=df['W'][i]
        self.q=df['PrintCount'][i]
        self.Order_NO=i
        self.two_side=df['Two Sided'][i]
        self.papersize=df["Paper size"][i]      
        self.papertype=df["Paper type"][i]
        self.area=self.h*self.w
        self.papercost = 0.1325808
        self.orientation_correcting()
        
    
        
class Input:
    def __init__(self,ItemN,Time,items,pro_Cap,Lateness,Earliness):
        self.N=ItemN
        self.T=Time
        self.H=118.9
        self.W=84.1
        self.BinCost = 60
        self.papercost= 0.1325808
        self.items=items
        self.proCap=pro_Cap
        Min_printing_cost=0
        for j,i in enumerate(items):
            Min_printing_cost+=(items[j].area/(self.H*self.W))*items[j].papercost*items[j].q

        self.MinBinNo=math.ceil((sum([a.area for a in self.items.values()]))/(self.H*self.W))
        self.Min_printing_cost=Min_printing_cost+self.BinCost*self.MinBinNo
        self.Lateness  = Lateness
        self.Earliness = Earliness