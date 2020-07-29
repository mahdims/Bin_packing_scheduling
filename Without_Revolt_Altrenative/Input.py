import pandas as pd
import math
class Item:
    def __init__(self,i,df):
        self.ID=i
        self.h=df['H'][i]
        self.w=df['W'][i]
        self.d=df['D'][i]
        self.q=df['Tiraj'][i]
        self.Order_NO=df['Order NO'][i]
        self.two_side=df['Two Sided'][i]
        self.papersize=df["Paper size"][i]      
        self.papertype=df["Paper type"][i]
        self.revenue=float(df["Selling price"][i].split()[0])
        self.area=self.h*self.w
        
        #earliness limit currently is equal for every item but can be load from data file
        self.e=2
        # lateness limit 
        self.l=1  
        # can be define diffrently by reading the paper type and determine its cost 
        self.papercost=0.1325808

        
        self.orientation_correcting()
        
    def orientation_correcting(self):
        if self.papersize=="A2":
            temp=self.h
            self.h=self.w
            self.w=temp
        if self.papersize=="A4":
            temp=self.h
            self.h=self.w
            self.w=temp
        if self.papersize=="A6":
            temp=self.h
            self.h=self.w
            self.w=temp
        
        
class Input:
    def __init__(self,ItemN,Time):
        # sheet 3 all orders is one sided in sheet 2 we have two sided order 
        df = pd.read_excel('Data2.xlsx', sheetname='Sheet2')   
        self.N=ItemN
        self.T=Time
        self.proCap=70000
        self.H=130
        self.W=90
        items={}
        Min_printing_cost=0
        for i in range (self.N):
            items[i]=Item(i,df)
            Min_printing_cost+=(items[i].area/(self.H*self.W))*items[i].papercost*items[i].q
            
        self.items=items
        self.MinBinNo=math.ceil((sum([a.area for a in self.items.values()]))/(self.H*self.W))
        self.Min_printing_cost=Min_printing_cost+60*self.MinBinNo

        
