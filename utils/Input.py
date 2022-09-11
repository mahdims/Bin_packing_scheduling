import pandas as pd
import numpy as np
import math
import copy


class Item:
    df = []

    def __init__(self, j, i, df):
        self.name = "%d" % i
        self.ID = j
        self.h = df['H'][i]
        self.w = df['W'][i]
        self.q = df['PrintCount'][i]
        self.d = 0
        self.Order_NO = i
        self.two_side = df['Two Sided'][i]
        self.papersize = df["Paper size"][i]
        self.papertype = df["Paper type"][i]
        self.area = self.h*self.w
        # earliness limit currently is equal for every item but can be load from data file
        self.e = 2
        # lateness limit 
        self.l = 1
        # can be defined differently by reading the paper type and determine its cost
        self.papercost = 0.1325808
        self.orientation_correcting()
        
    def orientation_correcting(self):
        if self.papersize == "A2":
            temp = self.h
            self.h = self.w
            self.w = temp
        if self.papersize == "A4":
            temp = self.h
            self.h = self.w
            self.w = temp
        if self.papersize == "A6":
            temp = self.h
            self.h = self.w
            self.w = temp
            
    def duplicate(self, k, q):
        NewItem = copy.deepcopy(self)
        NewItem.name = "%d_%d" % (self.ID, k)
        NewItem.q = q
        return NewItem

    def __repr__(self):
        return f"{self.ID}: {self.name}: {self.q}"


class Input:
    def __init__(self, ItemN, Time):
        Quantities = pd.read_excel('Traj.xlsx', sheetname='Quantities')  
        self.Possible_Quantities = np.array(Quantities)
        self.items = {}
        self.proCap = 0
        self.MinBinNo = 0
        self.Min_printing_cost = 0
        self.Lateness = 0
        self.Earliness = 0
        self.N = ItemN
        self.T = Time
        self.H = 118.9
        self.W = 84.1
        self.BinCost = 60
        self.papercost = 0.1325808

    def randomly_select(self, df, proCap, Lateness, Earliness):
        order_inx = np.random.choice(range(64929), self.N, replace=False)
        Item.df = df
        items = {}
        Min_printing_cost = 0
        for j,i in enumerate(order_inx):
            items[j] = Item(j, i, df)
            items[j].d = np.random.randint(0, self.T)
            Min_printing_cost += (items[j].area/(self.H*self.W))*items[j].papercost*items[j].q
        self.items = items
        MaxPrint = max([it.q for it in self.items.values()])
        self.proCap = int((MaxPrint/2)*proCap+(1-proCap)*MaxPrint)
        self.MinBinNo = math.ceil((sum([a.area for a in self.items.values()]))/(self.H*self.W))
        self.Min_printing_cost = Min_printing_cost+self.BinCost*self.MinBinNo
        self.Lateness = Lateness
        self.Earliness = Earliness

