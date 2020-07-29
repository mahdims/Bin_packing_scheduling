# -*- coding: utf-8 -*-
"""
Created on Wed Jan 16 21:29:57 2019

@author: User
"""
import math

def lower_bound(Data):
    Paper_C=Data.papercost
    Bin_area=Data.H*Data.W
    MinBin=Data.MinBinNo
    Total_Bin_Cost=60*MinBin
    Total_printing_Cost=0
    for it in Data.items.values():
        split_no=math.floor(Bin_area/it.area)
        printing_Q = math.ceil( it.q/float(split_no) )
        Total_printing_Cost +=  Paper_C *   printing_Q   
    
    return Total_Bin_Cost + Total_printing_Cost
        
        
        