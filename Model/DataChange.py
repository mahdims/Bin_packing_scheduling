# -*- coding: utf-8 -*-
"""
Created on Fri Mar 12 20:47:34 2010

@author: generic
"""
import cPickle as Pick 

def read_object(FileName,folder):
    if folder=="Input":
        path = "G:\My Drive\Side Projects\Cutting stock (APA)\Code\Model\%s\%s" %(folder,FileName)
    else:
        path = "G:\My Drive\Side Projects\Cutting stock (APA)\Code\Model\%s\%s_ModelSol" %(folder,FileName)
    
    with open(path , 'rb') as input:
        obj= Pick.load(input)
    input.close()
    return  obj 

def save_object(obj, filename):
    folder="Input"
    path = "G:\My Drive\Side Projects\Cutting stock (APA)\Code\Model\%s\%s" %(folder,FileName)
    with open(path, 'wb') as output:  # Overwrites any existing file.
        Pick.dump(obj, output, Pick.HIGHEST_PROTOCOL)

N = 5
T=1

for rep in range(10):
    FileName="Data_%d_%d_%d" %(N,T,rep)
    Data= read_object(FileName,"Input")
    Min_printing_cost=0
    for i in range(Data.N):
        Min_printing_cost+=(Data.items[i].area/(Data.H*Data.W))*Data.items[i].papercost*Data.items[i].q
        Data.Min_printing_cost=Min_printing_cost+Data.BinCost*Data.MinBinNo        
        #Data.items[i].Order_NO=Data.items[i].ID
        #Data.items[i].ID = i
    save_object(Data, FileName)
        