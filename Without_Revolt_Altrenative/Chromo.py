# -*- coding: utf-8 -*-
"""
Created on Wed Sep 26 13:20:26 2018

@author: mahdi
"""

from copy import copy
import numpy as np


def g(item,day,Data):
    if day<item.d-item.e:
        return Data.T
    elif day <=item.d:
        return item.d-day
    elif day<= item.d+item.l:
        return day-item.d
    elif day> item.d+item.l:
        return Data.T
    
class day:
    def __init__(self,Data,t):
        self.t=t
        self.cap=Data.proCap
        self.gamma=2/self.cap
        self.bin2print=[]
    def addBin(self,Bin):
        self.bin2print.append(Bin)
        self.cap-=Bin.quantity
        self.gamma=2/self.cap
        
class stack:
    def __init__(self,level_h,item):
        self.current_h=item.h
        self.w=item.w
        self.items=[item]
        self.level_h=level_h
        self.remaining_h=self.level_h-self.current_h
    def add(self,item):
        if self.remaining_h>=item.h and self.w>=item.w:
            self.items.append(item)
            self.current_h+=item.h
            self.remaining_h-=item.h
        else:
            print("Item did not fit in stack")
     
            
class level:
    level_list=[]
    def __init__(self,h,w):
        self.items_list=[]   
        self.remain_space=w*h
        self.h=h
        self.w=w
        self.stacks=[]
        self.remain_width=w
        level.level_list.append(self)
        
    def add(self,item):
        assign=0
        self.stacks=sorted(self.stacks,key=lambda x:x.remaining_h,reverse=False)
        for s in self.stacks:
            if s.remaining_h>=item.h and s.w==item.w:
                s.add(item)
                self.items_list.append(item)
                self.remain_space=self.remain_space-item.h*item.w
                assign=1
                break
        if assign==0 and item.w<=self.remain_width:
            self.stacks.append(stack(self.h,item))
            self.remain_width-=item.w
            self.items_list.append(item)
            self.remain_space=self.remain_space-item.h*item.w
            assign=1
        
        return (assign)
        
        
    def item_remove(self,item):
        ISlevel_eliminate=0
        
        self.remain_space-=item.area
        for st in self.stacks:
            if item in st.items:
                st.current_h-=item.h
                st.remaining_h+=item.h
                if st.current_h==0:
                    self.stacks.remove(st)
                    self.items_list.remove(item)
                    if len(self.stacks)==0:
                        ISlevel_eliminate=1
                    else:
                        self.remain_width-=st.w
                        self.h=max(s.current_h for s in self.stacks) 
                else:        
                    self.h=max(s.current_h for s in self.stacks)        
        
        return ISlevel_eliminate
        
    
    @classmethod
    def sort_levels(cls):
        cls.level_list=sorted(cls.level_list,key=lambda x:x.remain_space,reverse=False)
    
    @classmethod
    def reset(cls):
        cls.level_list=[]
        
class Bin:
    Bin_list=[]
    def __init__(self,Binlevel,Data,Revolta):
        self.h=Data.H
        self.w=Data.W/2*(2-Revolta)
        self.quantity=0
        self.remain_h=self.h
        self.levels=[]
        self.items=[]
        self.printing_weight=None
        self.printing_perferance=None
        self.utility=0
        self.Revolta=Revolta
        self.remain_space=self.h*self.w
        self.add(Binlevel)
        
        Bin.Bin_list.append(self)
    
    def add(self, Binlevel):
        
        self.levels.append(Binlevel)
        self.remain_h-=Binlevel.h
        self.items+=Binlevel.items_list
        self.remain_space-=sum([i.h*i.w for i in Binlevel.items_list])
        self.utility=1-(self.remain_space/(self.h*self.w))        
        
        
    def level_remove(self,Binlevel):
        self.levels.remove(Binlevel)  
        self.remain_h+=Binlevel.h
        for it in Binlevel.items_list:
            self.item_remove(it)
        
    def item_remove(self,it):
        self.items.remove(it)  
        self.remain_space+=it.h*it.w
        self.utility=1-(self.remain_space/(self.h*self.w))    
        # where the item was which level
        for le in self.levels:
            if it in le.items_list:
                if le.item_remove(it):
                    self.levels.remove(le)
                
            
    @classmethod
    def reset(cls):
        cls.Bin_list=[]        

class Chromo:


    def __init__(self, value ):
        self.value = value
        self.Fitness_Value = None
        self.early_lateness=None
        self.total_cost=None
        self.total_revenue=None
        self.Layout=None
        self.days=[]
        self.Bins=None
        self.Unassigned_items=[]
        

    def Fintess_Calc(self, Data,Frominitpop,pop) :
        # Objective weight
        weight=0.5
        
        self.Evaluation(Data)
        # Update the choromosome value    
        self.value=np.array(self.value)
        for i,b in enumerate(self.Bins): 
            # we will find the items index in bin b and change the self.value in these position to i 
            self.value[[it.ID for it in b.items]]=i
            
        self.value=list(self.value)
        # Calculate unit paper cost by paper type (self.papertype)
        papercost=Data.items[0].papercost
        
        LT=0
        for t,day in enumerate(self.days):
            for b in day.bin2print:
                for it in b.items:
                    if t<it.d-it.e:
                        LT+=(it.d-it.e)-t
                    elif t<=it.d+it.l:
                        LT+=0
                    elif t>it.d+it.l:
                        LT+=t-(it.d+it.l)
                        
                    
        
        
        
        #### Calulate the bin quantity ####
        
        for b in self.Bins:
            if b.Revolta==0:
                b.quantity=max([i.q for i in b.items])
            else:
                b.quantity=max([i.q/2 for i in b.items])
                    
                
            
        Printing_cost=sum([60+b.quantity*papercost for b in self.Bins])
        Total_revenue=sum([i.revenue for b in self.Bins for i in b.items])
        
           
        self.early_lateness=LT
        self.total_cost=Printing_cost
        self.total_revenue=Total_revenue
        #Total_profit=Total_revenue-Printing_cost
        #Max_profit=Total_revenue-Data.Min_printing_cost
        
        
        
        if Frominitpop:    
            self.Fitness_Value=weight*(Printing_cost/Data.Min_printing_cost)+ (1-weight)*(float(LT)/(Data.T))
        else:
            previous_LT=[p.early_lateness for p in pop]
            Range_LT=max(previous_LT)+1
            
            previous_TC=   [p.total_cost for p in pop]         
            Range_TC=max(previous_TC)+1
            
            #self.Fitness_Value=weight*(Printing_cost/Range_TC)+ (1-weight)*(float(LT)/Range_LT)
            self.Fitness_Value=weight*(Printing_cost/Data.Min_printing_cost)+ (1-weight)*(float(LT)/(Data.T))
        
        return(self)        

    def Evaluation(self,Data):
        
        Num_Bin=max(self.value)
        for b in range(Num_Bin+1):
            if b in self.value:
                Item_in_Bin=np.array(Data.items.values())[  np.where(np.array(self.value)==b) ]
                self.Finite_Best_Fit(Data,Item_in_Bin)
            
        self.Bins=Bin.Bin_list
        self.Corrective_procedure(Data)
        
        self.Scheduling_routine(Data)
        level.reset()
        Bin.reset()
        
        
    def Finite_Best_Fit(self,Data,Item_in_Bin):
        
        Revolta=0
        Bin_effective_width=Data.W
        # check if the bin will have two sided items or not      
        for it in Item_in_Bin:
            if it.two_side==1:
                Revolta=1
                Bin_effective_width=Data.W/2
                break
       
        # Sort item in one bin by their hights
        sorted_index=np.argsort([i.h for i in Item_in_Bin])[::-1]
        Item_in_Bin=Item_in_Bin[sorted_index]
        # Create the first level by the highest item
        current_level=level(Item_in_Bin[0].h,Bin_effective_width)
        current_level.add(Item_in_Bin[0])
        # Delete the item form item list
        Item_in_Bin=np.delete(Item_in_Bin,0)
        while len(Item_in_Bin)!=0 : # continue the loop until all items assigned to levels
            level.sort_levels() # sort levels by remainning space
            assign=0 # indicator for first item in item list if assigned or not
            # first we check to see if any of existed levels can accommodate the item or not            
            for l in level.level_list:       
                if l.add(Item_in_Bin[0]):
                    Item_in_Bin=np.delete(Item_in_Bin,0)
                    assign=1
                    break
            if assign==0: # if non of existed levels can accommodate the item create a new one
                current_level=level(Item_in_Bin[0].h,Bin_effective_width)
                current_level.add(Item_in_Bin[0])
                Item_in_Bin=np.delete(Item_in_Bin,0)
         
        level.sort_levels()
        Current_Bin=Bin(level.level_list[0],Data,Revolta)
        del level.level_list[0]
         
        for l in level.level_list:
            if l.h<=Current_Bin.remain_h:
                Current_Bin.add(l)
            else:
                self.Unassigned_items=self.Unassigned_items+l.items_list
        level.reset()
   
   
    def Corrective_procedure(self,Data):

        self.Bins=sorted(self.Bins,key=lambda x:x.remain_space,reverse=False)
        
        unassign_item_list=copy(self.Unassigned_items)
        for it in unassign_item_list:
            assign=0
            for b in self.Bins:
                if b.remain_space>=it.h*it.w and assign==0:
                    
                    for le in b.levels:
                        if le.add(it):
                            self.Unassigned_items.remove(it)
                            b.items.append(it)
                            assign=1
                            break
                    if assign==0 and b.remain_h>=it.h:
                        current_level=level(it.h,b.w)
                        current_level.add(it)
                        b.add(current_level)
                        self.Unassigned_items.remove(it)
                        assign=1
                        break
            if assign==0:
                # lets create a new bin to accomodate item "it"
                
                Bin_effective_width=Data.W
                Revolta=0
                if it.two_side==1:
                     Bin_effective_width=Data.W/2
                     Revolta=1
                
                current_level=level(it.h,Bin_effective_width)
                current_level.add(it)
                self.Bins.append(Bin(current_level,Data,Revolta))
                self.Unassigned_items.remove(it)
                assign=1
        
        if  len(self.Unassigned_items)!=0:
            print("Baaaang!!!")
        

        ###### Reduce number of bins if possible #####
        #if   np.average([b.utility for b in self.Bins])<=0.60:
        self.Bin_reduction(Data)
        
     
        
        
    def Random_bin_no_change(self,Data):
        Value=self.value      # for sake of safty not neccesary line   
        Current_BinNO=max(self.value)
        # randomly select the nuber of the bins
        BinNO=np.random.randint(Data.MinBinNo,2*Data.MinBinNo)

        Item_in_Bin=[]
        for b in range(Current_BinNO+1):
            if b in self.value:
                Item_in_Bin.append(np.array(Data.items.values())[  np.where(np.array(self.value)==b) ])
        
        if len(Item_in_Bin)==BinNO: # if the current number of bins is already BinNO we do not need run the following lines
            
            return self.value
            
            
        elif len(Item_in_Bin)>BinNO:
            # we decide on which bin tokeep and which bin to omit by the number of items they already have
            Item_in_Bin=sorted(Item_in_Bin,key=lambda x:len(x),reverse=True)
            Bins2keep=Item_in_Bin[:BinNO]
            Bins2omit=Item_in_Bin[BinNO:]
            
            #randomly transfer the items in Bins2omit to the bins we want to keep
            Items2reallocate=np.concatenate(tuple(Bins2omit))
            for it in Items2reallocate: 
                binindex=np.random.randint(BinNO)
                Bins2keep[binindex]=np.append(Bins2keep[binindex],it)
            
            
            #Update the self.value according to above changes
            Value=np.array(self.value)
            for i,b in enumerate(Bins2keep): 
                Value[[it.ID for it in b]]=i
            Value=list(Value)    
        
                
        return Value
            

        
    def Item_rearranging(self,Data,Bins2change,OtherBins) :
        # We first try to move levels between bins and then if some items still remains we try to move items
        
        for i in Bins2change:
            levels2move=copy(self.Bins[i].levels)
            
            #####  rearranging levels  #####
            for le in levels2move:
                
                for b in OtherBins:
                    
                    if le.h<=self.Bins[b].remain_h and le.w==self.Bins[b].w: # check the width to make sure that the bin and level are revolta or not 
                        self.Bins[i].level_remove(le)
                        self.Bins[b].add(le)
                        break
            
            ### rearranging the items   (Not ok with two sided)          
            if len (self.Bins[i].items) !=0:
                items2reassign=copy(self.Bins[i].items)
                for it in items2reassign:
                    assign=0 # check if item it is assigned somewhere
                    for b in OtherBins:
                        where2try=np.where(np.array([le.remain_space for le in self.Bins[b].levels])>=it.area)[0]
                        if len(where2try)!=0:
                            for j in where2try:
                                if self.Bins[b].levels[j].add(it):
                                    self.Bins[b].items.append(it)
                                    self.Bins[i].item_remove(it)
                                    assign=1
                                    break
                                
                        
                         # Existing levels can not accomudate item.
                         # But bin may have space for a new level
                        if assign==0 and it.h<=self.Bins[b].remain_h :
                            if self.Bins[b].Revolta==0 and it.two_side==1:
                                Do_nothing=0
                            elif it.w<=self.Bins[b].w:
                                current_level=level(it.h,self.Bins[b].w)
                                current_level.add(it)
                                self.Bins[b].add(current_level)
                                self.Bins[i].item_remove(it)
                                assign=1
                                
                        if assign==1: 
                            break        
        
    def Bin_reduction(self,Data):
        
        # find the bins with utility less than 50% to eliminate
        Bins2change=np.where(np.array([b.utility for b in self.Bins])<=0.6)[0]
        # the bins we want to add itmes to :
        OtherBins=np.array(list(set(range(len(self.Bins)))-set(Bins2change)),dtype=int)
        
        while len(OtherBins)<Data.MinBinNo :
            # If number of the bins we want to keep is less than the minimum number needed then we do the following:
            swapindex=np.random.randint(len(Bins2change))
            OtherBins=np.append(OtherBins,Bins2change[swapindex])            
            Bins2change=np.delete(Bins2change,swapindex)
            
        self.Item_rearranging(Data,Bins2change,OtherBins)
                        
        # finally we omit the bins with no items
        Bin_list=copy(self.Bins)        
        for b in Bin_list:
            if len(b.items)==0:
                self.Bins.remove(b)
        
        Low_Utility_Bins=np.where(np.array([b.utility for b in self.Bins])<=0.6)[0]
        if len(Low_Utility_Bins)>=2:        
            Bins2change=Low_Utility_Bins[1:]        
            OtherBins=[Low_Utility_Bins[0]]
            self.Item_rearranging(Data,Bins2change,OtherBins)
        
        Bin_list=copy(self.Bins)        
        for b in Bin_list:
            if len(b.items)==0:
                self.Bins.remove(b)
    def Scheduling_routine(self,Data):
        
        
  
        
        for t in range(Data.T): 
            self.days.append(day(Data,t))
            
            
        Bins=copy(self.Bins)
        while len(Bins)!=0:
            Qi=[]
            for b in Bins:
                
                sortedIndex=[]
                Printing_weight=[]
                Printing_perferance=[]
                for j,it in enumerate(b.items):
                    GIT=[g(it,t,Data) for t in range(Data.T)]
                    sortedIndex.append(np.argsort(GIT))
                    Printing_weight.append(GIT[sortedIndex[j][1]]-GIT[sortedIndex[j][0]])
                    Printing_perferance.append(sortedIndex[j][0])
                
                
                b.printing_weight=np.mean(Printing_weight)
                # the printing perference of the bin decides when a bin should printed 
                # we obtain bin printing preference equal to the printing perfernce of the most urgent itme in the bin 
                b.printing_perferance=Printing_perferance[np.argmax(Printing_weight)]
                
                Qi.append( b.printing_weight )
                    
            Bin2assign=np.argmax(Qi)
            self.days[Bins[Bin2assign].printing_perferance].addBin(Bins[Bin2assign])
            del Bins[Bin2assign]
            
            
            
    
