from copy import copy
import itertools
from Draw_the_Bin import Draw_the_Bin
import numpy as np
import math
import random
from scipy import stats


def One_Bin_Balance(Data,IB_threshold,Item_in_Bin,Revolta):
    #### This function ### CHANGE ###    
    
    indicator=0
    Unassigns=[]

    if Item_in_Bin.inbalance > IB_threshold:
        Item_in_Bin.Qsort()

        itemsList=copy(Item_in_Bin.list)
        divisorset=[a for a in Data.Possible_Quantities if a >= Item_in_Bin.qmode and a<Item_in_Bin[0].q]
        possible_bins=[]
        for divisor in divisorset:
            items=ItemSet(Data,itemsList)
            for it in itemsList:
                bigest_item_splits=divide(it,divisor[0])  # split the large quantity item
                if items.bin_remain_area >= (len(bigest_item_splits)-1)*it.area and len(bigest_item_splits) !=1 :
                    items.add(bigest_item_splits) # add splits to items
                    items.remove(it) # remove the orginal item
                    items.Qsort() # sort items according to their quantity
                    
                    # check if really splits can fit in one bin
                    (Unassigns,Current_Bin )=Chromo.Finite_Best_Fit(Data,np.array(items.list),Revolta) 
                    
                    if all([it.name == '%s'%it.ID for it in Unassigns]): # check function  to see if there is unassign form the same ID
                        possible_bins.append( (items.inbalance,items.MaxQ(),Current_Bin,Unassigns) )   
                        
                        if items.MaxQ() ==  Item_in_Bin.qmode:
                            Current_Bin.add2list()      
                            indicator=1
                            return (indicator,Unassigns)    
                            
                                           
                    else :
                        break
                    
                
                else:
                    break
                
        if len(possible_bins) != 0 :
            Best_Bin = min(possible_bins, key = lambda x: x[1])
            Best_Bin[2].add2list()
            indicator=1
            Unassigns = Best_Bin[3]

        #(Unassigns, Current_Bin )=Chromo.Finite_Best_Fit(Data,np.array(items.list),Revolta)       
        #Current_Bin.add2list() 
        #indicator=1
    
    return (indicator,Unassigns)
     
def split_item_set(Data, IM_goal ,item_set):
    ### this function ### CHANGE ###
    IM=item_set.inbalance
    item_set.Qsort()
    # check if we need splitting at all
    if IM!=0:
        i=len(item_set)/2
        
        # select part a (if it is the big quntities or small ones)
        if inbalance_measure(item_set.list[:i])<inbalance_measure(item_set.list[i:]):
            part_a=item_set.list[:i]
            part_b=item_set.list[i:]
        else:
            item_set.list=item_set.list[::-1]
            part_a=item_set.list[:-i]
            i=len(part_a)
            part_b=item_set.list[i:]
            
        
        IM_A=inbalance_measure(part_a)
        IM_B=inbalance_measure(part_b)
        Measure0=sum([abs(IM_A-IM_goal),abs(IM_B-IM_goal)])
        # add item to part a while the total IM diffrence with IM goal is improving
        Bpart=copy(part_b)[1:]
        IM_A=inbalance_measure(part_a+[item_set.list[i]])
        IM_B=inbalance_measure(Bpart)
        while sum([abs(IM_A-IM_goal),abs(IM_B-IM_goal)]) <= Measure0:
            
            Measure0=sum([abs(IM_A-IM_goal),abs(IM_B-IM_goal)])
            
            part_a.append(item_set.list[i])
            part_b.remove(item_set.list[i])
            if (i+1)==len(item_set):
                break
            
            i+=1
            Bpart.remove(item_set.list[i])
            IM_A=inbalance_measure(part_a+[item_set.list[i]])
            IM_B=inbalance_measure(Bpart)
                 
            
            
        part_a=ItemSet(Data , part_a)
           
        if len(part_b)==0:
            part_b=[]
        else:
            part_b=ItemSet(Data , part_b) 
            
        return (part_a, part_b)
    else:
         return (item_set, [])
         
def divide(Item,divisor):
    BigNo=Item.q
    m=int(math.floor(BigNo/divisor))
    out=m*[divisor]
    remain=BigNo-m*divisor
    if remain!=0:
        out+=[BigNo-m*divisor]
    out1=[]
    for k,q in enumerate(out):
        out1.append(Item.duplicate(k,q))
    return out1
    
def Quantity_similarity(Item_Q,Bin_Q): ### CHANGE ###
    indicator=0
    if  Item_Q <= (1+0.2)*Bin_Q:
        indicator=1
    return indicator
    
    
def inbalance_measure(items):
    if len(items)==0: return 0
    lst=[ it.q for it in items ]
    return np.std(lst)/np.mean(lst)
    
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

class ItemSet:
    def __init__(self,Data,List):
        self.list=copy(List)
        self.qmode=stats.mode([it.q for it in List])[0][0]
        self.inbalance=inbalance_measure( List ) 
        self.bin_remain_area=Data.W*Data.H-sum([it.w*it.h for it in self.list])  

    def __getitem__(self, position):
        return self.list[position]
        
    def __len__(self):  ### CHANGE ###
         return len(self.list)   
         
    def MaxQ(self):
        return max([it.q for it in self.list])
        
    def Qsort(self):
        self.list=sorted(self.list,key=lambda it:it.q,reverse=True)

    def remove(self,item):
        self.list.remove(item)
        if len(self.list)!=0:
            self.qmode=stats.mode([i.q for i in self.list])[0][0]
            self.inbalance=inbalance_measure( self.list ) 
            self.bin_remain_area+= item.w*item.h
        else:
            self.qmode=0
            self.inbalance==0
            self.bin_remain_area+= item.w*item.h
            
    def add(self,items):
        for i in items :
            self.list.append(i)
            self.bin_remain_area-= i.w*i.h
        self.qmode=stats.mode([it.q for it in self.list])[0][0]
        self.inbalance=inbalance_measure( self.list ) 
        
        
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
    def width_increase(self):
        self.w=self.w*2
        self.remain_width+=self.w/2
        self.remain_space+=(self.w/2)*self.h
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
                self.items_list.remove(item)
                st.current_h-=item.h
                st.remaining_h+=item.h
                if st.current_h==0:
                    self.stacks.remove(st)

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
        self.Two_sided_item=None
        self.remain_space=self.h*self.w
        self.add(Binlevel)
    def add2list(self):
        Bin.Bin_list.append(self)
    
    def add(self, Binlevel):
        
        self.levels.append(Binlevel)
        self.remain_h-=Binlevel.h
        self.items+=Binlevel.items_list
        self.remain_space-=sum([i.h*i.w for i in Binlevel.items_list])
        self.utility=1-(self.remain_space/(self.h*self.w))    
        self.quantity=max([it.q for it in Binlevel.items_list])*(1/2*(2-self.Revolta))
        
        
    def level_remove(self,Binlevel):
        self.levels.remove(Binlevel)  
        self.remain_h+=Binlevel.h
        for it in Binlevel.items_list:
            self.item_remove(it)
        
    def item_remove(self,it):
        self.items.remove(it)  
        self.remain_space+=it.h*it.w
        self.utility=1-(self.remain_space/(self.h*self.w))   
        if it.q==self.quantity*(1+self.Revolta):
            self.quantity=max([i.q for i in self.items])*(1/2*(2-self.Revolta))
        # where the item was? \\which level
        for le in self.levels:
            if it in le.items_list:
                if le.item_remove(it):
                    self.levels.remove(le)
                    self.remain_h+=le.h
                
            
    @classmethod
    def reset(cls):
        cls.Bin_list=[]        

class Chromo:


    def __init__(self, value,Revolting ):
        self.value = value
        self.Revolting=Revolting
        self.Fitness_Value = None
        self.early_lateness=None
        self.total_cost=None
        self.total_revenue=None
        self.Layout=None
        self.days=[]
        self.Bins=None
        self.Unassigned_items=[]
     
    def Random_bin_no_change(self,Data): # it is a mutation operator
        Value=self.value      # for sake of safty not neccesary line   
        Current_BinNO=max(self.value)
        # randomly select the number of the bins
        BinNO=np.random.randint(Data.MinBinNo,2*Data.MinBinNo)

        Item_in_Bin=[]
        for b in range(Current_BinNO+1):
            if b in self.value:
                Item_in_Bin.append(np.array(Data.items.values())[  np.where(np.array(self.value)==b) ])
        
        if len(Item_in_Bin)==BinNO: # if the current number of bins is already BinNO we do not need run the following lines
            
            return self.value
            
            
        elif len(Item_in_Bin)>BinNO:
            # we decide on which bin to keep and which bin to omit by the number of items they already have
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
        else :
            NOBins2add = BinNO - len(Item_in_Bin)
            NewBins=[len(Item_in_Bin) + a for a in range(NOBins2add) ]

            for NewB in NewBins:
                for it,b in enumerate(Value):
                    if random.random() <= 0.4:
                       Value[it]= NewB

            
            
                
        return Value
        
    def Fitness_Calc(self, Data,pop,iterationNO,Maxit) :
        # Objective weight
        weight = 1 # just printing cost
        #weight = 0 # just lateness and earliness cost
        #weight=0.5
                
        self.Evaluation(Data,iterationNO,Maxit)
        #### Update the choromosome value and revolting #### 
        self.value=np.array(self.value)
        Revolting=[]
        for i,b in enumerate(self.Bins): 
            # we will find the items index in bin b and change the self.value in these position to i 
            self.value[[it.ID for it in b.items]]=i
            Revolting.append(b.Revolta)
            
        self.value=list(self.value)
        self.Revolting=Revolting
        
        # Obtain unit paper cost by paper type (self.papertype)
        papercost=Data.items[0].papercost
        
        
        # Calculate lateness and earliness
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
                        
        
        ### See if the bin include tow sided items AND
        ### Calulate the bin quantity ####
        
        for b in self.Bins:
            
            Two_sided_item=any([it.two_side for it in b.items])
            b.Two_sided_item=Two_sided_item
            
            if b.Revolta==0:
                b.quantity=max([i.q for i in b.items])
            else:
                b.quantity=max([i.q/2 for i in b.items])
                    
                
        ### Kalep cost is 60 \\ Bin without revolting that have two sided items dubel the cost of kalep    
        Printing_cost=sum([60*(1 + 1*(1-b.Revolta)*b.Two_sided_item )+b.quantity*papercost for b in self.Bins])
        #Total_revenue=sum([i.revenue for b in self.Bins for i in b.items])
        
           
        self.early_lateness=LT
        self.total_cost=Printing_cost
        #self.total_revenue=Total_revenue
        #Total_profit=Total_revenue-Printing_cost
        #Max_profit=Total_revenue-Data.Min_printing_cost
        
        
        
        if iterationNO == 1:    
            self.Fitness_Value=weight*(Printing_cost/Data.Min_printing_cost)+ (1-weight)*(float(LT)/(Data.T))
        else:
            previous_LT = [p.early_lateness for p in pop]
            Range_LT=max(previous_LT)+1
            
            previous_TC = [p.total_cost for p in pop]         
            Range_TC=max(previous_TC)+1
            
            #self.Fitness_Value=weight*(Printing_cost/Range_TC)+ (1-weight)*(float(LT)/Range_LT)
            self.Fitness_Value=weight*(Printing_cost/Data.Min_printing_cost)+ (1-weight)*(float(LT)/(Data.T))
        
        
        ### CHANGE ###]
        # add penalty cost if splitted items are not in same bins #
        for i,b in enumerate(self.Bins):
            itemset=[bi.items for bi in self.Bins[:i]+self.Bins[i+1:]]
            itemset=list(itertools.chain.from_iterable(itemset))
            IDset=map(lambda x:x.ID,itemset)
            condition=[item.ID in IDset for item in b.items]
            if any(condition):
                self.Fitness_Value+=50000*self.Fitness_Value                
                break
                    
                    
        return(self)        

    def Evaluation(self,Data,iterationNO,Maxit):
        Unassigned_items=[]
        Num_Bin=max(self.value)

        for b in range(Num_Bin+1):
            Item_in_Bin=np.array(Data.items.values())[  np.where(np.array(self.value)==b) ]
            Unassigns=self.Create_balanced_bins(Data,Unassigned_items,Item_in_Bin,self.Revolting[b])  
            Unassigned_items+=Unassigns
        
        
        self.Bins=Bin.Bin_list
        Bin.reset()
        
        self.Corrective_procedure(Data,Unassigned_items)
        
        #### Reduce number of bins if possible ####
        # just use it in even iterations 
        
        if  np.average([b.utility for b in self.Bins])<=0.9 and float(iterationNO/float(1)).is_integer() :
            self.Bin_reduction(Data)
                  
        self.Scheduling_routine(Data)
        level.reset()
        
        
    def Create_balanced_bins(self, Data,Unassigns,Item_in_Bin,Revolta) : ## This function is ### CHANGE ###
        IM_goal=0.35    
        Item_in_Bin=ItemSet(Data,list(Item_in_Bin))
        
        (indicator,Unassigns)=One_Bin_Balance(Data,IM_goal,Item_in_Bin,Revolta)


        if indicator==0:
            (Unassigns,Current_Bin )=Chromo.Finite_Best_Fit(Data,np.array(Item_in_Bin.list),Revolta) 
            Current_Bin.add2list()

                    
        return Unassigns                   

    @classmethod    
    def Finite_Best_Fit(cls,Data,Item_in_Bin,Revolta):
        Unassigned_items=[]
        level.reset()        
        
        # check if the bin will have two sided items or not        
        Bin_effective_width=(1-Revolta)*Data.W+Revolta*Data.W/2
       
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
                Unassigned_items=Unassigned_items+l.items_list
        level.reset()
        
        return (Unassigned_items,Current_Bin )
   
    def Corrective_procedure(self,Data,Unassigned_items):
        
        self.Bins=sorted(self.Bins,key=lambda x:x.remain_space,reverse=False)
        
        unassign_item_list=copy(Unassigned_items)
        for it in unassign_item_list:
            assign=0
            for b in self.Bins:
                if b.remain_space>=it.h*it.w and assign==0:
                    
                    for le in b.levels:
                        if le.add(it):
                            Unassigned_items.remove(it)
                            b.items.append(it)
                            assign=1
                            break
                    if assign==0 and b.remain_h>=it.h:
                        current_level=level(it.h,b.w)
                        current_level.add(it)
                        b.add(current_level)
                        Unassigned_items.remove(it)
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
                Unassigned_items.remove(it)
                assign=1
        
        if  len(Unassigned_items)!=0:
            print("Baaaang!!!")
        
       
    def Item_rearranging(self,Data,Bins2change,OtherBins) :
        ### CHANGE ###        
        
        ### We first try to move levels between bins and then if some items 
        #still remains we try to move items
        for i in Bins2change:
            levels2move = copy( self.Bins[i].levels )
            
            ###  Rearranging levels  ###
            for le in levels2move:
                Q_of_level = max([it.q for it in le.items_list])
                
                for b in OtherBins:
                    
                    if le.h <= self.Bins[b].remain_h and le.w <= self.Bins[b].w and Quantity_similarity(Q_of_level,self.Bins[b].quantity): # check the width to make sure that the bin and level are revolta or not 
                                                
                        self.Bins[i].level_remove(le)
                        if le.w<self.Bins[b].w: 
                            le.width_increase()
                        self.Bins[b].add(le)
                        break
            
            ### rearranging the items             
            if len (self.Bins[i].items) != 0:
                items2reassign = copy(self.Bins[i].items)
                for it in items2reassign:
                    assign = 0 # check if item it is assigned somewhere
                    for b in OtherBins:
                        if Quantity_similarity(it.q,self.Bins[b].quantity): break 
                        where2try = np.where(np.array([le.remain_space for le in self.Bins[b].levels])>=it.area)[0]
                        if len(where2try)!=0 and assign == 0:
                            for j in where2try:
                                if self.Bins[b].levels[j].add(it):
                                    self.Bins[b].items.append(it)
                                    self.Bins[i].item_remove(it)
                                    assign = 1
                                    break
                                
                        
                         # Existing levels can not accomudate item.
                         # But bin may have space for a new level
                        if assign == 0 and it.h <= self.Bins[b].remain_h :
                             if it.w <= self.Bins[b].w:
                                current_level = level(it.h,self.Bins[b].w)
                                current_level.add(it)
                                self.Bins[b].add(current_level)
                                self.Bins[i].item_remove(it)
                                assign = 1
                                
                        if assign == 1: 
                            break        
        
    
    def Bin_reduction(self,Data):
        ### CHANGE ###
        # find the bins with utility less than 50% to eliminate
        Min_utility_level = 0.5
        Bins2change = np.where(np.array([b.utility for b in self.Bins])<=Min_utility_level)[0]
        # the bins we want to add itmes to :
        for b in   Bins2change:      
            OtherBins = np.array(list(set(range(len(self.Bins)))-set(Bins2change)),dtype=int)
            self.Item_rearranging(Data,Bins2change,OtherBins)
        
        '''        
        while len(OtherBins)<Data.MinBinNo :
            # If number of the bins we want to keep is less than the minimum number needed then we do the following:
            swapindex=np.random.randint(len(Bins2change))
            OtherBins=np.append(OtherBins,Bins2change[swapindex])            
            Bins2change=np.delete(Bins2change,swapindex)
            
        '''
        # finally we omit the bins with no items
        Bin_list=copy(self.Bins)        
        for b in Bin_list:
            if len(b.items)==0:
                self.Bins.remove(b)
        ## (New) Mergeing low utility bins together
        Low_Utility_Bins=np.where(np.array([b.utility for b in self.Bins])<=Min_utility_level)[0]
        counter=0
        while len(Low_Utility_Bins)>=2 and counter<=5:  
            counter+=1
            Low_Utility_Bins=sorted(Low_Utility_Bins,key=lambda x:self.Bins[x].quantity,reverse=False)
            OtherBins=[Low_Utility_Bins[0]]
            Bins2change=Low_Utility_Bins[1:]        
            if counter>=2:
                OtherBins=[Low_Utility_Bins[-1]]
                Bins2change=Low_Utility_Bins[:-1] 
            self.Item_rearranging(Data,Bins2change,OtherBins)
        
            Bin_list=copy(self.Bins)        
            for b in Bin_list:
                if len(b.items)==0:
                    self.Bins.remove(b)
                    
            Low_Utility_Bins=np.where(np.array([b.utility for b in self.Bins])<=Min_utility_level)[0]
                    
                    
                
    def Scheduling_routine(self,Data):
        
        #Check if I consider production capacity....  I did !! 
        
        
        for t in range(Data.T): 
            self.days.append(day(Data,t))
  
        Bins=copy(self.Bins)

        if Data.T==1:
            for b in Bins:
                self.days[0].addBin(b)
            return
            
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
                b.printing_preference=Printing_perferance[np.argmax(Printing_weight)]
                
                Qi.append( b.printing_weight )
                    
            Bin2assign=np.argmax(Qi)
            self.days[Bins[Bin2assign].printing_preference].addBin(Bins[Bin2assign])
            del Bins[Bin2assign]
            
            
            
    
