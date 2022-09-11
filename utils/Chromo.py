from copy import copy
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import numpy as np
import math
import random
from scipy import stats


def divsorset_gen(Data, Items):

    Div_Set = {}
    biggest = []
    max_ratio = 0
    for it in Items:
        RealM = int(math.floor(Data.H*Data.W/it.area))
        if math.ceil(it.q/float(RealM)) >= max_ratio:
            max_ratio = math.ceil(it.q/float(RealM))
        Div_Set[it.ID] = [math.ceil(it.q/float(m)) for m in range(2, RealM+1)]
        biggest += Div_Set[it.ID]
        
    # Set=[div for div in  Div_Set[minix] if all([div in a for a in Div_Set.values() ])]
    biggest = list(set(biggest))
    biggest.sort()
    biggest = biggest[biggest.index(max_ratio):]
    return biggest


def RGA_divsorset_gen(Data, Items):
    Div_Set = []
    max_ratio = 0  # the smallest possible PQ
    for it in Items:
        RealM = int(math.floor(Data.H / it.h))
        Div_Set += [math.ceil(it.q / float(m)) for m in range(2, RealM + 1)]

        if math.ceil(it.q/float(RealM)) >= max_ratio:
            max_ratio = math.ceil(it.q/float(RealM))

    Div_Set = list(set(Div_Set))
    Div_Set.sort()
    Div_Set = Div_Set[Div_Set.index(max_ratio):]

    return Div_Set

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


def RGA_divide(Data, Item, divisor):
    BigNo = Item.q
    m = int(math.ceil(BigNo / float(divisor)))
    RealM = int(math.floor(Data.H / Item.h))
    m = min(m, RealM)
    if float(BigNo) / m == BigNo / m:
        Fianl_Qs = [BigNo / m for _ in range(m)]
    else:
        Fianl_Qs = [BigNo / m for _ in range(m)]
        remain = BigNo - sum(Fianl_Qs)
        Fianl_Qs[-1] += remain

    q = max(Fianl_Qs)
    m = len(Fianl_Qs)
    Splits_item = Item.duplicate(m, q)
    Splits_item.h = m * Item.h
    Splits_item.area = Splits_item.h * Splits_item.w

    return m, [Splits_item]


def One_Bin_Balance(Data,RGA_flag, IB_threshold, Item_in_Bin, Revolta):
    indicator = 0
    Unassigns = []
    survive_mode = 1
    IB_threshold= 0
    if Item_in_Bin.inbalance >= IB_threshold :
        Item_in_Bin.Qsort()
        itemsList = copy(Item_in_Bin.list)

        if RGA_flag:
            divisorset = RGA_divsorset_gen(Data, itemsList)
        else:
            divisorset = divsorset_gen(Data, itemsList)
        
        if len(divisorset) == 0 and Item_in_Bin.inbalance != 0:
            survive_mode = 1
            divisorset = [[Item_in_Bin[-1].q], [Item_in_Bin[0].q/2]]
            
        possible_bins = []

        if len(itemsList) == Data.N and Revolta==1:
            stop = 0

        for divisor in divisorset:
            items = ItemSet(Data, itemsList)
            for it in itemsList:

                if RGA_flag:
                    m, bigest_item_splits = RGA_divide(Data, it, divisor)  # split the large quantity item
                else:
                    bigest_item_splits = divide(it, divisor)
                    m = len(bigest_item_splits)

                if items.bin_remain_area >= (m-1)*it.area and m != 1:
                    items.add(bigest_item_splits) # add splits to items
                    items.remove(it) # remove the orginal item
                    items.Qsort() # sort items according to their quantity

                    # TODO: Have an alternative FBF that gives priority to splited items
                    Unassigns, Current_Bin = Chromo.Finite_Best_Fit(Data, items.list, Revolta)

                    # check to see if there is unassign item split !
                    if all([it.name == str(it.ID) for it in Unassigns]) :
                        possible_bins.append((items.inbalance, items.MaxQ(), Current_Bin, Unassigns))
                        
                        if items.MaxQ() == Item_in_Bin.qmode and not survive_mode:
                            Current_Bin.add2list()      
                            indicator = 1
                            return indicator, Unassigns
                    else:
                        break
                        #continue
                else:
                    break
                
        if len(possible_bins) != 0:

            Best_Bin = min(possible_bins, key=lambda x: (x[1], x[0]))

            Best_Bin[2].add2list()
            indicator = 1
            Unassigns = Best_Bin[3]
        #(Unassigns, Current_Bin )=Chromo.Finite_Best_Fit(Data,np.array(items.list),Revolta)       
        #Current_Bin.add2list() 
        #indicator=1
    
    return indicator, Unassigns


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


def Quantity_similarity(Item_Q, Bin_Q):
    indicator = 0
    if Item_Q <= (1+0.2)*Bin_Q:
        indicator = 1
    return indicator


def inbalance_measure(items):
    if len(items) == 0:
        return 0
    lst = [it.q for it in items]
    return np.std(lst)/np.mean(lst)


class day:

    def __init__(self, Data, t):
        self.t = t
        self.cap = Data.proCap
        self.gamma = 2/(max(0, self.cap)+0.0001)
        self.bin2print = []

    def addBin(self, Bin):
        self.bin2print.append(Bin)
        self.cap -= Bin.quantity
        self.gamma = 2/(max(0, self.cap) + 0.0001)


class ItemSet:
    def __init__(self, Data, List):
        self.list = copy(List)
        self.qmode = stats.mode([it.q for it in List])[0][0]
        self.inbalance = inbalance_measure(List)
        self.bin_remain_area = Data.W*Data.H-sum([it.w*it.h for it in self.list])

    def __getitem__(self, position):
        return self.list[position]
        
    def __len__(self):
         return len(self.list)   
         
    def MaxQ(self):
        return max([it.q for it in self.list])
        
    def Qsort(self):
        self.list = sorted(self.list, key=lambda it: (it.q, it.h), reverse=True)

    def remove(self, item):
        self.list.remove(item)
        if len(self.list) != 0:
            self.qmode = stats.mode([i.q for i in self.list])[0][0]
            self.inbalance = inbalance_measure( self.list )
            self.bin_remain_area += item.w*item.h
        else:
            self.qmode = 0
            self.inbalance = 0
            self.bin_remain_area += item.w*item.h
            
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
        if assign==0 and item.w<=self.remain_width and self.h>=item.h:
            self.stacks.append(stack(self.h,item))
            self.remain_width-=item.w
            self.items_list.append(item)
            self.remain_space=self.remain_space-item.h*item.w
            assign=1
        
        return (assign)

    def item_remove(self,item):
        ISlevel_eliminate = 0
        
        self.remain_space -= item.area
        for st in self.stacks:
            if item in st.items:
                self.items_list.remove(item)
                st.items.remove(item)
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

    Bin_list = []

    def __init__(self, Binlevel, Data, Revolta, id=0):
        self.id = id
        self.h = Data.H
        self.w = Data.W/2*(2-Revolta)
        self.quantity = 0
        self.remain_h = self.h
        self.levels = []
        self.items = []
        self.printing_weight = None
        self.printing_perferance = None
        self.P_day = 0
        self.utility =0
        self.Revolta = Revolta
        self.Two_sided_item = None
        self.remain_space = self.h*self.w
        self.add(Binlevel)
        self.calc_quantity()

    def __repr__(self):
        st = ""
        for it in self.items:
            st = st + ","+str(it.ID)
        return st+f" | Q: {self.quantity}"

    def add2list(self):
        self.id = len(Bin.Bin_list)
        Bin.Bin_list.append(self)
    
    def add(self, Binlevel):
        
        self.levels.append(Binlevel)
        self.remain_h -= Binlevel.h
        self.items += Binlevel.items_list
        self.remain_space-=sum([i.h*i.w for i in Binlevel.items_list])
        self.utility=1-(self.remain_space/(self.h*self.w))
        self.calc_quantity()

    def calc_quantity(self):
        if len(self.items) == 0:
            self.quantity = 0
        else:
            self.quantity = max([it.q for it in self.items])*(1/(1+self.Revolta))

    def level_remove(self, Binlevel):
        self.levels.remove(Binlevel)  
        self.remain_h += Binlevel.h
        for it in Binlevel.items_list:
            self.item_remove(it)
        
    def item_remove(self,it):
        self.items.remove(it)  
        self.remain_space+=it.h*it.w
        self.utility=1-(self.remain_space/(self.h*self.w))   

        # where the item was? \\which level
        for le in self.levels:
            if it in le.items_list:
                # print([itm.ID for itm in le.items_list])
                if le.item_remove(it):
                    self.levels.remove(le)
                    self.remain_h += le.h

        self.calc_quantity()

    def Draw(self):
        # Create figure and axes
        fig, ax = plt.subplots(1)
        W = self.w * (1 + self.Revolta)
        ax.set_xlim([0, W])
        ax.set_ylim([0, self.h])

        plt.gca().set_aspect('equal', adjustable='box')
        level_StartY = 0
        for le in self.levels:
            rect = patches.Rectangle((0, level_StartY), le.w, le.h, linewidth=3, edgecolor='b', facecolor='none')
            ax.add_patch(rect)

            stack_StartX = 0
            for st in le.stacks:
                rect = patches.Rectangle((stack_StartX, level_StartY), st.w, st.level_h, linewidth=3, edgecolor='g',
                                         facecolor='none')
                ax.add_patch(rect)

                item_StartY = level_StartY
                for it in st.items:
                    rect = patches.Rectangle((stack_StartX, item_StartY), it.w, it.h, linewidth=1, edgecolor='r',
                                             facecolor='none')
                    ax.add_patch(rect)
                    if it.name == str(it.Order_NO):
                        it.name = it.ID
                    ax.text(stack_StartX + it.w / 2, item_StartY + it.h / 2, it.name, horizontalalignment='center',
                            verticalalignment='center')
                    item_StartY += it.h

                item_StartY = level_StartY
                stack_StartX += st.w

            level_StartY += le.h

        try:
            plt.show()
        except UserWarning:
            print(f"We couldn't show the Bins here but we save the bin in Bin_{self.id}.png file.")
            plt.savefig(f'Bin_{self.id}.png', bbox_inches='tight')

    @classmethod
    def reset(cls):
        cls.Bin_list = []


class Chromo:

    def __init__(self, Pars, iterationNO, value, Revolting):

        self.RGA_flag = Pars.RGA_flag
        self.IM = Pars.IM
        self.iterationNO = iterationNO
        self.Maxit = Pars.Maxit
        self.BR_rep = Pars.BR_rep
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

    def __repr__(self):

        return "Fit: %s TC: %s TL: %s" %(round(self.Fitness_Value,2), round(self.total_cost,1), self.early_lateness)
     
    def Random_bin_no_change(self, Data): # it is a mutation operator
        Value = self.value      # for sake of safty not neccesary line
        Current_BinNO = max(self.value)
        # randomly select the number of the bins
        BinNO = np.random.randint(Data.MinBinNo, 2*Data.MinBinNo)
        Item_in_Bin = []
        for b in range(Current_BinNO+1):
            if b in self.value:
                Item_in_Bin.append(np.array(list(Data.items.values()))[np.where(np.array(self.value) == b)])
        
        if len(Item_in_Bin) == BinNO: # if the current number of bins is already BinNO we do not need run the following lines
            
            return self.value

        elif len(Item_in_Bin) > BinNO:
            # we decide on which bin to keep and which bin to omit by the number of items they already have
            Item_in_Bin = sorted(Item_in_Bin, key=lambda x: len(x), reverse=True)
            Bins2keep = Item_in_Bin[:BinNO]
            Bins2omit = Item_in_Bin[BinNO:]
            
            #randomly transfer the items in Bins2omit to the bins we want to keep
            Items2reallocate = np.concatenate(tuple(Bins2omit))
            for it in Items2reallocate: 
                binindex = np.random.randint(BinNO)
                Bins2keep[binindex] = np.append(Bins2keep[binindex], it)

            #Update the self.value according to above changes
            Value = np.array(self.value)
            for i, b in enumerate(Bins2keep):
                Value[[it.ID for it in b]] = i
            Value = list(Value)
        else:
            NOBins2add = BinNO - len(Item_in_Bin)
            NewBins = [len(Item_in_Bin) + a for a in range(NOBins2add)]

            for NewB in NewBins:
                for it, b in enumerate(Value):
                    if random.random() <= 0.4:
                       Value[it] = NewB

        return Value

    def Fitness_measure(self, N, it, pop):
        lambda_min = 2
        lambda_max = 2 * N
        LT_weight = lambda_min + it * (lambda_max - lambda_min) / self.Maxit

        previous_LT = [p.early_lateness for p in pop]
        Max_LT = max(previous_LT) + 1

        previous_TC = [p.total_cost for p in pop]
        Max_TC = max(previous_TC) + 1

        self.Fitness_Value = (self.total_cost/Max_TC) + LT_weight * (float(self.early_lateness)/Max_LT)

    def Fitness_Calc(self, Data, pop):

        self.Evaluation(Data)
        # Update the choromosome value and revolting
        self.value=np.array(self.value)
        Revolting=[]
        for i,b in enumerate(self.Bins): 
            # we will find the items index in bin b and change the self.value in these position to i 
            self.value[[it.ID for it in b.items]]=i
            Revolting.append(b.Revolta)
            
        self.value = list(self.value)
        self.Revolting = Revolting
        papercost=Data.items[0].papercost

        # Calculate lateness and earliness
        LT = 0
        for t, day in enumerate(self.days):
            for b in day.bin2print:
                for it in b.items:
                    if t<it.d-it.e:
                        LT+=(it.d-it.e)-t
                    elif t<=it.d+it.l:
                        LT+=0
                    elif t>it.d+it.l:
                        LT+=t-(it.d+it.l)

        # Check the bin include tow sided items AND
        # Calculate the bin quantity ####
        for b in self.Bins:
            
            Two_sided_item=any([it.two_side for it in b.items])
            b.Two_sided_item=Two_sided_item

            b.calc_quantity()

        self.early_lateness = LT
        self.total_cost = sum([60*(1 + 1*(1-b.Revolta)*b.Two_sided_item )+b.quantity*papercost for b in self.Bins])

        # add penalty cost if split items are not in same bins O(n2)
        # Time consuming check that we avoid ! Since it is ensured by

        return self

    def Evaluation(self, Data):
        Unassigned_items = []
        Num_Bin = max(self.value)

        for b in range(Num_Bin+1):
            Item_in_Bin = np.array(list(Data.items.values()))[np.where(np.array(self.value) == b)]
            Unassigns = self.Create_balanced_bins(Data, Item_in_Bin, self.Revolting[b])
            Unassigned_items += Unassigns

        self.Bins = Bin.Bin_list
        Bin.reset()
        self.Corrective_procedure(Data, Unassigned_items)
        
        #### Reduce number of bins if possible ####
        # just use it in even iterations 
        
        if float(self.iterationNO/float(self.BR_rep)).is_integer():
            self.Bin_reduction(Data)
                  
        self.Scheduling_routine(Data)
        level.reset()

    def Create_balanced_bins(self, Data, Item_in_Bin, Revolta):

        Item_in_Bin = ItemSet(Data, list(Item_in_Bin))

        (indicator, Unassigned) = One_Bin_Balance(Data, self.RGA_flag, self.IM, Item_in_Bin, Revolta)

        if indicator == 0:
            (Unassigned, Current_Bin) = Chromo.Finite_Best_Fit(Data, Item_in_Bin.list, Revolta)
            Current_Bin.add2list()

        return Unassigned

    @classmethod    
    def Finite_Best_Fit(cls, Data, Item_in_Bin, Revolta):
        Unassigned_items=[]
        level.reset()        
        
        # check if the bin will have two sided items or not        
        Bin_effective_width=(1-Revolta)*Data.W+Revolta*Data.W/2
       
        # Sort item in one bin by their hights
        # sorted_index=np.argsort([i.h for i in Item_in_Bin])[::-1]
        # Item_in_Bin=Item_in_Bin[sorted_index]
        Item_in_Bin = sorted(Item_in_Bin, key=lambda x: (x.h, x.name != str(x.ID)), reverse=True)

        # Create the first level by the highest item
        current_level=level(Item_in_Bin[0].h,Bin_effective_width)
        current_level.add(Item_in_Bin[0])
        # Delete the item form item list
        del Item_in_Bin[0]
        while len(Item_in_Bin)!=0 : # continue the loop until all items assigned to levels
            level.sort_levels() # sort levels by remainning space
            assign=0 # indicator for first item in item list if assigned or not
            # first we check to see if any of existed levels can accommodate the item or not            
            for l in level.level_list:       
                if l.add(Item_in_Bin[0]):
                    del Item_in_Bin[0]
                    assign=1
                    break
            if assign==0: # if non of existed levels can accommodate the item create a new one
                current_level=level(Item_in_Bin[0].h,Bin_effective_width)
                current_level.add(Item_in_Bin[0])
                del Item_in_Bin[0]
         
        level.level_list = sorted(level.level_list, key=lambda x:x.h, reverse=True)
        Current_Bin = Bin(level.level_list[0], Data, Revolta)
        
        del level.level_list[0]
         
        for l in level.level_list:
            if l.h <= Current_Bin.remain_h:
                Current_Bin.add(l)
            else:
                Unassigned_items = Unassigned_items+l.items_list
        level.reset()
        
        return Unassigned_items, Current_Bin
   
    def Corrective_procedure(self,Data,Unassigned_items):
        
        self.Bins=sorted(self.Bins,key=lambda x:x.remain_space,reverse=False)
        
        unassign_item_list=copy(Unassigned_items)
        for it in unassign_item_list:
            assign=0
            for b in self.Bins:

                if b.remain_space >= it.h*it.w and assign == 0:
                    
                    for le in b.levels:
                        if le.add(it):
                            Unassigned_items.remove(it)
                            b.items.append(it)
                            assign = 1
                            break
                    if assign == 0 and b.remain_h >= it.h:
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
                self.Bins.append(Bin(current_level,Data,Revolta, id=len(self.Bins)))
                Unassigned_items.remove(it)
                assign=1
        
        if  len(Unassigned_items)!=0:
            print("Baaaang!!!")

    def Item_rearranging(self, Data, Bins2change, OtherBins):
        # We first try to move levels between bins and then if some items
        # still remains we try to move items
        for i in Bins2change:
            levels2move = copy(self.Bins[i].levels)
            
            ###  Rearranging levels  ###
            for le in levels2move:
                Q_of_level = max([it.q for it in le.items_list])
                
                for b in OtherBins:
                    
                    if le.h <= self.Bins[b].remain_h and le.w <= self.Bins[b].w \
                            and Quantity_similarity(Q_of_level, self.Bins[b].quantity):
                                                
                        self.Bins[i].level_remove(le)
                        if le.w<self.Bins[b].w: 
                            le.width_increase()
                        self.Bins[b].add(le)
                        break
            # rearranging the items
            if len (self.Bins[i].items) != 0:
                items2reassign = copy(self.Bins[i].items)
                for it in items2reassign:
                    assign = 0  # check if item it is assigned somewhere
                    for b in OtherBins:

                        if not Quantity_similarity(it.q, self.Bins[b].quantity):
                            continue

                        where2try = np.where(np.array([le.remain_space for le in self.Bins[b].levels]) >= it.area)[0]
                        if len(where2try) != 0 and assign == 0:
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
        # find the bins with utility less than 50% to eliminate
        Min_utility_level = 0.5
        Bins2change = np.where(np.array([b.utility for b in self.Bins]) <= Min_utility_level)[0]
        # the bins we want to add itmes to :
        OtherBins = np.array(list(set(range(len(self.Bins)))-set(Bins2change)), dtype=int)
        self.Item_rearranging(Data, Bins2change, OtherBins)
        
        '''        
        while len(OtherBins)<Data.MinBinNo :
            # If number of the bins we want to keep is less than the minimum number needed then we do the following:
            swapindex=np.random.randint(len(Bins2change))
            OtherBins=np.append(OtherBins,Bins2change[swapindex])            
            Bins2change=np.delete(Bins2change,swapindex)
            
        '''
        # finally we omit the bins with no items
        Bin_list = copy(self.Bins)
        for b in Bin_list:
            if len(b.items) == 0:
                self.Bins.remove(b)
        # (New) Mergeing low utility bins together
        Low_Utility_Bins = np.where(np.array([b.utility for b in self.Bins]) <= Min_utility_level)[0]
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
                if set([ite.ID for ite in b.items]) == set([0,4,3]):
                    stop = 0
                sortedIndex=[]
                Printing_weight=[]
                Printing_perferance=[]
                GIT = []
                for j,it in enumerate(b.items):
                    GIT.append([g(it,t,Data) for t in range(Data.T)])
                    sortedIndex.append(np.argsort(GIT[j]))
                    Printing_weight.append(GIT[j][sortedIndex[j][1]]-GIT[j][sortedIndex[j][0]])
                    Printing_perferance.append(sortedIndex[j][0])
                
                
                b.printing_weight=np.mean(Printing_weight)

                # the printing preference of the bin decides when a bin should printed
                # we obtain bin printing preference equal to the printing perfernce of the most urgent item in the bin
                most_urgent_value = max(Printing_weight)
                indices = [i for i, x in enumerate(Printing_weight) if x == most_urgent_value]
                if len(indices) == 1:
                    b.printing_preference = Printing_perferance[indices[0]]
                else:
                    good_days = set([a for a in sortedIndex[indices[0]] if GIT[0][a] != Data.T])
                    for i in indices[1:]:
                        good_days = good_days.intersection(set([a for a in sortedIndex[i] if GIT[i][a] != Data.T]))

                    if len(good_days) == 0:
                        b.printing_preference = Printing_perferance[indices[0]]
                    else:
                        b.printing_preference = list(good_days)[0]

                Qi.append(b.printing_weight)
                    
            Bin2assign=np.argmax(Qi)
            self.days[Bins[Bin2assign].printing_preference].addBin(Bins[Bin2assign])
            del Bins[Bin2assign]


def g(item, day, Data):
    if day < item.d-item.e:
        return Data.T
    elif day <= item.d:
        return item.d-day
    elif day <= item.d+item.l:
        return day-item.d
    elif day > item.d+item.l:
        return Data.T