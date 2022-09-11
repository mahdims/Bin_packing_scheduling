# -*- coding: utf-8 -*-
"""
Created on Mon Oct 01 01:50:02 2018

@author: mahdi
"""

class stock:
    def __init__(self,level_h,item):
        self.current_h=item.h
        self.w=item.w
        self.items=[item]
        self.level_h=level_h
        self.remaining_h=self.level_h-self.current_h
    def add(self,item):
        if self.remaining_h>=item.h and self.w<=item.w:
            self.items.append(item)
            self.current_h+=item.h
            self.remaining_h-=item.h
        else:
            print("Item did not fit in stock")
class level:
    level_list=[]
    def __init__(self,h,w):
        self.items_list=[]   
        self.remain_space=w*h
        self.h=h
        self.w=w
        self.stocks=[]
        self.remain_width=w
        level.level_list.append(self)
        
    def add(self,item):
        assign=0
        self.stocks=sorted(self.stocks,key=lambda x:x.remaining_h,reverse=True)
        for s in self.stocks:
            if s.remaining_h>=item.h and s.w<=item.w:
                s.add(item)
                self.items_list.append(item)
                self.remain_space=self.remain_space-item.h*item.w
                assign=1
                break
        if assign==0 and item.w<=self.remain_width:
            self.stocks.append(stock(self.h,item))
            self.remain_width-=item.w
            self.items_list.append(item)
            self.remain_space=self.remain_space-item.h*item.w
            assign=1
        return (assign)
    
    @classmethod
    def sort_levels(cls):
        cls.level_list=sorted(cls.level_list,key=lambda x:x.remain_space,reverse=True)
    
    @classmethod
    def reset(cls):
        cls.level_list=[]
        
class Bin:
    Bin_list=[]
    def __init__(self,level):
        self.h=Data.H
        self.w=Data.W
        self.remain_h=self.h-level.h
        self.levels=[level]
    
    def add(self, level):
        self.levels.append(level)
        self.remain_h-=level.h