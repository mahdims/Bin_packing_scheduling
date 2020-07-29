# -*- coding: utf-8 -*-
"""
Created on Wed Oct 03 13:02:53 2018

@author: mahdi
"""

import matplotlib.pyplot as plt
import matplotlib.patches as patches

def Draw_the_Bin(Data, Bin):
    
    # Create figure and axes
    fig,ax = plt.subplots(1)
    ax.set_xlim([0,Data.W])
    ax.set_ylim([0,Data.H])
    
    
    plt.gca().set_aspect('equal', adjustable='box')
    level_StartY=0
    for le in Bin.levels   :
        rect = patches.Rectangle((0,level_StartY),le.w,le.h,linewidth=3,edgecolor='b',facecolor='none')
        ax.add_patch(rect)        
        
        
        stack_StartX=0
        for st in le.stacks:
            rect = patches.Rectangle((stack_StartX,level_StartY),st.w,st.level_h,linewidth=1,edgecolor='r',facecolor='none')
            ax.add_patch(rect)            
            
            
            item_StartY=level_StartY
            for it in st.items:
                rect = patches.Rectangle((stack_StartX,item_StartY),it.w,it.h,linewidth=1,edgecolor='g',facecolor='y')
                ax.add_patch(rect)  
                ax.text(stack_StartX+it.w/2, item_StartY+it.h/2 ,it.ID,horizontalalignment='center',verticalalignment='center')
                item_StartY+=it.h
            
            item_StartY=level_StartY
            stack_StartX+=st.w
        
        
        
        level_StartY+=le.h
   
   
   

    plt.show()