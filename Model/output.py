import os
import pickle as Pick
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import pandas as pd
from Input import Input
from Draw_the_Bin import Draw_the_Bin

class Output():
    def __init__(self,y,o,al,be,ga,x,Q,revolt,objval,low,low2,runt,tn,tp):
        self.y=y
        self.o=o
        self.alpha=al
        self.beta=be
        self.gamma=ga
        self.x=x
        self.BinQ=Q
        self.revolt=revolt
        self.objval=objval
        self.lowerbound=low
        self.lowerbound2=low2
        self.runtime=runt
        self.tn=tn
        self.tp=tp


def read_object(filename, folder):
    WD = os.getcwd()
    if folder == "Input":
        path = WD + "/%s/%s" % (folder, filename)
    else:
        path = WD + "/%s/%s_ModelSol" % (folder, filename)

    with open(path, 'rb') as input:
        obj = Pick.load(input, encoding="latin1")
    input.close()
    return obj


def check_if_there(FName):
    WD = os.getcwd()
    Output_path = WD + "/Output/"
    FName = FName + "_ModelSol"
    flag = os.path.isfile(Output_path + FName)
    return flag


def Draw(Data, bin, stripe, levels):
    # Create figure and axes
    fig, ax = plt.subplots(1)
    W = self.w * (1 + self.Revolta)
    ax.set_xlim([0, W])
    ax.set_ylim([0, self.h])

    plt.gca().set_aspect('equal', adjustable='box')
    level_StartY = 0
    for le in levels:
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

    plt.show()


def solution_display(Data, Result):
    print("############### Results ################")

    print(Data.H, Data.W)
    for (i,j) ,m in Result.o:
        print("Item %s: M= %s Q= %s h = %s w = %s Assign to = %s" %
              (j,int(m), int(Data.items[j].q/m), int(m * Data.items[j].h), Data.items[j].w ,i) )

        #print("Printing Quqntity: %s" %b.quantity)
        #print("Items quantity: ", [it.q for it in b.items] )
        #print("Revolting Bin= %d" % b.Revolta)
    # print(Result.o)
    print(Result.alpha)
    print(Result.beta)
    print(Result.gamma)
    print(Result.BinQ)
    print(Result.revolt)
    print(Result.x )


if __name__ == "__main__":
    To_excel = []
    Ns = [5]
    for T in [3]:
        for N in Ns:  # range(15, 17):
            for TW in ['WS']:
                for PC in ["PL"]:
                    for rep in [3]: #range(1,2):
                        FileName = 'Data_%d_%d_%s_%s_%d' % (T, N, TW, PC, rep)
                        # FileName += "_NoSplit"
                        if not check_if_there(FileName):
                            print("%s is not there" % FileName)
                            To_excel.append([N, T, TW, PC, rep, "", "", ""])
                            continue

                        Data = read_object(FileName, "Input")
                        Result = read_object(FileName, "Output")
                        solution_display(Data, Result)
                        print("%s : %s %s %s" %(FileName, Result.lowerbound, Result.objval, Result.runtime))
                        To_excel.append([N, T, TW, PC, rep, round(Result.lowerbound, 1), round(Result.objval, 1), round(Result.runtime, 1)])

    #To_excel = pd.DataFrame(To_excel, columns=['N', 'T', 'TW', 'PC', 'index', 'LB', 'UB', 'Time'])
    #To_excel.to_csv("L_Model_13_15t.csv")
