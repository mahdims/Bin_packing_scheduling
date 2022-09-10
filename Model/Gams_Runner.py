from gams import *
import cPickle as Pick
from output import Output
import os



def save_object(obj, filename):
    with open("C:\Users\generic\Dropbox\APA printing\Code\Model\Output\%s_ModelSol" %FileName , 'wb') as out:  # Overwrites any existing file.
        Pick.dump(obj, out, Pick.HIGHEST_PROTOCOL)
    out.close()


def check_if_there(FName):
    Output_path = "C:\Users\generic\Dropbox\APA printing\Code\Model\Output"
    FName = FName + "_ModelSol"
    flag = os.path.isfile(Output_path + FName)
    return flag

ws = GamsWorkspace(working_directory="C:\Users\generic\Dropbox\APA printing\Code\Model")
Objs=[]
Ns =[13,15]#, 7, 10, 13, 15]
for T in [1, 3]:
    for N in Ns:  # range(15, 17):
        for TW in ['WL', 'WS']:
            for PC in ["PL", "PS"]:
                for rep in range(5):
                    FileName = 'Data_%d_%d_%s_%s_%d' %(T,N,TW,PC,rep)

                    if check_if_there(FileName):
                        print("I solved %s before" %FileName)
                        continue

                    db = ws.add_database_from_gdx("C:\Users\generic\Dropbox\APA printing\Code\Model\Input\%s.gdx" %FileName)
                    opt = GamsOptions(ws)
                    opt.defines["gdxincname"] = db.name
                    # opt.all_model_types = "xpress"
                    t1 = ws.add_job_from_file("Spliting_model-Correction.gms")
                    try:
                        t1.run(opt, databases = db)
                        # print("We just solve problem %s" %FileName )
                        Values={}
                        for name in ['alpha','beta','gamma','o','y','x','tn','tp','BQ','r','P','runtime','lowerbound','lowerbound2','objval']:
                            Values[name] =  t1.out_db[name]

                        y = [int(a.level) for a in Values['y'] if a.level!=0]
                        o = [((int(a.keys[0]),int(a.keys[1])),a.level) for a in Values['o'] if a.level!=0]
                        alpha=[ (int(a.keys[0]),int(a.keys[1])) for a in Values['alpha'] if a.level==1]
                        beta = [ (int(a.keys[0]),int(a.keys[1])) for a in Values['beta'] if a.level==1]
                        gamma= [ (int(a.keys[0]),int(a.keys[1])) for a in Values['gamma'] if a.level==1]
                        x = [ (int(a.keys[0]),int(a.keys[1])) for a in Values['x'] if a.level==1]
                        tn    = [a.level for a in Values['tn'] ]
                        tp    = [a.level for a in Values['tp'] ]
                        Printing_Q =[a.level for a in Values['P']]
                        Revolt = [a.level for a in Values['r'] ]
                        objval = t1.out_db["objval"].first_record().value
                        lowerbound = t1.out_db["lowerbound"].first_record().value
                        lowerbound2 = t1.out_db["lowerbound2"].first_record().value
                        runtime= t1.out_db["runtime"].first_record().value

                        out=Output(y,o,alpha,beta,gamma,x,Printing_Q,Revolt,objval,lowerbound,lowerbound2,runtime,tn,tp)
                        save_object(out,FileName)
                        print("%s %s %s %s" %(FileName, lowerbound, objval, runtime))
                        Objs.append(objval)
                    except:
                        print("Error in solving the problem Infeasibility Ckeck")


