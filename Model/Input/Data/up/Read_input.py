import cPickle as Pick 
import os
from Input import Input


# The Data we read from the files are an object from class Input
#The items in self.items are objects of class Item.
# after importing the input objcet from file by cPickle we can call all the atributes 

       

# this function read the file 
def read_object(FileName,In_path=None):
    if In_path is None:
        path = "%s\%s" %(os.getcwd(),FileName)
    else:
        path = "%s\%s" %(In_path,FileName)
    
    with open(path , 'rb') as input:
        obj= Pick.load(input)
    return  obj 

FileName = "ACap2"
Data= read_object(FileName)
print(Data.N, Data.T, Data.proCap)