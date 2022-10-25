import Store  # i don't believe this causes problems as it's defined in ToolAnalysis
#import numpy as np
#import tensorflow as tf

def Initialise(pyinit):
    print("i've initialized")
    return 1

def Execute():
    print("i'm executing")
    print(type(Store))
    Store.SetStoreVariable('a',6)
    return 1

def Finalise():
    print("i've finalised")
    #a=Store.GetStoreInt('a')
    #print(a)
    return 1

