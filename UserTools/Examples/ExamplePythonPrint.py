#####my python print script
import Store
import numpy as np
import matplotlib
matplotlib.use('Qt5Agg')
import matplotlib.pyplot as plt

def Initialise(pyinit):
    return 1

def Finalise():
    return 1

def Execute():
    # old syntax
    #a=Store.GetInt('a')
    #b=Store.GetDouble('b')
    #c=Store.GetString('c')

    # new syntax: Store.GetStoreVariable('StoreName','varName')
    a=Store.GetStoreVariable('Config','a') 
    b=Store.GetStoreVariable('Config','b')
    c=Store.GetStoreVariable('Config','c')
    print(a)
    print(b)
    print(c)

    #Store.SetStoreVariable('CStore','a',69)
    #num = Store.GetStoreVariable('CStore','a')
    #print(num)

    x_arr = np.array([1.0,2.0,3.0,4.0,5.0])

    print("Creating figure...")
    fig, ax = plt.subplots()
#    fig = plt.figure()
    print("Setting chart titles...")
    ax.set_title("test plot")
    ax.set_xlabel("x")
    ax.set_ylabel("y")

    print("Plotting...")
    ax.plot(x_arr,x_arr)
    print("Showing...")
#    plt.show()
    return 1
