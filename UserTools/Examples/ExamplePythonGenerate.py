#### my python value Generator
import Store

def Initialise(pyinit):
    return 1

def Finalise():
    return 1

def Execute():
    print('using python to change values')
    # old syntax
    #Store.SetInt('a',6)
    #Store.SetDouble('b',8.0)
    #Store.SetString('c','hello')

    f = Store.GetStoreVariable('Config','Factor')
    print("Factor: ",f)
    # new syntax: Store.SetStoreVariable('StoreName','varName',newValue)
    #Store.SetStoreVariable('Config','a',66+f)
    Store.SetStoreVariable('Config','b',4.20)
    Store.SetStoreVariable('Config','c','bye')
    return 1

