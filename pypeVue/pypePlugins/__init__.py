# __init__.py -- The purpose of this file is fourfold: (1) Its
# existence makes this directory represent a module.  (2) It adds this
# directory to python's import search list, so that other programs in
# this directory can easily import the FunctionList class. (3) It
# gives an example of a `tell` function. (4) It creates the
# FunctionList class, with method `registrar` and with class variables
# which provide pypeStubs and plugins with links to functions.

import os, sys, pypePlugins, inspect, importlib
print ('\nIn pypePlugins, __init__.py was called\n')
sys.path.insert(1, os.path.abspath(os.path.join('..', 'pypePlugins')))

def tell():                     # Example of a tell() function ...
    '''A tell() function returns a list or tuple of functions.  Each
    program in this directory should contain a tell() function, by
    which pypeVue discovers user functions to use as addons or
    overrides of its own functions.  Functions not listed in a tell
    will be ignored by pypeVue and not used as addons or overrides.
    If an unqualified user function name is in several tell's, the one
    from the highest-lexing module will be used.  Note, registrar does
    not call this example.    '''
    def whatFunc(layout): # Functions can be module level or local, etc
        return None
    return (someFunc, whatFunc)
# For user functions, pypeVue provides one argument, a layout.
def someFunc(layout): return otherstuff(layout, 4, 3)

class FunctionList:
    # The next lines create variables that embed correspondences
    # between function names in fNames and functions listed in fList.
    fNames = [] # pypeStubs.tell() will load this with function names
    fList  = []  # Will be list of functions corresponding to names
    fDict  = {}  # Will be dictionary of correspondences from zip(fNames,fList)
    fTotal = []  # Will be overall list of triples, (name, func, module)
    # next are for same accesses to user-function plugins.
    uNames = [];    uList = [];    uDict = {};    uTotal  = [] 
    
    def registrar():
        # pypeStubs will give us lists of base-level functions + names
        fs = pypePlugins.pypeStubs.tell()
        loaded = ['pypeStubs']
        reg, ureg = {}, {}
        for f in fs:
            reg[f.__name__] = f
        plugDir = 'pypePlugins'
        d = os.path.join('.', plugDir)
        pfiles = sorted(fn[:-3] for fn in os.listdir(d) if fn.endswith(".py"))    
        for pfi in pfiles:
            if pfi == '__init__' or pfi == 'pypeStubs':  continue
            toImp = f'{plugDir}.{pfi}'
            m = importlib.import_module(toImp) # m is a module
            mkeys = [obj for obj, pred in inspect.getmembers(m)]
            # If the module contains a `tell` object, try calling it.
            if 'tell' in mkeys:
                try: 
                    for f in m.tell(): # Add functions from tell() into reg{}
                        name = f.__name__ # Get the function name
                        if name in FunctionList.fNames:
                            reg[name] = f
                            fTotal.append((name, f, pfi))
                        else:
                            ureg[name] = f
                            uTotal.append((name, f, pfi))
                    loaded.append(pfi)
                except AttributeError:
                    print (f"Calling `tell` for {pfi} failed")
        print (f"Obtained functions from {', '.join(loaded)}")
        # Note, FunctionList.fNames was set up by pypeStubs.tell()
        FunctionList.fDict  = reg
        FunctionList.fList  = [reg[x] for x in FunctionList.fNames]
        FunctionList.uDict  = ureg
        FunctionList.uNames = sorted(ureg.keys())
        FunctionList.uList  = [ureg[x] for x in FunctionList.uNames]
        