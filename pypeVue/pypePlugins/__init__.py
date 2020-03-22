#!/usr/bin/env python3
'''__init__.py -- Purposes of this file include: (1) [Existence of
   __init__.py makes this directory represent a module. [Need to see
   if this claim is true]] (2) Via sys.path.insert(), it adds this
   directory to python's import search list, so that other modules can
   import FunctionList without needing the directory name. (3) It
   creates the FunctionList class, with method `registrar` and with
   class variables that link to functions. (4) It describes how `tell`
   functions work.

'''

import os.path, sys
#print ('\nIn pypePlugins, __init__.py was called')
sys.path.insert(1, os.path.abspath(os.path.join('..', 'pypePlugins')))

class FunctionList:
    # The next lines initialize dicts for correspondences between
    # functions and function names.
    fNames = [] # names of base-level functions
    fDict  = {} # dictionary with fDict[name] = function with given name
    fTotal = [] # Raw list of function spec triples, (name, func, module)
    # Next, the same things but for user-function plugins:
    uNames = [];    uDict = {};    uTotal = [] 
    
    def registrar():
        import pypePlugins.baseFuncs, inspect, importlib
        fReg, uReg = {}, {}
        # baseFuncs will give us a complete list of base-level functions
        fs = pypePlugins.baseFuncs.tell()
        for f in fs: # Make canonical list of fixed function names
            fReg[f.__name__] = f
        fNames = sorted(fReg.keys())
        loaded = ['baseFuncs']
        plugDir = 'pypePlugins'
        ppdir = os.path.join('.', plugDir)
        pfiles = sorted(fn[:-3] for fn in os.listdir(ppdir) if fn.endswith(".py"))
        for pfi in pfiles:
            if pfi == '__init__' or pfi == 'baseFuncs':  continue
            toImp = f'{plugDir}.{pfi}'
            m = importlib.import_module(toImp) # m is a module
            mkeys = [obj for obj, pred in inspect.getmembers(m)]
            # If the module contains a `tell` object, try calling it.
            if 'tell' in mkeys:
                try: 
                    for f in m.tell(): # Add functions from tell() into fReg{}
                        name = f.__name__ # Get the function name
                        if name in FunctionList.fNames:
                            fReg[name] = f
                        else:
                            uReg[name] = f
                    loaded.append(pfi)
                except AttributeError:
                    print (f"Calling `tell` for {pfi} failed")
        print (f"Registrar got functions from {', '.join(loaded)}")
        uNames = sorted(uReg.keys())
        
        FunctionList.fNames = fNames
        FunctionList.fDict  = fReg        
        FunctionList.uNames = uNames
        FunctionList.uDict  = uReg

        # Set class variables for all plugin functions
        for n in fNames: setattr(FunctionList, n, fReg[n])
        for n in uNames: setattr(FunctionList, n, uReg[n])

def tell():                     # Example of a tell() function ...
    '''A tell() function returns a list or tuple of functions.  Each
    program in this directory should contain a tell() function, by
    which pypeVue discovers functions to use as addons or overrides of
    its own fixed functions.  Functions not listed in a tell will be
    ignored by pypeVue and not used as addons or overrides.  If an
    unqualified user function name is in several tell's, the one from
    the highest-lexing module will be used.  Note, registrar does not
    call this example.    '''
    def whatFunc(layout): # Functions can be module level or local, etc
        return None
    return (someFunc, whatFunc)
# For user functions, pypeVue provides one argument, a layout.
def someFunc(layout): return otherstuff(layout, 4, 3)
