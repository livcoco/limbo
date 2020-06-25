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
#import os.path, sys
#print ('\nIn pypePlugins, __init__.py was called')
#sys.path.insert(1, os.path.abspath(os.path.join('..', 'pypePlugins')))

class FunctionList:
    # The next lines initialize dicts for correspondences between
    # functions and function names.
    fNames = [] # names of base-level functions
    fDict  = {} # dictionary with fDict[name] = function with given name
    fTotal = [] # Raw list of function spec triples, (name, func, module)
    # Next, the same things but for user-function plugins:
    uNames = [];    uDict = {};    uTotal = [] 
    
    def registrar(pll):
        '''Load plugins, either from baseFuncs (if plugins list pll is empty)
        or from files listed in pll.  Eg, if pll is "abc,def,," then
        registrar will get plugins from files pypePlugins/abc.py and
        pypePlugins/def.py.  A plugin mentioned in multiple files will
        be taken from the last-mentioned file. '''
        import pypevue.pypePlugins.baseFuncs, inspect, importlib
        ref = FunctionList
        #print (f'Registrar pll = {pll}')
        if pll=='':
            ref.fDict, ref.uDict = {}, {}
            # baseFuncs will give us a complete list of base-level functions
            fs = pypevue.pypePlugins.baseFuncs.tell()
            for f in fs: # Make canonical list of fixed function names
                ref.fDict[f.__name__] = f
            ref.fNames = sorted(ref.fDict.keys())
            ref.loaded = ['baseFuncs']
        else:
            plugDir = 'pypePlugins'
            finn = [fn for fn in pll.split(',') if fn != '']            
            #print (f'Registrar pll = {pll}, finn = {finn}')
            for pfi in finn:
                try:
                    toImp = f'{plugDir}.{pfi}'
                    m = importlib.import_module(toImp) # m is a module in the pypePlugins directory
                except ModuleNotFoundError as e:
                    m = importlib.import_module(pfi) # m is a module somewhere else
                mkeys = [obj for obj, pred in inspect.getmembers(m)]
                # If the module contains a `tell` object, try calling it.
                if 'tell' in mkeys:
                    try: 
                        for f in m.tell(): # Add functions from tell() into ref.fDict{}
                            name = f.__name__ # Get the function name
                            if name in ref.fNames:
                                ref.fDict[name] = f
                            else:
                                ref.uDict[name] = f
                        ref.loaded.append(pfi)
                    except AttributeError:
                        print (f"Calling `tell` for {pfi} failed")
            print (f"Registrar got functions from {', '.join(ref.loaded)}")

        # Set class variables for all plugin functions
        ref.uNames = sorted(ref.uDict.keys())
        for n in ref.fNames: setattr(ref, n, ref.fDict[n])
        for n in ref.uNames: setattr(ref, n, ref.uDict[n])

    def clear():
        # in order to make changes to a dome which uses plugin(s) while still
        # running the same process it is necessary to clear out the old plugin(s)
        ref = FunctionList
        ref.fNames = [] # names of base-level functions
        ref.fDict  = {} # dictionary with fDict[name] = function with given name
        ref.fTotal = [] # Raw list of function spec triples, (name, func, module)
        # Next, the same things but for user-function plugins:
        ref.uNames = [];    ref.uDict = {};    ref.uTotal = [] 
    
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
