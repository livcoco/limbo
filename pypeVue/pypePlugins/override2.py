# jiw March 2020
'''Test-program for pypeVue plugins setup'''

from pypeSim import Point, Post, Cylinder, Layout
from pypePlugins import FunctionList as ref
from sys import _getframe
import os.path
def saywhat():
    fname = _getframe(1).f_code.co_name
    print (f'{fname:15} @ {os.path.basename(__file__)}')

def addEdge(v,w, layout):      saywhat()
def addEdges(v,w, layout):     saywhat()
def arithmetic(line, xTrace):  saywhat()
def autoAdder(fout):           saywhat()
def scriptCyl(ss, preCyl):     saywhat()
def scriptPost(ss, prePost):
    saywhat()
    print (f'len(ref.fNames) = {len(ref.fNames)},  len(ref.fList) = {len(ref.fList)},  len(ref.fDict) = {len(ref.fDict)}')

def tell():
    def whatFunc(layout): # Functions can be module level or local, etc
        return None
    return (addEdge, addEdges, arithmetic, autoAdder, scriptCyl,scriptPost, someFunc, whatFunc)
# For user functions, pypeVue provides one argument, a layout.
def someFunc(layout): return otherstuff(layout, 4, 3)