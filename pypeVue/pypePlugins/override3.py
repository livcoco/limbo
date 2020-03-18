# jiw March 2020
'''Test-program for pypeVue plugins setup'''

from pypeSim import Point, Post, Cylinder, Layout
from pypePlugins import FunctionList as ref
from sys import _getframe
import os.path
def saywhat(fr): print (f'{fr.f_code.co_name:15} @ {os.path.basename(__file__)}')

def addEdge(v,w, layout):      saywhat(_getframe())
def addEdges(v,w, layout):     saywhat(_getframe())
def arithmetic(line, xTrace):  saywhat(_getframe())
def autoAdder(fout):           saywhat(_getframe())
def scriptCyl(ss, preCyl):     saywhat(_getframe())
def scriptPost(ss, prePost):
    saywhat(_getframe())
    print (f'len(ref.fNames) = {len(ref.fNames)},  len(ref.fList) = {len(ref.fList)},  len(ref.fDict) = {len(ref.fDict)}')

def dotell():                   # An example error to test error sidestep
    return (addEdge, addEdges, arithmetic, autoAdder, scriptCyl,scriptPost)
