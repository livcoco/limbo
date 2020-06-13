#/usr/bin/python3
# jiw March 2020
'''Test-program for pypeVue plugins setup'''

from pypevue.pypeVue import Point, Post, Cylinder, Layout
from pypevue.pypePlugins import FunctionList as ref
from sys import _getframe
import os.path
def saywhat():
    fname = _getframe(1).f_code.co_name
    print (f'{fname:15} @ {os.path.basename(__file__)}')

def autoAdder(fout):
    saywhat()
    print (f'len(ref.fNames) = {len(ref.fNames):2},  len(ref.fDict) = {len(ref.fDict)}\nlen(ref.uNames) = {len(ref.uNames):2},  len(ref.uDict) = {len(ref.uDict)}')

def tell():
    def whatFunc(layout): # Functions can be module level or local, etc
        saywhat()
        return None
    return (autoAdder, someFunc, whatFunc)
# For user functions, pypeVue provides one argument, a layout.
def someFunc(layout):
    saywhat()
    return otherstuff(layout, 4, 3)
