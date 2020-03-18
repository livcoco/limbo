# jiw March 2020
'''Test-program for pypeVue plugins setup'''

from pypeSim import Point, Post, Cylinder, Layout
from pypePlugins import FunctionList
from sys import _getframe
import os.path
def saywhat():
    fname = _getframe(1).f_code.co_name
    print (f'{fname:15} @ {os.path.basename(__file__)}')

def postTop(p, OP):            saywhat()    
def runScript(scripts):        saywhat()
def scriptCyl(ss, preCyl):     saywhat()
def scriptPost(ss, prePost):   saywhat()
def thickLet(thix):
    saywhat()
    ref = FunctionList
    print (f'len(ref.fNames) = {len(ref.fNames)},  len(ref.fList) = {len(ref.fList)},  len(ref.fDict) = {len(ref.fDict)}')

def tell():
    return (postTop, runScript, scriptCyl, scriptPost, thickLet)
