#!/usr/bin/env python3

# jiw March 2020
'''Create pypeVue plugins stubs'''

import os
from sys import argv, exit, exc_info, stderr, _getframe
from datetime import datetime
from math import sqrt, pi, cos, sin, asin, atan2
from pypeSim import Point, Post, Cylinder, Layout
from pypePlugins import FunctionList

def tell():
    '''This tell() has a principal role in initializing FunctionList.  It
    creates a list of stubs of functions.  Most or all of these stubs
    should be overridden by functions that appear in modules like
    baseFuncs.py, to implement core features of pypeVue.  Users also
    can override these functions to change how core features work, and
    can also define "user functions" for use in scripts, where they
    can be invoked by U codes in =L sections or by name within =A
    sections.

    If a function name is in several tell's, the one taken last will
    be used.  Process order: Files in pypePlugins/ are taken in
    alphanumerical order.  Next, if pypeVue gets a command line or
    script "plugin=filelist" parameter, that filelist is taken in order.

    '''
    
    FunctionList.fNames = ('addEdge', 'addEdges', 'arithmetic',
                           'autoAdder', 'generatePosts', 'installParams', 'levelAt',
                           'loadScriptFile', 'postTop', 'runScript', 'scriptCyl',
                           'scriptPost', 'thickLet', 'writeCylinders', 'writeLabels',
                           'writePosts')
    return (addEdge, addEdges, arithmetic, autoAdder, generatePosts,
            installParams, levelAt, loadScriptFile, postTop, runScript,
            scriptCyl, scriptPost, thickLet, writeCylinders, writeLabels,
            writePosts)

def saywhat():
    fname = _getframe(1).f_code.co_name
    print (f'{fname:15} @ {os.path.basename(__file__)}')

def addEdge(v,w, layout):      saywhat()
def addEdges(v,w, layout):     saywhat()
def arithmetic(line, xTrace):  saywhat()
def autoAdder(fout):           saywhat()
def generatePosts(code, numberTexts):    saywhat()
def installParams(script):     saywhat()
def levelAt(lev, p):           saywhat()
def loadScriptFile(fiName):    saywhat()
def postTop(p, OP):            saywhat()    
def runScript(scripts):        saywhat()
def scriptCyl(ss, preCyl):     saywhat()
def scriptPost(ss, prePost):   saywhat()
def thickLet(thix):            saywhat()
def writeCylinders(fout, clo, chi, listIt):    saywhat()
def writeLabels(fout):         saywhat()  
def writePosts(fout):          saywhat()

if __name__ == '__main__':
    from datetime import datetime
    date = datetime.today().strftime('%Y-%m-%d  %H:%M:%S')
    print (f'\n\n\n\nThis is {__name__} from {__file__} at {date}\n')
