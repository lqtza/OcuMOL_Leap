#!/usr/bin/env python

import os
import re
from setuptools import find_packages, setup
from setuptools.command.develop import develop
from setuptools.command.install import install

thisScriptDir = os.path.dirname(os.path.realpath(__file__))

def EnsureNewline(line):
    line = line.rstrip()
    line = line + '\n'
    return line

def ReplaceBtwnLines(fPath, startLine, endLine, txtBlock):
    repRe = re.compile(startLine + '.*' + endLine, re.DOTALL)
    replTxt = EnsureNewline(startLine) + EnsureNewline(txtBlock) + endLine
    
    try:
        with open(fPath, 'r') as f:
            oldTxt = f.read()
    except IOError:
        return False
    if repRe.search(oldTxt):
        newTxt = repRe.sub(replTxt, oldTxt)
        with open(fPath, 'w') as f:
            f.write(newTxt)
        return True
    else:
        return False
    
def AppendBlockToFile(fPath, startLine, endLine, txtBlock):
    oldTxt = ''
    newTxt = EnsureNewline(startLine) + EnsureNewline(txtBlock) + endLine
    try:
        with open(fPath, 'r') as f:
            oldTxt = f.read()
            oldTxt = EnsureNewline(oldTxt) + '\n'
    except IOError:
        pass
    with open(fPath, 'w') as f:
        f.write(oldTxt)
        f.write(newTxt)
    
class CustomSetupCommand:
    '''
    customized setup command base class
    meant to be used in a subclass that also inherits either setuptools.command.install.install or .develop
    '''
    leapSDKEnvVar = 'LEAPSDK_DIR'
    
    def run(self):
        self._writeLeapConfig()
        self._writePymolrc()
        
    def _writeLeapConfig(self):
        '''
        write config file to set path to leap motion packages
        '''
        if self.leapSDKEnvVar not in os.environ:
            raise EnvironmentError('%s%s' % ('You need to set the %s environment variable before installing ocumol, ' % self.leapSDKEnvVar,
                                             'e.g. you could run `LEAPSDK_DIR=/usr/local/LeapSDK pip install ocumol`'))
        leapSDKDir = os.environ.get(self.leapSDKEnvVar)
        with open('ocumol/leapConfig.py', 'w') as f:
            f.write("leapPath = '%s'" % os.path.join(leapSDKDir, 'lib'))
         
    def _writePymolrc(self):
        pymolrcPath = os.path.expanduser('~/.pymolrc.py')
        ocumolStartBumper = '#'*4 + 'START_OCUMOL_PLUGIN' + '#'*4
        ocumolEndBumper = '#'*4 + 'END_OCUMOL_PLUGIN' + '#'*4
        with open(os.path.join(thisScriptDir, 'pymolrc'), 'r') as f:
            ocumolPluginTxt = f.read()
        if not ReplaceBtwnLines(pymolrcPath, ocumolStartBumper, ocumolEndBumper, ocumolPluginTxt):
            AppendBlockToFile(pymolrcPath, ocumolStartBumper, ocumolEndBumper, ocumolPluginTxt)

class CustomDevelopCommand(CustomSetupCommand, develop):
    def run(self):
        CustomSetupCommand.run(self)
        develop.run(self)

class CustomInstallCommand(CustomSetupCommand, install):
    def run(self):
        CustomSetupCommand.run(self)
        install.run(self)

setup(
    author = 'Max Klein, Jeliazko Jeliazkov, Henry Lessen, Mariusz Matyszewski',
    cmdclass = {'develop': CustomDevelopCommand,'install': CustomInstallCommand},
    description = 'adds VR support to various molecular viewers',
    license = 'apache',
    name = "ocumol",
    packages=find_packages(),
    url='https://github.com/lqtza/OcuMOL_Leap'
    )
