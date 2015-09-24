import os
import re
from setuptools import find_packages, setup
from setuptools.command.install import install

thisScriptDir = os.path.dirname(os.path.realpath(__file__))

def EnsureNewline(line):
    line = line.rstrip()
    line = line + '\n'
    return line

def ReplaceBtwnLines(fPath, startLine, endLine, txtBlock):
    repRe = re.compile(startLine + '.*' + endLine, re.DOTALL)
#     startLine, endLine = EnsureNewline(startLine), EnsureNewline(endLine)
    
    try:
        with open(fPath, 'r') as f:
            oldTxt = f.read()
    except IOError:
        return False
    if repRe.search(oldTxt):
        newTxt = repRe.sub(startLine + txtBlock + endLine, oldTxt)
        with open(fPath, 'w') as f:
            f.write(newTxt)
        return True
    else:
        return False
    
def AppendBlockToFile(fPath, txtBlock):
    oldTxt = ''
    try:
        with open(fPath, 'r') as f:
            oldTxt = f.read()
    except IOError:
        pass
    with open(fPath, 'w') as f:
        f.write(oldTxt)
        f.write('\n')
        f.write(txtBlock)
    
class CustomInstallCommand(install):
    '''
    customized install command
    '''
    def run(self):
        self._writeLeapConfig()
        self._writePymolrc()
        
    def _writeLeapConfig(self):
        '''
        write config file to set path to leap motion packages
        '''
        if 'LEAPSDK_DIR' not in os.environ:
            raise KeyError
        leapsdk_dir = os.environ.get('LEAPSDK_DIR')
        with open('ocumol/leap_config.py', 'w') as f:
            f.write("leap_path = '%s'" % os.path.join(leapsdk_dir, 'lib'))
        install.run(self)
        
    def _writePymolrc(self):
        pymolrcPath = os.path.expanduser('~/.pymolrc')
        ocumolStartBumper = '#'*4 + 'START_OCUMOL_PLUGIN' + '#'*4
        ocumolEndBumper = '#'*4 + 'END_OCUMOL_PLUGIN' + '#'*4
        with open(os.path.join(thisScriptDir, 'pymolrc'), 'r') as f:
            ocumolPluginTxt = f.read()
        if not ReplaceBtwnLines(pymolrcPath, ocumolStartBumper, ocumolEndBumper, ocumolPluginTxt):
            AppendBlockToFile(pymolrcPath, ocumolPluginTxt)

setup(
    author = 'Max Klein, Jeliazko Jeliazkov, Henry Lessen, Mariusz Matyszewski',
    cmdclass = {'install': CustomInstallCommand},
    description = 'adds VR support to various molecular viewers',
    license = 'apache',
    name = "ocumol",
    packages=find_packages(),
    url='https://github.com/lqtza/OcuMOL_Leap'
    )