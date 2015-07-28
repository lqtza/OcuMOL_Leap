import os
from setuptools import find_packages, setup
from setuptools.command.install import install

class CustomInstallCommand(install):
    '''
    install command customized to set path to leap motion packages
    '''
    def run(self):
        leap_lib_dir = os.environ.get('LEAP_LIB_DIR')
        with open('ocumol/leap_config.py', 'w') as f:
            f.write('leap_path = %s' % leap_lib_dir)
        install.run(self)

setup(
    author = 'Max Klein, Jeliazko Jeliazkov, Henry Lessen, Mariusz Matyszewski',
    cmdclass = {'install': CustomInstallCommand},
    description = 'adds VR support to various molecular viewers',
    license = 'apache',
    name = "ocumol",
    packages=find_packages(),
    url='https://github.com/lqtza/OcuMOL_Leap'
    )