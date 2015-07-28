import os
from setuptools import find_packages, setup
from setuptools.command.install import install

class CustomInstallCommand(install):
    '''
    install command customized to set path to leap motion packages
    '''
    def run(self):
        leapsdk_dir = os.environ.get('LEAPSDK_DIR')
        with open('ocumol/leap_config.py', 'w') as f:
            f.write("leap_path = '%s'" % os.path.join(leapsdk_dir, 'lib'))
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