import sys

from ocumol.leapConfig import leapPath
if leapPath not in sys.path:
    sys.path.append(leapPath)

#from ocumol.src.hands.leap_only import PymolListener
# from ocumol.src.pymol.pymolHmd import PymolHmd, pymolHmdScript