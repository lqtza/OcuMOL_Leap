import sys

from ocumol.leap_config import leap_path
if leap_path not in sys.path:
    sys.path.append(leap_path)

#from ocumol.src.hands.leap_only import PymolListener
# from ocumol.src.pymol.pymolHmd import PymolHmd, pymolHmdScript