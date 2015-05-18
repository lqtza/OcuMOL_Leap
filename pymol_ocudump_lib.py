from pymol import cmd
import time
import sys
sys.path.append("/Users/jjeliazkov/Downloads/ocudump/build")
from src.cython.ocudump_cython import Ocudump
#sys.path.append("/Users/jjeliazkov/Downloads/python-rift/")

#initialize
o = Ocudump()

prevrf = [0,0,0]
 
while True:
    #Rift support
    o.getPose()
    print o.pose
    currf = o.pose

    diff = [currf[0]-prevrf[0],currf[1]-prevrf[1],currf[2]-prevrf[2]]

    cmd.turn('x',diff[0]*100)
    cmd.turn('y',diff[1]*100)
    cmd.turn('z',diff[2]*100)
    
    time.sleep(0.6)

    prevrf = currf
