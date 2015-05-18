from pymol import cmd
from rift import PyRift
import sys
import time
#sys.path.append("/Users/jjeliazkov/Downloads/OcuMOL_Leap/")
#sys.path.append("/Users/jjeliazkov/Downloads/python-rift/")
sys.path.append("/Users/lqtza/OcuMOL_Leap/")
import transformations as trans


 #initialize
circom=1
rft=PyRift() #breaks macpymol
prevrf = [0,0,0]
currf = [0,0,0,0]
 
while True:
    #Rift support
    rft.poll()
    currf = rft.rotation
    currf = trans.euler_from_quaternion([currf[3],currf[0],currf[1],currf[2]])

    diff = [currf[0]-prevrf[0],currf[1]-prevrf[1],currf[2]-prevrf[2]]

    cmd.turn('x',diff[0]*100)
    cmd.turn('y',diff[1]*100)
    cmd.turn('z',diff[2]*100)

    prevrf = currf
