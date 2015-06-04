from pymol import cmd
from pymol import util
import time
import sys
sys.path.append("/Users/lqtza/Hacks/ocudump/build")
from src.cython.ocudump import Ocudump

# oculus tracking data refresh rate, in Hz
trackingRefresh = 60

#initialize
o = Ocudump()
#load pdb
cmd.fetch('3ceg',async=0)
#color each chain differently
cmd.hide('everything','all')
cmd.show('cartoon','all')
#cmd.show('surface','all')
#cmd.set('transparency',0.5)
util.cbc()

#move slab away to give a comfortable viewing area
cmd.clip('move',10)

#set stereo mode
cmd.set('stereo_mode',3)
cmd.stereo()

#get rid of gui
cmd.set('internal_gui',0)

#set resolution to HD for rift
cmd.viewport(1920,1080)
 
#full screen
cmd.full_screen('on')

prevrf = [0,0,0]
 
while True:
    #Rift support
    o.getPose()
    #print o.pose [pitch (x), yaw (z), roll (y)]
    currf = o.pose

    diff = [currf[0]-prevrf[0],currf[1]-prevrf[1],currf[2]-prevrf[2]]
    #diff = currf

    cmd.turn('x',diff[0]*75)
    cmd.turn('y',diff[1]*75)
    cmd.turn('z',diff[2]*75)
    
    time.sleep(1/float(trackingRefresh))

    prevrf = currf
