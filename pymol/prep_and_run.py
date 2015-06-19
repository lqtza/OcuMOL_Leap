from pymol import cmd
from pymol import util
import time
import numpy as np
import sys
sys.path.append("/Users/lqtza/Hacks/ocudump/build/src/cython")
from ocudump import Ocudump

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

def set_origin_at_camera():
    '''
    places PyMol's origin (center of rotation) at the camera position, allowing for more natural interaction with PyMol via the head rotations in the Rift.
    Should probably be called every time the camera makes a translational move
    explanation:
        what we want: the displacement vector of the origin of cameraspace (which is also the camera's position) relative to the origin of modelspace, in modelspace coordinates.
                      If we feed this vector into origin(position=x), the origin will be placed appropriately
        how to get it: cameraspace and modelspace are two simple 3D cartesian spaces. PyMol's origin of rotation is guaranteed to be a cooincident point in the two spaces.
                       Based on the info in get_view, we can figure out first the rotation matrix (model->camera). Once the rotation has been applied to modelspace (or it's inverse
                       to cameraspace), the axes of the two spaces will be aligned. The displacement between the two spaces' origins can then be calculated via a simple 
                       vector difference of the location of the origin in each space (get_view[9:12] for cameraspace, get_view[12:15] for modelspace)
    '''
    view = np.array(cmd.get_view())
    
    # simple version
#     rot = view[0:9].reshape((3,3))
#     camera = view[9:12]
#     model = view[12:15]
#     cmd.origin(position=model - camera.dot(rot.T))

    # faster (?) version
    cmd.origin(position=view[12:15] - view[9:12].dot(view[0:9].reshape((3,3)).T))

set_origin_at_camera()
while True:
    #Rift support
    o.getPose()
    #print o.pose [pitch (x), yaw (z), roll (y)]
    currf = o.pose

    rot_diff = [currf[0]-prevrf[0],currf[1]-prevrf[1],currf[2]-prevrf[2]]
    pos_diff = [curff[3], currf[4], curff[5]]
    #diff = currf

    cmd.turn('x',-rot_diff[0]*25)
    cmd.turn('y',-rot_diff[1]*25)
    cmd.turn('z',-rot_diff[2]*25)

    cmd.translate(pos_diff)
    
    time.sleep(1/float(trackingRefresh))

    prevrf = currf
