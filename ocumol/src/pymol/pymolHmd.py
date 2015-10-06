from pymol import cmd
from pymol import util
import time
import numpy as np
import os
import threading

pymolHmdScript = os.path.realpath(__file__)

from ocudump import Ocudump, OcudumpDebug
from ocumol.src.helper.transformations import euler_matrix

class PymolHmd(threading.Thread):
    def __init__(self, debugMode=False, editMolecule=False, naturalRotation=True, pdb=''):
        threading.Thread.__init__(self)
        
        # flag that indicates if the instance is in active control of a pymol session
        self.running = False
        
        # flag that indicates if we should run in debug mode (indifferent to presence of Rift)
        self.debugMode = debugMode
        
        # pdb to load, if any
        self.pdb = pdb
        
        # set pymol editing mode on/off
        self.editMolecule = editMolecule

        # define how view rotation is updated
        self.naturalRotation = naturalRotation
        
        # whether to run in fullscreen mode or not
        self.fullscreen = False

        # set previous head pitch, yaw, roll
        self.previousPose = [0,0,0]

        # set previous head x,y,z
        self.previousXYZ = [0,0,0]

        # oculus tracking data refresh rate, in Hz
        self.trackingRefresh = 30

        # set scaling factor for translation
        self.translationScaling = 1

        # set scaling factor for rotation
        self.rotationScaling = 25
    
    def initCamera(self):
        self.ocudump.getPose()

        # this makes the camera surprisingly jumpy at first
        # max, do we need it?
        # self.basePose = self.ocudump.pose
    
    def initOcudump(self):
        if self.debugMode:
            self.ocudump = OcudumpDebug()
        else:
            self.ocudump = Ocudump()
        return self.ocudump.init()
    
    def initPymol(self):
        # move slab away to give a comfortable viewing area
        cmd.clip('move',10)

        # set stereo mode
        cmd.set('stereo_mode',3)
        cmd.stereo()

        # get rid of gui
        cmd.set('internal_gui',0)

        # set resolution to HD for rift
        cmd.viewport(1920,1080)

        # full screen?
        if self.fullscreen: 
            cmd.full_screen('on')

        # load pdb
        if len(self.pdb) == 4:
            cmd.fetch(self.pdb,async=0)
        
        # color each chain differently
        if self.editMolecule:
            cmd.hide('everything','all')
            cmd.show('cartoon','all')
            util.cbc()

        # set origin at camera
        self.setOriginAtCamera()

    def setOriginAtCamera(self):
        view = np.array(cmd.get_view())

        # concise version
        cmd.origin(position=view[12:15] - view[9:12].dot(view[0:9].reshape((3,3)).T))

        # simple version
        #     rot = view[0:9].reshape((3,3))
        #     camera = view[9:12]
        #     model = view[12:15]
        #     cmd.origin(position=model - camera.dot(rot.T))


    def visualize(self):
        # load pdb into pymol and set up view
        self.initPymol()
        
        # again, do we need this function?
        self.initCamera()
        
        # listen ...
        while True:

            # Rift support
            self.ocudump.getPose()

            #print ocudump.pose [pitch (x), yaw (z), roll (y)]
            pose = self.ocudump.pose

            # calculate head rotation
            x_rot = pose[0] - self.previousPose[0]
            y_rot = pose[1] - self.previousPose[1]
            z_rot = pose[2] - self.previousPose[2]

            # update view
            cmd.turn('x',x_rot*self.rotationScaling)
            cmd.turn('y',y_rot*self.rotationScaling)
            cmd.turn('z',z_rot*self.rotationScaling)

            # rift head tracking
            if self.ocudump.positionTracked:
                try:
                    # move head if sensor detects head motion
                    # not tested
                    cmd.move('x', pose[4] - self.prev_xyz[4])
                    cmd.move('y', pose[5] - self.prev_xyz[5])
                    cmd.move('z', pose[6] - self.prev_xyz[6])
                except AttributeError:
                    self.prev_xzy = pose[3:6]

            # update camera for "natural rotation"
            if self.naturalRotation:
                # this breaks the code... not sure why
                self.setOriginAtCamera()
                pass

            # record previous pose for accurate rotation diff
            self.previousPose = pose

            time.sleep(1/float(self.trackingRefresh))

    def run(self):
        self.running = True
        if not self.initOcudump():
            return False
        self.visualize()

if __name__ == '__main__':
    # do not run if imported as module
    pcls=PymolHmd()
    # using the start command from the inherited thread class
    pcls.start()
