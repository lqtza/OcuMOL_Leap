from pymol import cmd
from pymol import util
import time
import numpy as np
import os
import re
import threading

pymolHmdScript = os.path.realpath(__file__)

from ocudump import Ocudump, OcudumpDebug
from ocumol.src.helper.transformations import euler_matrix

class PymolHmd(threading.Thread):
    poseCoordDict = {'xrot':0, 'yrot':1, 'zrot':2, 
                     'x':3, 'y':4, 'z':5}
    
    def __init__(self, debugMode=False, editMolecule=False, naturalRotation=True, pdb=''):
        # set up the Thread stuff
        threading.Thread.__init__(self)
        
        # pre-initialize the ocudumpDebug, in case we want to set animations
        self.ocudumpDebug = OcudumpDebug()
        
        # initialize previousOrientation (pitch, yaw, roll) and previousPostion (x, y, z)
        self.previousOrientation = [0,0,0]
        self.previousPosition = [0,0,0]
        
        # internal flags
        self.fullscreen = False # whether to run in fullscreen mode or not
        self.running = False    # indicates if this instance's run method has been called
        
        # internal options
        self.rotationScaling = 25   # set scaling factor for rotation
        self.trackingRefresh = 30   # oculus tracking data refresh rate, in Hz
        self.translationScaling = 1 # set scaling factor for translation
        
        # user flags
        self.debugMode = debugMode              # should we run in debug mode (ie indifferent to presence of Rift)?
        self.editMolecule = editMolecule        # set pymol editing mode on/off
        self.naturalRotation = naturalRotation  # define how view rotation is updated
        
        # user options
        self.pdb = pdb  # pdb to load, if any
    
    def initAnimateElem(self, poseCoordName, minim, maxim, period):
        self.ocudumpDebug.initAnimateElem(self.poseCoordDict[poseCoordName], minim, maxim, period)
    
    def initCamera(self):
        self.ocudump.getPose()

        # this makes the camera surprisingly jumpy at first
        # max, do we need it?
        # self.basePose = self.ocudump.pose
    
    def initOcudump(self):
        if self.debugMode:
            self.ocudump = self.ocudumpDebug
        else:
            self.ocudump = Ocudump()
        return self.ocudump.init()
    
    def initPymol(self):

        # set stereo mode
        cmd.set('stereo_mode', 3)
        cmd.stereo()

        # get rid of gui
        cmd.set('internal_gui', 0)

        # set resolution to HD for rift
        cmd.viewport(1920, 1080);

        # full screen?
        if self.fullscreen: 
            cmd.full_screen('on')

        # load pdb
        if len(self.pdb)==4 and self.pdb==re.match('[a-zA-Z0-9]*', self.pdb).group(0):
            cmd.fetch(self.pdb, async=0)
        
        cmd.hide('everything', 'all')
        cmd.show('cartoon', 'all')
        
        # color each chain differently
        if self.editMolecule:
            util.cbc()
        
        # set origin at camera
        if self.naturalRotation:
            self.setOriginAtCamera()
    
        # move slab away to give a comfortable viewing area
        cmd.clip('move', 10)
        
        # jiggle the camera to deal with blank screen error
        cmd.move('z', .0001)
    
    def reinitPymol(self, pdb):
        if len(self.pdb)==4 and self.pdb==re.match('[a-zA-Z0-9]*', self.pdb).group(0):
            cmd.delete(self.pdb)
        if len(pdb)==4 and pdb==re.match('[a-zA-Z0-9]*', pdb).group(0):
            self.pdb = pdb
            cmd.fetch(self.pdb, async=0)
        
        cmd.hide('everything', 'all')
        cmd.show('cartoon', 'all')
        
        # color each chain differently
        if self.editMolecule:
            util.cbc()
        
        # set origin at camera
        if self.naturalRotation:
            self.setOriginAtCamera()
        
        # move slab away to give a comfortable viewing area
        cmd.clip('move', 10)
        
    def setOriginAtCamera(self):
        view = np.array(cmd.get_view())

        # concise version
        cmd.origin(position=view[12:15] - view[9:12].dot(view[0:9].reshape((3,3)).T))

        # simple version
        #     rot = view[0:9].reshape((3,3))
        #     camera = view[9:12]
        #     model = view[12:15]
        #     cmd.origin(position=model - camera.dot(rot.T))

    def setOriginAtMolecule(self):
        view = np.array(cmd.get_view())
        # TODO: replace this with something that sets the origin at the center of mass
        cmd.origin(position=view[9:12])
    
    def visualize(self):
        # load pdb into pymol and set up view
        self.initPymol()
        
        # TODO: check if this is necessary
        self.initCamera()
        
        # listen ...
        while True:
            # Rift support
            self.ocudump.getPose()

            #print ocudump.pose [pitch (x), yaw (z), roll (y)]
            pose = self.ocudump.pose

            # calculate head rotation
            xRot = pose[0] - self.previousOrientation[0]
            yRot = pose[1] - self.previousOrientation[1]
            zRot = pose[2] - self.previousOrientation[2]

            # update view
            cmd.turn('x', xRot*self.rotationScaling)
            cmd.turn('y', yRot*self.rotationScaling)
            cmd.turn('z', zRot*self.rotationScaling)

            # rift head tracking
            if self.ocudump.positionTracked:
                # move head if sensor detects head motion
                # not tested
                cmd.move('x', pose[3] - self.previousPosition[0])
                cmd.move('y', pose[4] - self.previousPosition[1])
                cmd.move('z', pose[5] - self.previousPosition[2])

            # update camera for "natural rotation"
            if self.naturalRotation:
                self.setOriginAtCamera()

            # record previous pose for accurate rotation diff
            self.previousOrientation = pose[0:3]

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
