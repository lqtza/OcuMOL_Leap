from pymol import cmd
from pymol import util
import time
import numpy as np
import os
import threading

pymolViewerThreadedScript = os.path.realpath(__file__)

from ocudump import Ocudump

class PyMOLViewerThreaded(threading.Thread):
    
    def __init__(self):
        threading.Thread.__init__(self)
        # oculus tracking data refresh rate, in Hz
        self.tracking_refresh = 60
        self.prev_frame_pos = [0,0,0]
        self.ocudump = Ocudump()
        self.init_pymol("3ceg")
        #self.visualize()
    
    def init_pymol(self,pdb='3ceg'):
        #load pdb
        cmd.fetch(pdb,async=0)
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
        #cmd.full_screen('on')

        #set origin at camera
        self.set_origin_at_camera()

    def set_origin_at_camera(self):
        view = np.array(cmd.get_view())
        # faster (?) version
        cmd.origin(position=view[12:15] - view[9:12].dot(view[0:9].reshape((3,3)).T))
            
        # simple version
        #     rot = view[0:9].reshape((3,3))
        #     camera = view[9:12]
        #     model = view[12:15]
        #     cmd.origin(position=model - camera.dot(rot.T))
        

    def visualize(self):
        while True:
            #Rift support
            self.ocudump.getPose()
            #print ocudump.pose [pitch (x), yaw (z), roll (y)]
            currf = self.ocudump.pose
            prevrf = self.prev_frame_pos

            rot_diff = [currf[0]-prevrf[0],currf[1]-prevrf[1],currf[2]-prevrf[2]]
            pos_diff = [currf[3], currf[4], currf[5]]
            #diff = currf

            cmd.turn('x',-rot_diff[0]*25)
            cmd.turn('y',-rot_diff[1]*25)
            cmd.turn('z',-rot_diff[2]*25)

            #cmd.translate(pos_diff)
            
            time.sleep(1/float(self.tracking_refresh))

            self.prev_frame_pos = currf

    def run(self):
        self.visualize()

#if __name__ == '__main__':
    #do not run if imported as module
# pcls=PyMOLViewerThreaded()
# pcls.start()
