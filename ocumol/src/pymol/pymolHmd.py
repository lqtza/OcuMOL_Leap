from pymol import cmd
from pymol import util
import time
import numpy as np
import os
import threading

pymolHmdScript = os.path.realpath(__file__)

from ocudump import Ocudump
from ocumol.tranformations import euler_matrix

class PymolHmd(threading.Thread):
    
    def __init__(self, naturalRotation=True):
        threading.Thread.__init__(self)
        # oculus tracking data refresh rate, in Hz
        
        self.init_pymol("3ceg")
        self.naturalRotation = naturalRotation
        self.ocudump = Ocudump()
        self.prev_frame_pos = [0,0,0]
        self.rotation_scaling = 25
        self.tracking_refresh = 60
        self.translation_scaling = 1
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
    
    def init_camera(self):
        self.ocudump.getPose()
        self.base_pose = self.ocudump.pose
    
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
        self.init_camera()
        while True:
            #Rift support
            self.ocudump.getPose()
            #print ocudump.pose [pitch (x), yaw (z), roll (y)]
            pose = self.ocudump.pose
            #prevrf = self.prev_frame_pos

            x_rot = pose[0] - self.base_pose[0]
            y_rot = pose[1] - self.base_pose[1]
            z_rot = pose[2] - self.base_pose[2]
            
            cmd.turn('x',x_rot*self.rotation_scaling)
            cmd.turn('y',y_rot*self.rotation_scaling)
            cmd.turn('z',z_rot*self.rotation_scaling)
            
            if self.ocudump.positionTracked:
                try:
                    cmd.move('x', pose[4] - self.prev_xyz[0])
                    cmd.move('y', pose[5] - self.prev_xyz[1])
                    cmd.move('z', pose[6] - self.prev_xyz[2])
                except AttributeError:
                    self.prev_xzy = pose[3:6]
            
                
#             rot_diff = [currf[0]-prevrf[0],currf[1]-prevrf[1],currf[2]-prevrf[2]]
#             pos_diff = [currf[3], currf[4], currf[5]]
            #diff = currf

#             cmd.turn('x',-rot_diff[0]*25)
#             cmd.turn('y',-rot_diff[1]*25)
#             cmd.turn('z',-rot_diff[2]*25)
            
            

            #cmd.translate(pos_diff)
            
            if self.naturalRotation:
                self.set_origin_at_camera()
            
            time.sleep(1/float(self.tracking_refresh))

#             self.prev_frame_pos = currf

    def run(self):
        self.visualize()

if __name__ == '__main__':
    #do not run if imported as module
    pcls=PymolHmd()
    pcls.start()
