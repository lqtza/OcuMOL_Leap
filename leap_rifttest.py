import sys
import math
from pymol import cmd
from rift import PyRift

import transformations as tran
 
sys.path.append("/Users/mariusz/code/LeapSDK/lib")
 
import Leap
from Leap import Matrix, Vector, SwipeGesture
 
class PymolListener(Leap.Listener):
    def __init__(self, *args, **kwargs):
        super(PymolListener, self).__init__(*args, **kwargs)
 
        self.prev_frame = None
 
        self.controller = Leap.Controller()
        self.controller.add_listener(self)
	self.controller.enable_gesture(Leap.Gesture.TYPE_SWIPE)
	self.controller.config.set("Gesture.Swipe.MinLenght",10.0)
	self.controller.config.set("Gesture.Swipe.MinVelocity",200)
	self.controller.config.save()
	self.rf = PyRift()
	self.prevrf = [0,0,0]
	self.currf = [0,0,0]
 
    def __del__(self):
        self.controller.remove_listener(self)
 
        super(PymolListener, self).__del__()
 
    def update_view(self, frame):
        if not self.prev_frame:
            return
 
        view = list(cmd.get_view())

	'''if len(frame.gestures())>=1:
	    print "point"'''

	for gest in frame.gestures():
	    if gest.type is Leap.Gesture.TYPE_SWIPE:
		swipe=Leap.SwipeGesture(gest)
		start = swipe.start_position
		end= swipe.position
		mag=math.sqrt((start.x-end.x)**2+(start.y-end.y)**2+(start.z-end.z)**2)
		#cmd.turn('x',mag/100)
		dire=swipe.direction
		cmd.rotate([dire.y,dire.x,dire.z],mag/50)

	self.rf.poll()
	self.currf = self.rf.rotation
	self.currf = tran.euler_from_quaternion([self.currf[3],self.currf[0],self.currf[1],self.currf[2]])

	diff = [self.currf[0]-self.prevrf[0], self.currf[1]-self.prevrf[1],self.currf[2]-self.prevrf[2]]

	#print diff

	cmd.turn('x',diff[0]*100)
	cmd.turn('y',diff[1]*100)
	cmd.turn('z',diff[2]*100)

	self.prevrf = self.currf
 
        '''if frame.rotation_probability(self.prev_frame) > 0.1:
            m = frame.rotation_matrix(self.prev_frame)
            m *= Matrix(Vector(*view[0:3]),
                        Vector(*view[3:6]),
                        Vector(*view[6:9]))
            view[:9] = m.to_array_3x3()'''
 
        '''if frame.scale_probability(self.prev_frame) > 0.1:
            s = frame.scale_factor(self.prev_frame)
            delta_z = math.log(s) * 100.0
            view[11] += delta_z
            view[15] -= delta_z
            view[16] -= delta_z'''
 
#        cmd.set_view(view)
 
    def on_frame(self, controller):
        frame = controller.frame()
 
        self.update_view(frame)
 
        self.prev_frame = frame
 
listener = PymolListener()
