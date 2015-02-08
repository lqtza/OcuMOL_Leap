import sys
import math
#from pymol import cmd
 
#sys.path.append("/Users/mariusz/code/LeapSDK/lib")
sys.path.append("/Users/lqtza/Downloads/LeapDeveloperKit/LeapSDK/lib")
 
import Leap
from Leap import Matrix, Vector, SwipeGesture
 
class PymolListener(Leap.Listener):
    def __init__(self, *args, **kwargs):
        super(PymolListener, self).__init__(*args, **kwargs)
 
        self.prev_frame = None
        self.do_rotation = False
 
        self.controller = Leap.Controller()
        self.controller.add_listener(self)
 
    def __del__(self):
        self.controller.remove_listener(self)
 
        super(PymolListener, self).__del__()
 
    def on_init(self, controller):
        print "Initialized"

    def on_connect(self, controller):
        print "Connected"
        
        # Enable gestures
        self.controller.enable_gesture(Leap.Gesture.TYPE_SWIPE)

        # Configure gestures
        self.controller.config.set("Gesture.Swipe.MinLenght",10.0)
        self.controller.config.set("Gesture.Swipe.MinVelocity",200)
        self.controller.config.save()

    def on_disconnect(self, controller):
        print "Disconnected"

    def on_exit(self, controller):
        print "Exited"

    def on_frame(self, controller):
        frame = controller.frame()
        print self.do_rotation

        # Two hands and open hand on the leftmost should allow for rotation
        if len(frame.hands) == 2 and frame.hands.leftmost.sphere_radius > 75:
            self.do_rotation = True

        # Two hands and closed hand on the leftmost should allow for translation
        elif len(frame.hands) == 2 and frame.hands.leftmost.sphere_radius < 40:
            self.do_rotation = False
 
        self.update_view(frame,self.do_rotation)
        self.prev_frame = frame
 
    def update_view(self, frame, do_rotation):
        if not self.prev_frame:
            return
 
        #view = list(cmd.get_view())

	'''if len(frame.gestures())>=1:
	    print "point"'''

        #for gest in frame.gestures():
	        #if gest.type is Leap.Gesture.TYPE_SWIPE:
		    #swipe=Leap.SwipeGesture(gest)
		    #start = swipe.start_position
		    #end= swipe.position
		    #mag=math.sqrt((start.x-end.x)**2+(start.y-end.y)**2+(start.z-end.z)**2)
		    #cmd.turn('x',mag/100)
		    #dire=swipe.direction
		    #cmd.rotate([dire.y,dire.x,dire.z],mag/50)
 
        if frame.hands.rightmost.rotation_probability(self.prev_frame) > 0.1 and do_rotation == True:
            rotation_about_x = frame.hands.rightmost.rotation_angle(self.prev_frame,Vector.x_axis)
            rotation_about_y = frame.hands.rightmost.rotation_angle(self.prev_frame,Vector.y_axis)
            rotation_about_z = frame.hands.rightmost.rotation_angle(self.prev_frame,Vector.z_axis)
            cmd.rotate('x'math.,rotation_about_x)
            #m = frame.hands.rightmost.rotation_matrix(self.prev_frame)
            #print m
            #m *= Matrix(Vector(*view[0:3]),
            #            Vector(*view[3:6]),
            #            Vector(*view[6:9]))
            #view[:9] = m.to_array_3x3()
 
        '''if frame.scale_probability(self.prev_frame) > 0.1:
            s = frame.scale_factor(self.prev_frame)
            delta_z = math.log(s) * 100.0
            view[11] += delta_z
            view[15] -= delta_z
            view[16] -= delta_z'''
 
#        cmd.set_view(view)
 
listener = PymolListener()

try:
    sys.stdin.readline()
except KeyboardInterrupt:
    pass
