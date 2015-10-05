import sys
import math
import time
from pymol import cmd
import os
import Leap
from Leap import Matrix, Vector, CircleGesture

class PymolListener(Leap.Listener):
    def __init__(self, *args, **kwargs):
        super(PymolListener, self).__init__(*args, **kwargs)

        self.prev_frame = None
        self.view_do_rotation = False
        self.view_do_translation = False
        self.mode = 'view' # this should be binary edit or view

        self.controller = Leap.Controller()
        self.controller.add_listener(self)
        self.controller.set_policy(Leap.Controller.POLICY_BACKGROUND_FRAMES)

        self.circom=1

    def __del__(self):
        self.controller.remove_listener(self)
        super(PymolListener, self).__del__()

    # inheriting from Leap.Listener and extending
    def on_init(self, controller):
        print "Initialized"

    # inheriting from Leap.Listener and extending
    def on_connect(self, controller):
        print "Connected"

        # Enable gestures
        self.controller.enable_gesture(Leap.Gesture.TYPE_CIRCLE)
        self.controller.enable_gesture(Leap.Gesture.TYPE_SWIPE)
        self.controller.config.set("Gesture.Swipe.MinVelocity",500)
        self.controller.config.save()

    # inheriting from Leap.Listener and extending
    def on_disconnect(self, controller):
        print "Disconnected"

    # inheriting from Leap.Listener and extending
    def on_exit(self, controller):
        print "Exited"

    # inheriting from Leap.Listener and extending
    def on_frame(self, controller):
        frame = controller.frame()
        #print self.view_do_rotation

        if self.mode == 'view':
            # Two hands and open hand on the leftmost should allow for rotation
            if len(frame.hands) == 2 and frame.hands.leftmost.sphere_radius > 75:
                self.view_do_rotation = True

            # Two hands and closed hand on the leftmost should allow for translation
            elif len(frame.hands) == 2 and frame.hands.leftmost.sphere_radius < 40:
                self.view_do_translation = True

    	    else:
                self.view_do_rotation = False
                self.view_do_translation = False

        self.update_view(frame, self.view_do_rotation, self.view_do_translation)
        self.prev_frame = frame

    def update_view(self, frame, do_rotation, do_translation):
        if not self.prev_frame:
            return

        # check what mode to set, also make directional in future
        # switch modes, TODO: deprecate?
        if len(frame.hands) == 2:
            for gest in frame.gestures():
                if gest.type is Leap.Gesture.TYPE_SWIPE:
                    if Leap.SwipeGesture(gest).direction.y > 0.5 and gest.duration_seconds > 0.15:
                        time.sleep(0.3)
                        if self.mode == 'view':
                            self.mode = 'edit'
                            cmd.bg_color("white")
                        else:
                            self.mode = 'view'
                            cmd.bg_color("black")
                        do_rotation = False
                        do_translation = False
                        print 'Changing mode to: ' + self.mode
                        time.sleep(0.6)
                        break # why?

        for gest in frame.gestures():
            if gest.type is Leap.Gesture.TYPE_CIRCLE:
                circle=Leap.CircleGesture(gest)

		if math.floor(circle.progress)>=1:
            self.circom=0

        if self.circom==0 and len(frame.gestures())==0:
            self.circom=1

	    if len(frame.hands)==1:
            cmd.center("all")
	    elif len(frame.hands)==2:
            cmd.orient("all")

        if frame.hands.rightmost.rotation_probability(self.prev_frame) > 0.1 and do_rotation == True:
            #print 'rotating'
            rotation_about_x = frame.hands.rightmost.rotation_angle(self.prev_frame,Vector.x_axis)
            rotation_about_y = frame.hands.rightmost.rotation_angle(self.prev_frame,Vector.y_axis)
            rotation_about_z = frame.hands.rightmost.rotation_angle(self.prev_frame,Vector.z_axis)
	    #print rotation_about_x, rotation_about_y, rotation_about_z
            cmd.rotate('x',rotation_about_x*100)
            cmd.rotate('y',rotation_about_y*100)
            cmd.rotate('z',rotation_about_z*100)
            #m = frame.hands.rightmost.rotation_matrix(self.prev_frame)
            #print m
            #m *= Matrix(Vector(*view[0:3]),
            #            Vector(*view[3:6]),
            #            Vector(*view[6:9]))
            #view[:9] = m.to_array_3x3()

        elif frame.hands.rightmost.translation_probability(self.prev_frame) > 0.1 and do_translation == True:
            translation = frame.hands.rightmost.translation(self.prev_frame)
            #print translation.to_float_array()
            cmd.translate(translation.to_float_array())

        elif self.mode == 'edit' and len(frame.hands) == 1:
            if frame.hands[0].is_right:
                translation = frame.hands[0].translation(self.prev_frame)
                cmd.translate(translation.to_float_array(),'EcZapA.clean')
            elif frame.hands[0].is_left:
                translation = frame.hands[0].translation(self.prev_frame)
                cmd.translate(translation.to_float_array(),'EcFtsZ_AB.clean')

        '''view = list(cmd.get_view())

        if frame.scale_probability(self.prev_frame) > 0.1 and len(frame.hands)==1:
            s = frame.scale_factor(self.prev_frame)
            delta_z = math.log(s) * 100.0
            view[11] += delta_z
            view[15] -= delta_z
            view[16] -= delta_z
	    cmd.set_view(view)'''

if __name__ == '__main__':
    listener = PymolListener()
