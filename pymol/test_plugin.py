"""
    A test plugin for PyMOL
"""
from Tkinter import *
import Pmw
import pymol

def __init__(self):
    """
    Add OcuMOL Leap to the PyMOL Plugins menu.
    """

    self.menuBar.addmenuitem('Plugin', 'command', 'Launch OcuMOL Leap',
                                label='OcuMOL Leap',
                                command= lambda s=self: OcuMOLLeap(s))

### OcuMOLLeap code ###
class OcuMOLLeap:
    def __init__(self,app):
        self.parent = app.root

        self.dialog = Pmw.Dialog(self.parent,
                                    buttons = ("Test Button1","Test Button2", "Test Button3"),
                                    title = 'OcuMOL Leap',
                                    command = self.execute)

    def execute(self,result):
        if result:
            print 'You clicked on, ' + result
        else:
            self.quit()

    def quit(self):
        self.dialog.destroy()
