"""
    A test plugin for PyMOL
"""
import Tkinter
from Tkinter import *
import Pmw
import pymol
import sys

#lazy way to add Rift Mover to PYTHONPATH
sys.path.append("/Users/lqtza/Hacks/OcuMOL_Leap/pymol/")
from oo_prep_and_run import PyMOLViewer

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
                                    buttons = ("Run Rift Only","Run Leap Only", "Run Both"),
                                    title = 'OcuMOL Leap',
                                    command = self.execute)
        self.dialog.withdraw()
        Pmw.setbusycursorattributes(self.dialog.component('hull'))
        w = Tkinter.Label(self.dialog.interior(),
                            text = 'PyMOL Oculus Rift Viewer + Leap Motion Mover\nJeliazko R. Jeliazkov, Max Klein, Henry Lessen, and Mariusz Matyszewski, 2015.\nhttps://github.com/lqtza/OcuMOL_Leap',
                            background = 'black',
                            foreground = 'white',
                            #pady = 20,
                            )
        w.pack(expand = 1, fill = 'both', padx = 4, pady = 4)
        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill='both',expand=1,padx=10,pady=10)

        #Create Oculus Rift Page
        page = self.notebook.add('Rift Visualizer')

        #Create Leap Motion Page
        page = self.notebook.add('Leap Mover')

        #Create About Page
        page = self.notebook.add('About')
        group = Pmw.Group(page, tag_text='About OcuMOL Leap')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)

        about_text = """Magic text goes here."""
        interior_frame = Frame(group.interior())
        bar = Scrollbar(interior_frame,)
        text_holder = Text(interior_frame, yscrollcommand=bar.set,background='#ddddff',font="Times 14")
        bar.config(command=text_holder.yview)
        text_holder.insert(END,about_text)
        text_holder.pack(side=LEFT,expand="yes",fill="both")
        bar.pack(side=LEFT,expand="yes",fill="y")
        interior_frame.pack(expand="yes",fill="both")

        #self.notebook.setnaturalsize()

        self.dialog.show()

    def execute(self,result):
        if result:
            print 'You clicked on, ' + result
            if result == 'Run Rift Only':
                # Rift needs to be connected for this to run.
                # PyMOL will crash if Rift is off or not connected.
                test_viewer = PyMOLViewer()
            elif result == 'Run Leap Only':
                print 'Inner if, should create an object of the LeapMover class.'
            elif result == 'Run Both':
                print 'Inner if, should create a Viewer and Mover object.'
        else:
            self.quit()

    def quit(self):
        self.dialog.destroy()
