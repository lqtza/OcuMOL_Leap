"""
    A test plugin for PyMOL
"""
import Tkinter
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

        #self.notebook.setnaturalsize()

        self.dialog.show()

    def execute(self,result):
        if result:
            print 'You clicked on, ' + result
        else:
            self.quit()

    def quit(self):
        self.dialog.destroy()
