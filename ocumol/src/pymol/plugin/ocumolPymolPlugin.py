"""
    An OcuMOL Leap plugin for PyMOL. To install go to the Plugin menu in PyMOL
    and under the Install New Plugin tab click Choose file... and select this file.
    Set the environmental variable, OCUMOLPATH, to point to the repository directory:
    export OCUMOLPATH=/Users/lqtza/Hacks/OcuMOL_Leap
"""
import Tkinter
from Tkinter import *
import Pmw
import pymol
import sys
import os
from pymol import cmd

from ocumol.src.pymol.pymolHmd import PymolHmd#, PymolListener

#hard coded paths... just don't
# sys.path.append("/Users/tel/git/OcuMOL_Leap/pymol")
# sys.path.append("/Users/tel/git/OcuMOL_Leap/hands")
#sys.path.append("/Users/lqtza/Hacks/LeapDeveloperKit/LeapSDK/lib")

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
        group = Pmw.Group(page, tag_text = 'Oculus Rift Visualizer')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)

        self.radiobuttons_rot = Pmw.RadioSelect(group.interior(),buttontype = 'radiobutton',orient = 'vertical',labelpos = 'w',)
        self.radiobuttons_rot.add('Natural Rotation')
        self.radiobuttons_rot.add('Molecule Rotation')
        self.radiobuttons_rot.setvalue('Natural Rotation')
        self.radiobuttons_rot.pack(fill='x',padx=4,pady=1)

        self.radiobuttons_mol = Pmw.RadioSelect(group.interior(),buttontype = 'radiobutton',orient = 'vertical',labelpos = 'w',label_text = 'Edit Molecule?',)
        self.radiobuttons_mol.add('Yes')
        self.radiobuttons_mol.add('No')
        self.radiobuttons_mol.setvalue('No')
        self.radiobuttons_mol.pack(fill='x',padx=4,pady=1)

        self.pdb_text = Pmw.EntryField(group.interior(),labelpos='w',label_text='Load PDB ID:',command=None,)
        self.pdb_text.pack(fill='x', padx=4, pady=1)

        self.stereo_shift_text = Pmw.EntryField(group.interior(),labelpos='w',label_text='Stereo shift:',value='5.24',validate = {'validator' : 'real'},modifiedcommand=self.changed,)
        self.stereo_shift_text.pack(fill='x', padx=4, pady=1)

        self.stereo_angle_text = Pmw.EntryField(group.interior(),labelpos='w',label_text='Stereo angle:',value='2.0',validate = {'validator' : 'real'},modifiedcommand=self.changed,)
        self.stereo_angle_text.pack(fill='x', padx=4, pady=1)

        #Create Leap Motion Page
        page = self.notebook.add('Leap Mover')
        group = Pmw.Group(page, tag_text = 'Leap Motion Mover')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)
        leap_opt = Pmw.OptionMenu(group.interior(),
                                    labelpos = 'w',
                                    label_text = 'Mover Mode',
                                    items = ('Move','Edit',),
                                    initialitem = 'Move',
                                    )
        leap_opt.pack(fill='x',expand=1,padx=4,pady=1)

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

        #create placeholders for listeners
        self.hmd=0
        self.hand=0

        #self.notebook.setnaturalsize()

        self.dialog.show()

    def execute(self,result):
        if result:
            print 'You clicked on, ' + result
            if result == 'Run Rift Only':
                # PyMOL will crash if Rift is off or not connected.
                # how do I make these run in the background?
                if self.radiobuttons_rot.getvalue() == 'Natural Rotation':
                    if self.radiobuttons_mol.getvalue() == 'Yes':
                        self.hmd = PymolHmd(naturalRotation=True,editMolecule=True,pdb=self.pdb_text.getvalue())
                    else:
                        self.hmd = PymolHmd(naturalRotation=True,editMolecule=False,pdb=self.pdb_text.getvalue())
                else:
                    if self.radiobuttons_mol.getvalue() == 'Yes':
                        self.hmd = PymolHmd(naturalRotation=False,editMolecule=True,pdb=self.pdb_text.getvalue())
                    else:
                        self.hmd = PymolHmd(naturalRotation=False,editMolecule=False,pdb=self.pdb_text.getvalue())

                self.hmd.Run()

            elif result == 'Run Leap Only':
                # Leap Motion needed for this... name convention is poor.
                # Currently doesn't work... exits on init.
                self.hand = PymolListener()

            elif result == 'Run Both':
                print 'Inner if, should create a Viewer and Mover object.'
        else:
            self.quit()

    def changed(self):
        # change stereo settings
        cmd.set('stereo_shift',self.stereo_shift_text.getvalue())
        cmd.set('stereo_angle',self.stereo_angle_text.getvalue())

    def quit(self):
        self.dialog.destroy()
