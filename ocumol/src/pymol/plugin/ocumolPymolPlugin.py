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

def __init__(self):
    """
    Add OcuMOL Leap to the PyMOL Plugins menu.
    """

    self.menuBar.addmenuitem('Plugin', 'command', 'Launch OcuMOL Leap',
                                label='OcuMOL Leap',
                                command= lambda s=self: OcuMOLLeapPlugin(s))

### OcuMOLLeap code ###
class OcuMOLLeapPlugin:
    def __init__(self,app):
        self.parent = app.root

        self.dialog = Pmw.Dialog(self.parent,
                                buttons = ( "Run Rift Only",
                                            "Run Leap Only",
                                            "Run Both",
                                            "About" ),
                                title = 'OcuMOL Leap',
                                command = self.execute)

        self.dialog.withdraw()

        Pmw.setbusycursorattributes(self.dialog.component('hull'))

        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill='both',expand=1,padx=5,pady=5)

        # Create Oculus Rift Page
        page = self.notebook.add('Rift Visualizer')
        riftGroup = Pmw.Group(page, tag_text = 'Oculus Rift Visualizer')
        riftGroup.pack(fill = 'both', expand = 1, padx = 5, pady = 5)

        # Radio buttons to select for rotation method
        self.rotationRadio = Pmw.RadioSelect(riftGroup.interior(),
                                                orient = 'horizontal',
                                                labelpos = 'w',
                                                label_text = 'Rotation Method:',
                                                frame_relief = 'ridge',
                                                frame_borderwidth = 2,)

        self.rotationRadio.add('Natural Rotation')
        self.rotationRadio.add('Molecule Rotation')

        self.rotationRadio.setvalue('Natural Rotation')
        self.rotationRadio.pack(padx=1,pady=1)

        # Radio buttons to select for view editing at initialization
        self.moleculeRadio = Pmw.RadioSelect(riftGroup.interior(),
                                                orient = 'horizontal',
                                                labelpos = 'w',
                                                label_text = 'Edit View:',
                                                frame_relief = 'ridge',
                                                frame_borderwidth = 2,)
        self.moleculeRadio.add('Yes')
        self.moleculeRadio.add('No')

        self.moleculeRadio.setvalue('No')
        self.moleculeRadio.pack(padx=1,pady=1)

        # Text to input pdb id, for testing
        self.pdbText = Pmw.EntryField(riftGroup.interior(),
                                        labelpos='w',
                                        label_text='Load PDB ID:',
                                        command=None,)

        self.pdbText.pack(padx=1, pady=1)

        # Text to input stereo shift parameter, see wiki for more info
        # TODO make this a slider
        self.stereoShiftText = Pmw.EntryField(riftGroup.interior(),
                                                labelpos='w',
                                                label_text='Stereo shift:',
                                                value='1.0',
                                                validate = {'validator' : 'real'},
                                                modifiedcommand=self.changed,)

        self.stereoShiftText.pack(padx=1, pady=1)

        # Text to input stereo angle parameter, see wiki for more info
        # TODO make this a slider
        self.stereoAngleText = Pmw.EntryField(riftGroup.interior(),
                                                labelpos='w',
                                                label_text='Stereo angle:',
                                                value='1.0',
                                                validate = {'validator' : 'real'},
                                                modifiedcommand=self.changed,)

        self.stereoAngleText.pack(padx=1, pady=1)

        # Create Leap Motion Page
        page = self.notebook.add('Leap Mover')
        group = Pmw.Group(page, tag_text = 'Leap Motion Mover')
        group.pack(fill = 'both', expand = 1, padx = 5, pady = 5)

        self.leapOpt = Pmw.OptionMenu(group.interior(),
                                    labelpos = 'w',
                                    label_text = 'Mover Mode',
                                    items = ('Move','Edit',),
                                    initialitem = 'Move',)

        self.leapOpt.pack(padx=1,pady=1)

        # Create About Pop-up
        Pmw.aboutversion('1.0')
        Pmw.aboutcopyright(
        'Apache License\n' +
        'Version 2.0, January 2004'
        )
        Pmw.aboutcontact(
        'PyMOL Oculus Rift Viewer and Leap Motion Mover\n' +
        'Max Klein, Jeliazko Jeliazkov, Henry Lessen, and Mariusz Matyszewski, 2015.\n' +
        'https://github.com/lqtza/OcuMOL_Leap'
        ) # note github link cannot be copied for some reason...
        self.about = Pmw.AboutDialog(self.parent,applicationname="OcuMOL Leap")
        self.about.withdraw()

        # create placeholders for listeners
        self.hmd=0
        self.hand=0

        self.notebook.setnaturalsize()

        self.dialog.show()

    def execute(self,result):
        if result:
            if result == 'Run Rift Only':
                # PyMOL will crash if Rift is off or not connected.
                if self.rotationRadio.getvalue() == 'Natural Rotation':
                    if self.moleculeRadio.getvalue() == 'Yes':
                        self.hmd = PymolHmd(naturalRotation=True,
                                            editMolecule=True,
                                            pdb=self.pdbText.getvalue())
                        self.hmd.start() # this is an inherited class function
                    else:
                        self.hmd = PymolHmd(naturalRotation=True,
                                            editMolecule=False,
                                            pdb=self.pdbText.getvalue())
                        self.hmd.start() # this is an inherited class function
                else:
                    if self.moleculeRadio.getvalue() == 'Yes':
                        self.hmd = PymolHmd(naturalRotation=False,
                                            editMolecule=True,
                                            pdb=self.pdbText.getvalue())
                        self.hmd.start() # this is an inherited class function
                    else:
                        self.hmd = PymolHmd(naturalRotation=False,
                                            editMolecule=False,
                                            pdb=self.pdbText.getvalue())
                        self.hmd.start() # this is an inherited class function

                # set initial stereo properties
                #cmd.set('stereo_shift',self.stereoShiftText.getvalue())
                #cmd.set('stereo_angle',self.stereoAngleText.getvalue())

            elif result == 'Run Leap Only':
                # Leap Motion needed for this... name convention is poor.
                # Currently doesn't work... exits on init.
                self.hand = PymolListener()

            elif result == 'Run Both':
                print 'Inner if, should create a Viewer and Mover object.'
            elif result == 'About':
                self.about.show()
        else:
            self.quit()

    def changed(self):
        # change stereo settings
        cmd.set('stereo_shift',self.stereoShiftText.getvalue())
        cmd.set('stereo_angle',self.stereoAngleText.getvalue())

    def quit(self):
        self.dialog.destroy()
