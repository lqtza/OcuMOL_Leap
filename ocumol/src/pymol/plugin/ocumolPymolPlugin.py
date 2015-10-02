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
                                            "Run Both" ),
                                title = 'OcuMOL Leap',
                                command = self.execute)

        self.dialog.withdraw()

        Pmw.setbusycursorattributes(self.dialog.component('hull'))

        authorList = Tkinter.Label(self.dialog.interior(),
                            text = """PyMOL Oculus Rift Viewer and Leap Motion Mover\n
                                        Max Klein, Jeliazko Jeliazkov, Henry Lessen,
                                        and Mariusz Matyszewski, 2015.\n
                                        https://github.com/lqtza/OcuMOL_Leap""",
                            background = 'black',
                            foreground = 'white',
                            #pady = 20,)

        authoList.pack(expand = 1, fill = 'both', padx = 4, pady = 4)

        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill='both',expand=1,padx=10,pady=10)

        # Create Oculus Rift Page
        page = self.notebook.add('Rift Visualizer')
        group = Pmw.Group(page, tag_text = 'Oculus Rift Visualizer')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)

        # Radio buttons to select for rotation method
        self.rotationRadio = Pmw.RadioSelect(group.interior(),
                                                buttontype = 'radiobutton',
                                                orient = 'vertical',
                                                labelpos = 'w',)

        self.rotationRadio.add('Natural Rotation')
        self.rotationRadio.add('Molecule Rotation')

        self.rotationRadio.setvalue('Natural Rotation')
        self.rotationRadio.pack(fill='x',padx=4,pady=1)

        # Radio buttons to select for view editing at initialization
        self.moleculeRadio = Pmw.RadioSelect(group.interior(),
                                                buttontype = 'radiobutton',
                                                orient = 'vertical',
                                                labelpos = 'w',
                                                label_text = 'Edit Molecule?',)

        self.moleculeRadio.add('Yes')
        self.moleculeRadio.add('No')

        self.moleculeRadio.setvalue('No')
        self.moleculeRadio.pack(fill='x',padx=4,pady=1)

        # Text to input pdb id, for testing
        self.pdbText = Pmw.EntryField(group.interior(),
                                        labelpos='w',
                                        label_text='Load PDB ID:',
                                        command=None,)

        self.pdbText.pack(fill='x', padx=4, pady=1)

        # Text to input stereo shift parameter, see wiki for more info
        self.stereoShiftText = Pmw.EntryField(group.interior(),
                                                labelpos='w',
                                                label_text='Stereo shift:',
                                                value='5.24',
                                                validate = {'validator' : 'real'},
                                                modifiedcommand=self.changed,)

        self.stereoShiftText.pack(fill='x', padx=4, pady=1)

        # Text to input stereo angle parameter, see wiki for more info
        self.stereoAngleText = Pmw.EntryField(group.interior(),
                                                labelpos='w',
                                                label_text='Stereo angle:',
                                                value='2.0',
                                                validate = {'validator' : 'real'},
                                                modifiedcommand=self.changed,)

        self.stereoAngleText.pack(fill='x', padx=4, pady=1)

        # Create Leap Motion Page
        page = self.notebook.add('Leap Mover')
        group = Pmw.Group(page, tag_text = 'Leap Motion Mover')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)

        self.leapOpt = Pmw.OptionMenu(group.interior(),
                                    labelpos = 'w',
                                    label_text = 'Mover Mode',
                                    items = ('Move','Edit',),
                                    initialitem = 'Move',)

        self.leapOpt.pack(fill='x',expand=1,padx=4,pady=1)

        # Create About Page
        page = self.notebook.add('About')
        group = Pmw.Group(page, tag_text='About OcuMOL Leap')
        group.pack(fill = 'both', expand = 1, padx = 10, pady = 5)

        aboutText = """Magic text goes here."""
        interiorFrame = Frame(group.interior())
        bar = Scrollbar(interior_frame,)
        textHolder = Text(interior_frame, yscrollcommand=bar.set,background='#ddddff',font="Times 14")
        bar.config(command=text_holder.yview)
        textHolder.insert(END,aboutText)
        textHolder.pack(side=LEFT,expand="yes",fill="both")
        bar.pack(side=LEFT,expand="yes",fill="y")
        interiorFrame.pack(expand="yes",fill="both")

        # create placeholders for listeners
        self.hmd=0
        self.hand=0

        self.notebook.setnaturalsize()

        self.dialog.show()

    def execute(self,result):
        if result:
            print 'You clicked on, ' + result
            if result == 'Run Rift Only':
                # PyMOL will crash if Rift is off or not connected.
                if self.rotationRadio.getvalue() == 'Natural Rotation':
                    if self.moleculeRadio.getvalue() == 'Yes':
                        self.hmd = PymolHmd(naturalRotation=True,
                                            editMolecule=True,
                                            pdb=self.pdbText.getvalue())
                    else:
                        self.hmd = PymolHmd(naturalRotation=True,
                                            editMolecule=False,
                                            pdb=self.pdbText.getvalue())
                else:
                    if self.moleculeRadio.getvalue() == 'Yes':
                        self.hmd = PymolHmd(naturalRotation=False,
                                            editMolecule=True,
                                            pdb=self.pdbText.getvalue())
                    else:
                        self.hmd = PymolHmd(naturalRotation=False,
                                            editMolecule=False,
                                            pdb=self.pdbText.getvalue())

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
        cmd.set('stereo_shift',self.stereoShiftText.getvalue())
        cmd.set('stereo_angle',self.stereoAngleText.getvalue())

    def quit(self):
        self.dialog.destroy()
