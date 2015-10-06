"""
    An OcuMOL Leap plugin for PyMOL. To install go to the Plugin menu in PyMOL
    and under the Install New Plugin tab click Choose file... and select this file.
"""
import os, sys
import Pmw
import re
import Tkinter
from Tkinter import *
import urllib2

import pymol
from pymol import cmd

from ocumol.src.pymol.pymolHmd import PymolHmd
from ocumol.src.hands.leap_only import PymolListener

def __init__(self):
    """
    Add OcuMOL Leap to the PyMOL Plugins menu.
    """

    self.menuBar.addmenuitem('Plugin', 'command', 'Launch OcuMOL Leap',
                             command = lambda s=self: OcuMOLLeapPlugin(s),
                             label='OcuMOL Leap',)

def RadioMultipleClosure(obj, buttonAttrNameDict):
    def Callback(buttonName, val):
        obj.__setattr__(buttonAttrNameDict[buttonName], val)
#         if obj.running:
#             obj.InitPymol()
    return Callback

def RadioSingleClosure(obj, attrName, targetButtonName):
    def Callback(buttonName):
        if buttonName==targetButtonName:
            obj.__setattr__(attrName, True)
            if obj.running:
                obj.initPymol()
    return Callback

def RotationClosure(obj):
    def Callback(buttonName):
        if buttonName=='Natural rotation':
            obj.__setattr__('naturalRotation', True)
        else:
            obj.__setattr__('naturalRotation', False)
            if obj.running:
                obj.setOriginAtMolecule()
    return Callback

def ValidatePDBClosure(obj):
    def ValidatePDB(txt):
        if len(txt)==4 and txt==re.match('[a-zA-Z0-9]*', txt).group(0):
            if obj.running:
                obj.reinitPymol(pdb=txt)
            else:
                obj.pdb = txt
            return Pmw.OK
        else:
            return Pmw.PARTIAL
    return ValidatePDB

### OcuMOLLeap code ###
class OcuMOLLeapPlugin:
    def __init__(self, app):
        self.app = app
        
        # create placeholders for listeners
        self.hand=0
        
        # set up hmd and callbacks
        self.hmd = PymolHmd()
        
        # initialize the main ocumol megawidget
        self.initNotebook()

        # initialize the subwidgets
        self.initRiftPage()
        self.initLeapPage()
        self.initRiftDebugPage()
        self.initAboutPopup()
        
        # finish up the initialization of the main megawidget
        self.initNotebookFinalize()
    
    def initAboutPopup(self):
        # Create About Pop-up
        Pmw.aboutversion('1.0')
        Pmw.aboutcopyright(
            'Apache License\n' +
            'Version 2.0, January 2004')
        Pmw.aboutcontact(
            'PyMOL Oculus Rift Viewer and Leap Motion Mover\n' +
            'Max Klein, Jeliazko Jeliazkov, Henry Lessen, and Mariusz Matyszewski, 2015.\n' +
            'https://github.com/lqtza/OcuMOL_Leap') # note github link cannot be copied for some reason...
        self.about = Pmw.AboutDialog(self.parent,applicationname="OcuMOL Leap")
        self.about.withdraw()
    
    def initLeapPage(self):
        # Create Leap Motion Page
        page = self.notebook.add('Leap Mover')
        group = Pmw.Group(page, tag_text='Leap Motion Mover')
        group.pack(fill='both', expand=1, padx=5, pady=5)

        self.leapOpt = Pmw.OptionMenu(group.interior(),
                                      initialitem='Move',
                                      items=('Move','Edit'),
                                      labelpos='w',
                                      label_text='Mover Mode')
        self.leapOpt.pack(padx=1,pady=1)
    
    def initNotebook(self):
        # set up the main ocumol megawidget (which is a notebook)
        self.parent = self.app.root
        self.dialog = Pmw.Dialog(self.parent,
                                 buttons=("Run Rift Only",
                                          "Run Leap Only",
                                          "Run Both",
                                          "About" ),
                                 command=self.execute,
                                 title='OcuMOL Leap')
        self.dialog.withdraw()
        Pmw.setbusycursorattributes(self.dialog.component('hull'))
        self.notebook = Pmw.NoteBook(self.dialog.interior())
        self.notebook.pack(fill='both', expand=1, padx=5, pady=5)
    
    def initNotebookFinalize(self):
        # finish setting up the notebook (which is the main ocumol megawidget)
        self.notebook.setnaturalsize()
        self.dialog.show()
    
    def initRiftPage(self):
        # set up some callbacks
        DebugModeCallback = RadioMultipleClosure(self.hmd, {'Debug mode':'debugMode'})
        MoleculeCallback = RadioSingleClosure(self.hmd, 'editMolecule', 'Yes')
        RotationCallback = RotationClosure(self.hmd)
        ValidatePDB = ValidatePDBClosure(self.hmd)
        
        # Create Oculus Rift Page
        page = self.notebook.add('Rift Visualizer')
        riftGroup = Pmw.Group(page, tag_text='Oculus Rift Visualizer')
        riftGroup.pack(fill='both', expand=1, padx=5, pady=5)

        # Radio buttons to select for rotation method
        self.rotationRadio = Pmw.RadioSelect(riftGroup.interior(),
                                             command=RotationCallback,
                                             frame_borderwidth=2,
                                             frame_relief='ridge',
                                             labelpos='w',
                                             label_text='Rotation method:',
                                             orient='horizontal')

        self.rotationRadio.add('Natural rotation')
        self.rotationRadio.add('Molecule rotation')

        self.rotationRadio.setvalue('Natural rotation')
        self.rotationRadio.pack(padx=1,pady=1)

        # Radio buttons to select for view editing at initialization
        self.moleculeRadio = Pmw.RadioSelect(riftGroup.interior(),
                                             command=MoleculeCallback,
                                             frame_borderwidth=2,
                                             frame_relief='ridge',
                                             labelpos='w',
                                             label_text='Edit view:',
                                             orient='horizontal')
        self.moleculeRadio.add('Yes')
        self.moleculeRadio.add('No')

        self.moleculeRadio.setvalue('No')
        self.moleculeRadio.pack(padx=1,pady=1)

        # Text to input pdb id, for testing
        self.pdbText = Pmw.EntryField(riftGroup.interior(),
                                      labelpos='w',
                                      label_text='Load PDB ID:',
                                      validate=ValidatePDB)
        self.pdbText.pack(padx=1, pady=1)

        # Text to input stereo shift parameter, see wiki for more info
        # TODO make this a slider
        self.stereoShiftText = Pmw.EntryField(riftGroup.interior(),
                                              labelpos='w',
                                              label_text='Stereo shift:',
                                              modifiedcommand=self.changed,
                                              validate={'validator' : 'real'},
                                              value='1.0')
        self.stereoShiftText.pack(padx=1, pady=1)

        # Text to input stereo angle parameter, see wiki for more info
        # TODO make this a slider
        self.stereoAngleText = Pmw.EntryField(riftGroup.interior(),
                                              labelpos='w',
                                              label_text='Stereo angle:',
                                              modifiedcommand=self.changed,
                                              value='1.0',
                                              validate={'validator' : 'real'})
        self.stereoAngleText.pack(padx=1, pady=1)
        
        # Check button to run the HMD in debug mode, which will keep going even if the Rift is not attached/detected
        self.debugModeCheck = Pmw.RadioSelect(riftGroup.interior(),
                                              buttontype='checkbutton',
                                              command=DebugModeCallback,
                                              frame_borderwidth=2,
                                              frame_relief='ridge',
                                              labelpos='w',
                                              label_text='',
                                              orient='horizontal')
        self.debugModeCheck.add('Debug mode')

#         self.moleculeRadio.setvalue('Debug mode')
        self.debugModeCheck.pack(padx=1,pady=1)

    def initRiftDebugPage(self):
        # Create Rift Debug Page
        riftDebugPage = self.notebook.add('Rift Debug')
        riftDebugGroup = Pmw.Group(riftDebugPage, tag_text='Rift Debug')
        riftDebugGroup.pack(fill='both', expand=1, padx=5, pady=5)

#         self.stereoAngleText = Pmw.EntryField(riftGroup.interior(),
#                                               labelpos='w',
#                                               label_text='Stereo angle:',
#                                               modifiedcommand=self.changed,
#                                               value='1.0',
#                                               validate={'validator' : 'real'})

    def execute(self, result):
        if result:
            if result=='Run Rift Only':
                # set initial stereo properties
                cmd.set('stereo_shift',self.stereoShiftText.getvalue())
                cmd.set('stereo_angle',self.stereoAngleText.getvalue())
                
                self.hmd.start()
                
                # we've already run the rift, so let's run the hand support
                # if that's what has been called for
                if result=='Run Both':
                    self.hand = PymolListener()

            elif result=='Run Leap Only':
                # Leap Motion needed for this... name convention is poor.
                # Currently doesn't work... exits on init.
                self.hand = PymolListener()

            elif result=='About':
                self.about.show()
        else:
            self.quit()

    def changed(self):
        # change stereo settings
        cmd.set('stereo_shift',self.stereoShiftText.getvalue())
        cmd.set('stereo_angle',self.stereoAngleText.getvalue())

    def quit(self):
        self.dialog.destroy()
