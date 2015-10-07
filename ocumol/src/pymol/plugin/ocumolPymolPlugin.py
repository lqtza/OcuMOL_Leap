"""
    An OcuMOL Leap plugin for PyMOL. To install go to the Plugin menu in PyMOL
    and under the Install New Plugin tab click Choose file... and select this file.
"""
import locale
import os, sys
import Pmw
import re
import Tkinter
from Tkinter import *

import pymol
from pymol import cmd

from ocumol.src.helper.helper import PrintException
from ocumol.src.hands.leap_only import PymolListener
from ocumol.src.pymol.pymolHmd import PymolHmd

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

def AnimateElemInnerLoop(pluginObj, poseCoordName):
    args = [poseCoordName]
    for argN in pluginObj.animateElemColumnStrings:
        key = '%s_%s' % (poseCoordName, argN)
        if pluginObj.animateElemFields[key].getvalue()=='':
            return False
        convFunc = int if argN=='period' else float
        val = convFunc(pluginObj.animateElemFields[key].getvalue())
        args.append(val)
    pluginObj.hmd.initAnimateElement(*args)
    return True

def AnimateElemClosure(argName, pluginObj):
    def AnimateElemSubclosure(poseCoordName):
        def AnimateElemValidateFloat(txt):
            if txt=='':
                return Pmw.OK
            try:
                locale.atof(txt)
                AnimateElemInnerLoop(pluginObj, poseCoordName)
                return Pmw.OK
            except ValueError:
                return Pmw.PARTIAL

        def AnimateElemValidateInt(txt):
            if txt=='':
                return Pmw.OK
            try:
                locale.atoi(txt)
                AnimateElemInnerLoop(pluginObj, poseCoordName)
                return Pmw.OK
            except ValueError:
                return Pmw.PARTIAL

        if argName=='period':
            return AnimateElemValidateInt
        else:
            return AnimateElemValidateFloat
    return AnimateElemSubclosure

### OcuMOLLeap code ###
class OcuMOLLeapPlugin:
    def __init__(self, app):
        # the t_e block here ensures that we actually get an error message when stuff goes wrong
        try:
            self.app = app

            # create placeholders for listeners
            self.hand=0

            # set up hmd and callbacks
            self.hmd = PymolHmd()

            # initialize the main ocumol megawidget
            self.initNotebook()

            # initialize the subwidgets
            self.initRiftPage()
            self.initRiftDebugPage()
            self.initLeapPage()
            self.initAboutPopup()

            # finish up the initialization of the main megawidget
            self.initNotebookFinalize()

        except:
            PrintException()

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
        self.rotationRadio.pack(padx=1, pady=1)

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
        self.debugModeCheck.pack(padx=1,pady=1)

    def initRiftDebugPage(self):
        # Create Rift Debug Page
        riftDebugPage = self.notebook.add('Rift Debug')
        riftDebugGroup = Pmw.Group(riftDebugPage, tag_text='Rift Debug')
        riftDebugGroup.pack(fill='both', expand=1, padx=5, pady=5)

        frame = Tkinter.Frame(riftDebugGroup.interior())
        frame.pack(fill='both', expand=1)

        animateLabel = Tkinter.Label(frame, text='Debug animation controls:')
        animateLabel.grid(column=0, row=0, sticky='nw')

        self.animateElemColumnStrings = ['min', 'max', 'period']
        columnLabels = {}
        for i in range(3):
            columnLabels[i] = Tkinter.Label(frame, text=self.animateElemColumnStrings[i])
            columnLabels[i].grid(column=i+1, row=1, sticky='nw')

        self.animateElemRowStrings = ['xRot', 'yRot', 'zRot', 'x', 'y', 'z']
        rowLabels = {}
        for i in range(6):
            rowLabels[i] = Tkinter.Label(frame, text=self.animateElemRowStrings[i])
            rowLabels[i].grid(column=0, row=i+2, sticky='nw')

        animateElemClosures = {'min':AnimateElemClosure(argName='min', pluginObj=self),
                               'max':AnimateElemClosure(argName='max', pluginObj=self),
                               'period':AnimateElemClosure(argName='period', pluginObj=self)}

        self.animateElemFields = {}
        for i,rs in enumerate(self.animateElemRowStrings):
            for j,cs in enumerate(self.animateElemColumnStrings):
                key = '%s_%s' % (rs, cs)
                AnimateElemValidate = animateElemClosures[cs](poseCoordName=rs)
                self.animateElemFields[key] = Pmw.EntryField(frame,
                                                             validate=AnimateElemValidate)
                self.animateElemFields[key].grid(column=j+1, row=i+2, stick='nw')

            # all of the fields in a given row must exist before validation can work, so we defer setting values until here
            for cs in self.animateElemColumnStrings:
                key = '%s_%s' % (rs, cs)
                val = '1' if cs=='period' else '0.0'
                self.animateElemFields[key].setentry(val)

        # for some reason, pmw checkbuttons don't work with grid :(
#         buttons = {}
#         for i in range(6):
#             buttons = Pmw.RadioSelect(frame,
#                                       buttontype='checkbutton',
#                                       orient='horizontal')
#             buttons[i].grid(column=3, row=i, sticky='w')

        animateExplanation = ['You can use the above panel to set independent ',
                              'animations on each of the 6 degrees of freedom ',
                              'OcuMol reads from the Rift.\n\n',
                              'min: sets the low end of the range of the ',
                              'animation.\n',
                              'max: sets the high end of the range of the ',
                              'animation.\n',
                              '(note: min & max have units of radians for the ',
                              'rotations, and standard PyMol unit for the ',
                              'translations)\n',
                              'period: sets the time it takes the animation to ',
                              'cycle from min to max and back again.\n',
                              '(note: period is measured in terms of count of ',
                              'tracking refreshes, which currently defaults to ',
                              'a rate of %d per second)\n' % self.hmd.trackingRefresh]
        animateExplanation = ''.join(animateExplanation)
        animateExplanationLabel = Tkinter.Label(frame, justify=Tkinter.LEFT, text=animateExplanation, wraplength=500)
        animateExplanationLabel.grid(column=0, row=8, columnspan=4, sticky='nw')

        frame.grid_rowconfigure(8, weight=1)
        frame.grid_columnconfigure(3, weight=1)

    def execute(self, result):
        if result:
            if result=='Run Rift Only':
                # set initial stereo properties
                cmd.set('stereo_shift', self.stereoShiftText.getvalue())
                cmd.set('stereo_angle', self.stereoAngleText.getvalue())

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
