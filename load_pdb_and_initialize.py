""" A simple script which takes in a single argument (pdb) and fetches it, if
    possible, to the local dir. The script then optimizes the view for OcuMOL Leap.
"""
import sys

from pymol import cmd

def set_stage(pdb_id):
    #load pdb
    load_pdb(pdb_id)
    #color each chain differently
    cmd.hide('everything','all')
    cmd.show('cartoon','all')
    cmd.show('surface','all')
    cmd.set('transparency',0.5)
    util.cbc()

    #move slab away to give a comfortable viewing area
    cmd.clip('move',10)

    #set stereo mode
    cmd.set('stereo_mode',3)
    cmd.stereo()

    #get rid of gui
    cmd.set('internal_gui',0)

    #set resolution to HD for rift
    cmd.viewport(1920,1080)
    return

def load_pdb(pdb_id):
    if len(pdb_id) != 4:
        print 'Please input a PDB ID of length 4.'
        sys.exit(1)
    else:
        cmd.fetch(pdb_id,async=0)

    return

cmd.extend("set_stage",set_stage)
