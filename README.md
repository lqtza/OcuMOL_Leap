# OcuMOL_Leap
Incorporating the Oculus Rift and Leap Motion into PyMOL.

**Note PyMOL must use system Python (freeware PyMOL should work fine, but only one of the licensed version will do).**

## Oculus Rift Only
1) Clone the repository from https://github.com/telamonian/ocudump/. Follow the directions. This should generate a cython
    library from the Oculus SDK. Your versions of the Oculus SDK and Runtime should be identical (this has been tested with 0.5
    on Mac OS X 10.10). You must run the 'make ocudump_cython' command to generate the cython library.

2) Assuming you have PyMOL (we used the latest Schrodinger system Python version), running 'spawn prep_and_run.py' script from
    within PyMOL should work. If it doesn't, check to make sure your PYTHONPATH includes the 'ocudump' library. 
    It can be set within the prep_and_run.py script.

## Oculus Rift and Leap Motion
This will be coming soon. There's some issues to iron out.
