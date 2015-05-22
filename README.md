# OcuMOL_Leap
Incorporating the Oculus Rift and Leap Motion into PyMOL.

**Note PyMOL must use system Python (freeware PyMOL should work fine, but only the `syspython` versions of licensed PyMOL will work).**

## Oculus Rift Only
1. Clone the ocudump repository from https://github.com/telamonian/ocudump/. Follow the directions. This should generate a cython library (called `ocudump.so`) from the Oculus SDK. Your versions of the Oculus SDK and Runtime should be identical (this has been tested with 0.5 on Mac OS X 10.10).

2. Clone this repository (OcuMOL_Leap). Open the `prep_and_run.py` script in an editor. Near the top of the script there is a line which reads:  
`sys.path.append("/Users/lqtza/Hacks/ocudump/build/src/cython")`  
Change this path so it points to the directory in which you built your `ocudump.so` file. If you followed the instructions from the ocudump readme, this would be:  
`sys.path.append("<path-to-your-ocudump-reppository>/build/src/cython")`

3. Assuming you have PyMOL (we used the latest Schrodinger system Python version, `MacPyMOL-v1.7.6.0-syspython.dmg`), running 'spawn prep_and_run.py' script from within PyMOL should work. If it doesn't, check to make sure your PYTHONPATH includes the 'ocudump' library. It can be set within the prep_and_run.py script.

## Oculus Rift and Leap Motion
This will be coming soon. There's some issues to iron out.
