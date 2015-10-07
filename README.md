# OcuMOL_Leap
Incorporating the Oculus Rift and Leap Motion into PyMOL.

**Note PyMOL must use system Python (freeware PyMOL should work fine, but only the `syspython` versions of licensed PyMOL will work).**

## Oculus Rift Only
### Hardware Setup
- Install latest Oculus Runtime. Make sure that the version matches your version of the Oculus SDK.
- On OSX
  - After plugging the Rift in, set your displays to mirroring mode in `System Preferences->Displays`.
- On Windows
  - After plugging the Rift in, set your Rift to Extended Display Mode.
  - Open your monitor settings in Control Panel and make sure that it's set to `Extend these displays` rather than `Duplicate these displays`.
    - You may also need to rotate the display corresponding to your Rift.
  - Set the Rift as your primary display.
    - If you don't do this, when PyMOL goes full screen it may end up on the wrong display.

### Software Setup
- Clone the ocudump repository from https://github.com/telamonian/ocudump/. Follow the ocdump Build and Install directions.     - This installs a python module called ocudump that helps OcuMol interfaces with the Oculus Rift.
  - Your versions of the Oculus SDK and Runtime should be identical (this has been tested with SDK v0.5 on Mac OS X 10.10, and with SDK v0.6 on Windows 7).
- Download the latest version of the Leap SDK from https://developer.leapmotion.com/downloads. Untar it, and move the `LeapSDK` directory within to somewhere convenient.
  - For example, I placed my copy of the LeapSDK at `/usr/local/LeapSDK`
- Clone this repository (OcuMOL_Leap), and `cd` to it.
- Run `LEAPSDK_DIR=<path-to-your-leapSDK> pip install -e .`
  - Be sure to replace `<path-to-your-leapSDK>` with the actual path to *your* copy of the Leap SDK.
  - If you are running with the `sudo` command, you can pass environment variables by `-E` (at least on OS X), so you would run: `LEAPSDK_DIR=<path-to-your-leapSDK> sudo -E pip install -e .`
  - The above command will get `pip` to install OcuMol in development mode, meaning that it will create a kind of soft link between your python module directory and the OcuMol directory.
  - Eventually OcuMol will also be available directly through pypi.
- The OcuMol plugin should now be available from the `Plugin` menu the next time you start PyMol.
  - OcuMol tells PyMol where to look for the plugin by modifying your `~/.pymolrc` PyMol config file (or creating it if doesn't exist). 

**Note: The Python script will crash if the Oculus Rift is not connected.**

## Oculus Rift and Leap Motion
This will be coming soon. There's some issues to iron out.
