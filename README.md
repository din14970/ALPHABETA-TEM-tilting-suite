# ALPHABETA TEM tilting suite

Python software to calculate alpha and beta tilt angles easily of a double tilt stage in a TEM (transmission electron microscope). Most suited for old TEMs and when one has to tilt to many different zone axes or 2-beam / Weak beam conditions.

My motivation for making this was the extremely frustrating experience of trying to set up correct
diffraction conditions on a piece of highly deformed steel (FCC), with dense dislocation network and
high density of deformation twins. While conventional wisdom dictates that one needs to only follow
kikuchi lines to arrive at zone axes, the iterative process made this extremely complex and time consuming.
In addition, many times after having spent time arriving at a new zone, I could discover it was outside of the tilt range.

Therefore, my idea was to be able to simply calculate the required tilts, so I could tilt in image mode.
Hereby I present you my first attempt, implemented in python.

The engine behind it are some classes and methods in the temCalcs.py file.
The main classes:

1) TEM: stores detector and stage objects (even though a TEM usually has only one stage, the option is left open to have multtipe)
To access a particular detecter, stage and chrystal objects must pass through the TEM.

2) Stage: stages contain Crystal objects (see later). It keeps track of its tilt range and current alpha and beta
setting. It can calculate which alpha and beta it should tilt to given a vector expressed in Stage coordinates.

3) Detector: a detector mainly stores various calibration dictionaries: in which direction do the alpha and beta
tilt axes lie at different magnifications or camera lenghts. Imaging, Diffraction, STEM are supported. Soon I will
generalize this so that a user can add any mode he or she pleases. Size calibration (real size/pixel) can be 
stored as well but so far has not really proven usefull.

4) Structure: a structure is the prototype of a crystal. The 7 crystal classes are supported, but no distinction
is made between space groups - sorry, I work on FCC materials and it seemed too much of a hastle. A structure
does not have an orientation and can exist without the TEM. It stores crystallographic data and can be used to 
perform basic crystallographic calculations - length of a lattice vector (real or recyprocal), angle between 
lattice vectors, conversion to cartesian coordinates in a cartesian coordinate system stuck to the structure,
converting a real lattice vector to the recyprocal lattice vector pointing in the same direction, etc.

5) Crystal: a Crystal is basically a piece of material with a particular structure stuck to a stage.
A crystal can perform all the functions of its derivative structure, but it also stores an orientation matrix.
This matrix determines how the cartesian coordinate system stuck to the crystal is rotated with respect to
the coordinate system stuck to the stage. The way to calculate this orientation matrix is very easy:
 
 a) Tilt your piece of material to any zone axis and index it correctly
 
 b) Change the alpha and beta of the stage object so that it reflects reality
 
 c) Call the calcOrient function on the crystal object, check in the code for the arguments to be provided.
 From thereon out your crystal is correctly anchored to the stage.
 
 Note that it is important that the detectors are accurately callibrated for this to work. Look at the jupyter notebook
 for more details
 
 The crystal is the workhorse class. With it you can calculate what alpha and beta you need to tilt to another zone
 or to a two-beam (implemented but not properly tested) or weak beam condition (not implemented yet). It does this by making reference to functions from its higher ups: the stage, the TEM and the detector.
 
 There are some additional methods just in the file. Rotation matrixes can be defined according to a rotation around X, Y or Z
 or simply any rotation matrix defined by an arbitrary axis and angle.
 
 There are some functions in there to try and get an integer representation of a vector (example (0.4, 0.6, 0.2) -> (2, 3, 1))
 
 Also there are some functions in there to make a stereographic projection figure.
 
 Hope with the things I have provided you can at least already get some use out of this.
 
 temCalcs is dependent on the following packages:
 
 numpy, math, matplotlib, re,  fractions, scipy, pyQT5
 
 Since version 2, a user interface is available that will be sufficient for most purposes. It is not meant to be idiot proof, so I'm sure there are tons of buggs to be found. The user interface is implemented using pyQT5 in the GUI.py file. The program can be run by calling "python GUI.py" in the command line. The user interface is rather self explanatory and an upcomming paper will go more into detail on the program.
