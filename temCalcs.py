"""These scripts serve as the backbone of TEM and crystallography calculations. They mainly serve for tilting help in the TEM."""

import numpy as np
import math as ma
import pandas as pd
import matplotlib.pyplot as plt
import pickle
import re
from fractions import Fraction
import math
from scipy.optimize import fsolve

class Constants(object):
    def __init__(self):
        # Planck's constant  
        self.h=6.62607004e-34 #m^2*kg/s
        # Elementary charge
        self.e=1.60217662e-19 #couloms/electron
        # Pi
        self.pi = np.pi
        # speed of light in vacuum
        self.c = 299792458 #m/s
        # electron rest mass
        self.m0 = 9.10938356e-31 #kg
    
    def gamma(self, v):
        """relativistic corrector"""
        return 1/np.sqrt(1-v**2/self.c**2)

#Make sure all constants are calable simply from c        
c=Constants()

#At the moment not all that in favor of a structureManager. So all available materials/structures are handled globally in the structures dictionary:
structures = {}

def addStructure(name = "structure"+str(len(structures)), **kwargs): ##doesn't add another structure if the name already exists but creates a new structure object
    structures[name] = Structure(**kwargs)
    return structures[name]
    
def removeStructure(name):
    del structures[name]
    
def getStructure(name):
    return structures[name]

def changeStructureName(oldname, newname):
    structures[newname]=structures[oldname]
    del structures[oldname]

#Slightly Easier ways of defining vectors
def miller(h, k, l):
    return np.array([h, k, l])

def vector(h, k, l):
    return np.array([h, k, l])
    
def sizeAngle(mag, direc):
    """Only valid for detector vectors. Returns X, Y, Z=0 given a magnitude and angle(degrees) with the x-axis"""
    return mag*vector(np.cos(direc/360*2*np.pi),np.sin(direc/360*2*np.pi),0)
    
#Another fundamental list of objects is the number of TEMS. From tems, stages, detectors and crystals are defined.
microscopes = {}

def addMicroscope(name = "microscope"+str(len(microscopes)), **kwargs):
    microscopes[name] = TEM(**kwargs)
    return microscopes[name]
    
def removeMicroscope(name):
    del microscopes[name]
    
def getMicroscope(name):
    return microscopes[name]

def changeMicroscopeName(oldname, newname):
    microscopes[newname]=microscopes[oldname]
    del microscopes[oldname]

def saveSession(filname = "session.txt"):
    pass

def loadSession(filname = "session.txt"):
    pass    

"""Rotation matrixes"""
def XR(d, rnd=10):
    d=float(d)/360*2*np.pi
    mat=np.array([[1, 0, 0],[0, round(np.cos(d), rnd), round(-np.sin(d), rnd)], [0, round(np.sin(d), rnd), round(np.cos(d), rnd)]])
    return mat
    
def YR(d, rnd=10):
    d=float(d)/360*2*np.pi
    mat=np.array([[round(np.cos(d), rnd), 0, round(np.sin(d), rnd)], [0, 1, 0], [round(-np.sin(d), rnd), 0, round(np.cos(d), rnd)]])
    return mat

def ZR(d, rnd=10):
    d=float(d)/360*2*np.pi
    mat=np.array([[round(np.cos(d), rnd), round(-np.sin(d), rnd), 0], [round(np.sin(d), rnd), round(np.cos(d), rnd), 0], [0, 0, 1]])
    return mat

def axisAngle(axis, d, rnd=10):
    d=float(d)/360*2*np.pi
    axis = np.array(normMatrix(axis).T)[0,:]
    vx = axis[0]
    vy = axis[1]
    vz = axis[2]
    mat = np.array([[vx**2 + (vy**2+vz**2)*np.cos(d), vx*vy*(1-np.cos(d))-vz*np.sin(d), vx*vz*(1-np.cos(d))+vy*np.sin(d)], [vx*vy*(1-np.cos(d))+vz*np.sin(d), vy**2 + (vx**2+vz**2)*np.cos(d),  vy*vz*(1-np.cos(d))-vx*np.sin(d)], [vx*vz*(1-np.cos(d))-vy*np.sin(d), vy*vz*(1-np.cos(d))+vx*np.sin(d), vz**2 + (vx**2+vy**2)*np.cos(d)]])
    return mat

def invmat():
    return ZR(180)

class TEM(object):
    def __init__(self, kv=300):
        #Overpotential in volts
        self.v = kv*1e3
        self.wavelength = self.calcLambda(self.v)
        self.ewaldR = 1/self.wavelength
        
        #detectors is a dictionary name - detector object
        self.detectors = {}
        #Maybe should leave option open of having more than one stage
        self.stages = {}
        
    def addStage(self, name = "", **kwargs):
        if name == "":
            name = "stage"+str(len(self.stages))
        self.stages[name] = Stage(self, **kwargs)
        return self.stages[name]
    
    def getStage(self, name):
        return self.stages[name]
        
    def removeStage(self, name):
        del self.stages[name]
        
    def changeStageName(self, oldname, newname):
        self.stages[newname]=self.stages[oldname]
        del self.stages[oldname]
        
    def addDetector(self, name="", **kwargs):
        if name == "":
            name = "detector"+str(len(self.detectors))
        self.detectors[name] = Detector(self, **kwargs)
        return self.detectors[name]
    
    def getDetector(self, name):
        return self.detectors[name]
        
    def removeDetector(self, name):
        del self.detectors[name]
        
    def changeDetectorName(self, oldname, newname):
        self.detectors[newname]=self.detectors[oldname]
        del self.detectors[oldname]
    
    def setKv(self, kv):
        #Overpotential in volts
        self.v = kv*1e3
        self.wavelength = self.calcLambda(self.v)
        self.ewaldR = 1/self.wavelength
    
    def getKv(self):
        return self.v/1e3
    
    def getEwaldR(self, units = "m"):
        if units == "m":
            return self.ewaldR
        elif units == "nm":
            return self.ewaldR/1e9
        else:
            print("Those units are not recognized at the moment.")
            return none
        
    def getLambda(self):
        return self.wavelength
        
    def calcLambda(self, v):
        """Calculates the relativistic electron wavelength from the voltage"""
        return c.h/np.sqrt(2*c.m0*c.e*v*(1+c.e*v/(2*c.m0*c.c**2)))


class Stage(object):
    """The stage of the TEM can tilt with alpha and beta and stores its own state."""
    def __init__(self, TEM, alpha =0, beta = 0):
        #the stage needs to know which TEM it belongs to
        self.TEM = TEM
        
        self.alpha = alpha
        self.beta = beta
        
        self.alow = -30
        self.atop = 30
        self.blow = -20
        self.btop = 20
        
        #list of crystals (this determines how the crystal is oriented with respect to the stage)
        self.crystals = {}
        
        #it should be possible to store interesting orientations as name - [alpha, beta]
        self.orientations = {}
    
    def setAlphaRange(self, *args):
        self.alow, self.atop = args
        
    def setBetaRange(self, *args):
        self.blow, self.btop = args
    
    def getTEM(self):
        return self.TEM
    
    def addCrystal(self, structurename, name=""):
        if name=="":
            name = "crystal"+str(len(self.crystals))
        self.crystals[name]=Crystal(self, structurename)
        return self.crystals[name]
        
    def changeCrystalName(self, oldname, newname):
        self.crystals[newname]=self.crystals[oldname]
        del self.crystals[oldname]
        
    def removeCrystal(self, name):
        del self.crystals[name]
        
    def getCrystal(self, name):
        return self.crystals[name]
    
    def getAlpha(self):
        return self.alpha
        
    def getBeta(self):
        return self.beta
    
    def setAlpha(self, alpha):
        self.alpha=alpha
    
    def tiltTo(self, args):
        alpha, beta = args
        self.alpha = alpha
        self.beta = beta
        
    def tiltBy(self, args):
        alpha, beta = args
        self.alpha += alpha
        self.beta += beta
    
    def setBeta(self, beta):
        self.beta=beta
    
    def stageToAbs(self, vecs):
        """Turn the absolute coordinates into stage coordinates. This is identical to doing the active rotation. [010]->[0, cos(alpha), sin(alpha)] in stage coordinates with only alpha rotation."""
        alpha=float(self.alpha)
        beta=float(self.beta)
        
        #first rotate around xf, then rotate around yf. Because yf rotates with xf we have to apply yr first then xr
        rot = np.dot(XR(alpha), YR(beta))
        return np.dot(np.array(rot), np.array(vecs))
        
    def absToStage(self, vecs):
        """Turn stage coordinates into absolute coordinates. This is identical to doing the reverse rotation. [010]->[0, cos(alpha), -sin(alpha)] in real with only alpha rotation."""
        alpha=float(self.alpha)
        beta=float(self.beta)
        
        #The inverse of the absolute to stage
        rot = np.dot(YR(-beta), XR(-alpha))
        return np.dot(np.array(rot), np.array(vecs))
    
    def inrange(self, a, b):
        a=a/(2*np.pi)*360
        b=b/(2*np.pi)*360
        atop = self.atop
        alow = self.alow
        btop = self.btop
        blow = self.blow
        if a<=atop and a>=alow and b<=btop and b>=blow:
            return True
        else:
            return False
    
    def inrangeIncr(self, a, b):
        return self.inrange(a + self.alpha/360*(2*np.pi), b + self.beta/360*(2*np.pi))
    
    def calcAlphaBeta(self, tozall, verbose=False, rnd = 10):
        """A general rotation matrix is not possible with just alpha and beta. It can only be arranged that a certain vector is set to the optical axis.
        The way to do this was calculated on paper and with sympy. Toz is the vector that needs to be turned on the Z-axis"""
        return toAllCols(self.calcAlphaBetaVec, tozall, rnd = rnd, verbose=verbose) 
    
    def calcAlphaBetaVec(self, toz, verbose=False, rnd = 10):
        """Calc the alpha and beta one needs to go to for stage vector toz to be moved to the z-axis"""
        
        #normalize
        toz = toz/np.linalg.norm(toz)
            
        c = toz[0]#+1e-15
        d = toz[1]#+1e-15
        f = toz[2]#+1e-15

        length = np.sqrt(c**2+d**2+f**2)

        #there is a 180 degree uncertainty on the answer. If the angle is negative (fourth quadrant) it could be meant for the second quadrant.
        beta = np.arctan(-c/f)
        beta2 = 0
        if beta<0:
            beta2 = np.pi + beta
        else:
            beta2 = -np.pi + beta
        
        #there is an exception: for all vectors in the y-z plane (c=0), beta can be any value. So it is best to set it to 0
        if round(c,6)==0:
            beta=0
            beta2=np.pi
        
        #this results in two possibilities for the solution of the second equation
        denom1 = f*np.cos(beta) - c*np.sin(beta)
        denom2 = f*np.cos(beta2) - c*np.sin(beta2)

        alpha1a = np.arctan(d/denom1)
        alpha1b = 0
        if alpha1a<0:
            alpha1b = np.pi + alpha1a
        else:
            alpha1b = -np.pi + alpha1a

        alpha2a = np.arctan(d/denom2)
        alpha2b = 0
        if alpha2a<0:
            alpha2b = np.pi + alpha2a
        else:
            alpha2b = -np.pi + alpha2a

        #the third equation serves as a test - this needs to be 0 for a correct solution
        def test(a, b):
            return (np.round(-c*np.sin(b)*np.cos(a) + f*np.cos(a)*np.cos(b) + d*np.sin(a) - length, 5) == 0 and self.inrange(a, b))

        possible = []
        #test all the possibilities
        if test(alpha1a, beta):
            if verbose:
                print("Option 1 is possible")
            possible.append([np.round(alpha1a/(2*np.pi)*360, rnd), np.round(beta/(2*np.pi)*360, rnd)])
        if test(alpha1b, beta):
            if verbose:
                print("Option 2 is possible")
            possible.append([round(alpha1b/(2*np.pi)*360, rnd), round(beta/(2*np.pi)*360, rnd)])
        if test(alpha2a, beta2):
            if verbose:
                print("Option 3 is possible")
            possible.append([round(alpha2a/(2*np.pi)*360, rnd), round(beta2/(2*np.pi)*360, rnd)])
        if test(alpha2b, beta2):
            if verbose:
                print("Option 4 is possible")
            possible.append([round(alpha2b/(2*np.pi)*360, rnd), round(beta2/(2*np.pi)*360, rnd)])
        
        #it is possible that nothing is possible
        if not possible:
            possible.append([np.nan, np.nan])
        
        #Only return the first option
        return np.array(possible[0]) 
        
            
class Detector(object):
    """A detector object stores it's own rotation with respect to the x-axis of the holder."""
    def __init__(self, TEM, diffrot=None, imrot=None, stemrot=None, diffpixcal=None, impixcal=None, stempixcal=None):
        #the detector needs to know which TEM it belongs to
        self.TEM = TEM
        
        #the rotation of the diffraction pattern
        self.diffrot = diffrot
        #This is the rotation of the image on any detector except HAADF
        self.imrot = imrot
        #this is the rotation of the HAADF image in STEM or the rotation of the Ronchigram on the TV or GIF CCD
        self.stemrot = stemrot
        
        #image mode size calibration
        self.impixcal = impixcal
        #diffraction mode size calibration
        self.diffpixcal = diffpixcal
        #Stem size calibration
        self.stempixcal = stempixcal
    
    def getTEM(self):
        return self.TEM
        
    def setCalibration(self, filename, mode = "diffraction", type = "rotation"):
        """This function sets the correct calibration dictionary in the right variables"""
        calib = open(filename).read()
        lin = re.compile("([1-9][0-9]*) *\t*: *\t*(\-?[0-9]+\.?[0-9]*) *\t*;")
        lst = re.findall(lin, calib)
        diction = dict([tuple(map(float, i)) for i in lst])
        try:
            if mode == 0 or mode=="diffraction" or mode == "diff" or mode == "d":
                if type ==0 or type == "rotation" or type == "r" or type =="rot":
                    self.diffrot = diction
                if type ==1 or type == "size" or type == "s" or type =="scale":
                    self.diffpixcal = diction
            if mode == 1 or mode=="imaging" or mode == "img" or mode=="i":
                if type ==0 or type == "rotation" or type == "r" or type =="rot":
                    self.imrot = diction
                if type ==1 or type == "size" or type == "s" or type =="scale":
                    self.impixcal = diction
            if mode == 2 or mode=="STEM" or mode == "stem" or mode=="s":
                if type ==0 or type == "rotation" or type == "r" or type =="rot":
                    self.stemrot = diction
                if type ==1 or type == "size" or type == "s" or type =="scale":
                    self.stempixcal = diction
            
        except:
            print("This is not a valid calibration mode or type.")
    
    def getRot(self, mode, setting):
        """This function gets the rotation angle of the detector. The rotation angle is how the absolute X and Y appear rotated on the detector. This must be tabulated and supplied."""
        try:
            if mode == 0 or mode=="diffraction" or mode == "diff" or mode == "d":
                return self.diffrot[setting]
            if mode == 1 or mode=="imaging" or mode == "img" or mode=="i":
                return self.imrot[setting]
            if mode == 2 or mode=="STEM" or mode == "stem" or mode=="s":
                return self.stemrot[setting]
        except:
            print("This detector seems to not have a calibration for this setting")
            return 0

    def getRealSize(self, val, mode, setting, binning=1):
        """Rescale a vector or scalar quantity by the calibration"""
        try:
            if mode == 0 or mode=="diffraction" or mode == "diff" or mode == "d":
                return val*self.diffpixcal[setting]*binning
            if mode == 1 or mode=="imaging" or mode == "img" or mode=="i":
                return val*self.impixcal[setting]*binning
            if mode == 2 or mode=="STEM" or mode == "stem" or mode=="s":
                return val*self.stempixcal[setting]*binning
        except:
            print("This detector seems to not have a calibration for these settings")
            return 0
    
    def detectorToAbs(self, vec, mode, setting):
        """Detector coordinates are transformed to absolute coordinates."""
        ang = self.getRot(mode, setting)
        #since ZR(ang) maps [100] -> [cos sin 0] which is the active rotation, we desire here the passive rotation -> to map [cos sin 0] (how X looks on the detector) -> [100] (the vector it should represent)
        #This involves doing the inverse opperation
        rm = ZR(ang).T
        return np.dot(rm, vec)
        
    def absToDetector(self, vec, mode, setting):
        """Absolute coordinates are transformed to detector coordinates."""
        ang = self.getRot(mode, setting)
        #ZR(ang) maps [100] -> [cos sin 0] This is precisely how the absolute X axis looks on the screen.
        rm = ZR(ang)
        return np.dot(rm, vec)
        
    def plotRealAxes(self, mode, setting):
        xaxis = np.array([1, 0, 0]);
        yaxis = np.array([0, 1, 0]);
        #Does the Z-axis go up the column or down the colunn ? This is to be determined from the way the tilts work. It determines which vector y is.
        
        rotx = self.absToDetector(xaxis, mode, setting)
        roty = self.absToDetector(yaxis, mode, setting)
        
        fig, ax = plt.subplots(1)
        
        #plot the x-axis as seen on the screen
        ax.arrow(0, 0, rotx[0], rotx[1], color = "red")
        #plot the y-axis as seen on the screen
        ax.arrow(0, 0, roty[0], roty[1], color = "green")
        
        sc = 1.1
        wd = 0.005
        #plot the x-axis as seen on the screen
        ax.arrow(0, 0, sc*xaxis[0], sc*xaxis[1], width = wd, color = "black")
        #plot the y-axis as seen on the screen
        ax.arrow(0, 0, sc*yaxis[0], sc*yaxis[1], width = wd, color = "black")
        
        #Label the detector axes
        sc2 = 1.2
        ax.text(sc2*xaxis[0], sc2*xaxis[1], "X", color = "black")
        ax.text(sc2*yaxis[0], sc2*yaxis[1], "Y", color = "black")

        plt.gca().set_aspect('equal', adjustable='box')
        
        #Label the tilt axes
        sc = 1.1
        ax.text(sc*rotx[0], sc*rotx[1], r"$\alpha$", color = "red")
        ax.text(sc*roty[0], sc*roty[1], r"$\beta$", color = "green")
        
        #turn off labels of axes
        ax.set_yticklabels([])
        ax.set_xticklabels([])
        
        sc2 = 1.3
        ax.set_xlim([-sc2, sc2])
        #this reverses the y-axis
        ax.set_ylim([sc2, -sc2])
        
        plt.show()
    
class Structure(object):
    """The structure class contains all the crystallographic data and calculations. The Crystal class can inherit from this class."""
    def __init__(self, a=1, b=1, c=1, alpha=90, beta=90, gamma=90):
        self.a=a
        self.b=b
        self.c=c
        self.alpha=alpha
        self.beta=beta
        self.gamma=gamma
        
        self.typ = self.getCrystalClass()
        
        #RM takes miller indices and transforms them to coordinates in the cartesian system stuck to the crystal
        self.RM = self.getRealMatrix()
        #invRM takes coordinates defined in the cartesian system stuck to the crystal and maps them to miller indices
        self.invRM = np.linalg.inv(self.RM)
        
        #RRM takes miller indices in recyprocal space and maps them to the cartesian system stuck to the crystal
        self.RRM = self.getRecypMatrix()
        #invRRM takes coordinates defined in the cartesian system stuck to the crystal and maps them to recyprocal space miller indices
        self.invRRM = np.linalg.inv(self.RRM)
        
        #direct metric tensor
        self.G = self.getGmatrix()
        #the recyprocal metric tensor is simply the inverse
        self.recypG = np.linalg.inv(self.G)
    
    def changeCrystallography(self, a, b, c, alpha, beta, gamma):
        self.__init__(a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma)
    
    def getCrystallography(self):
        return self.a, self.b, self.c, self.alpha, self.beta, self.gamma
        
    def millerToCartesian(self, vec, typ = "real"):
        """This function returns the coordinates in the cartesian coordinate system stuck to the crystal given a set of miller indices as columns or an array. Standard it will be assumed miller indices as defined by the real coordinate system, but recyp is also valid and must be supplied as typ = 'recyp' as argument"""
        if typ=="real":
            return np.dot(np.array(self.RM), np.array(vec))
        elif typ=="recyp":
            return np.dot(np.array(self.RRM), np.array(vec))
        else:
            return None
    
    def cartesianToMiller(self, vec, typ = "real"):
        """This function returns the miller indices given coordinates in the cartesian system stuck to the crystal. Standard it will be assumed miller indices as defined by the real coordinate system, but recyp is also valid and must be supplied as typ = 'recyp' as argument"""
        if typ=="real":
            return np.dot(np.array(self.invRM), np.array(vec))
        elif typ=="recyp":
            return np.dot(np.array(self.invRRM), np.array(vec))
        else:
            return None
            
    def recypToReal(self, vec):
        """get the coordinates of the real coordinates corresponding to a certain recyprocal space vector"""
        #first convert to cartesian
        car = self.millerToCartesian(vec, typ="recyp")
        #convert the cartesian to real space miller
        return self.cartesianToMiller(car)
    
    def realToRecyp(self, vec):
        """get the coordinates of the recyprocal coordinates corresponding to a certain real space vector"""
        #first convert to cartesian
        car = self.millerToCartesian(vec)
        #convert the cartesian to recyp space miller
        return self.cartesianToMiller(car, typ="recyp")
        
    def getCrystalClass(self):
        typ = ""
        if self.a==self.b and self.b==self.c and self.alpha==90 and self.beta==90 and self.gamma==90:
            typ = "cubic"
        elif self.a==self.b and self.alpha==90 and self.beta==90 and self.gamma==90:
            typ = "tetragonal"
        elif self.alpha==90 and self.beta==90 and self.gamma==90:
            typ = "orthorhombic"
        elif self.a==self.b and self.gamma==120 and self.alpha==90 and self.beta==90:
            typ = "hexagonal"
        elif self.a==self.b and self.b==self.c and self.alpha==self.beta and self.beta==self.gamma:
            typ = "rhombohedral"
        elif self.gamma == 90 and self.alpha==90:
            typ = "monoclinic"
        else:
            typ = "triclinic"
        return typ
        
    def getRealMatrix(self, rnd=10):
        """returns in the general triclinic case the real space vectors in terms of cartesian axes
        The cartesian axes stuck to the crystal is defined as: a is on x-axis, a-b plane is on x-y plane"""
        a=float(self.a)
        b=float(self.b)
        c=float(self.c)
        alpha=float(self.alpha)/360*2*np.pi
        beta=float(self.beta)/360*2*np.pi
        gamma=float(self.gamma)/360*2*np.pi
        av=np.array([a, 0, 0])
        bv=np.array([round(b*np.cos(gamma), rnd), round(b*np.sin(gamma), rnd), 0])
        c1=round(c*np.cos(beta), rnd)
        c2=round(c*(np.cos(alpha)-np.cos(beta)*np.cos(gamma)), rnd)
        c3=round(np.sqrt(c**2-c1**2-c2**2), rnd)
        cv=np.array([c1, c2, c3])
        matR=np.array([av, bv, cv]).T
        return matR
        
    def getRecypMatrix(self, rnd=10):
        """Matrix containing the basis vectors of the recyprocal lattice. Based on the real matrix."""
        a1 = self.RM[:, 0]
        a2 = self.RM[:, 1]
        a3 = self.RM[:, 2]
        vol = np.dot(a1, np.cross(a2, a3))
        astar = np.cross(a2, a3)/vol
        bstar = np.cross(a3, a1)/vol
        cstar = np.cross(a1, a2)/vol
        return np.array([astar, bstar, cstar]).T
        
    def getGmatrix(self, rnd=10):
        """metric tensor. This is the same as the RW.T * RW see definition on p 79 marc de graef"""
        a=float(self.a)
        b=float(self.b)
        c=float(self.c)
        alpha=float(self.alpha)/360*2*np.pi
        beta=float(self.beta)/360*2*np.pi
        gamma=float(self.gamma)/360*2*np.pi
        mat=np.array([[a**2, round(a*b*np.cos(gamma), rnd), round(a*c*np.cos(beta), rnd)],[round(b*a*np.cos(gamma), rnd), b**2, round(b*c*np.cos(alpha), rnd)],[round(c*a*np.cos(beta), rnd), round(c*b*np.cos(alpha), rnd), c**2]])
        return mat
        
    def getSym(self, vec, checksame = False, typ = "real", unc = 1e-9):
        """Returns all the possible permutations of miller indices/vectors as columns in a matrix. Only returns the vectors that are crystallographically identical if checksame is set True. Real and recyprocal vectors are both valid"""
        vec=np.array(vec) #make sure vec is an array. This way a list is also accepted.
        tmpmat = np.matrix([vec,-vec]).T #-vec and vec can already be entered as columns of the permutation matrix
        for i in range(3): #To make the permutations, the elements must be swapped.
            val1 = i
            val2 = (i+1)%3
            val3 = (i+2)%3
            vn = []
            vn.append(np.array([vec[val1], vec[val2], vec[val3]])) #depending on i, the values are switched. 8 extra vectors per permutations must possibly be added: the one only with switched numbers.
            vn.append(np.array([-vec[val1], vec[val2], vec[val3]])) #the one with the first element negative
            vn.append(np.array([vec[val1], -vec[val2], vec[val3]])) #the one with the second element negative
            vn.append(np.array([vec[val1], vec[val2], -vec[val3]])) #the one with the third element negative
            for j in vn: #all are checked to see whether they already exist in the matrix
                if not isExist(tmpmat, j): #if they don't they get added
                    tmpmat = np.c_[tmpmat, j]
                if not isExist(tmpmat, -j):
                    tmpmat = np.c_[tmpmat, -j]
                    
        if checksame and self.typ!="cubic":
            #in case we only want those vectors that are crystallographically the same length. If the matrix is cubic we know we don't have to eliminate anything.
            #tst is the length of the supplied vector
            tst = self.getVectorLength(vec, typ=typ)
            #others is the list of lengths of "equivalent" vectors
            others = self.getVectorLength(tmpmat.T, typ=typ)
            #get all the columns from tempmat where the difference between the length of the supplied vector and the equivalent vectors is negligible
            tmpmat2 = tmpmat[:, abs(others-tst)<unc]
            tmpmat = tmpmat2
        
        return tmpmat
        
    def getVectorDot(self, vec1, vec2, typ="real"):
        """Get the dot product between two miller indices. The vectors can be defined in real space or recyprocal space but must both be defined in the same space. Otherwise the normal dot product is applicable. A list of COLUMN vectors can also be supplied or 1 array"""
        if typ=="real":
            return np.array(np.dot(vec1.T,np.dot(self.G, vec2)))
        elif typ=="recyp":
            return np.array(np.dot(vec1.T,np.dot(self.recypG, vec2)))
        else:
            return None
    
    def getVectorLength(self, vec, typ="real"):
        """The length of a certain real space vector or plane normal depending on the type,  A list of ROW vectors can also be supplied, an array of 1 dimension will be returned"""
        return np.array(np.matrix(np.sqrt(self.getVectorDot(vec, vec, typ=typ))).diagonal())
            
    def getVectorAngle(self, vec1, vec2, typ="real", units = "radians"):
        #! still some strange behavior when testing when vec1 or vec2 is a one dimensional array and the other is not. Otherwise it works perfectly. This has to do with the way the division happens with matrix/arrays. Fix later.
        """The angle between two vectors in real or recyprocal space,  A list of Column vectors can also be supplied, an array of 2 dimensions will be returned"""
        num= self.getVectorDot(vec1, vec2, typ=typ)
        denom = np.outer(self.getVectorLength(vec1, typ=typ), self.getVectorLength(vec2, typ=typ))
        angls=  np.arccos(np.divide(num, denom))
        if units =="radians":
            return angls
        if units =="degrees":
            return angls/(2*np.pi)*360
        else:
            print("Those units aren't valid.")
            return None

        
class Crystal(object):
    
    def __init__(self, stage, structurename):
        self.structure = getStructure(structurename) #!note the globally defined getStructure as opposed to self.getstructure
        self.stage = stage #the crystal needs to know which stage it belongs to
        
        #orient matrix takes vectors from the stage coordinate system to the cartesian axes stuck to the crystal
        self.orient = np.identity(3)
        #the inverse of the orientation matrix does not need to be defined -> orthogonal so .T is the inverse already. 
        #self.orient.T takes vectors from the cartesian system stuck on the crystal to the stage coordinates
    
    def changeCrystallography(self, *args):
        pass #can't change the crystallography of the structure from the crystal so it is overwritten here
    
    def getCrystallography(self):
        return self.structure.getCrystallography()
    
    def getStructure(self):
        return self.structure
    
    def cartesianToStage(self, vec):
        """Coordinates from the cartesian system stuck to the crystal  are converted to screen coordinates."""
        return np.dot(self.orient, vec)
    
    def stageToCartesian(self, vec):
        """Coordinates from the stage system  are converted to cartesian (stuck to crystal) coordinates."""
        return np.dot(self.orient.T, vec)
        
    def millerToDetector(self, vec, detector, mode, setting, typ = "real", verbose = False):
        """Miller indices converted to screen coordinates. Miller indices can be real or recyp."""
        #convert miller to cartesian
        car = self.millerToCartesian(vec, typ=typ)
        #convert cartesian to stage
        stg = self.cartesianToStage(car)
        #convert stage to detector
        abscords = self.stage.stageToAbs(stg) #change stage coordinates to absolute
        dec = self.stage.getTEM().getDetector(detector) #find the detector
        deccoords = dec.absToDetector(abscords, mode, setting)
        if verbose:
            print("Miller (%s)" %(typ))
            print(vec)
            print("From miller indices to Cartesian coordinates:")
            print(car)
            print("From crystal cartesian to stage:")
            print(stg)
            print("From stage to absolute:")
            print(abscords)
            print("From absolute to detector:")
            print(deccoords)
        return deccoords #change absolute coordinates to detector
        
        
    def detectorToMiller(self, vec, detector, mode, setting, typ="real", verbose=False):
        """screen coordinates converted to Miller indices . Miller indices can be real or recyp."""
        #convert detector to stage
        dec = self.stage.getTEM().getDetector(detector) #find the detector
        realcords = dec.detectorToAbs(vec, mode, setting) #change to absolute coordinates
        stagecoords = self.stage.absToStage(realcords) #change absolute to stage coordinates
        #convert stage to cartesian
        car = self.stageToCartesian(stagecoords)
        #convert cartesian to miller
        cryscoords = self.cartesianToMiller(car, typ=typ)
        if verbose:
            print("Detector coordinates:")
            print(vec)
            print("From detector to absolute coordinates:")
            print(realcords)
            print("From absolute coordinates to stage coordinates:")
            print(stagecoords)
            print("From stage coordinates to cartesian coordinates (crystal):")
            print(car)
            print("From cartesian coordinates (crystal) to miller (%s):" %(typ))
            print(cryscoords)
        return cryscoords
    
    def getZone(self, typ="real", verbose=False, integer = False):
        realcords = miller(0,0,1)
        stagecoords = self.stage.absToStage(realcords)
        car = self.stageToCartesian(stagecoords)
        cryscoords = self.cartesianToMiller(car, typ=typ)
        if verbose:
            print("Absolute coordinates:")
            print(realcords)
            print("Stage coordinates:")
            print(stagecoords)
            print("Cartesian coordinates (crystal):")
            print(car)
            print("Miller (%s):" %(typ))
            print(cryscoords)
        if integer:
            return integerRep(cryscoords)
        else:
            return cryscoords
    
    def calcOrient(self, za, ref, ang, detector, mode, setting, acur = 1e-9):
        """The crystal has a certain orientation with respect to the stage. The orientation is most easily found when studying the zone axis (real space) || screen z-axis and a visible reflection on the detector (recyprocal space) defined by an angle from the detector x-axis.
        The orientation must be simply a rotation matrix and hence defines how the cartesian system stuck to the crystal is rotated compared to the stage coordinate system."""
        #first check that za (real space) and ref (recyprocal space) are indeed perpendicular. This follows the normal h*u + k*v + l*w = 0 relationship valid for any crystal system.
        if abs(np.dot(za, ref))<acur:
            #turn angle from degrees to radians
            ang = ang/360*2*np.pi
            
            #calculate the cartesian equivalents of the vectors
            zaC = self.millerToCartesian(za)
            refC = self.millerToCartesian(ref, typ = "recyp")
            #normalize the vectors
            zaC = zaC/np.linalg.norm(zaC)
            refC = refC/np.linalg.norm(refC)
            depC = np.cross(zaC, refC)
            #the vectors of the crystal to be transformed
            mat1 = np.array([zaC, refC, depC]).T
            
            #the matrix of corresponding detector vectors
            c1 = np.array([0,0,1])
            c2 = np.array([np.cos(ang), np.sin(ang), 0])
            c3 = np.array([np.cos(ang+np.pi/2), np.sin(ang+np.pi/2), 0])
            mat2 = np.array([c1, c2, c3]).T
            
            #these must be converted to stage coordinates.
            dec = self.stage.getTEM().getDetector(detector) #find the detector
            realcords = dec.detectorToAbs(mat2, mode, setting) #change to absolute coordinates
            stagecoords = self.stage.absToStage(realcords)
            
            
            #the rotation matrix needs to turn mat 1 (cartesian vectors stuck to crystal) into stagecoords (stage vectors). Therefore
            ormat = np.dot(stagecoords, np.linalg.inv(mat1))
            self.setOrient(ormat)
            #multiplying by ormat goes from crystal cartesian vector to stage coordinates, ormat.T (inverse) goes from stage to cartesian.
            return ormat
        else:
            print("ZA vector and reflection vector are not perpendicular")
            return np.identity(3)
    
    def setOrient(self, ormat):
        self.orient=ormat
    
    def calcSampleTilt(self, g, n=1, units = "degrees"):
        """Sample tilt necessary to put reflection ng into two-beam condition assuming you start from zone axis"""
        #get the length of vector g in nm-1
        length = self.getVectorLength(g, typ="recyp")[0,0]
        #get the length of the radius of the ewald sphere
        #convert to nm-1
        #K0 = 507.93 for kv 300
        K0 = self.stage.getTEM().getEwaldR(units = "nm")
        ans = np.pi/2 - np.arccos(n*length/2/K0)
        if units == "degrees":
            return ans/(2*np.pi)*360
        elif units == "mrad":
            return ans*1000
        elif units == "radians":
            return ans
        else:
            print("Units not recognized")
            return 0
    
    def calcBeamTilt(self, g, n, k, units = "degrees"):
        """This calculates how much the beam needs to be tilted when beam ng is active in order for beam kg to get on the optical axis"""
        #get the length of vector g in nm-1
        length = self.getVectorLength(g, typ="recyp")[0,0]
        #get the length of the radius of the ewald sphere
        #convert to nm-1
        #K0 = 507.93
        K0 = self.stage.getTEM().getEwaldR(units = "nm")
        
        yot = self.calcSampleTilt(g, n, units = "radians")
        thet = np.pi/2 - np.arccos(k*length*np.cos(yot)/K0)
        if units == "degrees":
            return thet/(2*np.pi)*360
        elif units == "mrad":
            return thet*1000
        elif units == "radians":
            return thet
        else:
            print("Units not recognized")
            return 0
        
    def calcNfromAngles(self, g, sa, ba):
        """This method calculates the active 'n' after a given sample angle tilt and beam angle tilt"""
        #get the length of vector g in nm-1
        length = self.getVectorLength(g, typ="recyp")[0,0]
        #get the length of the radius of the ewald sphere
        #convert to nm-1
        #K0 = 507.93
        K0 = self.stage.getTEM().getEwaldR(units = "nm")
        
        #takes sa and ba as radian quantitites!
        #this calculation was done on paper, sa = sample angle = phi, ba = beam angle = theta. This was purely geometric.
        ang = np.pi/2-sa-ba
        dst = 2*K0*np.cos(ang)
        return dst/length

    def calcNfromTilts(self, g, n, k):
        """This method calculates what the active 'n' value is after a sample tilt and beamtilt"""
        #n is the "reflection" tilted to in two beam condition
        #k is the "reflection brought to the optical axis"
        yot = self.calcSampleTilt(g, n, units = "radians")
        thet = self.calcBeamTilt(g, n, k, units = "radians")
        ans = self.calcNfromAngles(g, yot, thet)
        return ans

    def calcWB2beam(self, g, k, i):
        """This method calculates what the correct 2-beam tilt condition should be (n) if we want the ewald sphere to cut through ig after the beam tilt so that reflection k is centered"""
        #get the length of vector g in nm-1
        length = self.getVectorLength(g, typ="recyp")[0,0]
        #get the length of the radius of the ewald sphere
        #convert to nm-1
        #K0 = 507.93
        K0 = self.stage.getTEM().getEwaldR(units = "nm")
        
        #the angle between the beam and the reflections can be easily calculated 
        ang = np.arccos(i*length/2/K0)
        #need to solve the following equation numericaly:
        def eq(phi):
            return np.arccos(k*length*np.cos(phi)/K0)-ang-phi
        
        phi = fsolve(eq, 0)
        #from the known phi, the initial 2-beam can be determined
        lent = 2*K0*np.sin(phi)
        
        return (lent/length)[0]

    def calcSg(self, g, n, k, n2, verbose=False):
        #get the length of vector g in nm-1
        length = self.getVectorLength(g, typ="recyp")[0,0]
        #get the length of the radius of the ewald sphere
        #convert to nm-1
        #K0 = 507.93
        K0 = self.stage.getTEM().getEwaldR(units = "nm")
        
        #n2 is used for calculating sg, usually this will be k, the reflection used for imaging 
        yot = self.calcSampleTilt(g, n, units = "radians")
        thet = self.calcBeamTilt(g, n, k, units = "radians")
        #slope and intercept of the line
        slop = np.tan((np.pi-thet)/2)
        intc = n2*length*(np.sin(yot)-slop*np.cos(yot))
        #solve the intersection with the circle
        c = intc/K0 + slop*np.sin(thet)-np.cos(thet)
        a = -slop
        b = 1
        #v = np.arcsin(c/np.sqrt(a**2+b**2))-np.arcsin(a/np.sqrt(a**2+b**2))
        fp = np.arcsin(-c/np.sqrt(a**2+b**2))
        sp = np.arcsin(a/np.sqrt(a**2+b**2))
        v1 = fp - sp - np.pi
        v2 = -fp - sp
        
        #usually there will be two possible intersections, we want the one below the center of the ewald sphere
        if np.sin(v1)<0:
            v = v1
        elif np.sin(v2)<0:
            v = v2
        else:
            print("Something went terribly wrong.")
            v=0
        
        #print(fp/(2*np.pi)*360)
        #print(sp/(2*np.pi)*360)
        #print(ad/(2*np.pi)*360)
        
        #what are the points on the ewald sphere
        x = K0*(np.cos(v)+np.sin(thet))
        y = K0*(np.sin(v)+np.cos(thet))
        #choose the point where y is lower than the center of the ewald sphere
        #location of x and y of n2g
        x0 = n2*length*np.cos(yot)
        y0 = n2*length*np.sin(yot)
        
        if verbose:
            print("v: %s" %(v/np.pi/2*360))
            print("x_gn2: %s" %x0)
            print("y_gn2: %s" %y0)
            print("slope: %s" %slop)
            print("intercept: %s" %intc)
            print("x_int: %s" %x)
            print("y_int: %s" %y)
        
        #the distance:
        if y0<y: #the ewald sphere is above the n2g point - the point is outside - negative sg
            return round(-np.sqrt((x0-x)**2+(y0-y)**2), 10)
        else: #the reflection is inside - positive sg
            return round(np.sqrt((x0-x)**2+(y0-y)**2), 10)    
    
    def getAlphaBetaMiller(self, vec, typ = "real", verbose = False, **kwargs):
        """For miller indices or a list of miller indices calculate the alpha and beta one needs to go to have those vectors on the Z-axis"""
        car = self.millerToCartesian(vec, typ=typ)
        stg = self.cartesianToStage(car)
        if verbose:
            print("Cartesian:")
            print(car)
            print("Stage:")
            print(stg)
        return self.stage.calcAlphaBeta(stg, verbose=verbose, **kwargs)
    
    def getAlphaBetaWBeam(self, g, n=1, k=3, outzone = 8):
        pass
    
    def getAlphaBetac2Beam(self, g, n=1, outzone = 8, za = "current", rnd=5, verbose = False):
        """Calculate alpha and beta for getting g in a two-beam condition (centered). The beam tilt that is necessary is not considered.
        g must be a recyprocal space vector. You must start from zone axis!"""
        #####To implement later -> can give automatically 2 options for tilting 8 degrees away in either direction.
        #test if g is perpendicular to the zone axis
        if str(za)=="current":
            za = self.getZone(typ="real")
        #ELSE a ZA must be supplied! Otherwise error.
        #check if the zone axis and reflection are perpendicular
        if abs(np.round(np.dot(za, g), rnd))==0:
            #get the absolute coordinates of the zone axis
            zax = self.millerToAbs(za)
            #the z-axis must be rotated outzone. The rotation axis is the reflection.
            rotax = self.millerToAbs(g, typ="recyp")
            zdash = np.dot(axisAngle(rotax, outzone), zax) #the new z-axis that we should see when we excite te systematic row
            #the rotation axis  to tilt to the two beam condition is now the cross product of this and the reflection
            rotax2 = np.cross(rotax, zdash)
            #calculate the angle necessary to tilt to -g
            ang = self.calcSampleTilt(g, n = -n)
            #the zdash axis is rotated by this amount around rotax2
            zab = np.dot(axisAngle(rotax2, ang), zdash)
            #this axis must be converted to miller
            zabmil = self.absToMiller(zab, typ = "recyp")
            #This axis must then be brought to the z-axis
            if verbose:
                print("Absolute coordinates of the zone axis:")
                print(zax)
                print("Z axis rotated out of zone:")
                print(zdash)
                print("Z axis rotated to the two-beam:")
                print(zab)
                print("Vector to Z axis in miller indices")
                print(zabmil)
            #return the alpha and beta
            return self.getAlphaBetaMiller(zabmil, typ="recyp", verbose=verbose)
        else:
            print("This reflection and zone axis are not perpendicular")
            return np.array([0,0])
            
    
    def millerToAbs(self, vec, typ="real"):
        car = self.millerToCartesian(vec, typ=typ)
        stg = self.cartesianToStage(vec)
        return self.stage.stageToAbs(stg)
        
    def absToMiller(self, vec, typ="real"):
        stg = self.stage.absToStage(vec)
        car = self.stageToCartesian(stg)
        return self.cartesianToMiller(car, typ=typ)
    
    def millerToCartesian(self, vec, typ = "real"):
        return self.structure.millerToCartesian(vec, typ)
    
    def cartesianToMiller(self, vec, typ = "real"):
        return self.structure.cartesianToMiller(vec, typ)
            
    def recypToReal(self, vec):
        return self.structure.recypToReal(vec)
    
    def realToRecyp(self, vec):
        return self.structure.realToRecyp(vec)
        
    def getCrystalClass(self):
        return self.structure.getCrystalClass()
        
    def getRealMatrix(self, rnd=10):
        return self.structure.getRealMatrix()
        
    def getRecypMatrix(self, rnd=10):
        return self.structure.getRecypMatrix()
        
    def getGmatrix(self, rnd=10):
        return self.structure.getGmatrix()
        
    def getSym(self, vec, checksame = False, typ = "real", unc = 1e-9):
        return self.structure.getSym(vec, checksame, typ, unc)
        
    def getVectorDot(self, vec1, vec2, typ="real"):
        return self.structure.getVectorDot(vec1, vec2, typ)
    
    def getVectorLength(self, vec, typ="real"):
        return self.structure.getVectorLength(vec, typ)
            
    def getVectorAngle(self, vec1, vec2, typ="real"):
        return self.structure.getVectorAngle(vec1, vec2, typ)
        
    def plotCrystalonDetector(self, detector, mode, setting, vecincl = None, typ="real", plotaxes = True, verbose=False):
        d001 = self.getSym(miller(1, 0, 0), typ=typ)
        d111 = self.getSym(miller(1, 1, 1), typ=typ)
        d011 = self.getSym(miller(1, 1, 0), typ=typ)
        toplotd001 = self.millerToDetector(d001, detector, mode, setting, typ=typ, verbose=verbose)
        toplotd111 = self.millerToDetector(d111, detector, mode, setting, typ=typ, verbose=verbose)
        toplotd011 = self.millerToDetector(d011, detector, mode, setting, typ=typ, verbose=verbose)
        
        #make the actual plot
        fig, ax = plt.subplots(1)
        stereographicCanvas(ax)
        
        #plotStereographic(d001, d001, ax, verbose = verbose, s=100, c="b")
        #plotStereographic(d111, d111, ax, verbose = verbose, s=100, c="b")
        #plotStereographic(d011, d011, ax, verbose = verbose, s=100, c="g")
        
        if plotaxes:
            xaxis = np.array([1, 0, 0]);
            yaxis = np.array([0, 1, 0]);
            #Does the Z-axis go up the column or down the colunn ? This is to be determined from the way the tilts work. It determines which vector y is.
            
            dec = self.stage.getTEM().getDetector(detector) #find the detector
            rotx = dec.absToDetector(xaxis, mode, setting)
            roty = dec.absToDetector(yaxis, mode, setting)
            
            sc = 1.1
            wd = 0.005
            #plot the x-axis as seen on the screen
            ax.arrow(0, 0, sc*rotx[0], sc*rotx[1],  width = wd, color = "red")
            #plot the y-axis as seen on the screen
            ax.arrow(0, 0, sc*roty[0], sc*roty[1],  width = wd, color = "green")
            
            #Label the axes
            sc2 = 1.35
            ax.text(sc2*rotx[0], sc2*rotx[1], r"$\alpha$", color = "red", fontsize=14)
            ax.text(sc2*roty[0], sc2*roty[1], r"$\beta$", color = "green", fontsize=14)
        
        plotStereographic(toplotd001, d001, ax, verbose = verbose, s=100, c="r")
        plotStereographic(toplotd111, d111, ax, verbose = verbose, s=100, c="b")
        plotStereographic(toplotd011, d011, ax, verbose = verbose, s=100, c="g")
        
        plt.show()
        
def normMatrix(vec):
    """normalize a vector or columns of a matrix. Returns the matrix with normalized columns. If an array was passed a matrix with one column is returned."""
    #make vec a matrix in case it is one dimensional
    if vec.ndim==1:
        vec = np.matrix(vec).T
    #normalize the column vectors
    invl = 1/np.matrix(np.linalg.norm(vec, axis=0))
    lrp = np.repeat(invl, vec.shape[0], axis=0)
    return np.array(np.multiply(lrp, vec))

def calcStereo(vec, rnd=1e-9):
    """takes a matrix containing real vectors (x, y, z) in the columns and calculates the (xproj, yproj) of te stereographic projection."""
    #normalize all columns. Assume the sphere has unit radius.
    vec = np.matrix(normMatrix(vec))
    #the vectors lying above or on the plane - they must be projected using pole [00-1]
    whrpos = np.where(vec[2,:]>=-rnd)[1]
    vecpos = vec[:,whrpos]
    #projection scaling factor
    u = np.matrix(1/(vecpos[2,:]+1))
    u = np.repeat(u, 2, axis=0)
    vecposproj = np.multiply(vecpos[0:2, :], u)
    
    #the vectors lying below or on the plane - they must be projected using pole [001]                
    whrneg = np.where(vec[2,:]<=rnd)[1]
    vecneg = vec[:,whrneg]
    #projection scaling factor
    u = np.matrix(1/(vecneg[2,:]-1))
    u = np.repeat(u, 2, axis=0)
    vecnegproj = np.multiply(vecneg[0:2, :], u)
    return whrpos, vecposproj, whrneg, vecnegproj
    

def isExist(mat, vec):
    """Checks whether an array vec exists as a column in matrix mat"""
    #if vec is an array turn it into a 1 column matrix
    if type(vec)==np.ndarray:
        vec=np.matrix(vec).T
    #isit is a row vector that states whether vec is the same as a certain column in mat
    isit=np.all(mat==vec, axis=0)
    #if there is any true value in isit, it must mean the column exists in mat
    if True in isit:
        return True
    else:
        return False

def integerRep(vec, decimals = 9):
    """Make a vector or array of column vectors more presentable by giving their closest integer representation."""
    vec = np.array(vec)
    vec = np.round(vec, decimals = decimals)
    return getIntRepAny(vec)
    
def lcm(a, b):
    """Get least common multiple of 2 numbers"""
    a = int(a)
    b = int(b)
    return a*b/math.gcd(a,b)

def lcm3(lst):
    """Get least common multiple of 3 numbers"""
    return lcm(lst[0], lcm(lst[1], lst[2]))

def getDenominator(x):
    """Get the denominator of a decimal number"""
    return Fraction(str(x)).limit_denominator().denominator

def getIntRep(vect):
    """Get the integer representation of a vector of decimal numbers"""
    vect = vect/np.min(np.abs(vect[np.nonzero(vect)]))
    g=[]
    for i in vect.tolist():
        i = abs(i)
        if i!=0.0:
            g.append(getDenominator(i))
        else:
            g.append(1)
    return vect*lcm3(g)

def getIntRepAny(vect):
    """Get the integer representation of a vector or list of vectors of decimal numbers"""    
    if vect.ndim ==1:
        return getIntRep(vect)
        
    elif vect.ndim ==2:
        toret = np.zeros(vect.shape)
        for i in range(vect.shape[1]):
            vecy = vect[:,i]
            toret[:,i] = getIntRep(vecy)
        return toret   

    else:
        print("The data is not supplied in the right format")
        return vect

def toAllCols(func, lstvc, **kwargs):
    """Apply a function to all columns of a vector. This is nice to 'vectorize' functions that only take one dimensional arrays as input"""
    if lstvc.ndim==2:
        cols = np.size(lstvc, 1)
        ans = np.array([])
        for i in range(cols):
            act = lstvc[:,i]
            #add the first entry
            if not ans.size:
                ans = np.array(func(act, **kwargs))
            else:
                ans = np.vstack([ans, np.array(func(act, **kwargs))])
        #if the answer is only one column, return a one dimensional array
        if np.size(ans, 1)==1:
            return ans[:,0]
        else: #otherwise return a matrix
            return ans.T
    elif lstvc.ndim==1:
        return func(lstvc, **kwargs)
    else:
        print("This input format is not of the correct dimension")
        return lstvc
        
def stereographicCanvas(ax):
    """The canvas is set up for stereographic projection"""
    #plot the unit circle
    deg = np.linspace(0, 2*np.pi, 100)
    xcirc = np.cos(deg)
    ycirc = np.sin(deg)
    ax.plot(xcirc, ycirc, c="black")
    #plot the lines of the axes
    ax.plot(np.array([0, 0]), np.array([-1, 1]), c="black")
    ax.plot(np.array([-1, 1]), np.array([0, 0]), c="black")
    
    #axes settings
    xyscale = 1.5
    ax.set_xlim([-xyscale, xyscale])
    ax.set_ylim([xyscale, -xyscale])
    
    #Plot the detector X and Y axes
    xaxis = np.array([1, 0, 0]);
    yaxis = np.array([0, 1, 0]);
    
    sc = 1.1
    wd = 0.005
    #plot the x-axis as seen on the screen
    ax.arrow(0, 0, sc*xaxis[0], sc*xaxis[1], width = wd, color = "black")
    #plot the y-axis as seen on the screen
    ax.arrow(0, 0, sc*yaxis[0], sc*yaxis[1], width = wd, color = "black")
    
    #Label the axes
    sc2 = 1.35
    ax.text(sc2*xaxis[0], sc2*xaxis[1], "X", color = "black")
    ax.text(sc2*yaxis[0], sc2*yaxis[1], "Y", color = "black")

    plt.axis("off")
    plt.gca().set_aspect('equal', adjustable='box')
    plt.draw()

def plotStereographic(vecs, labs, ax, verbose=False,  **kwargs):
    """The annotated poles of certain vectors are plotted in the stereogram. Vecs should be the vectors after any transformation/rotation"""
    projpospos, projpos, projnegpos, projneg  = calcStereo(vecs)
    if verbose:
        print("The coordinates that will be plot:")
        print(projpos)
    labs = labs[:, projpospos]
    x = projpos[0, :].T
    y = projpos[1, :].T
    for i, xy in enumerate(zip(x, y)):
        ax.annotate("%s" %(labs[:,i].T), xy=xy, textcoords = "data")
    ax.scatter(x,y, **kwargs)       
        