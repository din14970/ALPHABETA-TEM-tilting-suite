"""These scripts serve as the backbone of TEM and crystallography calculations. They mainly serve for tilting help in the TEM."""

import numpy as np
import math as ma
import matplotlib.pyplot as plt
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
    if name is None or name == "":
        name = "structure"+str(len(structures))
    structures[name] = Structure(name, **kwargs)
    return structures[name]
    
def removeStructure(name):
    #we also have to get rid of all the potential crystals that have this structure
    #loop over the microscopes
    for i in list(microscopes.keys()):
        micr = getMicroscope(i)
        for j in list(micr.stages.keys()):
            stag = micr.getStage(j)
            for k in list(stag.crystals.keys()):
                crys = stag.getCrystal(k)
                if crys.structure.name == name:
                    stag.removeCrystal(crys.name)
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

def addMicroscope(name = None, **kwargs):
    if name is None or name == "":
        name = "microscope"+str(len(microscopes))
    microscopes[name] = TEM(name, **kwargs)
    return microscopes[name]
    
def removeMicroscope(name):
    del microscopes[name]
    
def getMicroscope(name):
    return microscopes[name]

def changeMicroscopeName(oldname, newname):
    microscopes[newname]=microscopes[oldname]
    del microscopes[oldname]

def saveSession(filname = "session.txt"):
    #Open file
    f = open(filname+".txt", "w")
    #Loop through all structures
    for i in list(structures.keys()):
        i = getStructure(i)
        a, b, c, alpha, beta, gamma = i.getCrystallography()
        str = "structure;%s;%s;%s;%s;%s;%s;%s;" %(i.name, a, b, c, alpha, beta, gamma)
        f.write(str+"\n")
    #Loop through all microscopes
    for i in list(microscopes.keys()):
        i = getMicroscope(i)
        str = "microscope;%s;%s;" %(i.name, i.getKv())
        f.write(str+"\n")
        #loop through all stages
        for j in list(i.stages.keys()):
            j = i.getStage(j)
            str = "stage;%s;%s;%s;%s;%s;%s;%s;%s;%s;%s;" %(i.name, j.name, j.getAlpha(), j.getBeta(), j.getAlphaMin(), j.getAlphaMax(), j.getBetaMin(), j.getBetaMax(), j.alpharev, j.betarev)
            f.write(str+"\n")
            #loop through all the crystals
            for k in list(j.crystals.keys()):
                k = j.getCrystal(k)
                str = "crystal;%s;%s;%s;%s;%s;%s;" %(i.name, j.name, k.name, k.getStructure().name, k.dumpOrient(), k.getComment())
                f.write(str+"\n")
        #loop through all detectors
        for j in list(i.detectors.keys()):
            j = i.getDetector(j)
            str = "detector;%s;%s;imaging;%s;diffraction;%s;stem;%s;" %(i.name, j.name, j.getCalFileName("imaging"), j.getCalFileName("diffraction"), j.getCalFileName("stem"))
            f.write(str+"\n")
            
    f.close()

def clearSession():
    """Delete all objects (microscopes and structures)"""
    microscopes.clear()
    structures.clear()
        
def loadSession(filname = "session.txt"):
    #Open file
    f = open(filname)
    text = f.readlines()
    
    for i in text:
        info = i.replace("\n","").split(";")
        if info[0].lower() == "structure":
            ##create a structure
            addStructure(name = info[1].strip(), a = float(info[2]), b = float(info[3]), c = float(info[4]), alpha = float(info[5]), beta = float(info[6]), gamma = float(info[7]))
        elif info[0].lower() == "microscope":
            addMicroscope(name = info[1].strip(), kv=float(info[2]))
        elif info[0].lower() == "stage":
            micr = getMicroscope(info[1].strip())
            micr.addStage(name = info[2].strip(), alpha = float(info[3]), beta = float(info[4]), alphamin = float(info[5]), alphamax = float(info[6]), betamin = float(info[7]), betamax = float(info[8]), alpharev = info[9]=="True", betarev = info[10]=="True")
        elif info[0].lower() == "detector":
            micr = getMicroscope(info[1].strip())
            ccd = micr.addDetector(name=info[2].strip())
            #check if the calibrations were ever made
            if info[4]: #imaging
                ccd.setCalibration(info[4], mode = "imaging", type = "r")
            if info[6]: #diffraction
                ccd.setCalibration(info[6], mode = "diffraction", type = "r")
            if info[8]: #stem
                ccd.setCalibration(info[8], mode = "stem", type = "r")
        elif info[0].lower() == "crystal":
            micr = getMicroscope(info[1].strip())
            stag = micr.getStage(info[2].strip())
            crl = stag.addCrystal(info[4].strip(), info[3])
            crl.reconstructOrient(info[5])
            crl.setComment(info[6])
        else:
            pass
        
    f.close()

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
    def __init__(self, name, kv=300):
        #name of the microscope
        self.name = name
        #Overpotential in volts
        self.v = kv*1e3
        self.wavelength = self.calcLambda(self.v)
        self.ewaldR = 1/self.wavelength
        
        #detectors is a dictionary name - detector object
        self.detectors = {}
        #Maybe should leave option open of having more than one stage
        self.stages = {}
    
    def setName(self, name):
        #also change this in the TEM dictionary
        microscopes[name] = microscopes.pop(self.name)
        self.name = name
        
    def addStage(self, name = "", **kwargs):
        if name == "":
            name = "stage"+str(len(self.stages))
        self.stages[name] = Stage(self, name = name, **kwargs)
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
        self.detectors[name] = Detector(self, name = name, **kwargs)
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
        kv = float(kv)
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
        elif units == "angstrom":
            return self.ewaldR/1e10
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
    def __init__(self, TEM, name = "", alpha =0, beta = 0, alphamin = -30, alphamax = 30, betamin = -20, betamax = 20, alpharev = False, betarev = False):
        #the stage needs to know which TEM it belongs to
        self.TEM = TEM
        
        self.name = name
        
        self.alpha = alpha
        self.beta = beta
        
        self.alow = alphamin
        self.atop = alphamax
        self.blow = betamin
        self.btop = betamax
        
        self.alpharev = alpharev
        self.betarev = betarev
        
        #list of crystals (this determines how the crystal is oriented with respect to the stage)
        self.crystals = {}
        
        #it should be possible to store interesting orientations as name - [alpha, beta]
        self.orientations = {}
    
    def setName(self, name):
        #also change this in the TEM dictionary
        self.getTEM().stages[name] = self.getTEM().stages.pop(self.name)
        self.name = name
    
    def getAlphaMin(self):
        return self.alow
        
    def getAlphaMax(self):
        return self.atop
        
    def getBetaMin(self):
        return self.blow
        
    def getBetaMax(self):
        return self.btop
    
    def setAlphaRange(self, *args):
        self.alow, self.atop = args
        
    def setBetaRange(self, *args):
        self.blow, self.btop = args
    
    def setRev(self, *args):
        self.alpharev, self.betarev = args
    
    def getTEM(self):
        return self.TEM
    
    def addCrystal(self, structurename, name=""):
        if name=="":
            name = "crystal"+str(len(self.crystals))
        self.crystals[name]=Crystal(name, self, structurename)
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
        
        xx = XR(alpha)
        
        yy = YR(beta)
        
        ##Make sure reversed or not is considered in definining the rotation matrixes
        #this correction made 18/03/12
        if self.alpharev:
            xx = XR(-alpha)
        
        if self.betarev:
            yy = YR(-beta)
        
        rot = np.dot(xx, yy)
        return np.dot(np.array(rot), np.array(vecs))
        
    def absToStage(self, vecs):
        """Turn stage coordinates into absolute coordinates. This is identical to doing the reverse rotation. [010]->[0, cos(alpha), -sin(alpha)] in real with only alpha rotation."""
        alpha=float(self.alpha)
        beta=float(self.beta)
        
        xx = XR(-alpha)
        yy = YR(-beta)
        
        ###This correciton made 18/03/12
        if self.alpharev:
            xx = XR(alpha)
        
        if self.betarev:
            yy = YR(beta)
        
        #The inverse of the absolute to stage
        rot = np.dot(yy, xx)
        return np.dot(np.array(rot), np.array(vecs))
    
    def inrange(self, a, b, units = "radians"):
        if units == "radians":
            a=a/(2*np.pi)*360
            b=b/(2*np.pi)*360
        if units == "degrees":
            pass #already good to compare
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
    
    def calcAlphaBeta(self, tozall, rnd = 10):
        """A general rotation matrix is not possible with just alpha and beta. It can only be arranged that a certain vector is set to the optical axis.
        The way to do this was calculated on paper and with sympy. Toz is the vector that needs to be turned on the Z-axis"""
        return toAllCols(self.calcAlphaBetaVec, tozall, rnd = rnd) 
    
    def calcAlphaBetaVec(self, toz, rnd = 10):
        """18/03/12 The new simplified script to calc the alpha and beta one needs to go to for stage vector toz to be moved to the z-axis"""
        
        #the whole of this script works with radians until the conversion at the end
        #normalize
        toz = toz/np.linalg.norm(toz)
            
        c = toz[0]#+1e-15
        d = toz[1]#+1e-15
        f = toz[2]#+1e-15

        length = np.sqrt(c**2+d**2+f**2)
        
        #print("[%s %s %s] %s" %(c, d, f, length))
        
        #we have to solve a particular equation, namely tilt(alpha, beta)*stage vec = abstoStage([001]) for alpha and beta
        def equation(x, sv): #this is a wrapper for the next equation because the third element is redundant and should not be returned
            alpha, beta = x
            y = miller(0, 0, 1)
            #y = self.absToStage(miller(0, 0, 1))
            if self.alpharev:
                alpha = -alpha
            if self.betarev:
                beta = -beta
            A = np.dot(XR(alpha), YR(beta))
            res = np.dot(A, sv) - y
            return res
            
        def tosolve(x, sv):
            res = equation(x, sv)
            return (res[0], res[1]) #apparently we can't use 3 things to evaluate, the third will be supplied as test
        
        alpha, beta = fsolve(tosolve, (0,0), args = (toz)) #anser is already in degrees
        
        #check the vector is in range
        def test(a, b):
            #is the third element also 0? - it may have mapped on -1, then the third element will be almost -2
            thirdelem = equation((alpha, beta), toz)[2]
            return (round(thirdelem,6)==0 and self.inrange(a, b, units = "degrees"))

        if test(alpha, beta):
            return np.array([np.round(alpha, rnd), np.round(beta, rnd)])
        else:
            return np.array([np.nan, np.nan])
    
    def calcAlphaBetaVecOld(self, toz, verbose=False, rnd = 10):
        """Calc the alpha and beta one needs to go to for stage vector toz to be moved to the z-axis. This old script investigates all angle possibilities between 0-360 degrees, whereas the new one simplifies it."""
        
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
    def __init__(self, TEM, name = "", diffrot=None, imrot=None, stemrot=None, diffpixcal=None, impixcal=None, stempixcal=None):
        #the detector needs to know which TEM it belongs to
        self.TEM = TEM
        
        #name
        self.name = name
        
        #the rotation of the diffraction pattern
        self.diffrot = diffrot
        self.diffrotfile = ""
        #This is the rotation of the image on any detector except HAADF
        self.imrot = imrot
        self.imrotfile = ""
        #this is the rotation of the HAADF image in STEM or the rotation of the Ronchigram on the TV or GIF CCD
        self.stemrot = stemrot
        self.stemrotfile = ""
        
        #image mode size calibration
        self.impixcal = impixcal
        #diffraction mode size calibration
        self.diffpixcal = diffpixcal
        #Stem size calibration
        self.stempixcal = stempixcal
    
    def getMags(self, mod):
        if mod == "Diffraction":
            return sorted(list(self.diffrot.keys()))
        elif mod == "STEM":
            return sorted(list(self.stemrot.keys()))
        elif mod == "Imaging":
            return sorted(list(self.imrot.keys()))
        else:
            return None
    
    def getTEM(self):
        return self.TEM
        
    def getName(self):
        return self.name
    
    def setName(self, name):
        #also change this in the TEM dictionary
        self.getTEM().detectors[name] = self.getTEM().detectors.pop(self.name)
        self.name = name
    
    def getCalibration(self, mode):
        if mode=="STEM":
            return self.stemrot
        elif mode == "Diffraction":
            return self.diffrot
        elif mode == "Imaging":
            return self.imrot
        else:
            return None
    
    def getCalFileName(self, mode):
        """Return the file names of the particular mode"""
        if mode == 0 or mode=="diffraction" or mode == "diff" or mode == "d":
            return self.diffrotfile
        if mode == 1 or mode=="imaging" or mode == "img" or mode=="i":
            return self.imrotfile
        if mode == 2 or mode=="STEM" or mode == "stem" or mode=="s":
            return self.stemrotfile
    
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
                    self.diffrotfile = filename
                if type ==1 or type == "size" or type == "s" or type =="scale":
                    self.diffpixcal = diction
            if mode == 1 or mode=="imaging" or mode == "img" or mode=="i":
                if type ==0 or type == "rotation" or type == "r" or type =="rot":
                    self.imrot = diction
                    self.imrotfile = filename
                if type ==1 or type == "size" or type == "s" or type =="scale":
                    self.impixcal = diction
            if mode == 2 or mode=="STEM" or mode == "stem" or mode=="s":
                if type ==0 or type == "rotation" or type == "r" or type =="rot":
                    self.stemrot = diction
                    self.stemrotfile = filename
                if type ==1 or type == "size" or type == "s" or type =="scale":
                    self.stempixcal = diction
            
        except:
            print("This is not a valid calibration mode or type.")
    
    def getRot(self, mode, setting):
        """This function gets the rotation angle of the detector. The rotation angle is how the absolute X and Y appear rotated on the detector. This must be tabulated and supplied."""
        
        if mode == 0 or mode=="diffraction" or mode =="Diffraction" or mode == "diff" or mode == "d":
            return self.diffrot[setting]
        elif mode == 1 or mode=="imaging" or mode == "Imaging" or mode == "img" or mode=="i":
            return self.imrot[setting]
        elif mode == 2 or mode=="STEM" or mode == "stem" or mode=="s":
            return self.stemrot[setting]
        else:
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
    def __init__(self, name, a=1, b=1, c=1, alpha=90, beta=90, gamma=90):
        
        self.name = name
        
        self.a=a
        self.b=b
        self.c=c
        self.alpha=alpha
        self.beta=beta
        self.gamma=gamma
        
        #self.typ = self.getCrystalClass()
        
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
    
    def setName(self, name):
        #also change this in the TEM dictionary
        structures[name] = structures.pop(self.name)
        self.name = name
        
    def changeCrystallography(self, name, a, b, c, alpha, beta, gamma):
        self.setName(name)
        self.__init__(name, a=a, b=b, c=c, alpha=alpha, beta=beta, gamma=gamma)
    
    def getCrystallography(self):
        return self.a, self.b, self.c, self.alpha, self.beta, self.gamma
        
    def millerToCartesian(self, vec, typ = "real"):
        """This function returns the coordinates in the cartesian coordinate system stuck to the crystal given a set of miller indices as columns or an array. Standard it will be assumed miller indices as defined by the real coordinate system, but recyp is also valid and must be supplied as typ = 'recyp' as argument"""
        if typ=="real":
            return np.dot(np.array(self.RM), np.array(vec))
        elif typ=="recyp" or typ=="recyprocal":
            return np.dot(np.array(self.RRM), np.array(vec))
        else:
            return None
    
    def cartesianToMiller(self, vec, typ = "real"):
        """This function returns the miller indices given coordinates in the cartesian system stuck to the crystal. Standard it will be assumed miller indices as defined by the real coordinate system, but recyp is also valid and must be supplied as typ = 'recyp' as argument"""
        if typ=="real":
            return np.dot(np.array(self.invRM), np.array(vec))
        elif typ=="recyp" or typ=="recyprocal":
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
    
    @staticmethod
    def getCrystalClass(a, b, c, alpha, beta, gamma):
        typ = ""
        if a==b and b==c and alpha==90 and beta==90 and gamma==90:
            typ = "cubic"
        elif a==b and alpha==90 and beta==90 and gamma==90:
            typ = "tetragonal"
        elif alpha==90 and beta==90 and gamma==90:
            typ = "orthorhombic"
        elif a==b and gamma==120 and alpha==90 and beta==90:
            typ = "hexagonal"
        elif a==b and b==c and alpha==beta and beta==gamma:
            typ = "rhombohedral"
        elif gamma == 90 and alpha==90:
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
        c1=round(c*np.cos(beta), rnd) #this can be solved from the dot product between c and a vectors
        #c2=round(c*(np.cos(alpha)-np.cos(beta)*np.cos(gamma)), rnd) #this can be solved from the dot product between c and b vectors
        c2 = (c*b*np.cos(alpha) - c1*bv[0])/(bv[1])
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
            
            vn.append(np.array([vec[val1], vec[val3], vec[val2]])) #depending on i, the values are switched. 8 extra vectors per permutations must possibly be added: the one only with switched numbers.
            vn.append(np.array([-vec[val1], vec[val3], vec[val2]])) #the one with the first element negative
            vn.append(np.array([vec[val1], -vec[val3], vec[val2]])) #the one with the second element negative
            vn.append(np.array([vec[val1], vec[val3], -vec[val2]])) #the one with the third element negative
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
            
    def getVectorAngle(self, vec1, vec2, typ="real", units="radians"):
        #! still some strange behavior when testing when vec1 or vec2 is a one dimensional array and the other is not. Otherwise it works perfectly. This has to do with the way the division happens with matrix/arrays. Fix later.
        """The angle between two vectors in real or recyprocal space,  A list of Column vectors can also be supplied, an array of 2 dimensions will be returned"""
        num= self.getVectorDot(vec1, vec2, typ=typ)
        denom = np.outer(self.getVectorLength(vec1, typ=typ), self.getVectorLength(vec2, typ=typ))
        angls=  np.arccos(np.divide(num, denom))
        if units =="radians":
            return angls
        elif units =="degrees":
            return angls/(2*np.pi)*360
        else:
            print("Those units aren't valid.")
            return None
            
            
    def getZoneRefPairs(self, l1, a1, l2, a2, maxind = 5, err = 0.5, err2 = 2, verbose = False):
        """This function helps transform lengths and directions measured in the diffraction pattern to reflections and zones.
        This is essentially a help for indexing.
        As input you provide the length and angle measured from 2 reflections in a diffraction pattern. Make sure angles are measured properly!
        Enter lengths in nm^-1 because that is usually the case in DM. Angles in degrees.
        Optional arguments are the largest valu of h, k or l to search, the +/- error on the length where it's a match (in nm^-1), and the +/- error on the angle between reflections, in degrees 
        Returned is (matching reflections to l1, the calculated lengths vector, matching reflections to l2, the calculated lengths vector, the angle between them, the calculated zone axis)"""
        
        #l1 and l2 are often times in nm^-1. Convert to angstrom^-1.
        l1 = l1/10
        l2 = l2/10
        err = err/10
        #Crucial thing is angle between the vectors to differentiate them
        a = abs(a2-a1)
        
        #first get a length match, only positive vectors need to be considered
        #construct a list of possible combined zones in which to search.
        tryvecs = getHKLlist(maxind = maxind).T
        #Make all the possible permutations of each unique type of indices
        longlist = np.array([[0,0,0]]).T
        for i in tryvecs:
            cps = np.array(self.getSym(i))
            longlist = np.c_[longlist, cps]
        
        #find the length of the vectors assuming they are recyprocal lattice points
        lv = self.getVectorLength(longlist, typ = "recyp")[0]
        
        #find where l1 and l2 match the length array. Then find the vectors matching the indices.
        indxes1 = np.where(np.logical_and(lv<l1+err, lv>l1-err))[0]
        vecs1 = longlist[:, indxes1]
        
        indxes2 = np.where(np.logical_and(lv<l2+err, lv>l2-err))[0]
        vecs2 = longlist[:, indxes2]
        
        #find angles between all the vectors that are ok in length
        angls = self.getVectorAngle(vecs1, vecs2, typ = "recyp", units = "degrees")
        #find indexes of those vectors where the angle between them are ok
        anglindx = np.where(np.logical_and(angls<a+err2, angls>a-err2))
        
        #find the vectors that match the good fit for the angle
        #rows or anglindx[0] matches vec1, columns or anglindx[1] matches vec2
        match1 = vecs1[:, anglindx[0]]
        match2 = vecs2[:, anglindx[1]]
        matchangls = angls[anglindx[0], anglindx[1]]
        matchl1 = self.getVectorLength(match1, typ = "recyp")
        matchl2 = self.getVectorLength(match2, typ = "recyp")
        
        zones = calcCross(match1, match2)
        
        if verbose:
            
            print("All testing vectors:")
            print(longlist)
            print("Lengths of the vectors:")
            print(lv)
            print("Matches to l1")
            print(indxes1)
            print(vecs1)
            print("Matches to l2")
            print(indxes2)
            print(vecs2)
            print("Angles between l1 and l2:")
            print(angls)
            
        #put into right format to output
        return match1.T.tolist(), matchl1[0].tolist(), match2.T.tolist(), matchl2[0].tolist(), matchangls.tolist(), zones.T.tolist()
        
        

        
class Crystal(object):
    
    def __init__(self, name, stage, structurename, comment = ""):
        self.structure = getStructure(structurename) #!note the globally defined getStructure as opposed to self.getstructure
        self.stage = stage #the crystal needs to know which stage it belongs to
        
        #orient matrix takes vectors from the stage coordinate system to the cartesian axes stuck to the crystal
        self.orient = np.identity(3)
        #the inverse of the orientation matrix does not need to be defined -> orthogonal so .T is the inverse already. 
        #self.orient.T takes vectors from the cartesian system stuck on the crystal to the stage coordinates
        self.name = name
        self.comment = comment
    
    def dumpOrient(self):
        """Returns a string representation of the orient matrix for saving purposes"""
        s = str(self.orient.flatten().tolist()).replace("[", "").replace("]", "") #make a list of the form "f, f, ...."
        return s
    
    def reconstructOrient(self, s):
        """Reconstructs and saves the orient matrix from a string."""
        self.orient = np.fromstring(s, dtype = np.float64, sep = ",").reshape(3, 3)
        return self.orient
    
    def getComment(self):
        return self.comment
        
    def setComment(self, cmnt):
        self.comment = cmnt
    
    def setName(self, name):
        #also change this in the TEM dictionary
        self.getStage().crystals[name] = self.getStage().crystals.pop(self.name)
        self.name = name
    
    def changeCrystallography(self, *args):
        pass #can't change the crystallography of the structure from the crystal so it is overwritten here
    
    def getCrystallography(self):
        return self.structure.getCrystallography()
    
    def getStage(self):
        return self.stage
    
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
            temp = integerRep(cryscoords)
            if max(np.absolute(temp))>10:
                return cryscoords/np.linalg.norm(cryscoords)
            else:
                return temp
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
        """Sample tilt in degrees necessary to put reflection ng into two-beam condition assuming you start from zone axis"""
        #get the length of vector g in angstrom^-1
        length = self.getVectorLength(g, typ="recyp")[0,0]
        #get the length of the radius of the ewald sphere
        #convert to angstrom^-1
        #K0 = 507.93 for kv 300
        K0 = self.stage.getTEM().getEwaldR(units = "angstrom")
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
        #get the length of vector g in angstrom^-1
        length = self.getVectorLength(g, typ="recyp")[0,0]
        #get the length of the radius of the ewald sphere
        #convert to angstrom^-1
        #K0 = 507.93
        K0 = self.stage.getTEM().getEwaldR(units = "angstrom")
        
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
        #get the length of vector g in angstrom^-1
        length = self.getVectorLength(g, typ="recyp")[0,0]
        #get the length of the radius of the ewald sphere
        #convert to angstrom^-1
        #K0 = 507.93
        K0 = self.stage.getTEM().getEwaldR(units = "angstrom")
        
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
        """This method calculates what the correct 2-beam tilt condition should be (n) if we want the ewald sphere to cut through i*g after the beam tilt so that reflection k is centered"""
        #get the length of vector g in angstrom^-1
        length = self.getVectorLength(g, typ="recyp")[0,0]
        #get the length of the radius of the ewald sphere
        #convert to angstrom^-1
        #K0 = 507.93
        K0 = self.stage.getTEM().getEwaldR(units = "angstrom")
        
        #let us assume the sample is tilted by phy and there is a beam tilt theta such that kg is on the optical axis and ig is excited. i is n in many books.
        #Here we calculate n0 -> if there is no beam tilt, i.e. before putting k on the optical axis, what is the excited beam to tilt to.
        #the angle between the beam and the row of reflections can be easily calculated 
        ang = np.arccos(i*length/2/K0)
        #The sample has a certain tilt phi. When kg is on the optical axis, then k*l*cos(phi) = K0*cos(ang+phi). To solve sample tilt, we must hence solve the following equation.
        def eq(phi):
            return np.arccos(k*length*np.cos(phi)/K0)-ang-phi
        
        phi = fsolve(eq, 0)
        #from the known phi, the initial 2-beam can be determined at 0 beam tilt. This is finding the distance between the point 0,0 and the intersection between the ewald sphere over 0 and a line going at an angle phi with the x-axis.
        lent = 2*K0*np.sin(phi)
        
        return (lent/length)[0]
    
    def calcSgSimple(self, g, k, i, units = "nm"):
        """Calculate Sg on the same basis as calcWB2beam but this only works if the imaging reflection is inside the ewald wphere and n>k"""
        length = self.getVectorLength(g, typ="recyp")[0,0]
        K0 = self.stage.getTEM().getEwaldR(units = "angstrom")
        ang = np.arccos(i*length/2/K0)
        def eq(phi):
            return np.arccos(k*length*np.cos(phi)/K0)-ang-phi
        phi = fsolve(eq, 0)
        sg = K0-(K0*np.sin(ang+phi) - k*length*np.sin(phi))
        #Testing the simpler formula in the book by kirk
        #print((i-1)*length**2/(2*K0)*10) tried on 19/03/18 and it matches well
        if units == "angstrom":
            return sg[0] #in angstrom-1
        if units == "nm":
            return sg[0]*10
    
    def calcSg(self, g, n, k, n2, verbose=False):
        #get the length of vector g in nm-1
        length = self.getVectorLength(g, typ="recyp")[0,0]
        #get the length of the radius of the ewald sphere
        #convert to nm-1
        #K0 = 507.93
        K0 = self.stage.getTEM().getEwaldR(units = "angstrom")
        
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
        return self.stage.calcAlphaBeta(stg, **kwargs)
    
    def isReachable(self, vec, typ = "real", verbose = False, **kwargs):
        """This function returns true or false whether the particular vector can be reached by the stage tilt"""
        res = self.getAlphaBetaMiller(vec, typ=typ, verbose = verbose, **kwargs)
        ###test now which columns are nan
        #sum yields the sum over the rows or nan if it contains nan. Then isnan yields whether it is a nan or not.
        #if it's a nan it should say False, it is not reachable
        return np.invert(np.isnan(np.sum(res, 0)))
    
    def nearestZone(self, curzone = None, numb = 1):
        """Return the numb nearest zones (real indices) to the provided zone. If there is no zone provided the current zone is calculated."""
        ###Numb is as yet not implemented
        if curzone is None:
            curzone = self.getZone() #the zone could be some irrational combination which makes an integer representation impossible

        uvwn = normMatrix(curzone) #in case curzone is not provided as a unit vector
        
        mx = np.array([])
        nm = 0
        #list of types of zones to test
        
        lst = [np.array([1,0,0]), np.array([1,1,0]), np.array([1,1,1]), np.array([1,2,0]), np.array([1, 1, 2]), np.array([1, 2, 2]), np.array([1,3,0]), np.array([1,3,1])]
        #lst = [np.array([1,0,0])]
        #find the ZA that is closest from a list
        for i in lst:
            mat = self.getSym(i)
            matn = normMatrix(mat)
            #test which ones are reachable and only take those
            isr = self.isReachable(matn, verbose = False)
            matn = matn[:, isr]
            #test that matn isn't empty; if all are nan then matn should be empty
            if matn.size>0:
                #make the dot product with hkl
                dp = np.dot(matn.T, uvwn)
                #print(dp)
                #find where the dot product is maximum - this vector is closest to the current zone
                indx =  np.argmax(dp)
                #print(indx)
                if dp[indx,0]>nm:
                    #check also that the found zone is not identical to curzone. Since they are both unit vectors, dp = 1
                    if round(dp[indx,0], 6)!=1:
                        #if it's larger than the already found dot product, then replace the storage
                        nm = dp[indx,0]
                        mx = integerRep(matn[:, indx])
                        #print(nm)
                        #print(mx)
        return mx
    
    def getAlphaBetaWBeam(self, g, n=3, k=1, rnd = 5, outzone = 8, za = "current", verbose = False):
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
            
            #calculate the n necessary to go to 2beam condition i = n
            n0 = self.calcWB2beam(g, k, n)
            #calculate the angle necessary to tilt to n0g
            
            ang = self.calcSampleTilt(g, n = n0)
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
    
    def getZoneAtAlphaBeta(self, alpha, beta, typ = "real", rnd = 5, verbose = False, integer = True):
        """It is easiest to change and change back alpha and beta of the stage and just call getZone"""
        #keep the current stage angles
        storalpha = self.stage.getAlpha()
        storbeta = self.stage.getBeta()
        #set the angles temporarily
        self.stage.setAlpha(alpha)
        self.stage.setBeta(beta)
        #get the zone
        res = self.getZone(typ=typ, verbose = verbose, integer = integer)
        #set the angles back
        self.stage.setAlpha(storalpha)
        self.stage.setBeta(storbeta)
        #return the results
        return res
        
    
    def getAlphaBetaTiltRound(self, g, outzone, za = "current", rnd = 5, verbose = False):
        """Calculate alpha and beta for getting to a zone outzone away from the zone axis but not a two-beam condition"""
        if str(za)=="current":
            za = self.getZone(typ="real")
        if abs(np.round(np.dot(za, g), rnd))==0:
            #get the absolute coordinates of the zone axis
            zax = self.millerToAbs(za)
            #the z-axis must be rotated outzone. The rotation axis is the reflection.
            rotax = self.millerToAbs(g, typ="recyp")
            zdash = np.dot(axisAngle(rotax, outzone), zax) #the new z-axis that we should see when we excite te systematic row
            #this axis must be converted to miller
            zabmil = self.absToMiller(zdash, typ = "real")
            #This axis must then be brought to the z-axis
            if verbose:
                print("Absolute coordinates of the zone axis:")
                print(zax)
                print("Z axis rotated out of zone:")
                print(zdash)
                print("Vector to Z axis in miller indices (real)")
                print(zabmil)
            #return the alpha and beta
            return self.getAlphaBetaMiller(zabmil, typ="real", verbose=verbose)
        else:
            print("This reflection and zone axis are not perpendicular")
            return np.array([0,0])
    
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
            
    def getVectorAngle(self, vec1, vec2, typ="real", units = "radians"):
        return self.structure.getVectorAngle(vec1, vec2, typ=typ, units = units)
        
    def plotCrystalonDetector(self, detector, mode, setting, vecincl = [], typ="real", plotaxes = True, verbose=False, plotbond = True, hkm = 3):
        #d001 = self.getSym(miller(1, 0, 0), typ=typ)
        #d111 = self.getSym(miller(1, 1, 1), typ=typ)
        #d011 = self.getSym(miller(1, 1, 0), typ=typ)
        #toplotd001 = self.millerToDetector(d001, detector, mode, setting, typ=typ, verbose=verbose)
        #toplotd111 = self.millerToDetector(d111, detector, mode, setting, typ=typ, verbose=verbose)
        #toplotd011 = self.millerToDetector(d011, detector, mode, setting, typ=typ, verbose=verbose)
        
        vecincl = getHKLlistNoMultipes(maxind = hkm).T.tolist()
        
        millers = []
        tp = []
        for i in vecincl:
            i = np.array(i)
            mil = self.getSym(i, typ = typ)
            millers.append(mil)
            toplot = self.millerToDetector(mil, detector, mode, setting, typ=typ, verbose=verbose)
            tp.append(toplot)
        
        #print(millers)
        #print(tp)
        
        #change the type of brackets based on the typ passed
        brc = []
        if typ == "real":
            brc = ["[", "]"]
        elif typ== "recyprocal":
            brc = ["(", ")"]
            
        #make the actual plot
        fig, ax = plt.subplots(figsize = (7, 6))
        stereographicCanvas(ax)
        
        #plotStereographic(d001, d001, ax, verbose = verbose, s=100, c="b")
        #plotStereographic(d111, d111, ax, verbose = verbose, s=100, c="b")
        #plotStereographic(d011, d011, ax, verbose = verbose, s=100, c="g")
        
        if plotbond: #plot the boundaries also
            detail = 10
            
            ##first construct arrays of angles
            amx=self.stage.getAlphaMax()
            amn = self.stage.getAlphaMin()
            bmx = self.stage.getBetaMax()
            bmn = self.stage.getBetaMin()
            
            aar = np.linspace(amn, amx, detail).tolist()
            bar = np.linspace(bmn, bmx, detail).tolist()
            
            #Turn the angles into millers and construct it in the right way as to have 4 lines
            l1 = [] #the amin line
            l2 = [] #the amax line
            for i in aar:
                l1.append(self.getZoneAtAlphaBeta(i, bmn, typ = typ, integer = False))
                l2.append(self.getZoneAtAlphaBeta(i, bmx, typ = typ, integer = False))
                
            l3 = [] #the amin line
            l4 = [] #the amax line
            for i in bar:
                l3.append(self.getZoneAtAlphaBeta(amn, i, typ = typ, integer = False))
                l4.append(self.getZoneAtAlphaBeta(amx, i, typ = typ, integer = False))
            
            l1 = np.array(l1).T
            l2 = np.array(l2).T
            l3 = np.array(l3).T
            l4 = np.array(l4).T
            
            #Turn the millers into detector vectors
            tpl1 = self.millerToDetector(l1, detector, mode, setting, typ=typ, verbose=verbose)
            tpl2 = self.millerToDetector(l2, detector, mode, setting, typ=typ, verbose=verbose)
            tpl3 = self.millerToDetector(l3, detector, mode, setting, typ=typ, verbose=verbose)
            tpl4 = self.millerToDetector(l4, detector, mode, setting, typ=typ, verbose=verbose)
            
            #plotstereographic but editted
            plotLineStereographic(tpl1, ax)
            plotLineStereographic(tpl2, ax)
            plotLineStereographic(tpl3, ax)
            plotLineStereographic(tpl4, ax)
        
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
        
        #plotStereographic(toplotd001, d001, ax, verbose = verbose, brc = brc, s=100, c="r")
        #plotStereographic(toplotd111, d111, ax, verbose = verbose, brc = brc, s=100, c="b")
        #plotStereographic(toplotd011, d011, ax, verbose = verbose, brc = brc, s=100, c="g")
        
        for i, j in enumerate(millers):
            plotStereographic(tp[i], j, ax, verbose = verbose, brc = brc)
        
        plt.tight_layout()
        
        return fig, ax
        #plt.show()

def getSym(vec, unc = 1e-9):
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
            
            vn.append(np.array([vec[val1], vec[val3], vec[val2]])) #depending on i, the values are switched. 8 extra vectors per permutations must possibly be added: the one only with switched numbers.
            vn.append(np.array([-vec[val1], vec[val3], vec[val2]])) #the one with the first element negative
            vn.append(np.array([vec[val1], -vec[val3], vec[val2]])) #the one with the second element negative
            vn.append(np.array([vec[val1], vec[val3], -vec[val2]])) #the one with the third element negative
            for j in vn: #all are checked to see whether they already exist in the matrix
                if not isExist(tmpmat, j): #if they don't they get added
                    tmpmat = np.c_[tmpmat, j]
                if not isExist(tmpmat, -j):
                    tmpmat = np.c_[tmpmat, -j]
        
        return tmpmat
        
def normMatrix(vec):
    """normalize a vector or columns of a matrix. Returns the matrix with normalized columns. If an array was passed a matrix with one column is returned."""
    #make vec a matrix in case it is one dimensional
    if vec.ndim==1:
        vec = np.matrix(vec).T
    #normalize the column vectors
    invl = 1/np.matrix(np.linalg.norm(vec, axis=0))
    lrp = np.repeat(invl, vec.shape[0], axis=0)
    return np.array(np.multiply(lrp, vec))

def calcCross(vecs1, vecs2, intrep = True):
    """outputs UVW between the different sets of hkl. They need to be of the same size and in format 3xN. Optionally, a non-integer representation can be returned."""
    top = np.multiply(vecs1[1, :], vecs2[2,:])-np.multiply(vecs1[2, :], vecs2[1,:])
    middle = np.multiply(vecs1[2, :], vecs2[0,:])-np.multiply(vecs1[0, :], vecs2[2,:])
    bottom = np.multiply(vecs1[0, :], vecs2[1,:])-np.multiply(vecs1[1, :], vecs2[0,:])
    res = np.array([top, middle, bottom])
    if intrep:
        return getIntRepAny(res)   
    else:
        return res

def getHKLlist(maxind = 5):
    """Get a list of all unique combinations of 3x[0-maxind] excluding [000]. This is returned as a 3xN array"""
    lst = []
    for i in np.arange(maxind+1):
        for j in np.arange(i, maxind+1):
            for k in np.arange(j, maxind+1):
                if i==0 and j==0 and k==0:
                    pass
                else:
                    lst.append([i, j, k])
    return np.array(lst).T

def getHKLlistNoMultipes(maxind = 5):
    lst = []
    for i in np.arange(maxind+1):
        for j in np.arange(i, maxind+1):
            for k in np.arange(j, maxind+1):
                if i==0 and j==0 and k==0:
                    pass
                else:
                    lst.append([i, j, k])
    #go through the array once more and remove all those where there is multiples
    nl = getIntRepAny(np.array(lst).T)
    unique_rows = np.unique(nl, axis=1)
    return unique_rows
    
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

def plotStereographic(vecs, labs, ax, verbose=False, brc = ["[", "]"],  **kwargs):
    """The annotated poles of certain vectors are plotted in the stereogram. Vecs should be the vectors after any transformation/rotation"""
    projpospos, projpos, projnegpos, projneg  = calcStereo(vecs)
    if verbose:
        print("The coordinates that will be plot:")
        print(projpos)
    labs = labs[:, projpospos]
    x = projpos[0, :].T
    y = projpos[1, :].T
    #make the size of the pieces inversely proportional to their multiplicity
    mult = vecs.shape[1]
    for i, xy in enumerate(zip(x, y)):
        ax.annotate("%s" %(str(labs[:,i].T).replace("[[", brc[0]).replace("]]", brc[1])).replace(".", "").replace("-0", "0"), xy=xy, textcoords = "data")
    ax.scatter([x],[y], s = 100*12/mult, **kwargs)

def plotLineStereographic(vecs, ax, **kwargs):
    """A line is plotted in the stereographic projection"""
    projpospos, projpos, projnegpos, projneg  = calcStereo(vecs)
    x = np.array(projpos[0, :])[0]
    y = np.array(projpos[1, :])[0]
    ax.plot(x,y, color = "black", **kwargs)           
        