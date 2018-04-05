#! /usr/bin/env python
# -*- coding: utf-8 -*-
#
import sys
from PyQt5.QtGui import *
from PyQt5.QtWidgets import *
from PyQt5.QtCore import *
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

import matplotlib.pyplot as plt

import numpy as np
import temCalcs as tc

sp = QSizePolicy()
sp.setRetainSizeWhenHidden(True)

class App(QMainWindow):
    
    def __init__(self):
        super().__init__()
        self.title = "TEMTilt v0.01 BETA"
        #self.left = 50
        #self.top = 50
        #self.width = 640
        #self.height = 400
        
        self.setWindowTitle(self.title)
        #self.setGeometry(self.left, self.top, self.width, self.height)
        self.setWindowIcon(QIcon('./Images/logo.png'))
        
        layout = QVBoxLayout()
        layout.setSpacing(0)
        marg = 0
        layout.setContentsMargins(marg, marg, marg, marg)
        ###########The loading and saving menu#########
        self.loadsavemenu = loadMenu(self)
        layout.addWidget(self.loadsavemenu)
        
        ###########The microscopes menu#################
        self.microscopemenu = microscopeMenu(self)
        layout.addWidget(self.microscopemenu)
        ###########The strucutres menu##################
        self.structuremenu = structuresMenu(self)
        layout.addWidget(self.structuremenu)
        ###The crystals menu#########
        self.crystalmenu = crystalMenu(self)
        layout.addWidget(self.crystalmenu)
        ####The calculation menu#########
        self.calcmenu = calculationMenu(self)
        layout.addWidget(self.calcmenu)
        
        ####The plotting menu#########
        self.plotmenu = plottingMenu(self)
        layout.addWidget(self.plotmenu)
        
        
        layout.setSizeConstraint(layout.SetFixedSize)
        
        wd = QWidget()
        wd.setLayout(layout)
        
        self.setCentralWidget(wd)
        
        #autoload
        self.loadsavemenu.autoload()
        
        #show the window
        self.show()
    
    def updateGUI(self):
        self.microscopemenu.updateAll()
        self.structuremenu.updateAll()
        self.crystalmenu.updateAll()
    
    def updateCrystals(self):
        self.crystalmenu.updateAll()
    
    def getCurrentTEM(self):
        return self.microscopemenu.currentTEM()
    
    def getCurrentDetector(self):
        return self.microscopemenu.currentDetector()
        
    def getCurrentStage(self):
        return self.microscopemenu.currentStage()
    
    def getAlpha(self):
        return self.microscopemenu.getAlpha()
    
    def getBeta(self):
        return self.microscopemenu.getBeta()
    
    def getMode(self):
        return self.microscopemenu.getMode()
    
    def getSetting(self):
        return self.microscopemenu.getSetting()
    
    def getTheta(self):
        return self.microscopemenu.getTheta()
    
    def getCurrentStructure(self):
        return self.structuremenu.currentStruc()
        
    def getCurrentCrystal(self):
        return self.crystalmenu.currentCrystal()
        
        
class loadMenu(QWidget):
    
    def __init__(self, caller):
        super().__init__()
        self.caller = caller
        
        self.mainbox = QGroupBox("Load/Save")
        mainlayout = QGridLayout()
        mainlayout.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        
        size = [24, 24]
        
        self.loadbutton = QPushButton("")
        self.loadbutton.setToolTip("Load configuration")
        self.loadbutton.setIcon(QIcon(".\Images\load.png"))
        self.loadbutton.clicked.connect(self.load)
        self.loadbutton.resize(size[0]+6, size[1]+6)
        self.loadbutton.setIconSize(QSize(size[0],size[1]))
        self.loadbutton.setSizePolicy(sp)
        mainlayout.addWidget(self.loadbutton, 1, 0)
        
        self.savebutton = QPushButton("")
        self.savebutton.setToolTip("Save configuration")
        self.savebutton.setIcon(QIcon(".\Images\save.png"))
        self.savebutton.clicked.connect(self.save)
        self.savebutton.resize(size[0]+6, size[1]+6)
        self.savebutton.setIconSize(QSize(size[0],size[1]))
        self.savebutton.setSizePolicy(sp)
        mainlayout.addWidget(self.savebutton, 1, 1)
        
        self.trashbutton = QPushButton("")
        self.trashbutton.setToolTip("Clear the session")
        self.trashbutton.setIcon(QIcon(".\Images\poubelle.png"))
        self.trashbutton.clicked.connect(self.trash)
        self.trashbutton.resize(size[0]+6, size[1]+6)
        self.trashbutton.setIconSize(QSize(size[0],size[1]))
        self.trashbutton.setSizePolicy(sp)
        mainlayout.addWidget(self.trashbutton, 1, 2)
        
        
        self.mainbox.setLayout(mainlayout)
        
        Layout = QVBoxLayout()
        Layout.addWidget(self.mainbox)
        #mainLayout.addWidget(self.buttonBox)
        self.setLayout(Layout)
        
        #show the window
        self.show()
    
    def trash(self):
        if self.checkEvent("Are you sure you want to clear everything in the session?"):
            tc.clearSession()
            self.caller.updateGUI()
    
    def checkEvent(self, msg = "Are you sure?"):
        reply = QMessageBox.question(self, ' ', msg, QMessageBox.Yes, QMessageBox.No)
        if reply == QMessageBox.Yes:
            return True
        else:
            return False
    
    def load(self):
        filename, okpres = QFileDialog.getOpenFileName(caption = "Load previous session", filter = "Text files (*txt)")
        if okpres:
            tc.clearSession()
            tc.loadSession(filename)
            self.caller.updateGUI()
    
    def autoload(self):
        #autoload the autoload.txt file if it exists
        try:
            tc.loadSession("autoload.txt")
            self.caller.updateGUI()
        except:
            pass
    
    def testload(self):
        tc.clearSession()
        tc.addStructure(name="Steel", a=3.66, b=3.66, c=3.66, alpha=90, beta=90, gamma=90)
        mic = tc.addMicroscope(name = "Jeol3000", kv=300)
        stag = mic.addStage(name = "Doubletilt", alpha = 20, beta = 10, alphamin = -29, alphamax = 31, betamin = -21, betamax = 22, alpharev = True, betarev = False)
        ccd = mic.addDetector(name="MSC")
        ccd.setCalibration("rotationCalibrationMSC-IMG.txt", mode = "imaging", type = "r")
        ccd.setCalibration("rotationCalibrationMSC-DIFF.txt", mode = "diffraction", type = "r")
        crl = stag.addCrystal("Steel", name = "grain1")
        crl.calcOrient(tc.miller(1, 1, 0), tc.miller(1, -1, 1), -30, "MSC", "d", 30)
        self.caller.updateGUI()
        
    def save(self):
        filname, okpress = QFileDialog.getSaveFileName(caption = "Save session", filter = "Text files (*txt)")
        #check that filname ends in txt
        if okpress:
            tc.saveSession(filname)

class microscopeMenu(QWidget):

    def __init__(self, caller):
        ###########The microscopes menu#################
        super().__init__()
        ##make a grouping box
        self.formInstrumentGroupBox = QGroupBox("Instrument settings")
        
        #caller is the app that calls it
        self.caller = caller
        
        #inside the grouping box have a grid layout
        temlayout = QGridLayout()
        temlayout.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        temlayout.setSpacing(10)
        
        #labels - never hide labels
        self.mcrlbl = QLabel("Microscope")
        self.mcrlbl.setToolTip("Transmission electron microscopes.")
        temlayout.addWidget(self.mcrlbl, 1, 0)
        self.stlbl = QLabel("Stage")
        self.stlbl.setToolTip("Double tilt stages on the current microscope.\nA microscope is prerequisite to adding and editing stages.")
        temlayout.addWidget(self.stlbl, 2, 0)
        self.dtlbl = QLabel("Detector")
        self.dtlbl.setToolTip("Detectors on the current microscope.\nA microscope is prerequisite to adding and editing detectors.")
        temlayout.addWidget(self.dtlbl, 3, 0)
        
        self.vtlbl = QLabel("Voltage")
        self.vtlbl.setToolTip("High tension of the current microscope.\nA microscope is prerequisite to changing these settings.")
        temlayout.addWidget(self.vtlbl, 1, 5)
        
        self.allbl = QLabel(u"\u03b1")
        self.allbl.setToolTip("\u03b1 tilt angle setting of the current stage.\nA stage is prerequisite to changing these settings.")
        temlayout.addWidget(self.allbl, 2, 5) #alpha
        self.btlbl = QLabel(u"\u03b2")
        self.btlbl.setToolTip("\u03b2 tilt angle setting of the current stage.\nA stage is prerequisite to changing these settings.")
        temlayout.addWidget(self.btlbl, 2, 7) #beta
        
        self.mdlbl = QLabel("Mode")
        self.mdlbl.setToolTip("Current imaging mode of the microscope.\nA detector is a prerequisite to changing these settings.")
        temlayout.addWidget(self.mdlbl, 3, 5)
        
        self.mglbl = QLabel("Mag/CL") #this one must be updated when mode is changed
        self.mglbl.setToolTip("Current magnification or camera length setting.\nA detector calibration is a prerequisite to changing these settings.")
        temlayout.addWidget(self.mglbl, 4, 5)
        
        self.ttlbl = QLabel(u"\u03b8")
        self.ttlbl.setToolTip("The angle the \u03b1-axis makes with the detector x-axis.\nA detector calibration is a prerequisite to changing these settings.")
        temlayout.addWidget(self.ttlbl, 4, 7)
        
        #add the microscopes dropdown list
        xbut = 50
        ybut = 50
        self.temlist = self.combobox(list(tc.microscopes.keys()), pos = [xbut, ybut], action = self.updateAll)
        temlayout.addWidget(self.temlist, 1, 1)
        
        #add microscope add, edit and delete buttons
        self.adtembut = self.button(logo=".\Images\plus.png", hint = "Add new microscope", position = [xbut+100, ybut], action = self.createMicroscope)
        self.edittembut = self.button(logo = ".\Images\edit.png", hint = "Edit microscope", position = [xbut+130, ybut], action = self.editMicroscope)
        self.deltembut = self.button(logo=".\Images\delete-icon.png", hint = "Delete microscope", position = [xbut+160, ybut], action = self.deleteMicroscope)
        
        temlayout.addWidget(self.adtembut, 1, 2)
        temlayout.addWidget(self.edittembut, 1, 3)
        temlayout.addWidget(self.deltembut, 1, 4)
        
        #add the kv box. Update done automatically.
        self.kvbox = QSpinBox(self)
        self.kvbox.setRange(1, 10000)
        self.kvbox.setSingleStep(10)
        #self.kvbox.move(xbut + 260, ybut)
        self.kvbox.setSuffix(" kV")
        self.kvbox.valueChanged[int].connect(self.updateVoltage)
        self.kvbox.setSizePolicy(sp)
        temlayout.addWidget(self.kvbox, 1, 6)
        
        ######The stage menu#################
        #add the stage dropdown list
        yoffset = 50
        self.stagelist = self.combobox([], pos = [xbut, ybut+yoffset], action = self.updateStages)
        temlayout.addWidget(self.stagelist, 2, 1)
        
        dec = 2
        #add boxes for alpha and beta
        self.alphabox = QDoubleSpinBox(self)
        self.alphabox.setDecimals(dec)
        self.alphabox.setSingleStep(10**(-dec))
        #maxima should only be updated once a stage actually exists
        self.alphabox.setMaximum(90)
        self.alphabox.setMinimum(-90)
        #self.alphabox.setValue(alpha)
        self.alphabox.setSuffix(u" \u00b0")
        self.alphabox.valueChanged.connect(self.updateAlpha)
        self.alphabox.setSizePolicy(sp)
        temlayout.addWidget(self.alphabox, 2, 6)
        
        self.betabox = QDoubleSpinBox(self)
        self.betabox.setDecimals(dec)
        self.betabox.setSingleStep(10**(-dec))
        self.betabox.setMaximum(90)
        self.betabox.setMinimum(-90)
        self.betabox.setSuffix(u" \u00b0")
        self.betabox.valueChanged.connect(self.updateBeta)
        self.betabox.setSizePolicy(sp)
        temlayout.addWidget(self.betabox, 2, 8)
        
        #add the stage add and delete buttongs
        self.adstagbut = self.button(logo=".\Images\plus.png", hint = "Add new stage", position = [xbut+100, ybut+yoffset], action = self.createStage)
        self.editstagbut = self.button(logo = ".\Images\edit.png", hint = "Edit stage", position = [xbut+130, ybut+yoffset], action = self.editStage)
        self.delstagbut = self.button(logo=".\Images\delete-icon.png", hint = "Delete stage", position = [xbut+160, ybut+yoffset], action = self.deleteStage)
        temlayout.addWidget(self.adstagbut, 2, 2)
        temlayout.addWidget(self.editstagbut, 2, 3)
        temlayout.addWidget(self.delstagbut, 2, 4)
        
        
        ######The detector menu#########
        self.detectorlist = self.combobox([], action = self.updateDetectors)
        temlayout.addWidget(self.detectorlist, 3, 1)
        
        self.addetbut = self.button(logo=".\Images\plus.png", hint = "Add new detector", action = self.createDetector)
        self.editdetbut = self.button(logo = ".\Images\edit.png", hint = "Edit detector",  action = self.editDetector)
        self.deldetbut = self.button(logo=".\Images\delete-icon.png", hint = "Delete detector",  action = self.deleteDetector)
        temlayout.addWidget(self.addetbut, 3, 2)
        temlayout.addWidget(self.editdetbut, 3, 3)
        temlayout.addWidget(self.deldetbut, 3, 4)
        
        #the mode dropdown button
        self.imgmode = self.combobox(["Diffraction", "Imaging", "STEM"], action = self.updateMagoptions)
        temlayout.addWidget(self.imgmode, 3, 6)
        
        #the load calibration button
        self.loadcalib = self.button(logo = ".\Images\impo.png", hint = "Load calibration file", action = self.addCalibration)
        temlayout.addWidget(self.loadcalib, 3, 8)
        
        #the mag/cl dropdown
        self.magcl = self.combobox([], action = self.updateTheta)
        temlayout.addWidget(self.magcl, 4, 6)
        
        #the theta view box. The angle must be added later when you update the field
        self.thetaview = QLineEdit(self)
        self.thetaview.setReadOnly(True)
        self.thetaview.setSizePolicy(sp)
        temlayout.addWidget(self.thetaview, 4, 8)
        
        ####set the layout of the microscope control box
        self.formInstrumentGroupBox.setLayout(temlayout)
        
        self.updatebuttons()
        
        mainLayout = QVBoxLayout()
        mainLayout.addWidget(self.formInstrumentGroupBox)
        #mainLayout.addWidget(self.buttonBox)
        self.setLayout(mainLayout)
        
        #show the window
        self.show()
    
    def getAlpha(self):
        return self.alphabox.value()
        
    def getBeta(self):
        return self.betabox.value()
        
    def getMode(self):
        return self.imgmode.currentText()
        
    def getSetting(self):
        return self.magcl.currentText()
        
    def getTheta(self):
        return self.thetaview.currentText()
    
    def button(self, text ="" , logo = "", hint = "", position = [0,0], size = [24, 24], action = None):
        if action is None:
            action = self.nothingHappens
        button = QPushButton(text, self)
        button.setToolTip(hint)
        #button.move(position[0], position[1])
        button.setIcon(QIcon(logo))
        button.clicked.connect(action)
        button.resize(size[0]+6, size[1]+6)
        button.setIconSize(QSize(size[0],size[1]))
        button.setSizePolicy(sp)
        return button
    
    def combobox(self, lst, pos = [0,0], action = None):
        """combobox creates a drop down list from a list of strings lst at position [x, y]. When one option is active an action is performed."""
        if action is None:
            action = self.nothingHappens
        combobx = QComboBox(self)
        combobx.addItems(lst)
        #combobx.move(pos[0], pos[1])
        combobx.activated[str].connect(action)
        combobx.setSizePolicy(sp)
        return combobx
    
    def nothingHappens(self):
        pass
        
    def createMicroscope(self):
        res, okPressed = microscopeDialog.getInfo()
        #print(res)
        #print(ok)
        #text, okPressed = QInputDialog.getText(self, "New microscope","Microscope name:", QLineEdit.Normal, "")
        #d, okPressed = QInputDialog.getDouble(self, "New microscope","Operating voltage (kV):", 300)
        if okPressed:
            #need to check if that name doesn't already exist or the TEM will overwrite items
            if res[0] not in tc.microscopes:
                newtem = tc.addMicroscope(name = res[0], kv = res[1])
                self.updateAll(focus = newtem.name)
            else:
                QMessageBox.warning(self, " ", "A microscope by this name already exists. Please choose another name.")
    
    def editMicroscope(self):
        res, okPressed = microscopeDialog.getInfo(windowtitle = "Edit microscope", name = self.currentTEM().name, voltage = self.currentTEM().getKv())
        if okPressed:
            #if the name is nothing or it already exists, keep old name
            if res[0]=="" or res[0] in tc.microscopes:
                res[0] = self.currentTEM().name
            self.currentTEM().setKv(res[1])
            self.currentTEM().setName(res[0]) #must update name last if you don't want to store current TEM in a temp variable
            self.updateAll(focus = res[0])
    
    def deleteMicroscope(self):
        nt = self.currentTEM().name
        check = self.checkEvent(msg = "Are you sure you want to delete %s?" %(nt))
        if check:
            tc.removeMicroscope(nt)
            self.updateAll()
    
    def createStage(self):
        res, okPressed = stageDialog.getInfo()
        if okPressed:
            #need to check if that name doesn't already exist or it will be overwritten
            if res[0] not in self.currentTEM().stages:
                newstage = self.currentTEM().addStage(name = res[0], alpha = res[3], beta = res[6], alphamin = res[1], alphamax = res[2], betamin = res[4], betamax = res[5], alpharev = res[7], betarev = res[8])
                self.updateStages(focus = newstage.name)
            else:
                QMessageBox.warning(self, " ", "A stage by this name already exists. Please choose another name.")
        
    def editStage(self):
        res, okPressed = stageDialog.getInfo(windowtitle = "Edit stage", name = self.currentStage().name, alpha = self.currentStage().getAlpha(), beta = self.currentStage().getBeta(), alphamin = self.currentStage().getAlphaMin(), alphamax = self.currentStage().getAlphaMax(), betamin = self.currentStage().getBetaMin(), betamax = self.currentStage().getBetaMax(), alpharev = self.currentStage().alpharev, betarev = self.currentStage().betarev )
        if okPressed:
            #if the name is nothing or already exists in the list keep the name
            if res[0]=="" or res[0] in self.currentTEM().stages:
                res[0] = self.currentStage().name
            self.currentStage().setAlphaRange(res[1], res[2])
            self.currentStage().setBetaRange(res[4], res[5])
            self.currentStage().setAlpha(res[3])
            self.currentStage().setBeta(res[6])
            self.currentStage().setRev(res[7], res[8])
            self.currentStage().setName(res[0]) #must change name last
            self.updateStages(focus = res[0])
    
    def deleteStage(self):
        nt = self.currentStage().name
        check = self.checkEvent(msg = "Are you sure you want to delete %s?" %(nt))
        if check:
            self.currentTEM().removeStage(nt)
            self.updateStages()
    
    def createDetector(self):
        #res, okPressed = detectorDialog.getInfo()
        res, okPressed = detectorDialog.getInfo()
        if okPressed:
            #need to check if that name doesn't already exist or it will be overwritten
            if res[0] not in self.currentTEM().detectors:
                newdet = self.currentTEM().addDetector(name = res[0])
                self.updateDetectors(focus = newdet.name)
            else:
                QMessageBox.warning(self, " ", "A detector by this name already exists. Please choose another name.")
        
    def editDetector(self):
        res, okPressed = detectorDialog.getInfo(windowtitle = "Edit stage", name = self.currentDetector().name)
        if okPressed:
            #if the name is nothing or already exists in the list keep the name
            if res[0]=="" or res[0] in self.currentTEM().detectors:
                res[0] = self.currentDetector().name
            self.currentDetector().setName(res[0]) #must change name last
            self.updateDetectors(focus = res[0])
    
    def deleteDetector(self):
        nt = self.currentDetector().name
        check = self.checkEvent(msg = "Are you sure you want to delete %s?" %(nt))
        if check:
            self.currentTEM().removeDetector(nt)
            self.updateDetectors()
    
    def updateVoltage(self, value):
        try:
            self.currentTEM().setKv(value)
        except:
            pass
    
    def updateAlpha(self):
        #if the alpha and beta box are visible it means a stage must be present
        currentstage = self.currentStage()
        currentstage.setAlpha(self.alphabox.value())
        
    def updateBeta(self):
        currentstage = self.currentStage()
        currentstage.setBeta(self.betabox.value())
    
    def addCalibration(self):
        filename, okpress = QFileDialog.getOpenFileName(caption = "Open calibration file", filter = "Text files (*txt)")
        if filename != "" and okpress:
            currentdetector = self.currentDetector()
            mod = self.imgmode.currentText()
            if mod == "Diffraction":
                mod = "diffraction"
            if mod == "STEM":
                mod = "stem"
            if mod == "Imaging":
                mod = "imaging"
            try:
                currentdetector.setCalibration(filename, mode = mod)
                self.updateMagoptions()
                self.updatebuttons()
            except:
                QMessageBox.warning(self, " ", "That does not appear to be a valid calibration file.")
        
    def updateMagoptions(self):
        mod = self.imgmode.currentText()
        self.magcl.clear()
        try:
            opts = self.currentDetector().getMags(mod)
            opts = list(map(int, opts))
            opts = list(map(str, opts))
            #update the Mag dropdown list
            #print(opts)
            self.magcl.addItems(opts)
        except:
            pass
        self.updateTheta()
        #we should also update the crystals
        self.caller.updateCrystals()
        
    def updateTheta(self):
        try:
            mod = self.imgmode.currentText()
            if mod == "Diffraction":
                mod = "diffraction"
            if mod == "STEM":
                mod = "stem"
            if mod == "Imaging":
                mod = "imaging"
            set = self.magcl.currentText()
            val = self.currentDetector().getRot(mod, float(set))
            self.thetaview.setText(str(val) + u" \u00b0")
        except: #if no mag data can be found
            self.thetaview.setText("")
    
    def currentTEM(self):
        if tc.microscopes:
            return tc.getMicroscope(self.temlist.currentText())
        else:
            return None
            
    def currentStage(self):
        if tc.microscopes:
            #check that stages aren't empty
            if self.currentTEM().stages:
                return self.currentTEM().getStage(self.stagelist.currentText())
            else:
                return None
        else:
            return None
            
    def currentDetector(self):
        if tc.microscopes:
            #check that detectors aren't empty
            if self.currentTEM().detectors:
                return self.currentTEM().getDetector(self.detectorlist.currentText())
            else:
                return None
        else:
            return None
    
    def checkEvent(self, msg = "Are you sure?"):
        reply = QMessageBox.question(self, ' ', msg, QMessageBox.Yes, QMessageBox.No)
        if reply == QMessageBox.Yes:
            return True
        else:
            return False
    
    def updateAll(self, **kwargs):
        self.updateTEMlist(**kwargs)
        self.updateStagelist()
        self.updateDetectorlist()
        self.updateMagoptions()
        self.updatebuttons()
        #the crystals are somewhat dependent on this so we call it here
        self.caller.updateCrystals()
        
    def updateStages(self, **kwargs):
        self.updateStagelist(**kwargs)
        self.updatebuttons()
        #the crystals are somewhat dependent on this so we call it here
        self.caller.updateCrystals()
        
    def updateDetectors(self, **kwargs):
        self.updateDetectorlist(**kwargs)
        self.updateMagoptions()
        self.updatebuttons()
        #the crystals are somewhat dependent on this so we call it here
        self.caller.updateCrystals()
    
    def updateTEMlist(self, focus = None):
        #save whichever name was active
        currentstrucname = ""
        if focus is None:
            try:
                currentstrucname = self.currentTEM().name
            except:
                pass
        else:
            currentstrucname = focus
        
        #update the TEM dropdown list
        self.temlist.clear()
        self.temlist.addItems(list(tc.microscopes.keys()))
        
        #set the TEM again that was active before if possible
        #first find the index of the current TEM
        try:
            index = self.temlist.findText(currentstrucname, Qt.MatchFixedString)
            if index>=0:
                self.temlist.setCurrentIndex(index)
        except:
            pass
        
    def updateStagelist(self, focus = None):
        #save whichever was active
        currentdetectorname = ""
        if focus is None:
            try:
                currentdetectorname = self.currentStage().name
            except:
                pass
        else:
            currentdetectorname = focus
            
        #update the stage dropdown list
        self.stagelist.clear()
        #the stages are only those stages that belong to the current TEM, if there is a current TEM
        try:
            self.stagelist.addItems(list(self.currentTEM().stages.keys()))
        except:
            #no TEM exists
            pass
        
        #set the stage again that was active before if possible
        #first find the index of the current stage
        try:
            index = self.stagelist.findText(currentdetectorname, Qt.MatchFixedString)
            if index>=0:
                self.stagelist.setCurrentIndex(index)
        except:
            pass
    
    def updateDetectorlist(self, focus = None):
        #save whichever was active
        currentdetectorname = ""
        if focus is None:
            try:
                currentdetectorname = self.currentDetector().name
            except:
                pass
        else:
            currentdetectorname = focus
            
        #update the stage dropdown list
        self.detectorlist.clear()
        #the stages are only those stages that belong to the current TEM, if there is a current TEM
        try:
            self.detectorlist.addItems(list(self.currentTEM().detectors.keys()))
        except:
            #no TEM exists
            pass
        
        #set the stage again that was active before if possible
        #first find the index of the current stage
        try:
            index = self.detectorlist.findText(currentdetectorname, Qt.MatchFixedString)
            if index>=0:
                self.detectorlist.setCurrentIndex(index)
        except:
            pass
    
    def updatebuttons(self):
        #Hide certain features if the list of microscopes is empty
        if not tc.microscopes:
            #microscope buttons
            self.edittembut.hide()
            self.deltembut.hide()
            self.kvbox.hide()
            #stage buttons
            self.stagelist.hide()
            self.adstagbut.hide()
            self.editstagbut.hide()
            self.delstagbut.hide()
            #the alpha and beta boxes are hidden
            self.alphabox.hide()
            self.betabox.hide()
            #the detector buttons
            self.detectorlist.hide()
            self.addetbut.hide()
            self.editdetbut.hide()
            self.deldetbut.hide()
            #detector calibration buttons
            self.loadcalib.hide()
            self.imgmode.hide()
            self.thetaview.hide()
            self.magcl.hide()
            
        else: #there is a microscope
            self.edittembut.show()
            self.deltembut.show()
            self.kvbox.show()
            #update the kV box
            self.kvbox.setValue(self.currentTEM().getKv())
            
            #stage buttons
            self.stagelist.show()
            self.adstagbut.show()
            
            #detector buttons
            self.detectorlist.show()
            self.addetbut.show()
            
            #check if there are stages
            if self.currentTEM().stages:
                self.editstagbut.show()
                self.delstagbut.show()
                
                #update alpha and beta boxes and their max and mins. First set maxes otherwise below zero doesn't work.
                self.alphabox.setMaximum(self.currentStage().getAlphaMax())
                self.alphabox.setMinimum(self.currentStage().getAlphaMin())
                self.alphabox.setValue(self.currentStage().getAlpha())
                
                self.betabox.setMaximum(self.currentStage().getBetaMax())
                self.betabox.setMinimum(self.currentStage().getBetaMin())
                self.betabox.setValue(self.currentStage().getBeta())
                
                self.alphabox.show()
                self.betabox.show()
            else:
                self.editstagbut.hide()
                self.delstagbut.hide()
                self.alphabox.hide()
                self.betabox.hide()
                
            #check if there are detectors
            if self.currentTEM().detectors:
                self.editdetbut.show()
                self.deldetbut.show()
                
                self.loadcalib.show()
                self.imgmode.show()
                #check if there is a relevant calibration file
                if self.currentDetector().getCalibration(self.imgmode.currentText()):
                    self.thetaview.show()
                    self.magcl.show()
                else:
                    self.thetaview.hide()
                    self.magcl.hide()
                
            
            else:
                self.editdetbut.hide()
                self.deldetbut.hide()
                self.loadcalib.hide()
                self.imgmode.hide()
                self.thetaview.hide()
                self.magcl.hide()
                
class structuresMenu(QWidget):
    
    def __init__(self, caller):
        ###########The structures menu#################
        super().__init__()
        ##make a grouping box
        self.formGroupBox = QGroupBox("Crystallography")
        
        #caller is the thing that calls this menu#########
        self.caller = caller
        
        #inside the grouping box have a grid layout
        layout = QGridLayout()
        layout.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        layout.setSpacing(10)
        
        #scrollable list of structures
        self.strclbl = QLabel("Structure")
        self.strclbl.setToolTip("Crystal structures serving as templates to Crystals/Grains objects.")
        layout.addWidget(self.strclbl, 1, 0)
        
        #structures
        self.struclist = self.combobox(list(tc.structures.keys()), action = self.doNothing)
        layout.addWidget(self.struclist, 1, 1)
        
        #buttons
        self.adstrucbut = self.button(logo=".\Images\plus.png", hint = "Add new structure",  action = self.createStructure)
        self.editstrucbut = self.button(logo = ".\Images\edit.png", hint = "Edit structure",  action = self.editStructure)
        self.delstrucbut = self.button(logo=".\Images\delete-icon.png", hint = "Delete structure",  action = self.deleteStructure)
        self.calcbut = self.button(logo=".\Images\calc.png", hint = "Crystallography calculator",  action = self.showCalculator)
        
        layout.addWidget(self.adstrucbut, 1, 2)
        layout.addWidget(self.editstrucbut, 1, 3)
        layout.addWidget(self.delstrucbut, 1, 4)
        layout.addWidget(self.calcbut, 1, 5)
        
        
        
        ####set the layout of the microscope control box
        self.formGroupBox.setLayout(layout)
        
        self.updateButtons()
        
        mainLayout = QVBoxLayout()
        mainLayout.addWidget(self.formGroupBox)
        #mainLayout.addWidget(self.buttonBox)
        self.setLayout(mainLayout)
        
        #show the window
        self.show()
    
    def checkEvent(self, msg = "Are you sure?"):
        reply = QMessageBox.question(self, ' ', msg, QMessageBox.Yes, QMessageBox.No)
        if reply == QMessageBox.Yes:
            return True
        else:
            return False
    
    def button(self, text ="" , logo = "", hint = "", position = [0,0], size = [24, 24], action = None):
        if action is None:
            action = self.nothingHappens
        button = QPushButton(text, self)
        button.setToolTip(hint)
        #button.move(position[0], position[1])
        button.setIcon(QIcon(logo))
        button.clicked.connect(action)
        button.resize(size[0]+6, size[1]+6)
        button.setIconSize(QSize(size[0],size[1]))
        button.setSizePolicy(sp)
        return button
    
    def combobox(self, lst, pos = [0,0], action = None):
        """combobox creates a drop down list from a list of strings lst at position [x, y]. When one option is active an action is performed."""
        if action is None:
            action = self.nothingHappens
        combobx = QComboBox(self)
        combobx.addItems(lst)
        #combobx.move(pos[0], pos[1])
        combobx.activated[str].connect(action)
        combobx.setSizePolicy(sp)
        return combobx
    
    def doNothing(self):
        pass
    
    def anglesArePossible(self, ang):
        """This function returns whether a list of 3 angles is possible or not. The sum of the angles needs to be smaller than 360 degrees. Also the sum of the two smallest angles needs to be larger than the largest angle. Angles are supplied in degrees."""
        ang = sorted(ang)
        cond1 = sum(ang)<360
        cond2 = ang[0] + ang[1]>ang[2]
        return (cond1 and cond2)
    
    def createStructure(self):
        res, okPressed = structureDialog.getInfo()
        if okPressed:
            #need to check if that name doesn't already exist or the structure will overwrite items
            if res[0] not in tc.structures:
                if self.anglesArePossible(res[4:7]):
                    newstruc = tc.addStructure(name = res[0], a=res[1], b=res[2], c=res[3], alpha=res[4], beta=res[5], gamma=res[6])
                    self.updateAll(focus = newstruc.name)
                else:
                    QMessageBox.warning(self, " ", "The crystal angles are invalid.")
            else:
                QMessageBox.warning(self, " ", "A structure by this name already exists. Please choose another name.")
    
    def editStructure(self):
        cur = self.currentStruc()
        res, okPressed = structureDialog.getInfo(windowtitle = "Edit structure", name = cur.name, a = cur.a, b = cur.b, c = cur.c, alpha = cur.alpha, beta = cur.beta, gamma = cur.gamma)
        if okPressed:
            if self.anglesArePossible(res[4:7]): #res[4, 5, 6]
                #if the name is nothing or it already exists, keep old name
                if res[0]=="" or res[0] in tc.structures:
                    res[0] = self.currentStruc().name
                self.currentStruc().changeCrystallography(*res)
                self.updateAll(focus = res[0])
            else:
                QMessageBox.warning(self, " ", "The crystal angles are invalid.")
    
    def deleteStructure(self):
        nt = self.currentStruc().name
        check = self.checkEvent(msg = "Are you sure you want to delete %s? All Crystals with this structure will also be deleted!" %(nt))
        if check:
            tc.removeStructure(nt)
            self.updateAll()
    
    def showCalculator(self):
        pass
    
    def currentStruc(self):
        try:
            return tc.getStructure(self.struclist.currentText())
        except:
            return None
        
    def updateAll(self, **kwargs):
        self.updateStructurelist(**kwargs)
        self.updateButtons()
        self.caller.updateCrystals() #also update the crystals
    
    def updateStructurelist(self, focus = None):
        #save whichever name was active
        currentstrucname = ""
        if focus is None:
            try:
                currentstrucname = self.currentStruc().name
            except:
                pass
        else:
            currentstrucname = focus
        
        #update the structure dropdown list
        self.struclist.clear()
        self.struclist.addItems(list(tc.structures.keys()))
        
        #set the structure again that was active before or the one that was newly added
        try:
            index = self.struclist.findText(currentstrucname, Qt.MatchFixedString)
            if index>=0:
                self.struclist.setCurrentIndex(index)
        except:
            pass
    
    def updateButtons(self):
        #if there are no structures
        if not tc.structures:
            self.editstrucbut.hide()
            self.delstrucbut.hide()
            self.calcbut.hide()
            
        else:
            self.editstrucbut.show()
            self.delstrucbut.show()
            self.calcbut.show()
    
class crystalMenu(QWidget):

    def __init__(self, caller):
        ###########The crystals menu#################
        super().__init__()
        
        self.caller = caller
        
        self.formGroupBox = QGroupBox("Crystals/Grain")
        
        
        #inside the grouping box have a grid layout
        layout = QGridLayout()
        layout.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        layout.setSpacing(10)
        
        #scrollable list of structures
        self.cryslabel = QLabel("Crystal")
        layout.addWidget(self.cryslabel, 1, 0)
        self.cryslabel.setToolTip("A microscope, stage, calibrated detector and structure are prerequisites to Crystal creation and editing.")
        
        
        #crystals - empty in the beginning
        self.crystallist = self.combobox(list([]), action = self.changeToolTip)
        layout.addWidget(self.crystallist, 1, 1)
        
        #buttons
        self.addcrystalbut = self.button(logo=".\Images\plus.png", hint = "Add new crystal",  action = self.addCrystal)
        self.editcrystalbut = self.button(logo = ".\Images\edit.png", hint = "Edit crystal",  action = self.editCrystal)
        self.delcrystalbut = self.button(logo=".\Images\delete-icon.png", hint = "Delete crystal",  action = self.deleteCrystal)
        
        layout.addWidget(self.addcrystalbut, 1, 2)
        layout.addWidget(self.editcrystalbut, 1, 3)
        layout.addWidget(self.delcrystalbut, 1, 4)
        
        ####set the layout of the control box
        self.formGroupBox.setLayout(layout)
        
        self.updateButtons()
        
        mainLayout = QVBoxLayout()
        mainLayout.addWidget(self.formGroupBox)
        #mainLayout.addWidget(self.buttonBox)
        self.setLayout(mainLayout)
        
        #show the window
        self.show()
        
    def addCrystal(self):
        res, okPressed = crystalDialog.getInfo(self.caller)
        if okPressed:
            stage = self.caller.getCurrentStage()
            struc = self.caller.getCurrentStructure()
            detc = self.caller.getCurrentDetector()
            mod = self.caller.getMode()
            set = float(self.caller.getSetting())
            #need to check if that name doesn't already exist or the structure will overwrite items
            if res[0] not in stage.crystals:
                name = res[0]
                u = res[1]
                v = res[2]
                w = res[3]
                h = res[4]
                k = res[5]
                l = res[6]
                thet = res[7]
                cmnt = res[8]
                if (u*h + v*k + w*l != 0) or (u==0 and v==0 and w==0) or (h==0 and k==0 and l==0):
                    QMessageBox.warning(self, " ", "Those indices are not valid! The reflection must lie in the zone (hu+kv+lw=0) and not all indices can be 0.")
                    
                else:
                    newcrystal = stage.addCrystal(struc.name, name = name)
                    newcrystal.calcOrient(tc.miller(u, v, w), tc.miller(h, k, l), thet, detc.name, mod, set)
                    newcrystal.setComment(cmnt)
                    self.updateAll(focus = newcrystal.name)
            else:
                QMessageBox.warning(self, " ", "A structure by this name already exists. Please choose another name.")
    
    def checkEvent(self, msg = "Are you sure?"):
        reply = QMessageBox.question(self, ' ', msg, QMessageBox.Yes, QMessageBox.No)
        if reply == QMessageBox.Yes:
            return True
        else:
            return False
    
    def button(self, text ="" , logo = "", hint = "", position = [0,0], size = [24, 24], action = None):
        if action is None:
            action = self.nothingHappens
        button = QPushButton(text, self)
        button.setToolTip(hint)
        #button.move(position[0], position[1])
        button.setIcon(QIcon(logo))
        button.clicked.connect(action)
        button.resize(size[0]+6, size[1]+6)
        button.setIconSize(QSize(size[0],size[1]))
        button.setSizePolicy(sp)
        return button
    
    def combobox(self, lst, pos = [0,0], action = None):
        """combobox creates a drop down list from a list of strings lst at position [x, y]. When one option is active an action is performed."""
        if action is None:
            action = self.nothingHappens
        combobx = QComboBox(self)
        combobx.addItems(lst)
        #combobx.move(pos[0], pos[1])
        combobx.activated[str].connect(action)
        combobx.setSizePolicy(sp)
        return combobx
    
    def doNothing(self):
        pass
        
    def changeToolTip(self):
        #must change the tooltip to the comment of the crystal
        try:
            self.crystallist.setToolTip(self.caller.getCurrentStage().getCrystal(self.crystallist.currentText()).getComment())
        except: #if no crystal is found
            self.crystallist.setToolTip("")
    
    def editCrystal(self):
        
        cur = self.currentCrystal()
        res, okPressed = crystalDialog.getInfo(self.caller, windowtitle = "Edit crystal", name = cur.name, comt = cur.getComment(), strucname = cur.getStructure().name)
        if okPressed:
            stage = self.caller.getCurrentStage()
            struc = self.caller.getCurrentStructure()
            detc = self.caller.getCurrentDetector()
            mod = self.caller.getMode()
            set = float(self.caller.getSetting())
            #need to check if that name doesn't already exist or the structure will overwrite items
            if res[0]=="" or res[0] in stage.crystals:
                res[0] = self.currentCrystal().name
            name = res[0]
            u = res[1]
            v = res[2]
            w = res[3]
            h = res[4]
            k = res[5]
            l = res[6]
            thet = res[7]
            cmnt = res[8]
            if (u*h + v*k + w*l != 0) or (u==0 and v==0 and w==0) or (h==0 and k==0 and l==0):
                QMessageBox.warning(self, " ", "Those indices are not valid! The reflection must lie in the zone (hu+kv+lw=0) and not all indices can be 0.")
    
            else:
                cur.calcOrient(tc.miller(u, v, w), tc.miller(h, k, l), thet, detc.name, mod, set)
                cur.setComment(cmnt)
                cur.setName(name)
                self.updateAll(focus = res[0])
                
    def deleteCrystal(self):
        nt = self.currentCrystal().name
        check = self.checkEvent(msg = "Are you sure you want to delete %s?" %(nt))
        if check:
            self.caller.getCurrentStage().removeCrystal(nt)
            self.updateAll()
    
    def currentCrystal(self):
        try:
            return self.caller.getCurrentStage().getCrystal(self.crystallist.currentText())
        except:
            return None
        
    def updateAll(self, **kwargs):
        self.updateCrystallist(**kwargs)
        self.updateButtons()
    
    def updateCrystallist(self, focus = None):
        #save whichever name was active
        currentCrystalname = ""
        if focus is None:
            try:
                currentCrystalname = self.currentCrystal().name
            except:
                pass
        else:
            currentCrystalname = focus
        
        #update the structure dropdown list
        self.crystallist.clear()
        self.crystallist.setToolTip("")
        try: #if there is a stage
            self.crystallist.addItems(list(self.caller.getCurrentStage().crystals.keys()))
        except: #there is no active stage or no stage at all
            pass
        
        #set the structure again that was active before or the one that was newly added
        try:
            index = self.crystallist.findText(currentCrystalname, Qt.MatchFixedString)
            if index>=0:
                self.crystallist.setCurrentIndex(index)
                #must change the tooltip to the comment of the crystal
                self.crystallist.setToolTip(self.caller.getCurrentStage().getCrystal(currentCrystalname).getComment())
        except:
            pass
    
    def updateButtons(self):
        #if there are no structures
        stage = self.caller.getCurrentStage()
        struc = self.caller.getCurrentStructure()
        detc = self.caller.getCurrentDetector()
        mod = self.caller.getMode()
        set = self.caller.getSetting()
        if stage is None or struc is None or detc is None or mod=="" or set=="":
            self.addcrystalbut.hide()
            self.editcrystalbut.hide()
            self.delcrystalbut.hide()
            
            
        else:
            
            self.addcrystalbut.show()
            #are there any crystals?
            if self.currentCrystal() is not None:
                self.editcrystalbut.show()
                self.delcrystalbut.show()
            else:
                self.editcrystalbut.hide()
                self.delcrystalbut.hide()

class calculationMenu(QWidget):
    
    def __init__(self, caller):
        ###########The menu for doing calculations#################
        super().__init__()
        
        self.caller = caller
        
        self.formGroupBox = QGroupBox("Calculation pane")
        layout = QHBoxLayout()
        layout.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        
        #mx is the maximum and minimum index of spinboxes
        mx = 20
        
        ##########Create a tabbed widget#############
        self.tabs = QTabWidget()
        
        ###the individual tabs
        self.tabzone = QWidget()
        self.tabzone.setToolTip("Calculate the \u03b1 and \u03b2 necessary for tilting a lattice vector to the optical axis.")
        
        tzlayout = QGridLayout()
        tzlayout.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        tzlayout.addWidget(QLabel("Space"), 1, 0)
        
        self.tzspacebox = QComboBox() #which space?
        self.tzspacebox.addItem("real")
        self.tzspacebox.addItem("recyprocal")
        self.tzspacebox.setSizePolicy(sp)
        tzlayout.addWidget(self.tzspacebox, 2, 0)
        
        
        self.tzuh = QLabel("u/h") #uvw/hkl boxes
        self.tzvk = QLabel("v/k")
        self.tzwl = QLabel("w/l")
        tzlayout.addWidget(self.tzuh, 1, 1)
        tzlayout.addWidget(self.tzvk, 1, 2)
        tzlayout.addWidget(self.tzwl, 1, 3)
        
        self.tzuhbox = QSpinBox()
        self.tzuhbox.setSizePolicy(sp)
        self.tzuhbox.setMaximum(mx)
        self.tzuhbox.setMinimum(-mx)
        tzlayout.addWidget(self.tzuhbox, 2, 1)
        
        self.tzvkbox = QSpinBox()
        self.tzvkbox.setSizePolicy(sp)
        self.tzvkbox.setMaximum(mx)
        self.tzvkbox.setMinimum(-mx)
        tzlayout.addWidget(self.tzvkbox, 2, 2)
        
        self.tzwlbox = QSpinBox()
        self.tzwlbox.setSizePolicy(sp)
        self.tzwlbox.setMaximum(mx)
        self.tzwlbox.setMinimum(-mx)
        tzlayout.addWidget(self.tzwlbox, 2, 3)
        
        self.tzcalcbut = QPushButton("") #calculate button
        self.tzcalcbut.setToolTip("Calculate")
        self.tzcalcbut.setIcon(QIcon(".\Images\calculate.png"))
        self.tzcalcbut.setSizePolicy(sp)
        self.tzcalcbut.clicked.connect(self.tzcalc)
        tzlayout.addWidget(self.tzcalcbut, 2, 4)
        
        tzlayout.addWidget(QLabel("[\u03b1,\u03b2]"), 1, 5) #alpha and beta display
        
        self.tzRES = []
        
        self.tzalbox = QLineEdit()
        self.tzalbox.setReadOnly(True)
        self.tzalbox.setSizePolicy(sp)
        tzlayout.addWidget(self.tzalbox, 2, 5)
        
        self.tzsetbut = QPushButton("Set")
        self.tzsetbut.setToolTip("Set the stage tilt to the calculated tilt.")
        self.tzsetbut.setSizePolicy(sp)
        self.tzsetbut.clicked.connect(self.tzset)
        tzlayout.addWidget(self.tzsetbut, 2, 6)
        
        self.tabzone.setLayout(tzlayout)
        
        ##################################################
        
        self.tab2beam = QWidget()
        self.tab2beam.setToolTip("Calculate the \u03b1 and \u03b2 necessary for tilting to -g Bragg condition.\nAfter tilting, tilt the beam so g is centered resulting in a g centered 2-beam.")
        
        t2blayout = QGridLayout()
        t2blayout.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        t2blayout.addWidget(QLabel("Reflection"), 2, 0)
        t2blayout.addWidget(QLabel("Zone"), 3, 0)
        
        self.t2buh = QLabel("h/u") #uvw/hkl boxes
        self.t2bvk = QLabel("k/v")
        self.t2bwl = QLabel("l/w")
        t2blayout.addWidget(self.t2buh, 1, 1)
        t2blayout.addWidget(self.t2bvk, 1, 2)
        t2blayout.addWidget(self.t2bwl, 1, 3)
        
        t2bmx = 20
        t2bdec = 3
        self.t2bhbox = QSpinBox()
        self.t2bhbox.setSingleStep(1)
        self.t2bhbox.setMinimum(-mx)
        self.t2bhbox.setMaximum(mx)
        self.t2bhbox.setSizePolicy(sp)
        t2blayout.addWidget(self.t2bhbox, 2, 1)
        
        self.t2bkbox = QSpinBox()
        self.t2bkbox.setSingleStep(1)
        self.t2bkbox.setMinimum(-mx)
        self.t2bkbox.setMaximum(mx)
        self.t2bkbox.setSizePolicy(sp)
        t2blayout.addWidget(self.t2bkbox, 2, 2)
        
        self.t2blbox = QSpinBox()
        self.t2blbox.setSingleStep(1)
        self.t2blbox.setMinimum(-mx)
        self.t2blbox.setMaximum(mx)
        self.t2blbox.setSizePolicy(sp)
        t2blayout.addWidget(self.t2blbox, 2, 3)
        
        self.t2bubox = QSpinBox()
        self.t2bubox.setSingleStep(1)
        self.t2bubox.setMinimum(-mx)
        self.t2bubox.setMaximum(mx)
        self.t2bubox.setSizePolicy(sp)
        t2blayout.addWidget(self.t2bubox, 3, 1)
        
        self.t2bvbox = QSpinBox()
        self.t2bvbox.setSingleStep(1)
        self.t2bvbox.setMinimum(-mx)
        self.t2bvbox.setMaximum(mx)
        self.t2bvbox.setSizePolicy(sp)
        t2blayout.addWidget(self.t2bvbox, 3, 2)
        
        self.t2bwbox = QSpinBox()
        self.t2bwbox.setSingleStep(1)
        self.t2bwbox.setMinimum(-mx)
        self.t2bwbox.setMaximum(mx)
        self.t2bwbox.setSizePolicy(sp)
        t2blayout.addWidget(self.t2bwbox, 3, 3)
        
        #Angle box
        t2blayout.addWidget(QLabel("Angle"), 1, 4)
        
        self.t2bangbox = QDoubleSpinBox()
        self.t2bangbox.setSuffix("\u00b0")
        self.t2bangbox.setDecimals(t2bdec)
        self.t2bangbox.setSingleStep(10**(-t2bdec))
        self.t2bangbox.setMinimum(-90)
        self.t2bangbox.setMaximum(90)
        self.t2bangbox.setValue(8)
        self.t2bangbox.setSizePolicy(sp)
        t2blayout.addWidget(self.t2bangbox, 2, 4)
        
        self.t2bsetbut = QPushButton("Current") #set the current zone
        self.t2bsetbut.setToolTip("Set the values to the current zone.\nYou must be in a rational index (max(u, v, w)<10) zone.")
        self.t2bsetbut.setSizePolicy(sp)
        self.t2bsetbut.clicked.connect(self.t2bset)
        t2blayout.addWidget(self.t2bsetbut, 3, 4)
        
        self.t2bcalcbut = QPushButton("") #calculate button
        self.t2bcalcbut.setToolTip("Calculate")
        self.t2bcalcbut.setIcon(QIcon(".\Images\calculate.png"))
        self.t2bcalcbut.setSizePolicy(sp)
        self.t2bcalcbut.clicked.connect(self.t2bcalc)
        t2blayout.addWidget(self.t2bcalcbut, 2, 5)
        
        t2blayout.addWidget(QLabel("[\u03b1,\u03b2]"), 1, 6) #alpha and beta display
        
        self.t2bRES = []
        self.t2balbox = QLineEdit()
        self.t2balbox.setReadOnly(True)
        self.t2balbox.setSizePolicy(sp)
        t2blayout.addWidget(self.t2balbox, 2, 6)
        
        self.t2bsetbut2 = QPushButton("Set") #set the current zone
        self.t2bsetbut2.setToolTip("Set tilt to the calculated values")
        self.t2bsetbut2.setSizePolicy(sp)
        self.t2bsetbut2.clicked.connect(self.t2bset2)
        t2blayout.addWidget(self.t2bsetbut2, 2, 7)
        
        self.tab2beam.setLayout(t2blayout)
        
        
        ###################################################
        self.tabtiltround = QWidget()
        
        self.tabtiltround.setToolTip("Calculate the \u03b1 and \u03b2 necessary for tilting around reflection (hkl) by a certain angle.")
        
        ttrlayout = QGridLayout()
        ttrlayout.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        ttrlayout.addWidget(QLabel("Reflection"), 2, 0)
        ttrlayout.addWidget(QLabel("Zone"), 3, 0)
        
        self.ttruh = QLabel("h/u") #uvw/hkl boxes
        self.ttrvk = QLabel("k/v")
        self.ttrwl = QLabel("l/w")
        ttrlayout.addWidget(self.ttruh, 1, 1)
        ttrlayout.addWidget(self.ttrvk, 1, 2)
        ttrlayout.addWidget(self.ttrwl, 1, 3)
        
        ttrmx = 20
        ttrdec = 3
        self.ttrhbox = QSpinBox()
        self.ttrhbox.setSingleStep(1)
        self.ttrhbox.setMinimum(-mx)
        self.ttrhbox.setMaximum(mx)
        self.ttrhbox.setSizePolicy(sp)
        ttrlayout.addWidget(self.ttrhbox, 2, 1)
        
        self.ttrkbox = QSpinBox()
        self.ttrkbox.setSingleStep(1)
        self.ttrkbox.setMinimum(-mx)
        self.ttrkbox.setMaximum(mx)
        self.ttrkbox.setSizePolicy(sp)
        ttrlayout.addWidget(self.ttrkbox, 2, 2)
        
        self.ttrlbox = QSpinBox()
        self.ttrlbox.setSingleStep(1)
        self.ttrlbox.setMinimum(-mx)
        self.ttrlbox.setMaximum(mx)
        self.ttrlbox.setSizePolicy(sp)
        ttrlayout.addWidget(self.ttrlbox, 2, 3)
        
        self.ttrubox = QSpinBox()
        self.ttrubox.setSingleStep(1)
        self.ttrubox.setMinimum(-mx)
        self.ttrubox.setMaximum(mx)
        self.ttrubox.setSizePolicy(sp)
        ttrlayout.addWidget(self.ttrubox, 3, 1)
        
        self.ttrvbox = QSpinBox()
        self.ttrvbox.setSingleStep(1)
        self.ttrvbox.setMinimum(-mx)
        self.ttrvbox.setMaximum(mx)
        self.ttrvbox.setSizePolicy(sp)
        ttrlayout.addWidget(self.ttrvbox, 3, 2)
        
        self.ttrwbox = QSpinBox()
        self.ttrwbox.setSingleStep(1)
        self.ttrwbox.setMinimum(-mx)
        self.ttrwbox.setMaximum(mx)
        self.ttrwbox.setSizePolicy(sp)
        ttrlayout.addWidget(self.ttrwbox, 3, 3)
        
        #Angle box
        ttrlayout.addWidget(QLabel("Angle"), 1, 4)
        
        self.ttrangbox = QDoubleSpinBox()
        self.ttrangbox.setSuffix("\u00b0")
        self.ttrangbox.setDecimals(ttrdec)
        self.ttrangbox.setSingleStep(10**(-ttrdec))
        self.ttrangbox.setMinimum(-90)
        self.ttrangbox.setMaximum(90)
        self.ttrangbox.setSizePolicy(sp)
        ttrlayout.addWidget(self.ttrangbox, 2, 4)
        
        self.ttrsetbut = QPushButton("Current") #set the current zone
        self.ttrsetbut.setToolTip("Set the values to the current zone.\nYou must be in a rational index (max(u, v, w)<10) zone.")
        self.ttrsetbut.setSizePolicy(sp)
        self.ttrsetbut.clicked.connect(self.ttrset)
        ttrlayout.addWidget(self.ttrsetbut, 3, 4)
        
        self.ttrcalcbut = QPushButton("") #calculate button
        self.ttrcalcbut.setToolTip("Calculate")
        self.ttrcalcbut.setIcon(QIcon(".\Images\calculate.png"))
        self.ttrcalcbut.setSizePolicy(sp)
        self.ttrcalcbut.clicked.connect(self.ttrcalc)
        ttrlayout.addWidget(self.ttrcalcbut, 2, 5)
        
        ttrlayout.addWidget(QLabel("[\u03b1,\u03b2]"), 1, 6) #alpha and beta display
        
        self.ttrRES = []
        self.ttralbox = QLineEdit()
        self.ttralbox.setReadOnly(True)
        self.ttralbox.setSizePolicy(sp)
        ttrlayout.addWidget(self.ttralbox, 2, 6)
        
        self.ttrsetbut2 = QPushButton("Set") #set the current zone
        self.ttrsetbut2.setToolTip("Set tilt to the calculated values")
        self.ttrsetbut2.setSizePolicy(sp)
        self.ttrsetbut2.clicked.connect(self.ttrset2)
        ttrlayout.addWidget(self.ttrsetbut2, 2, 7)
        
        self.tabtiltround.setLayout(ttrlayout)
        
        ##################################################
        
        self.tabweakbeam = QWidget()
        
        self.tabweakbeam.setToolTip("Calculate the \u03b1 and \u03b2 necessary for tilting around reflection (hkl) by a certain angle.")
        
        twblayout = QGridLayout()
        twblayout.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        twblayout.addWidget(QLabel("Reflection"), 2, 0)
        twblayout.addWidget(QLabel("Zone"), 3, 0)
        
        twbmx = 20
        twbdec = 3
        
        tlb = QLabel("N")
        tlb.setToolTip("The desired N in g(Ng) after tilting sample and beam.")
        twblayout.addWidget(tlb, 1,5)
        
        self.twbNbox = QDoubleSpinBox()
        self.twbNbox.setDecimals(twbdec)
        self.twbNbox.setSingleStep(1)
        self.twbNbox.setMinimum(1)
        self.twbNbox.setMaximum(mx)
        self.twbNbox.setSizePolicy(sp)
        twblayout.addWidget(self.twbNbox, 2, 5)
        
        tlb = QLabel("K")
        tlb.setToolTip("Reflection Kg will be centered by beam tilt and used for imaging.")
        twblayout.addWidget(tlb, 1,6)
        
        self.twbKbox = QSpinBox()
        self.twbKbox.setSizePolicy(sp)
        self.twbKbox.setMaximum(mx)
        self.twbKbox.setMinimum(1)
        self.twbKbox.setValue(1)
        twblayout.addWidget(self.twbKbox, 2, 6)
        
        self.twbuh = QLabel("h/u") #uvw/hkl boxes
        self.twbvk = QLabel("k/v")
        self.twbwl = QLabel("l/w")
        twblayout.addWidget(self.twbuh, 1, 1)
        twblayout.addWidget(self.twbvk, 1, 2)
        twblayout.addWidget(self.twbwl, 1, 3)
        
        
        self.twbhbox = QSpinBox()
        self.twbhbox.setSingleStep(1)
        self.twbhbox.setMinimum(-mx)
        self.twbhbox.setMaximum(mx)
        self.twbhbox.setSizePolicy(sp)
        twblayout.addWidget(self.twbhbox, 2, 1)
        
        self.twbkbox = QSpinBox()
        self.twbkbox.setSingleStep(1)
        self.twbkbox.setMinimum(-mx)
        self.twbkbox.setMaximum(mx)
        self.twbkbox.setSizePolicy(sp)
        twblayout.addWidget(self.twbkbox, 2, 2)
        
        self.twblbox = QSpinBox()
        self.twblbox.setSingleStep(1)
        self.twblbox.setMinimum(-mx)
        self.twblbox.setMaximum(mx)
        self.twblbox.setSizePolicy(sp)
        twblayout.addWidget(self.twblbox, 2, 3)
        
        self.twbubox = QSpinBox()
        self.twbubox.setSingleStep(1)
        self.twbubox.setMinimum(-mx)
        self.twbubox.setMaximum(mx)
        self.twbubox.setSizePolicy(sp)
        twblayout.addWidget(self.twbubox, 3, 1)
        
        self.twbvbox = QSpinBox()
        self.twbvbox.setSingleStep(1)
        self.twbvbox.setMinimum(-mx)
        self.twbvbox.setMaximum(mx)
        self.twbvbox.setSizePolicy(sp)
        twblayout.addWidget(self.twbvbox, 3, 2)
        
        self.twbwbox = QSpinBox()
        self.twbwbox.setSingleStep(1)
        self.twbwbox.setMinimum(-mx)
        self.twbwbox.setMaximum(mx)
        self.twbwbox.setSizePolicy(sp)
        twblayout.addWidget(self.twbwbox, 3, 3)
        
        #Angle box
        twblayout.addWidget(QLabel("Angle"), 1, 4)
        
        self.twbangbox = QDoubleSpinBox()
        self.twbangbox.setSuffix("\u00b0")
        self.twbangbox.setDecimals(twbdec)
        self.twbangbox.setSingleStep(10**(-twbdec))
        self.twbangbox.setMinimum(-90)
        self.twbangbox.setMaximum(90)
        self.twbangbox.setSizePolicy(sp)
        self.twbangbox.setValue(8)
        twblayout.addWidget(self.twbangbox, 2, 4)
        
        self.twbsetbut = QPushButton("Current") #set the current zone
        self.twbsetbut.setToolTip("Set the values to the current zone.\nYou must be in a rational index (max(u, v, w)<10) zone.")
        self.twbsetbut.setSizePolicy(sp)
        self.twbsetbut.clicked.connect(self.twbset)
        twblayout.addWidget(self.twbsetbut, 3, 4)
        
        self.twbcalcbut = QPushButton("") #calculate button
        self.twbcalcbut.setToolTip("Calculate")
        self.twbcalcbut.setIcon(QIcon(".\Images\calculate.png"))
        self.twbcalcbut.setSizePolicy(sp)
        self.twbcalcbut.clicked.connect(self.twbcalc)
        twblayout.addWidget(self.twbcalcbut, 2, 7)
        
        twblayout.addWidget(QLabel("[\u03b1,\u03b2]"), 2, 8) #alpha and beta display
        
        self.twbRES = []
        self.twbalbox = QLineEdit()
        self.twbalbox.setReadOnly(True)
        self.twbalbox.setSizePolicy(sp)
        twblayout.addWidget(self.twbalbox, 2, 9)
        
        self.twbsetbut2 = QPushButton("Set") #set the current zone
        self.twbsetbut2.setToolTip("Set tilt to the calculated values")
        self.twbsetbut2.setSizePolicy(sp)
        self.twbsetbut2.clicked.connect(self.twbset2)
        twblayout.addWidget(self.twbsetbut2, 2, 10)
        
        n0lab = QLabel("N0")
        n0lab.setToolTip("N0g is tilted to Bragg orientation before beam tilting to Kg")
        twblayout.addWidget(n0lab, 3, 8)
        
        self.twbn0box = QLineEdit()
        self.twbn0box.setReadOnly(True)
        self.twbn0box.setSizePolicy(sp)
        twblayout.addWidget(self.twbn0box, 3, 9)
        
        sglab = QLabel("Sg")
        sglab.setToolTip("Geometric excitation error")
        twblayout.addWidget(sglab, 4, 8)
        
        self.twbsgbox = QLineEdit()
        self.twbsgbox.setReadOnly(True)
        self.twbsgbox.setSizePolicy(sp)
        twblayout.addWidget(self.twbsgbox, 4, 9)
        
        self.tabweakbeam.setLayout(twblayout)
        
        ##################################################
        self.tabvecatangle = QWidget()
        self.tabvecatangle.setToolTip("Calculate the lattice vector aligned with the optical axis at \u03b1 and \u03b2.")
        
        tvalayout = QGridLayout()
        tvalayout.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        tvalayout.addWidget(QLabel("\u03b1"), 1, 0) #alpha
        
        dec = 2
        self.tvaalbox = QDoubleSpinBox()
        self.tvaalbox.setSuffix("\u00b0")
        self.tvaalbox.setDecimals(dec)
        self.tvaalbox.setSingleStep(10**(-dec))
        self.tvaalbox.setMinimum(-90)
        self.tvaalbox.setMaximum(90)
        self.tvaalbox.setSizePolicy(sp)
        tvalayout.addWidget(self.tvaalbox, 2, 0)
        
        tvalayout.addWidget(QLabel("\u03b2"), 1, 1) #beta
        
        self.tvabebox = QDoubleSpinBox()
        self.tvabebox.setSuffix("\u00b0")
        self.tvabebox.setDecimals(dec)
        self.tvabebox.setSingleStep(10**(-dec))
        self.tvabebox.setMinimum(-90)
        self.tvabebox.setMaximum(90)
        self.tvabebox.setSizePolicy(sp)
        tvalayout.addWidget(self.tvabebox, 2, 1)
        
        self.tvasetcur = QPushButton("Current")
        self.tvasetcur.setSizePolicy(sp)
        self.tvasetcur.clicked.connect(self.tvaSetCur)
        tvalayout.addWidget(self.tvasetcur, 2, 2)
        
        tvalayout.addWidget(QLabel("Space"), 1, 3)
        
        self.tvaspacebox = QComboBox() #which space?
        self.tvaspacebox.addItem("real")
        self.tvaspacebox.addItem("recyprocal")
        self.tvaspacebox.setSizePolicy(sp)
        tvalayout.addWidget(self.tvaspacebox, 2, 3)
        
        
        self.tvacalcbut = QPushButton("") #calculate button
        self.tvacalcbut.setToolTip("Calculate")
        self.tvacalcbut.setIcon(QIcon(".\Images\calculate.png"))
        self.tvacalcbut.setSizePolicy(sp)
        self.tvacalcbut.clicked.connect(self.tvacalc)
        tvalayout.addWidget(self.tvacalcbut, 2, 4)
        
        self.tvauh = QLabel("[u/h, v/k, w/l]") #uvw/hkl box
        tvalayout.addWidget(self.tvauh, 1, 5)
        
        self.tvaresbox = QLineEdit()
        self.tvaresbox.setToolTip("Lattice/recyprocal space vector in zone axis at \u03b1/\u03b2 tilt. An integer representation is given if one can be found.")
        self.tvaresbox.setReadOnly(True)
        self.tvaresbox.setSizePolicy(sp)
        tvalayout.addWidget(self.tvaresbox, 2, 5)
        
        self.tabvecatangle.setLayout(tvalayout)
        
        #######Add the tabs########
        self.tabs.addTab(self.tabzone, "Tab 1")
        self.tabs.addTab(self.tabvecatangle, "Tab 2")
        self.tabs.addTab(self.tabtiltround, "Tab 3")
        self.tabs.addTab(self.tab2beam, "Tab 4")
        self.tabs.addTab(self.tabweakbeam, "Tab 5")
        
        self.tabs.setTabText(0, "\u03b1 / \u03b2 to Zone")
        self.tabs.setTabText(1, "Zone at \u03b1 / \u03b2")
        self.tabs.setTabText(2, "\u03b8 around reflection")
        self.tabs.setTabText(3, "2-beam")
        self.tabs.setTabText(4, "Weak beam")
        
        
        ########Create the layout for the tabs########
        
        layout.addWidget(self.tabs)
        
        self.formGroupBox.setLayout(layout)
        
        mainLayout = QVBoxLayout()
        mainLayout.addWidget(self.formGroupBox)
        self.setLayout(mainLayout)
        
        #show the window
        self.show()
        
    def tzcalc(self):
        """Calculate the alpha and beta necessary to go to a certain zone and display it"""
        #check if the indexes are not all 0
        u = self.tzuhbox.value()
        v = self.tzvkbox.value()
        w = self.tzwlbox.value()
        typ = self.tzspacebox.currentText()
        if not (u==0 and v==0 and w==0):
            crys = self.caller.getCurrentCrystal()
            if crys is not None: #check that there is in fact a crystal
                vec = tc.miller(u, v, w)
                temp = crys.getAlphaBetaMiller(vec, typ = typ).T.tolist()
                if not np.isnan(sum(temp)): #check that the condition is reachable
                    self.tzRES = temp
                    self.tzalbox.setText("[%s\u00b0, %s\u00b0]" %(round(self.tzRES[0], 2), round(self.tzRES[1], 2)))
                else:
                    self.tzalbox.setText("Unreachable zone!")
                    self.tzRES = []
            else:
                self.tzalbox.setText("No active crystal!")
                self.tzRES = []
        else:
            self.tzalbox.setText("Invalid indices!")
            self.tzRES = []
        
        
    def tzset(self):
        """Set the current alpha and beta to the calculated one"""
        try:
            stag = self.caller.getCurrentStage()
            stag.setAlpha(self.tzRES[0])
            stag.setBeta(self.tzRES[1])
            self.caller.updateGUI()
        except: #the stage may have been deleted or the indices invalid (nothing calculated)
            pass
    
    def ttrset(self):
        """Set the zone to the currently calculated one"""
        try:
            crys = self.caller.getCurrentCrystal()
            zon = crys.getZone(typ="real", verbose=False, integer = True)
            if max(np.absolute(zon))<1: #some index exceeded 10 when trying to convert to ints
                self.ttrubox.setValue(0)
                self.ttrvbox.setValue(0)
                self.ttrwbox.setValue(0)
            else: #an integer zone
                self.ttrubox.setValue(zon[0])
                self.ttrvbox.setValue(zon[1])
                self.ttrwbox.setValue(zon[2])
        except: #there is no crystal for whatever reason
            self.ttrubox.setValue(0)
            self.ttrvbox.setValue(0)
            self.ttrwbox.setValue(0)
    
    def ttrset2(self):
        """Set the current alpha and beta to the calculated one"""
        try:
            stag = self.caller.getCurrentStage()
            stag.setAlpha(self.ttrRES[0])
            stag.setBeta(self.ttrRES[1])
            self.caller.updateGUI()
        except: #the stage may have been deleted or the indices invalid
            pass
    
    def ttrcalc(self):
        u = self.ttrubox.value()
        v = self.ttrvbox.value()
        w = self.ttrwbox.value()
        h = self.ttrhbox.value()
        k = self.ttrkbox.value()
        l = self.ttrlbox.value()
        ang = self.ttrangbox.value()
        
        if (u==0 and v==0 and w==0) or (h==0 and k==0 and l==0) or (u*h+k*v+l*w != 0): #invalid indexes
            self.ttralbox.setText("Zone/Ref. not perpendicular!")
            self.ttrRES = []
        else:
            crys = self.caller.getCurrentCrystal()
            if crys is not None: #check that there is in fact a crystal
                vec = tc.miller(u, v, w)
                g = tc.miller(h, k, l)
                temp = crys.getAlphaBetaTiltRound(g, ang, za = vec, rnd = 5, verbose = False).T.tolist()
                if not np.isnan(sum(temp)): #check that the condition is reachable
                    self.ttrRES = temp
                    self.ttralbox.setText("[%s\u00b0, %s\u00b0]" %(round(self.ttrRES[0], 2), round(self.ttrRES[1], 2)))
                else:
                    self.ttralbox.setText("Unreachable position!")
                    self.ttrRES = []
            else:
                self.ttralbox.setText("No active crystal!")
                self.ttrRES = []
        
    def tvacalc(self):
        """Calculate the lattice vector (real or recyprocal) that is at a certain alpha and beta"""
        al = self.tvaalbox.value()
        be = self.tvabebox.value()
        typ = self.tvaspacebox.currentText()
        stag = self.caller.getCurrentStage()
        if stag is not None: #check that there is a stage
            if stag.inrange(al, be, units = "degrees"):#check if you are in range
                crys = self.caller.getCurrentCrystal()
                if crys is not None: #check that there is in fact a crystal
                    temp = crys.getZoneAtAlphaBeta(al, be, typ = typ, rnd = 5, verbose = False, integer = True).T.tolist()
                    self.tvaresbox.setText("[%s %s %s]" %(round(temp[0], 3), round(temp[1], 3), round(temp[2], 3)))
                else:
                    self.tvaresbox.setText("No active crystal!")
            else:
                self.tvaresbox.setText("Tilt out of range!")
        else:
            self.tvaresbox.setText("No active stage!")
        
    def tvaSetCur(self):
        """Set the value of the alpha and beta box in the vec at alpha pane to the current tilt in the active stage."""
        try:
            stag = self.caller.getCurrentStage()
            al = stag.getAlpha()
            be = stag.getBeta()
            self.tvaalbox.setValue(al)
            self.tvabebox.setValue(be)
        except: #the stage may have been deleted
            pass
    
    def t2bset(self):
        """Set the zone to the currently calculated one"""
        try:
            crys = self.caller.getCurrentCrystal()
            zon = crys.getZone(typ="real", verbose=False, integer = True)
            if max(np.absolute(zon))<1: #some index exceeded 10 when trying to convert to ints
                self.t2bubox.setValue(0)
                self.t2bvbox.setValue(0)
                self.t2bwbox.setValue(0)
            else: #an integer zone
                self.t2bubox.setValue(zon[0])
                self.t2bvbox.setValue(zon[1])
                self.t2bwbox.setValue(zon[2])
        except: #there is no crystal for whatever reason
            self.t2bubox.setValue(0)
            self.t2bvbox.setValue(0)
            self.t2bwbox.setValue(0)
    
    def t2bset2(self):
        """Set the current alpha and beta to the calculated one"""
        try:
            stag = self.caller.getCurrentStage()
            stag.setAlpha(self.t2bRES[0])
            stag.setBeta(self.t2bRES[1])
            self.caller.updateGUI()
        except: #the stage may have been deleted or the indices invalid
            pass
    
    def t2bcalc(self):
        u = self.t2bubox.value()
        v = self.t2bvbox.value()
        w = self.t2bwbox.value()
        h = self.t2bhbox.value()
        k = self.t2bkbox.value()
        l = self.t2blbox.value()
        ang = self.t2bangbox.value()
        
        if (u==0 and v==0 and w==0) or (h==0 and k==0 and l==0) or (u*h+k*v+l*w != 0): #invalid indexes
            self.t2balbox.setText("Zone/Ref. not perpendicular!")
            self.t2bRES = []
        else:
            crys = self.caller.getCurrentCrystal()
            if crys is not None: #check that there is in fact a crystal
                vec = tc.miller(u, v, w)
                g = tc.miller(h, k, l)
                temp = crys.getAlphaBetac2Beam(g, n=1, outzone = ang, za = vec, rnd=5, verbose = False).T.tolist()
                if not np.isnan(sum(temp)): #check that the condition is reachable
                    self.t2bRES = temp
                    self.t2balbox.setText("[%s\u00b0, %s\u00b0]" %(round(self.t2bRES[0], 2), round(self.t2bRES[1], 2)))
                else:
                    self.t2balbox.setText("Unreachable position!")
                    self.t2bRES = []
            else:
                self.t2balbox.setText("No active crystal!")
                self.t2bRES = []
        
    def twbset(self):
        """Set the zone to the currently calculated one"""
        try:
            crys = self.caller.getCurrentCrystal()
            zon = crys.getZone(typ="real", verbose=False, integer = True)
            if max(np.absolute(zon))<1: #some index exceeded 10 when trying to convert to ints
                self.twbubox.setValue(0)
                self.twbvbox.setValue(0)
                self.twbwbox.setValue(0)
            else: #an integer zone
                self.twbubox.setValue(zon[0])
                self.twbvbox.setValue(zon[1])
                self.twbwbox.setValue(zon[2])
        except: #there is no crystal for whatever reason
            self.twbubox.setValue(0)
            self.twbvbox.setValue(0)
            self.twbwbox.setValue(0)
        
    def twbset2(self):
        try:
            stag = self.caller.getCurrentStage()
            stag.setAlpha(self.twbRES[0])
            stag.setBeta(self.twbRES[1])
            self.caller.updateGUI()
        except: #the stage may have been deleted or the indices invalid
            pass
        
    def twbcalc(self):
        u = self.twbubox.value()
        v = self.twbvbox.value()
        w = self.twbwbox.value()
        h = self.twbhbox.value()
        k = self.twbkbox.value()
        l = self.twblbox.value()
        ang = self.twbangbox.value()
        
        kref = self.twbKbox.value()
        nref = self.twbNbox.value()
        
        if (u==0 and v==0 and w==0) or (h==0 and k==0 and l==0) or (u*h+k*v+l*w != 0): #invalid indexes
            self.twbalbox.setText("Zone/Ref. not perpendicular!")
            self.twbRES = []
        else:
            #k should be smaller than n
            if kref<nref:
                crys = self.caller.getCurrentCrystal()
                if crys is not None: #check that there is in fact a crystal
                    vec = tc.miller(u, v, w)
                    g = tc.miller(h, k, l)
                    temp = crys.getAlphaBetaWBeam(g, n= nref, k = kref, outzone = ang, za = vec, rnd=5, verbose = False).T.tolist()
                    if not np.isnan(sum(temp)): #check that the condition is reachable
                        self.twbRES = temp
                        self.twbalbox.setText("[%s\u00b0, %s\u00b0]" %(round(self.twbRES[0], 2), round(self.twbRES[1], 2)))
                        n0 = crys.calcWB2beam(g, kref, nref)
                        self.twbn0box.setText(str(round(n0, 3)))
                        sg = crys.calcSgSimple(g, kref, nref)
                        self.twbsgbox.setText("%s 1/nm" %(round(sg, 3)))
                    else:
                        self.twbalbox.setText("Unreachable position!")
                        self.twbn0box.setText("")
                        self.twbsgbox.setText("")
                        self.twbRES = []
                else:
                    self.twbalbox.setText("No active crystal!")
                    self.twbn0box.setText("")
                    self.twbsgbox.setText("")
                    self.twbRES = []
            else:
                self.twbalbox.setText("K must be <N!")
                self.twbn0box.setText("")
                self.twbsgbox.setText("")
                self.twbRES = []
                
                
                
class plottingMenu(QWidget):
    
    def __init__(self, caller):
        ###########The menu for plotting#################
        super().__init__()
        
        self.caller = caller
        
        self.formGroupBox = QGroupBox("Stereographic plots")
        
        
        #inside the grouping box have a grid layout
        layout = QGridLayout()
        layout.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        layout.setSpacing(10)
        
        #scrollable list of structures
        self.cryslabel = QLabel("Space")
        layout.addWidget(self.cryslabel, 1, 0)
        self.cryslabel.setToolTip("Plot real or recyprocal space vectors.")
        
        self.hkllabel = QLabel("Max. HKL")
        layout.addWidget(self.hkllabel, 1, 2)
        self.hkllabel.setToolTip("Highest hkl index to include.")
        
        self.boundlabel = QLabel("Show bounds")
        layout.addWidget(self.boundlabel, 2, 0)
        self.boundlabel.setToolTip("Plot the tilting boundaries.")
        
        self.setslabel = QLabel("Show settings")
        layout.addWidget(self.setslabel, 2, 2)
        self.setslabel.setToolTip("Add a title to the plot with the microscope settings.")
        
        #The actual workings
        self.spacebox = QComboBox()
        self.spacebox.addItem("real")
        self.spacebox.addItem("recyprocal")
        layout.addWidget(self.spacebox, 1, 1)
        
        self.indexbox = QSpinBox()
        self.indexbox.setValue(1)
        self.indexbox.setMaximum(10)
        self.indexbox.setMinimum(1)
        self.indexbox.setSizePolicy(sp)
        layout.addWidget(self.indexbox, 1, 3)
        
        self.boundcheck = QCheckBox("")
        self.boundcheck.setChecked(True)
        layout.addWidget(self.boundcheck, 2, 1)
        
        self.settingcheck = QCheckBox("")
        self.settingcheck.setChecked(True)
        layout.addWidget(self.settingcheck, 2, 3)
        
        self.plotbut = QPushButton("Plot")
        self.plotbut.clicked.connect(self.makePlot)
        self.plotbut.setSizePolicy(sp)
        self.plotbut.setToolTip("Make a stereographic plot of the active Crystal on the active detector.\nA crystal and calibrated detector are prerequisites.")
        layout.addWidget(self.plotbut, 1, 4)
        
        self.commentbox = QLineEdit()
        self.commentbox.setReadOnly(True)
        layout.addWidget(self.commentbox, 2, 4)
        
        self.formGroupBox.setLayout(layout)
        
        self.updateButtons()
        
        mainLayout = QVBoxLayout()
        mainLayout.addWidget(self.formGroupBox)
        #mainLayout.addWidget(self.buttonBox)
        self.setLayout(mainLayout)
        
        #show the window
        self.show()
        
    
    def makePlot(self):
        try:
            crys = self.caller.getCurrentCrystal()
            dt = self.caller.getCurrentDetector().name
            mod = self.caller.getMode()
            set = float(self.caller.getSetting())
            typ = self.spacebox.currentText()
            plotbond = self.boundcheck.isChecked()
            plotax = True
            hkllim = self.indexbox.value()
            lgnd = self.settingcheck.isChecked()
            
            self.commentbox.setText("Making plots...")
            
            plotDialog.getInfo(crys = crys, dt = dt, mod=mod, set=set, typ = typ, plotax=plotax, plotbond=plotbond, hkllim = hkllim, lgnd = lgnd)
        
        except:
            self.commentbox.setText("Can't make a plot. You need a crystal and calibrated detector.")
        
    def updateButtons(self):
        pass
    
                
class indexWizard(QWidget):
    """The indexing wizard is separated so that it can be embedded both as a separate button for the crystallography part or when you add a new crystal"""
    def __init__(self, structure):
        super().__init__()
        
        #the measurements are stored as variables inside and are updated only when calculate is pressed. This is to prevent
        #that a calculation happens, after which angles and reflections are changed. When values are changed, the results box should be cleared
        self.l1 = 0
        self.theta1 = 0
        self.l2 = 0
        self.theta2 = 0
        
        
        ##make a grouping box
        self.formGroupBox = QGroupBox("Indexing Wizard")
        layout = QVBoxLayout()
        
        self.structure = structure
        
        ##the top part is a grid
        self.entertoolbox = QGroupBox("Parameters")
        self.enterTools = QGridLayout()
        self.enterTools.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        self.enterTools.setSpacing(10)
        self.entertoolbox.setLayout(self.enterTools)
        
        ##Display which structure is active
        self.enterTools.addWidget(QLabel("Structure"), 1, 0)
        strucbox = QLineEdit()
        strucbox.setText(structure.name)
        strucbox.setReadOnly(True)
        strucbox.setSizePolicy(sp)
        self.enterTools.addWidget(strucbox, 1, 1)
        
        #labels
        lab1 = QLabel("Length")
        lab1.setAlignment(Qt.AlignCenter)
        lab1.setSizePolicy(sp)
        
        lab2 = QLabel(u"\u03b8")
        lab2.setAlignment(Qt.AlignCenter)
        lab2.setSizePolicy(sp)
        
        lab3 = QLabel("Reflection 1")
        lab3.setSizePolicy(sp)
        
        lab4 = QLabel("Reflection 2")
        lab4.setSizePolicy(sp)
        
        lab5 = QLabel("Error")
        lab5.setSizePolicy(sp)
        
        lab6 = QLabel("Max. HKL index")
        lab6.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        lab6.setSizePolicy(sp)
        
        self.enterTools.addWidget(lab1,2, 1) 
        self.enterTools.addWidget(lab2,2, 2)
        self.enterTools.addWidget(lab3, 3, 0)
        self.enterTools.addWidget(lab4, 4, 0)
        self.enterTools.addWidget(lab5, 5, 0)
        self.enterTools.addWidget(lab6, 5, 3)
        
        #entering boxes
        dec = 2
        
        #the order of adding the widgets to entertools is important, it determines tab-ability
        self.l1box = QDoubleSpinBox()
        self.l1box.setValue(1)
        self.l1box.setDecimals(dec)
        self.l1box.setMaximum(100)
        self.l1box.setMinimum(0.01)
        self.l1box.setSingleStep(0.01)
        self.l1box.setSuffix(" 1/nm")
        self.l1box.setSizePolicy(sp)
        self.l1box.valueChanged.connect(self.anyValueChanged)
        self.enterTools.addWidget(self.l1box, 3, 1)
        
        self.theta1box = QDoubleSpinBox()
        self.theta1box.setValue(0)
        self.theta1box.setDecimals(dec)
        self.theta1box.setMaximum(180)
        self.theta1box.setMinimum(-180)
        self.theta1box.setSingleStep(0.01)
        self.theta1box.setSuffix(u" \u00b0")
        self.theta1box.setSizePolicy(sp)
        self.theta1box.valueChanged.connect(self.anyValueChanged)
        self.enterTools.addWidget(self.theta1box, 3, 2)
        
        
        self.l2box = QDoubleSpinBox()
        self.l2box.setValue(1)
        self.l2box.setDecimals(dec)
        self.l2box.setMaximum(100)
        self.l2box.setMinimum(0.01)
        self.l2box.setSingleStep(0.01)
        self.l2box.setSuffix(" 1/nm")
        self.l2box.setSizePolicy(sp)
        self.l2box.valueChanged.connect(self.anyValueChanged)
        self.enterTools.addWidget(self.l2box, 4, 1)
        
        self.theta2box = QDoubleSpinBox()
        self.theta2box.setValue(0)
        self.theta2box.setDecimals(dec)
        self.theta2box.setMaximum(180)
        self.theta2box.setMinimum(-180)
        self.theta2box.setSingleStep(0.01)
        self.theta2box.setSuffix(u" \u00b0")
        self.theta2box.setSizePolicy(sp)
        self.theta2box.valueChanged.connect(self.anyValueChanged)
        self.enterTools.addWidget(self.theta2box, 4, 2)
        
        ##Error boxxes
        self.errlbox = QDoubleSpinBox()
        self.errlbox.setValue(0.3)
        self.errlbox.setDecimals(dec)
        self.errlbox.setMaximum(1)
        self.errlbox.setMinimum(0.01)
        self.errlbox.setSingleStep(0.01)
        self.errlbox.setPrefix(u"\u00b1 ")
        self.errlbox.setSuffix(" 1/nm")
        self.errlbox.setSizePolicy(sp)
        self.errlbox.valueChanged.connect(self.anyValueChanged)
        self.enterTools.addWidget(self.errlbox, 5, 1)
        
        self.errthetabox = QDoubleSpinBox()
        self.errthetabox.setValue(2)
        self.errthetabox.setDecimals(dec)
        self.errthetabox.setMaximum(10)
        self.errthetabox.setMinimum(0.01)
        self.errthetabox.setSingleStep(0.01)
        self.errthetabox.setPrefix(u"\u00b1 ")
        self.errthetabox.setSuffix(u" \u00b0")
        self.errthetabox.setSizePolicy(sp)
        self.errthetabox.valueChanged.connect(self.anyValueChanged)
        self.enterTools.addWidget(self.errthetabox, 5, 2)
        
        
        #Calculate options and button
        self.indexbox = QSpinBox()
        self.indexbox.setValue(5)
        self.indexbox.setMaximum(10)
        self.indexbox.setMinimum(1)
        self.indexbox.setSizePolicy(sp)
        self.indexbox.valueChanged.connect(self.anyValueChanged)
        self.enterTools.addWidget(self.indexbox, 5, 4)
        
        self.calcbutton = QPushButton("Calculate", self)
        self.calcbutton.setToolTip("Find the possible indexations")
        self.calcbutton.clicked.connect(self.calculate)
        self.calcbutton.setSizePolicy(sp)
        #self.calcbutton.setSizePolicy(sp)
        self.enterTools.addWidget(self.calcbutton, 6, 1)
        
        lab7 = QLabel("Result")
        lab7.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.enterTools.addWidget(lab7, 6, 2)
        self.resultbox = QLineEdit()
        self.resultbox.setReadOnly(True)
        self.enterTools.addWidget(self.resultbox, 6, 3, 1, 2)
        
        self.helpbutton = QPushButton("", self)
        self.helpbutton.clicked.connect(self.showHelp)
        self.helpbutton.setToolTip("What is this?")
        self.helpbutton.resize(24+6, 24+6)
        self.helpbutton.setIconSize(QSize(24,24))
        self.helpbutton.setIcon(QIcon(".\Images\question.png"))
        self.helpbutton.setSizePolicy(sp)
        self.enterTools.addWidget(self.helpbutton, 1, 4)
        
        #we need a spacer to fill the hole on the right #actually no we don't, leftaligning of the grid will do the trick.
        #self.enterTools.addItem(QSpacerItem(100,1),1,4)
        
        
        ##the bottom part is a scrollable table list that is filled up when calculate is clicked
        self.data = []
        
        self.opts = QTableWidget(0, 6)
        self.colabels=["(hkl)\u2081", "Length (1/nm)", "(hkl)\u2082", "Length (1/nm)", "Angle (\u00b0)", "[uvw]"]
        self.opts.setHorizontalHeaderLabels(self.colabels)
        self.opts.resizeColumnsToContents()
        self.opts.resizeRowsToContents()
        self.opts.setFixedHeight(200)
        self.opts.setFixedWidth(self.opts.width()) #.horizontalHeader()
        self.opts.cellDoubleClicked[int, int].connect(self.passData)
        self.opts.cellClicked[int, int].connect(self.selecData)
        self.opts.setEditTriggers(QAbstractItemView.NoEditTriggers)
        #self.opts.clear()
        
        #self.opts.setColumnWidth(0, 50)
        #self.opts.setColumnWidth(1, 50)
        #self.opts.setColumnWidth(2, 50)
        #self.opts.setColumnWidth(3, 50)
        #self.opts.setColumnWidth(4, 50)
        #self.opts.setColumnWidth(5, 50)
        
        #self.scroll = QScrollArea()
        #scroll.setFixedHeight(200)
        layout.addWidget(self.entertoolbox)
        layout.addWidget(self.opts)
        
        self.formGroupBox.setLayout(layout)
        
        mainlayout = QHBoxLayout()
        mainlayout.addWidget(self.formGroupBox)
        self.setLayout(mainlayout)
    
    def anyValueChanged(self):
        self.resultbox.setText("")
    
    def updateValues(self):
        self.l1 = self.l1box.value()
        self.theta1 = self.theta1box.value()
        self.l2 = self.l2box.value()
        self.theta2 = self.theta2box.value()
    
    def getAngles(self):
        return self.theta1, self.theta2
    
    def passData(self):
        pass #haven't quite figured out yet how to emit and then catch this
        
    def selecData(self, r, c):
        self.opts.selectRow(r)
    
    def getActiveRow(self):
        if len(self.data)!=0: #there must be something returned, in the beginning data is completely empty
            if len(self.data[0])!=0: #there must be data when calculate has been clicked
                rw = self.opts.currentRow()
                if rw>=0: #it is -1 if nothing is selected
                    ref1, l1, ref2, l2, ang, za = self.data
                    return ref1[rw, :], l1[rw], ref2[rw, :], l2[rw], ang[rw], za[rw, :]
                else:
                    return None
            else:
                return None
        else:
            return None
    
    def showHelp(self):
        self.testIndexation()
    
    def testIndexation(self):
        self.l1box.setValue(5.6)
        self.l2box.setValue(4.9)
        self.theta1box.setValue(-22)
        self.theta2box.setValue(33)
        self.calculate()
    
    def calculate(self):
        vals = self.getEnteredValues()
        self.updateValues() #only once calculate is pressed are the internal values changed which are returned when requested
        self.data = self.structure.getZoneRefPairs(vals[0], vals[1], vals[2], vals[3], maxind = vals[4], err = vals[5], err2 = vals[6], verbose = False)
        
        ref1, l1, ref2, l2, ang, za = self.data
        self.data = list(self.data)
        self.data[0] = np.array(ref1)
        self.data[1] = np.array(l1)
        self.data[2] = np.array(ref2)
        self.data[3] = np.array(l2)
        self.data[4] = np.array(ang)
        self.data[5] = np.array(za)
        
        #update the visible table as well
        self.opts.clear() #first clear
        #change the size of the table
        self.opts.setRowCount(len(self.data[0]))
        self.opts.setHorizontalHeaderLabels(self.colabels)
        
        #check if the data is empty:
        if len(self.data[0])==0:
            self.resultbox.setText("No matches found!")
        else:
            self.resultbox.setText("Matches found!")
        
        tabfont = QFont()
        tabfont.setPointSize(8)
        tabfont.setFamily("Arial")
        
        #insert the data
        for i in range(len(self.data[0])):
            #reflection 1 - change square brackets to round ones
            n1 = QTableWidgetItem(str(self.data[0][i]).replace("]", ")").replace("[", "("))
            n1.setFont(tabfont)
            #length 1 - transform back to 1/nm!
            n2 = QTableWidgetItem(str(round(float(self.data[1][i])*10, 2)))
            n2.setFont(tabfont)
            #reflection 2 - change square brackets to round ones
            n3 = QTableWidgetItem(str(self.data[2][i]).replace("]", ")").replace("[", "("))
            n3.setFont(tabfont)
            #length 2 - transform back to 1/nm!
            n4 = QTableWidgetItem(str(round(float(self.data[3][i])*10, 2)))
            n4.setFont(tabfont)
            #angle
            n5 = QTableWidgetItem(str(round(float(self.data[4][i]), 2)))
            n5.setFont(tabfont)
            #zone axis
            n6 = QTableWidgetItem(str(self.data[5][i]))
            n6.setFont(tabfont)
            self.opts.setItem(i, 0, n1)
            self.opts.setItem(i, 1, n2)
            self.opts.setItem(i, 2, n3)
            self.opts.setItem(i, 3, n4)
            self.opts.setItem(i, 4, n5)
            self.opts.setItem(i, 5, n6)    
        
        self.opts.resizeColumnsToContents()
        self.opts.resizeRowsToContents()
        
    def getEnteredValues(self):
        l1 = self.l1box.value()
        l2 = self.l2box.value()
        theta1 = self.theta1box.value()
        theta2 = self.theta2box.value()
        er1 = self.errlbox.value()
        er2 = self.errthetabox.value()
        maxind = self.indexbox.value()
        return (l1, theta1, l2, theta2, maxind, er1, er2)
        
    def getTable(self):
        return self.data
    
    def getTableEntry(self):
        pass
        
class Dialog(QDialog):
    
    def __init__(self):
        super(Dialog, self).__init__()
        
        self.items = []
        
        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Ok | QDialogButtonBox.Cancel)
        self.buttonBox.accepted.connect(self.accept)
        self.buttonBox.rejected.connect(self.reject)
        
        self.setWindowIcon(QIcon('./Images/logo.png'))
        
    def getres(self):
        """Get results method. All the relevant input objects to a dialog are added to the self.item variable. This is looped over and the relevant value is put in a list. This list is returned."""
        values = []
        for i in self.items:
            if isinstance(i, QDoubleSpinBox) or isinstance(i, QSpinBox):
                values.append(i.value())
            elif isinstance(i, QLineEdit):
                values.append(i.text())
            elif isinstance(i, QCheckBox):
                values.append(i.isChecked())
            elif isinstance(i, QComboBox):
                values.append(i.currentText())
            else:
                values.append(None)
        return values
    
        
class microscopeDialog(Dialog):
 
    def __init__(self, windowtitle = "New microscope", **kwargs):
        super(microscopeDialog, self).__init__()
        
        self.createFormGroupBox(**kwargs)
 
        self.setWindowTitle(windowtitle)
        
        self.makeLayout()
 
    def createFormGroupBox(self, name = "", voltage = 300):
        self.formGroupBox = QGroupBox("Parameters")
        layout = QFormLayout()
        
        #name box
        nmbox = QLineEdit()
        nmbox.setMaxLength(20)
        nmbox.setText(name)
        layout.addRow(QLabel("Name"), nmbox)
        self.items.append(nmbox)
        
        #kv box
        kvbox = QSpinBox()
        kvbox.setRange(1, 10000)
        kvbox.setSingleStep(10)
        kvbox.setValue(voltage)
        kvbox.setSuffix(" kV")
        layout.addRow(QLabel("Voltage"), kvbox)
        self.items.append(kvbox)
        self.formGroupBox.setLayout(layout)
    
    def makeLayout(self):
        mainLayout = QVBoxLayout()
        mainLayout.addWidget(self.formGroupBox)
        mainLayout.addWidget(self.buttonBox)
        mainLayout.setSizeConstraint(mainLayout.SetFixedSize)
        self.setLayout(mainLayout)
    
    def getres(self):
        return super(microscopeDialog, self).getres()
    
    # static method to create dialog and return all inputs
    @staticmethod
    def getInfo(**kwargs):
        dialog = microscopeDialog(**kwargs)
        result = dialog.exec_()
        res = dialog.getres()
        return (res, result == QDialog.Accepted)

class structureDialog(Dialog):
 
    def __init__(self, windowtitle = "New structure", **kwargs):
        super(structureDialog, self).__init__()
        
        self.createFormGroupBox(**kwargs)
 
        self.setWindowTitle(windowtitle)
        
        self.makeLayout()
 
    def createFormGroupBox(self, name = "", a = 1, b = 1, c = 1, alpha = 90, beta = 90, gamma = 90):
        self.formGroupBox = QGroupBox("Parameters")
        layout = QFormLayout()
        
        #name box
        nmbox = QLineEdit()
        nmbox.setMaxLength(20)
        nmbox.setText(name)
        layout.addRow(QLabel("Name"), nmbox)
        self.items.append(nmbox)
        
        dec = 2
        #The crystal structure parameter boxes in angstrom and degrees
        self.abox = QDoubleSpinBox()
        self.abox.setRange(0.01, 10000.0)
        self.abox.setSingleStep(0.01)
        self.abox.setSuffix(u" \u212B")
        self.abox.valueChanged.connect(self.whatcrystal)
        self.abox.setDecimals(dec)
        self.abox.setValue(a)
        layout.addRow(QLabel("a"), self.abox)
        self.items.append(self.abox)
        
        self.bbox = QDoubleSpinBox()
        self.bbox.setRange(0.01, 10000.0)
        self.bbox.setSingleStep(0.01)
        self.bbox.setSuffix(u" \u212B")
        self.bbox.valueChanged.connect(self.whatcrystal)
        self.bbox.setValue(b)
        self.bbox.setDecimals(dec)
        layout.addRow(QLabel("b"), self.bbox)
        self.items.append(self.bbox)
        
        self.cbox = QDoubleSpinBox()
        self.cbox.setRange(0.01, 10000.0)
        self.cbox.setSingleStep(0.01)
        self.cbox.setSuffix(u" \u212B")
        self.cbox.setValue(c)
        self.cbox.setDecimals(dec)
        self.cbox.valueChanged.connect(self.whatcrystal)
        layout.addRow(QLabel("c"), self.cbox)
        self.items.append(self.cbox)
        
        self.alphabox = QDoubleSpinBox()
        self.alphabox.setRange(0.01, 179.99)
        self.alphabox.setSingleStep(0.01)
        self.alphabox.setValue(alpha)
        self.alphabox.setSuffix(u" \u00b0")
        self.alphabox.setDecimals(dec)
        self.alphabox.valueChanged.connect(self.whatcrystal)
        layout.addRow(QLabel(u"\u03b1"), self.alphabox)
        self.items.append(self.alphabox)
        
        self.betabox = QDoubleSpinBox()
        self.betabox.setRange(0.01, 179.99)
        self.betabox.setSingleStep(0.01)
        self.betabox.setValue(beta)
        self.betabox.setSuffix(u" \u00b0")
        self.betabox.setDecimals(dec)
        self.betabox.valueChanged.connect(self.whatcrystal)
        layout.addRow(QLabel(u"\u03b2"), self.betabox)
        self.items.append(self.betabox)
        
        self.gammabox = QDoubleSpinBox()
        self.gammabox.setRange(0.01, 179.99)
        self.gammabox.setSingleStep(0.01)
        self.gammabox.setValue(gamma)
        self.gammabox.setDecimals(dec)
        self.gammabox.setSuffix(u" \u00b0")
        self.gammabox.valueChanged.connect(self.whatcrystal)
        layout.addRow(QLabel(u"\u03b3"), self.gammabox)
        self.items.append(self.gammabox)
        
        #textbox that says the crystal system
        self.typbox = QLineEdit()
        self.typbox.setReadOnly(True)
        layout.addRow(QLabel("Crystal class"), self.typbox)
        self.whatcrystal()
        
        self.formGroupBox.setLayout(layout)
    
    def whatcrystal(self):
        try:
            typ = tc.Structure.getCrystalClass(self.abox.value(), self.bbox.value(), self.cbox.value(), self.alphabox.value(), self.betabox.value(), self.gammabox.value())
            self.typbox.setText(typ)
        except:
            pass
        
    def makeLayout(self):
        mainLayout = QVBoxLayout()
        mainLayout.addWidget(self.formGroupBox)
        mainLayout.addWidget(self.buttonBox)
        mainLayout.setSizeConstraint(mainLayout.SetFixedSize)
        self.setLayout(mainLayout)
    
    def getres(self):
        return super(structureDialog, self).getres()
    
    # static method to create dialog and return all inputs
    @staticmethod
    def getInfo(**kwargs):
        dialog = structureDialog(**kwargs)
        result = dialog.exec_()
        res = dialog.getres()
        return (res, result == QDialog.Accepted)        
        
class stageDialog(Dialog):
 
    def __init__(self, windowtitle = "New stage", **kwargs):
        super(stageDialog, self).__init__()
        
        self.createFormGroupBox(**kwargs)
 
        self.setWindowTitle(windowtitle)
        
        self.makeLayout()
 
    def createFormGroupBox(self, name = "", alpha = 0.0, alphamin = -30.0, alphamax= 30.0, beta = 0.0, betamin = -20.0, betamax = 20.0, alpharev = False, betarev = False):
        
        self.formGroupBox = QGroupBox("Parameters")
        layout = QGridLayout()
        layout.setSpacing(10)
        
        #labels
        layout.addWidget(QLabel("Name"), 1, 0)
        layout.addWidget(QLabel(u"\u03b1"), 3, 0)
        layout.addWidget(QLabel(u"\u03b2"), 4, 0)
        layout.addWidget(QLabel("Value"), 2, 1)
        layout.addWidget(QLabel("Min."), 2, 2)
        layout.addWidget(QLabel("Max."), 2, 3)
        #layout.addWidget(QLabel("Reversed?"), 2, 4)
        
        #self.formGroupBox = QGroupBox("Boundaries")
        
        #name box
        nmbox = QLineEdit()
        nmbox.setMaxLength(20)
        nmbox.setText(name)
        layout.addWidget(nmbox, 1, 1, 1, 3)
        self.items.append(nmbox)
        
        #help button
        qbutton = QPushButton("", self)
        qbutton.setToolTip("What is this?")
        qbutton.setIcon(QIcon(".\Images\question.png"))
        #qbutton.clicked.connect(action)
        qbutton.resize(24+6, 24+6)
        qbutton.setIconSize(QSize(24,24))
        qbutton.setSizePolicy(sp)
        layout.addWidget(qbutton, 1, 4)
        
        #angle boxes
        dec = 2
        
        #alphamin
        self.alphaminbox = QDoubleSpinBox()
        self.alphaminbox.setDecimals(dec)
        self.alphaminbox.setSingleStep(10**(-dec))
        self.alphaminbox.setMaximum(0)
        self.alphaminbox.setMinimum(-90+10**(-dec))
        self.alphaminbox.setSuffix(u" \u00b0")
        self.alphaminbox.setValue(alphamin)
        self.alphaminbox.valueChanged.connect(self.anglemaxupdate)
        layout.addWidget(self.alphaminbox, 3, 2)
        self.items.append(self.alphaminbox)
        
        #alphamax
        self.alphamaxbox = QDoubleSpinBox()
        self.alphamaxbox.setDecimals(dec)
        self.alphamaxbox.setSingleStep(10**(-dec))
        self.alphamaxbox.setMaximum(90-10**(-dec))
        self.alphamaxbox.setMinimum(0)
        self.alphamaxbox.setSuffix(u" \u00b0")
        self.alphamaxbox.setValue(alphamax)
        self.alphamaxbox.valueChanged.connect(self.anglemaxupdate)
        layout.addWidget(self.alphamaxbox, 3, 3)
        self.items.append(self.alphamaxbox)
        
        #Alpha
        self.alphabox = QDoubleSpinBox()
        self.alphabox.setDecimals(dec)
        self.alphabox.setSingleStep(10**(-dec))
        self.alphabox.setMaximum(self.alphamaxbox.value())
        self.alphabox.setMinimum(self.alphaminbox.value())
        self.alphabox.setSuffix(u" \u00b0")
        self.alphabox.setValue(alpha)
        layout.addWidget(self.alphabox, 3, 1)
        self.items.append(self.alphabox)
        
        #betamin
        self.betaminbox = QDoubleSpinBox()
        self.betaminbox.setDecimals(dec)
        self.betaminbox.setSingleStep(10**(-dec))
        self.betaminbox.setMaximum(-10**(-dec))
        self.betaminbox.setMinimum(-90+10**(-dec))
        self.betaminbox.setSuffix(u" \u00b0")
        self.betaminbox.setValue(betamin)
        self.betaminbox.valueChanged.connect(self.anglemaxupdate)
        layout.addWidget(self.betaminbox, 4, 2)
        self.items.append(self.betaminbox)
        
        #betamax
        self.betamaxbox = QDoubleSpinBox()
        self.betamaxbox.setDecimals(dec)
        self.betamaxbox.setSingleStep(10**(-dec))
        self.betamaxbox.setMaximum(90-10**(-dec))
        self.betamaxbox.setMinimum(10**(-dec))
        self.betamaxbox.setSuffix(u" \u00b0")
        self.betamaxbox.setValue(betamax)
        self.betamaxbox.valueChanged.connect(self.anglemaxupdate)
        layout.addWidget(self.betamaxbox, 4, 3)
        self.items.append(self.betamaxbox)
        
        #beta
        self.betabox = QDoubleSpinBox()
        self.betabox.setDecimals(dec)
        self.betabox.setSingleStep(10**(-dec))
        self.betabox.setMaximum(self.betamaxbox.value())
        self.betabox.setMinimum(self.betaminbox.value())
        self.betabox.setSuffix(u" \u00b0")
        self.betabox.setValue(beta)
        layout.addWidget(self.betabox, 4, 1)
        self.items.append(self.betabox)
        
        #reversed checkboxes
        acheck = QCheckBox("Reversed?")
        acheck.setChecked(alpharev)
        layout.addWidget(acheck, 3, 4)
        self.items.append(acheck)
        
        bcheck = QCheckBox("Reversed?")
        bcheck.setChecked(betarev)
        layout.addWidget(bcheck, 4, 4)
        self.items.append(bcheck)
        
        self.formGroupBox.setLayout(layout)
    
    def anglemaxupdate(self):
        self.alphabox.setMaximum(self.alphamaxbox.value())
        self.alphabox.setMinimum(self.alphaminbox.value())
        self.betabox.setMaximum(self.betamaxbox.value())
        self.betabox.setMinimum(self.betaminbox.value())
        
    
    def makeLayout(self):
        mainLayout = QVBoxLayout()
        mainLayout.addWidget(self.formGroupBox)
        mainLayout.addWidget(self.buttonBox)
        mainLayout.setSizeConstraint(mainLayout.SetFixedSize)
        self.setLayout(mainLayout)
    
    def getres(self):
        return super(stageDialog, self).getres()
    
    # static method to create dialog and return all inputs
    @staticmethod
    def getInfo(**kwargs):
        dialog = stageDialog(**kwargs)
        result = dialog.exec_()
        res = dialog.getres()
        return (res, result == QDialog.Accepted)        

class crystalDialog(Dialog):
    
    def __init__(self, caller, windowtitle = "New crystal",  **kwargs):
        super(crystalDialog, self).__init__()
        
        #caller is the mainwindow app. We need acces to the current settings.
        self.caller = caller
        
        self.createFormGroupBox(**kwargs)
 
        self.setWindowTitle(windowtitle)
        
        
        self.makeLayout()
    
    def createFormGroupBox(self, name = "", u = 0, v = 0, w = 0, h = 0, k = 0, l = 0, thet = 0, comt = "", strucname = None):
        #left pane layout
        self.lp = QGroupBox("")
        lpLayout = QVBoxLayout()
        
        lpLayout.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        lpLayout.setSizeConstraint(lpLayout.SetFixedSize)
        
        ###########################################################
        #a box with review information on the status of the microscope
        self.bx1 = QGroupBox("Review")
        bx1layout = QGridLayout()
        bx1layout.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        
        bx1layout.addWidget(QLabel("Structure:"), 1, 0)
        if strucname is None: #structure can be provided as argument (for already existing crystals). If there is none just take the currently active structure.
            self.struc = QLabel(self.caller.getCurrentStructure().name)
        else:
            self.struc = QLabel(strucname)
        bx1layout.addWidget(self.struc, 1, 1)
        bx1layout.addWidget(QLabel("Stage:"), 2, 0)
        bx1layout.addWidget(QLabel(self.caller.getCurrentStage().name), 2, 1)
        bx1layout.addWidget(QLabel("\u03b1 ="), 2, 2)
        bx1layout.addWidget(QLabel(str(round(self.caller.getAlpha(), 2))+" \u00b0"), 2, 3)
        bx1layout.addWidget(QLabel("\u03b2 ="), 2, 4)
        bx1layout.addWidget(QLabel(str(round(self.caller.getBeta(), 2))+" \u00b0"), 2, 5)
        bx1layout.addWidget(QLabel("Detector:"), 3, 0)
        bx1layout.addWidget(QLabel(self.caller.getCurrentDetector().name), 3, 1)
        bx1layout.addWidget(QLabel("Mode:"), 3, 2)
        bx1layout.addWidget(QLabel(self.caller.getMode()), 3, 3)
        bx1layout.addWidget(QLabel("Mag/CL:"), 3, 4)
        bx1layout.addWidget(QLabel(self.caller.getSetting()), 3, 5)
        #bx1layout.setSizeConstraint(bx1layout.SetFixedSize)
        self.bx1.setLayout(bx1layout)
        
        ############################################
        #a group box with new to enter information
        self.bx2 = QGroupBox("Parameters")
        bx2layout = QGridLayout()
        bx2layout.setAlignment(Qt.AlignLeft | Qt.AlignTop)
        
        #all the labels in this box
        bx2layout.addWidget(QLabel("Name"), 1, 0)
        bx2layout.addWidget(QLabel("Zone axis"), 3, 0)
        bx2layout.addWidget(QLabel("Reflection"), 5, 0)
        bx2layout.addWidget(QLabel("u"), 2, 1)
        bx2layout.addWidget(QLabel("v"), 2, 2)
        bx2layout.addWidget(QLabel("w"), 2, 3)
        bx2layout.addWidget(QLabel("h"), 4, 1)
        bx2layout.addWidget(QLabel("k"), 4, 2)
        bx2layout.addWidget(QLabel("l"), 4, 3)
        bx2layout.addWidget(QLabel("\u03b8"), 4, 4) 
        bx2layout.addWidget(QLabel("Comment:"), 6, 0)
        #bx2layout.setSizeConstraint(bx2layout.SetFixedSize)
        
        #All the entering boxes
        self.namebox = QLineEdit()
        self.namebox.setMaxLength(20)
        self.namebox.setText(name)
        self.namebox.setSizePolicy(sp)
        bx2layout.addWidget(self.namebox, 1, 1, 1, 3)
        self.items.append(self.namebox)
        
        mxindx = 20
        
        self.ubox = QSpinBox()
        self.ubox.setMaximum(mxindx)
        self.ubox.setMinimum(-mxindx)
        self.ubox.setValue(u)
        self.ubox.setSizePolicy(sp)
        bx2layout.addWidget(self.ubox, 3, 1)
        self.items.append(self.ubox)
        
        self.vbox = QSpinBox()
        self.vbox.setMaximum(mxindx)
        self.vbox.setMinimum(-mxindx)
        self.vbox.setValue(v)
        self.vbox.setSizePolicy(sp)
        bx2layout.addWidget(self.vbox, 3, 2)
        self.items.append(self.vbox)
        
        self.wbox = QSpinBox()
        self.wbox.setMaximum(mxindx)
        self.wbox.setMinimum(-mxindx)
        self.wbox.setValue(w)
        self.wbox.setSizePolicy(sp)
        bx2layout.addWidget(self.wbox, 3, 3)
        self.items.append(self.wbox)
        
        self.hbox = QSpinBox()
        self.hbox.setMaximum(mxindx)
        self.hbox.setMinimum(-mxindx)
        self.hbox.setValue(h)
        self.hbox.setSizePolicy(sp)
        bx2layout.addWidget(self.hbox, 5, 1)
        self.items.append(self.hbox)
        
        self.kbox = QSpinBox()
        self.kbox.setMaximum(mxindx)
        self.kbox.setMinimum(-mxindx)
        self.kbox.setValue(k)
        self.kbox.setSizePolicy(sp)
        bx2layout.addWidget(self.kbox, 5, 2)
        self.items.append(self.kbox)
        
        self.lbox = QSpinBox()
        self.lbox.setMaximum(mxindx)
        self.lbox.setMinimum(-mxindx)
        self.lbox.setValue(l)
        self.lbox.setSizePolicy(sp)
        bx2layout.addWidget(self.lbox, 5, 3)
        self.items.append(self.lbox)
        
        dec = 2
        
        self.thetabox = QDoubleSpinBox()
        self.thetabox.setDecimals(dec)
        self.thetabox.setSingleStep(10**(-dec))
        self.thetabox.setMaximum(180)
        self.thetabox.setMinimum(-180)
        self.thetabox.setSuffix(u" \u00b0")
        self.thetabox.setValue(thet)
        self.thetabox.setSizePolicy(sp)
        bx2layout.addWidget(self.thetabox, 5, 4)
        self.items.append(self.thetabox)
        
        self.commentbox = QLineEdit()
        self.commentbox.setText(comt)
        self.commentbox.setSizePolicy(sp)
        bx2layout.addWidget(self.commentbox, 6, 1, 1, 4)
        self.items.append(self.commentbox)
        
        #the indexing wizard button
        self.wizardbutton = QPushButton("", self)
        self.wizardbutton.clicked.connect(self.showWizard)
        self.wizardbutton.setToolTip("Show/Hide indexing wizard")
        self.wizardbutton.resize(24+6, 24+6)
        self.wizardbutton.setIconSize(QSize(24,24))
        self.wizardbutton.setIcon(QIcon(".\Images\wizard.png"))
        self.wizardbutton.setSizePolicy(sp)
        #( | Qt.AlignTop)
        bx2layout.addWidget(self.wizardbutton, 1, 4, Qt.AlignRight)
        
        
        self.bx2.setLayout(bx2layout)
        
        
        #####################################################
        #a box with testing
        self.bx3 = QGroupBox("Testing")
        bx3layout = QVBoxLayout()
        #bx3layout.setSizeConstraint(bx3layout.SetFixedSize)
        self.bx3.setLayout(bx3layout)
        
        self.testbut = QPushButton("Calculate a test")
        self.testbut.clicked.connect(self.calcTest)
        self.testbut.setToolTip("Calculate closeby tilt condition to test indexation.")
        self.testbut.setSizePolicy(sp)
        bx3layout.addWidget(self.testbut, 0, Qt.AlignCenter)
        
        self.testtextbox = QTextEdit()
        self.testtextbox.setReadOnly(True)
        bx3layout.addWidget(self.testtextbox)
        
        ################################################
        #Add all the boxes on the left pane
        lpLayout.addWidget(self.bx1)
        lpLayout.addWidget(self.bx2)
        lpLayout.addWidget(self.bx3)
        
        #also add the buttonbox (ok and cancel)!!!
        lpLayout.addWidget(self.buttonBox)
        
        self.lp.setLayout(lpLayout)
        
        #####Right pane
        self.rp = QGroupBox("")
        rpLayout = QVBoxLayout()
        #a box with the indexation wizard
        self.indxwiz1 = indexWizard(self.caller.getCurrentStructure())
        rpLayout.addWidget(self.indxwiz1)
        #self.indxwiz1.hide()
        
        #a box with the submit button
        self.sbmt = QPushButton("Set selected indexing")
        self.sbmt.clicked.connect(self.selectCalc)
        self.sbmt.setToolTip("Use selected indexing for the crystal")
        self.sbmt.setSizePolicy(sp)
        rpLayout.addWidget(self.sbmt, 0, Qt.AlignCenter)
        
        self.rp.setLayout(rpLayout)
        self.rp.hide()
    
    def showWizard(self):
        if self.rp.isVisible():
            self.rp.hide()
            self.resize(self.sizeHint())
        else:
            self.rp.show()
            self.resize(self.sizeHint())
   
    def selectCalc(self):
        data = self.indxwiz1.getActiveRow()
        angls = self.indxwiz1.getAngles()
        #set values of the fields - use the first reflection as reflection
        if data is not None:
            self.ubox.setValue(data[5][0])
            self.vbox.setValue(data[5][1])
            self.wbox.setValue(data[5][2])
            self.hbox.setValue(data[0][0])
            self.kbox.setValue(data[0][1])
            self.lbox.setValue(data[0][2])
            self.thetabox.setValue(angls[0])
        else:
            pass
        
    def calcTest(self):
        #first check that parameters are valid: no all zeros and dot(zone, reflection) is 0
        u = self.ubox.value()
        v = self.vbox.value()
        w = self.wbox.value()
        h = self.hbox.value()
        k = self.kbox.value()
        l = self.lbox.value()
        th = self.thetabox.value()
        #these conditions are invalid indices
        if (u*h + v*k + w*l != 0) or (u==0 and v==0 and w==0) or (h==0 and k==0 and l==0):
            self.testtextbox.setText("Those indices are not valid! The reflection must lie in the zone (hu+kv+lw=0) and not all indices can be 0.")
        #these conditions are valid but maybe nothing will be found. We create a temporary crystal.
        else:
            stage = self.caller.getCurrentStage()
            struc = self.caller.getCurrentStructure()
            det = self.caller.getCurrentDetector()
            mod = self.caller.getMode()
            set = float(self.caller.getSetting())
            #create a temporary crystal
            tempcrys = stage.addCrystal(struc.name)
            #Set the orientation of the crystal
            tempcrys.calcOrient(tc.miller(u, v, w), tc.miller(h, k, l), th, det.name, mod, set)
            #nearest zone
            zan = tempcrys.nearestZone(curzone = tc.miller(u, v, w))
            #if nothing is found
            if zan.size ==0:
                self.testtextbox.setText("No test could be computed, there appear to be no low index zones nearby.")
            else: #something is found
                za = str(zan)
                a, b = tempcrys.getAlphaBetaMiller(zan)
                self.testtextbox.setText("If this indexation is correct, you should find zone\n%s at \u03b1 = %s \u00b0 and \u03b2 = %s \u00b0.\n\nIf this indexation is not correct, try to invert the reflection (h = -h, k = -k, l = -l).\nIf this has been tried without success, the indexation may be wrong entirely, the settings of the microscope may be wrong and/or the calibration of the detector may be wrong." %(za, round(a,2), round(b, 2)))
            stage.removeCrystal(tempcrys.name)
    
    def makeLayout(self):
        mainLayout = QHBoxLayout()
        mainLayout.addWidget(self.lp)
        mainLayout.addWidget(self.rp)
        mainLayout.setSizeConstraint(mainLayout.SetFixedSize)
        self.setLayout(mainLayout)
    
    def getres(self):
        return super(crystalDialog, self).getres()
    
    # static method to create dialog and return all inputs
    @staticmethod
    def getInfo(caller, **kwargs):
        dialog = crystalDialog(caller, **kwargs)
        result = dialog.exec_()
        res = dialog.getres()
        return (res, result == QDialog.Accepted)   
        
class detectorDialog(Dialog):
 
    def __init__(self, windowtitle = "New detector", **kwargs):
        super(detectorDialog, self).__init__()
        
        self.createFormGroupBox(**kwargs)
 
        self.setWindowTitle(windowtitle)
        
        self.makeLayout()
 
    def createFormGroupBox(self, name = ""):
        
        self.formGroupBox = QGroupBox("Parameters")
        layout = QGridLayout()
        layout.setSpacing(10)
        
        #labels
        layout.addWidget(QLabel("Name"), 1, 0)
        
        #name box
        nmbox = QLineEdit()
        nmbox.setMaxLength(20)
        nmbox.setText(name)
        layout.addWidget(nmbox, 1, 1)
        self.items.append(nmbox)
        
        self.formGroupBox.setLayout(layout)
        
    
    def makeLayout(self):
        mainLayout = QVBoxLayout()
        mainLayout.addWidget(self.formGroupBox)
        mainLayout.addWidget(self.buttonBox)
        mainLayout.setSizeConstraint(mainLayout.SetFixedSize)
        self.setLayout(mainLayout)
    
    def getres(self):
        return super(detectorDialog, self).getres()
    
    # static method to create dialog and return all inputs
    @staticmethod
    def getInfo(**kwargs):
        dialog = detectorDialog(**kwargs)
        result = dialog.exec_()
        res = dialog.getres()
        return (res, result == QDialog.Accepted)       
        
class plotDialog(Dialog):
 
    def __init__(self, **kwargs):
        super(plotDialog, self).__init__()
        self.createFormGroupBox(**kwargs)
        
        self.buttonBox = QDialogButtonBox(QDialogButtonBox.Cancel)
        self.buttonBox.rejected.connect(self.reject)
        
        self.fig = 0
        self.ax = 0
        
        self.setWindowTitle("Stereographic plot")
        self.makeLayout()
    
    def saveplot(self):
        filename = QFileDialog.getSaveFileName(caption = "Save plot", filter = "PNG files (*png)")[0]
        plt.savefig(filename, dpi = 320)
    
    def createFormGroupBox(self, crys = None, dt = None, mod = None, set = None, typ="real", plotax = True, plotbond = True, hkllim = 1, lgnd = True):
        #dt is detector name
        vecincl = [] #need to connect this later to the maximum hkl index
        
        self.formGroupBox = QGroupBox("")
        layout = QVBoxLayout()
        
        self.fig, self.ax = crys.plotCrystalonDetector(dt, mod, set, vecincl = vecincl, typ=typ, plotaxes = plotax, plotbond = plotbond, hkm = hkllim, verbose=False)
        
        if lgnd: #set the title to the current things
            tit = "Crystal: %s; Space: %s\n\u03b1: %s\u00b0; \u03b2: %s\u00b0;\nDetector: %s; Mode: %s; Mag/CL: %s" %(crys.name, typ, crys.stage.getAlpha(), crys.stage.getBeta(), dt, mod, set)
            toleg = self.ax.plot([], [], color = "white")
            self.ax.legend(toleg, [tit], loc = 'best')
        
        #self.lim = 1
        #self.ax.set_xlim(-self.lim, self.lim)
        #self.ax.set_ylim(-self.lim, self.lim)
        
        self.canvas = FigureCanvas(self.fig)
        self.canvas.setSizePolicy(sp)
        self.canvas.draw()
        
        self.mpl_toolbar = NavigationToolbar(self.canvas, None)
        layout.addWidget(self.mpl_toolbar)
        layout.addWidget(self.canvas)
        self.formGroupBox.setLayout(layout)
        self.formGroupBox.setSizePolicy(sp)
        
    
    def makeLayout(self):
        mainLayout = QVBoxLayout()
        mainLayout.addWidget(self.formGroupBox)
        mainLayout.addWidget(self.buttonBox)
        mainLayout.setSizeConstraint(mainLayout.SetFixedSize)
        self.setLayout(mainLayout)
        
    # static method doesn't have to do much in this case
    @staticmethod
    def getInfo(**kwargs):
        dialog = plotDialog(**kwargs)
        result = dialog.exec_()
        
        
if __name__ == '__main__':
    # Create an PyQT5 application object.
    app = QApplication(sys.argv)
    # the widget is the main window which must be created
    ex = App()
    sys.exit(app.exec_())
    
 

