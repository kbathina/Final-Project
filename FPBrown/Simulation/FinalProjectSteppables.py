from PySteppables import *
import CompuCell
import sys
import random
import math as m
import numpy
import os

Init_max_div = 8.0
gl_initial_concentration = 15
gl_conversion_efficiency = .75

gl_death_threshold = .3
qui_death_rate = .2



OUTPUT_FILE_NAME = os.path.expanduser("~/Dropbox/CompuCell Stuff/Final Project/FPGreen/output.txt")
outputFile = None

def kill(steppable, cell, type, mcs):
    writeCell(steppable, type, mcs, cell)
    cell.targetVolume = 0
    cell.lambdaVolume = 100
    
def switchstate(steppable,cell, cellDict):
    if cellDict["Internal_Glucose_Storage"] >= cellDict["pro_threshold"]:
        cellDict["state"]  = "P"
        cellDict["gl_metabolic_rate"] = 1.2
        
    if cellDict["Internal_Glucose_Storage"] <= cellDict["qui_threshold"]:
        cellDict["gl_metabolic_rate"] = .7
        cell.type = 3  
        cellDict["state"]  = "Q"
        

def initOutputFile(fname = OUTPUT_FILE_NAME):
    
    global outputFile
    outputFile = open(fname, "w")
    
    s = "MCS"
    s += "; "
    s += "EVENT_ID"
    s += "; "
    s += "CELL_TYPE"
    s += "; "
    s += "CELL_ID"
    s += "; "
    s += "Parent_ID"
    s += "; "
    s += "Last_Stem_Cell_ID"
    s += "; "
    s += "P_STEM"
    s += "; "
    s += "MAX_DIVISIONS"
    s += "; "
    s += "CURRENT_DIVISIONS"
    s += "\n"
    
    outputFile.write(s)
    outputFile.flush()
    
    

def writeCell(steppable, eventId, mcs, cell):
    """ 
    Writes the state of a cell as a record in an output file.
    """
    global outputFile
    
    dict =  steppable.getDictionaryAttribute(cell)
   
    s = str(mcs)
    s += "; "
    s += str(eventId)
    s += "; "
    s += str(cell.type)
    s += "; "
    s += str(cell.id)
    s += "; "
    s += str(dict["Parent_ID"])
    s += "; "
    s += str(dict["Last_Stem_Cell_ID"])
    s += "; "
    s += str(dict["p_stem"])
    s += "; "
    s += str(dict["max_div"])
    s += "; "
    s += str(dict["cur_div"])
    s += "\n"
   
    outputFile.write(s)
    outputFile.flush()


from PySteppablesExamples import MitosisSteppableBase
            

class ConstraintInitializerSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def start(self):
        field=self.getConcentrationField("Glucose")
        field[:,:,:] = gl_initial_concentration

        initOutputFile()
        for cell in self.cellList:
            cellDict=self.getDictionaryAttribute(cell)
            cellDict["p_stem"]= .25
            cellDict["max_div"]= Init_max_div  
            cellDict["cur_div"]=0  
            cellDict["Parent_ID"] = cell.id
            cellDict["Last_Stem_Cell_ID"] = cell.id
            cellDict["state"] = None
            cellDict["Internal_Glucose_Storage"] = 8
            cellDict["qui_threshold"] = 2
            cellDict["pro_threshold"] = 4
            cellDict["gl_absorption_rate"] = .3
            cellDict["gl_growth_threshold"] = .25

            switchstate(self,cell,cellDict)            
         
            cell.targetVolume=25  
            cell.lambdaVolume=10.0 
            writeCell(self, "INIT", 0, cell)
        

class GrowthSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def step(self,mcs):
        
        GlucField=self.getConcentrationField("Glucose") 
        
        
        for cell in self.cellList:
            cellDict=self.getDictionaryAttribute(cell)
            if ((cell.type != 0) and (cell.type != 3)) \
            and (cell.targetVolume > 0):
                x = int(cell.xCOM)
                y = int(cell.yCOM)
                z = int(cell.zCOM)
                GlucValue = GlucField[x,y,z]
                GlucAbs = cellDict["gl_absorption_rate"] * GlucValue
                cellDict["Internal_Glucose_Storage"] += .1 * (GlucAbs - cellDict["gl_metabolic_rate"])
                #print cellDict["Internal_Glucose_Storage"] 
                
                switchstate(self,cell,cellDict)
                                
                GlucField[x,y,z]  = GlucValue - GlucAbs
                
                if (cell.targetVolume - cell.volume < 5)  and (cellDict["state"] == "P"):
                    cell.targetVolume += gl_conversion_efficiency * \
                    .9 * (GlucAbs - cellDict["gl_metabolic_rate"]) \
                    / (cellDict["gl_growth_threshold"] + .9 * (GlucAbs - cellDict["gl_metabolic_rate"]))
                    
                
            if cellDict["Internal_Glucose_Storage"] <= gl_death_threshold: 
                cell.type = 3
                writeCell(self, "Necrotic", mcs, cell)

            if cell.type == 3:
                cell.targetVolume -= qui_death_rate * cell.targetVolume
                if cell.targetVolume <= 1:
                    writeCell(self, "Nutrient-Death", mcs, cell)
                    cell.type = 0


                            
        

class MitosisSteppable(MitosisSteppableBase):
    def __init__(self,_simulator,_frequency=1):
        MitosisSteppableBase.__init__(self,_simulator, _frequency)
    
    def step(self,mcs):
        self.mcs = mcs
        cells_to_divide=[]
        for cell in self.cellList:
            cellDict=self.getDictionaryAttribute(cell)
            if (cell.volume>40) and ((cell.type == 2) or (cell.type == 1 and cellDict["state"]) == "P"):               
                cells_to_divide.append(cell)               
        for cell in cells_to_divide:
            self.divideCellRandomOrientation(cell)


    def updateAttributes(self):
        
        parentCell=self.mitosisSteppable.parentCell
        childCell=self.mitosisSteppable.childCell
        parentCell.targetVolume=25
        childCell.targetVolume=25
        childCell.lambdaVolume=parentCell.lambdaVolume
        
        pdict=self.getDictionaryAttribute(parentCell)
        cdict=self.getDictionaryAttribute(childCell)
        
        childCell.targetVolume=parentCell.targetVolume
        childCell.lambdaVolume=parentCell.lambdaVolume
        
        if parentCell.type == 1 or parentCell.type == 2:
        
            if random.random() < 0.01: 
                x = round(random.gauss(0.0,.5),2)
                if x >= 0:
                    pdict["max_div"] += x
                else:
                    pdict["max_div"] = pdict["max_div"] / (1 - x / Init_max_div)
                    
                writeCell(self, "MAXDIV-MUTATION", self.mcs, parentCell)
            
            if random.random() < 0.1:
                x = round(random.gauss(0.0,.05),2)
                if (pdict['p_stem'] < 1) and (pdict['p_stem'] > 0):
                    pdict["p_stem"] += x  
                    
                writeCell(self, "PSTEM-MUTATION", self.mcs, parentCell)
                
            if random.random() < 0.01:
                x = round(random.gauss(0.0,.05),2)
                pdict["qui_threshold"] += x  
                    
                writeCell(self, "Q-MUTATION", self.mcs, parentCell)
          
            if random.random() < 0.01:
                x = round(random.gauss(0.0,.05),2)
                pdict["pro_threshold"] += x  
                    
                writeCell(self, "P-MUTATION", self.mcs, parentCell)
                
            if random.random() < 0.01:
                x = round(random.gauss(0.0,.05),2)
                pdict["gl_absorption_rate"] += x  
                    
                writeCell(self, "ABSORPTION-MUTATION", self.mcs, parentCell)
          
            if random.random() < 0.01:
                x = round(random.gauss(0.0,.05),2)
                pdict["gl_growth_threshold"] += x  
                    
                writeCell(self, "GROWTH_THRESHOLD-MUTATION", self.mcs, parentCell)
         
        for key,value in pdict.items():
            cdict[key] = value
            
            
        if parentCell.type == 1:
            if random.random() < pdict['p_stem']:
                childCell.type=1
            else:
                cdict["cur_div"] += 1
                childCell.type = 2
                
        if parentCell.type == 2:
            pdict['cur_div'] += 1
            cdict['cur_div'] += 1
            childCell.type = 2
        
        cdict["Parent_ID"] = parentCell.id
        if parentCell.type == 1:
            cdict["Last_Stem_Cell_ID"] = parentCell.id   
        else:
            cdict["Last_Stem_Cell_ID"] = pdict["Last_Stem_Cell_ID"]  
            
        writeCell(self, "MITOSIS", self.mcs, childCell)    
        
        if parentCell.type == 2:  
            if (pdict['cur_div'] >= pdict['max_div']):
                kill(self, parentCell, "Limit-Death", self.mcs)
        
        

class DeathSteppable(SteppableBasePy):
    def __init__(self,_simulator,_frequency=1):
        SteppableBasePy.__init__(self,_simulator,_frequency)
    def step(self,mcs):
        pass
                

                        

        

class Graph(SteppableBasePy):
    def __init__(self,_simulator,_frequency=10):
        SteppableBasePy.__init__(self,_simulator,_frequency)
        
    def start(self):
        pass
        
    def step(self,mcs):
        pass

            
    def finish(self):
        pass
    
