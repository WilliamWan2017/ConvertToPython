# -*- coding: utf-8 -*-
"""
Created on Thu Aug 30 10:10:46 2018

@author: william
"""
from sympy import *
from sympy.parsing.sympy_parser import parse_expr 
#sympy.parsing.latex.parse_latex(s)[、、]
from sympy.abc import x,y,z,a,b,c,f,t,k,n
import json 
import Convert
from latex2sympy.process_latex import process_sympy
import FormatParseLatex
import sys 
def GetLcationName(location,ContinuourVariables,locations,Edge,Ext=[]):
    return location["boxName"]

def GetVar_ContsList(location,ContinuourVariables,locations,Edge,Ext=[]):
    return ','.join(ContinuourVariables)

def GetLocationInitList(location,ContinuourVariables,locations,Edge,Ext=[]):
    LocationNames=[location1["boxName"] for location1 in locations.values() if location1["boxName"] ]
    return '_FT,'.join(LocationNames)+"_FT"

def Getconvert_continuous_variable_to_simpy(location,ContinuourVariables,locations,Edge,Ext=[]):
    listResult=[]
    for variableName in ContinuourVariables:
        listResult.append('S.sympify(\''+variableName+'(t)\'): '+variableName+'')
        
    return ','.join(listResult)

def GetEdge_Guard(location,ContinuourVariables,locations,Edge,Ext=[]):
    if Ext[0]==1:
        strReturn='if ' +Edge["guard"]+":"
    else:
        strReturn='elif ' +Edge["guard"]+":"        
    return strReturn

def GetReset(location,ContinuourVariables,locations,Edge,Ext=[]):
    strReturn= Edge["reset"] 
    return strReturn

def GetLocationDestIndex(location,ContinuourVariables,locations,Edge,Ext=[]):
    for ik in locations.keys():
        if locations[ik]["boxName"]==Edge["strToLocation"]:
            return str(ik)
    return '0'

def GetCurrentlocationIndex(location,ContinuourVariables,locations,Edge,Ext=[]):
    for ik in locations.keys():
        if locations[ik]["boxName"]==location["boxName"]:
            return str(ik)
    return '0' 
def GetLocationDestName(location,ContinuourVariables,locations,Edge,Ext=[]):
    return Edge["strToLocation"]


def GetVar_Cont(location,ContinuourVariables,locations,Var_Cont,Ext=[]):
    return str(Var_Cont);
def GetVar_Cont_Format(location,ContinuourVariables,locations,Var_Cont,Ext=[]):
    
    listFormat=[]
    for v in ContinuourVariables:
        listFormat.append( '%7.4f')
    return ':'.join(listFormat);
def GetInvariant(location,ContinuourVariables,locations,Var_Cont,Ext=[]): 
    return location["invariant"]

def GetVar_contiCompute(BlockLines,HA,currentLocation,ContinuourVariables,locations):
    strResult=''
    for Var_Cont in ContinuourVariables :
        for strLine in BlockLines:             
            AnalysesLine=strLine.split('$')
            for i in range(len(AnalysesLine)):  
                if ('&' in AnalysesLine[i]):
                    func=switcherInFunction.get(AnalysesLine[i])
                    if func :                
                        strResult+=func(currentLocation,ContinuourVariables,locations,Var_Cont)
                    else:
                        func=Convert.switcher.get(AnalysesLine[i])
                        if func :                
                            strResult+=func(AnalysesLine,HA)
                else:
                    strResult+=AnalysesLine[i]
    return strResult

def GetinitCheckPoints(currentLocation,ContinuourVariables,locations,checkPoint,Ext=[]): 
    strResult=''  
    listVariables=[]
    listInistValus=[]
    CheckPointVariables=deleteDuplicatedElementFromList([currentLocation['checkPoints'][checkPointSeq]["VariableName"] for checkPointSeq in currentLocation['checkPoints'].keys() ])
    for VariableName in CheckPointVariables:
        listVariables.append('d'+VariableName) 
        listInistValus.append ('9999999' )
    if len(listVariables)    >0:
        strResult=','.join(listVariables)+'='+','.join(listInistValus)
    return strResult

def GetDelat_VariableList(currentLocation,ContinuourVariables,locations,checkPoint,Ext=[]): 
    strResult=''    
    CheckPointVariables=deleteDuplicatedElementFromList(['d'+currentLocation['checkPoints'][checkPointSeq]["VariableName"] for checkPointSeq in currentLocation['checkPoints'].keys() ])
    if len(CheckPointVariables)>0:
        strResult=','.join(CheckPointVariables)    
    else:
        strResult='0'
    return strResult
def GetCheckPointValue(currentLocation,ContinuourVariables,locations,checkPoint,Ext=[]):
    return checkPoint["Value"]
def GetVariableName(currentLocation,ContinuourVariables,locations,checkPoint,Ext=[]):
    return checkPoint["VariableName"]
def GetDelat_Variable(currentLocation,ContinuourVariables,locations,checkPoint,Ext=[]):
    return 'd'+checkPoint["VariableName"]
def GetOther_odes (currentLocation,ContinuourVariables,locations,checkPoint,Ext=[]):
    strResult=''
    listOthers=[]
    for Var_Cont in ContinuourVariables :
        if not Var_Cont== checkPoint["VariableName"]:
            listOthers.append(currentLocation["boxName"]+"_ode_"+Var_Cont)
    strResult=','.join(listOthers)    
    return strResult
switcherInFunction = {
        "&locationname": GetLcationName,
        "&Var_ContsList":GetVar_ContsList,
        "&Var_Cont":GetVar_Cont,
        "&locationInitList":GetLocationInitList,
        "&convert_continuous_variable_to_simpy":Getconvert_continuous_variable_to_simpy,  
        "&Edge_Guard":GetEdge_Guard,  
        "&Reset":GetReset,  
        "&locationDestIndex":GetLocationDestIndex,  
        "&locationDestName":GetLocationDestName,     
        "&Var_Cont_Format":GetVar_Cont_Format,
        "&Invariant":GetInvariant,
        "&VariableName":GetVariableName,
        "&CheckPointValue":GetCheckPointValue,
        "&Other_odes":GetOther_odes,
        "&Delat_VariableList":GetDelat_VariableList,
        "&Delat_Variable":GetDelat_Variable,
        "&initCheckPoints":GetinitCheckPoints,
        "&CurrentlocationIndex":GetCurrentlocationIndex
    } 

   
            
            
def GetLocationFunction(BlockLines,HA,Locations,ContinuourVariables):
    strResult=''  
    for iL in range(len(Locations.keys())-1): 
        print(str(Locations))
        print(iL)
        location=Locations[iL]
        subBlockLines=[]
        isBlockBegin=False
        CurrentBlockName='' 
        for strLine in BlockLines:      
            if "&Block Begin" in strLine:
                if not isBlockBegin:
                    isBlockBegin=True
                    CurrentBlockName=FormatParseLatex.formatBlockName(strLine)
                else:                    
                    subBlockLines.append(strLine)
            elif "&Block End "+CurrentBlockName in strLine:      
                strResult+=(formatSubBlock(CurrentBlockName,subBlockLines, HA,location,ContinuourVariables,Locations))
                isBlockBegin=False;
                subBlockLines=[]
            else:
                if not isBlockBegin:  
                    AnalysesLine=strLine.split('$')
                    for i in range(len(AnalysesLine)):  
                        if ('&' in AnalysesLine[i]):
                            func=switcherInFunction.get(AnalysesLine[i])
                            if func :                
                                strResult+=func(location,ContinuourVariables,Locations,None)
                            else:
                                func=Convert.switcher.get(AnalysesLine[i])
                                if func :                
                                    strResult+=func(AnalysesLine,HA)
                        else:
                            strResult+=AnalysesLine[i]
                else:
                    subBlockLines.append(strLine)
    return strResult

def GetEdges(BlockLines,HA,currentLocation,ContinuourVariables,locations):
    strResult=''
    iEdge=0 
    for EdgeName in HA["lines"].keys():
        Edge= HA["lines"][EdgeName]       
        if Edge["strFromLocation"]==currentLocation["boxName"] and not Edge["strToLocation"]==currentLocation["boxName"]:
            iEdge= iEdge+1
            for strLine in BlockLines:             
                AnalysesLine=strLine.split('$')
                for i in range(len(AnalysesLine)):  
                    if ('&' in AnalysesLine[i]):
                        func=switcherInFunction.get(AnalysesLine[i])
                        if func :                
                            strResult+=func(currentLocation,ContinuourVariables,locations,Edge,[iEdge])
                        else:
                            func=Convert.switcher.get(AnalysesLine[i])
                            if func :                
                                strResult+=func(AnalysesLine,HA)
                    else:
                        strResult+=AnalysesLine[i]
    return strResult
     
 
         
def formatSubBlock(CurrentBlockName, BlockLines,currentHA,currentLocation,ContinuourVariables,locations):
    if (CurrentBlockName== "Edge"):
        return GetEdges(BlockLines,currentHA,currentLocation,ContinuourVariables,locations)
    elif CurrentBlockName=='Var_contiCompute':         
        return GetVar_contiCompute(BlockLines,currentHA,currentLocation,ContinuourVariables,locations)
    elif CurrentBlockName=='CheckPoints':         
        return GetCheckPoints(BlockLines,currentHA,currentLocation,ContinuourVariables,locations)
      
    return ''  

def GetCheckPoints(BlockLines,HA,currentLocation,ContinuourVariables,locations):
    strResult=''
    for checkPointSeq in sorted(currentLocation["checkPoints"].keys()):
        checkPoint= currentLocation["checkPoints"][checkPointSeq]
        for strLine in BlockLines:             
            AnalysesLine=strLine.split('$')         
            for i in range(len(AnalysesLine)):  
                if ('&' in AnalysesLine[i]):
                    func=switcherInFunction.get(AnalysesLine[i])
                    if func :                
                        strResult+=func(currentLocation,ContinuourVariables,locations,checkPoint)
                    else:
                        func=Convert.switcher.get(AnalysesLine[i])
                        if func :                
                            strResult+=func(AnalysesLine,HA)
                else:
                    strResult+=AnalysesLine[i]
    return strResult
     
            
            
            
def deleteDuplicatedElementFromList(list):
        resultList = []
        for item in list:
                if not item in resultList:
                        resultList.append(item)
        return resultList