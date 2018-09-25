from sympy import *
from sympy.parsing.sympy_parser import parse_expr 
#sympy.parsing.latex.parse_latex(s)[、、]
from sympy.abc import x,y,z,a,b,c,f,t,k,n
import json 
from latex2sympy.process_latex import process_sympy
import FormatParseLatex
import ConvertFunction
import sys
ContinuourVariables=[]
OrderContinuourVariables=[]
def GetLocations(currentHA):
    Locations={}
    iLocation=1
    isHaveEndLocation=False
    iLocationCount=len(currentHA["boxes"].keys())
    locationkeys=sorted(currentHA["boxes"].keys())
    for LocationName in locationkeys:
        if currentHA["boxes"][LocationName]['isInitial']:
            Locations[0]=currentHA["boxes"][LocationName]
        elif currentHA["boxes"][LocationName]['isEnd']:            
            Locations[iLocationCount-1]=currentHA["boxes"][LocationName]
            isHaveEndLocation=True
        else:            
            Locations[iLocation]=currentHA["boxes"][LocationName]
            iLocation+=1
    if isHaveEndLocation==False:
        Locations[iLocationCount]={"boxName":"","isInitial":False}
    return Locations
    
def convertF(JsonFile='robot.json',OutputFile='GerenateResult.py',HAName='HA1',ModelName='',ttol=10**-2, iterations=1000,TempleteFile='CodeTemplete.txt',until=0.07,vtol=0 ):
    with open (JsonFile, 'r') as fh:        
        currentProject=json.load(fh) 
    
    if not (ModelName in currentProject["Models"].keys()):
        ModelName=next(iter(currentProject["Models"].keys()));
    if not (HAName in currentProject["Models"][ModelName]["HAs"].keys()):
        HAName=next(iter(currentProject["Models"][ModelName]["HAs"].keys()))
    currentHA=currentProject["Models"][ModelName]["HAs"][HAName]
    Locations=GetLocations(currentHA)
    print(currentHA)
    BlockLines=[];
    with open(TempleteFile) as ft:
        codeLines=ft.readlines()
    SaveCodes=[]
    isBlockBegin=False
    CurrentBlockName='' 
    for i in range(len(codeLines)):
        print(codeLines[i])
        if "&Block Begin" in codeLines[i]:
            if not isBlockBegin:
                isBlockBegin=True
                CurrentBlockName=FormatParseLatex.formatBlockName(codeLines[i])
            else:                
                BlockLines.append(codeLines[i])
        elif "&Block End "+CurrentBlockName in codeLines[i]:
            SaveCodes.append(formatBlock(CurrentBlockName,BlockLines,currentHA,ttol , [iterations,vtol,until] , Locations))
            isBlockBegin=False;
            BlockLines=[]
        else:
            if not isBlockBegin:  
                if '$' in codeLines[i]:
                    AnalysesLine=codeLines[i].split("$")
                    SaveCodes.append(formatCode(AnalysesLine,currentHA,Locations,[until]))
                else:
                    SaveCodes.append(codeLines[i])   
            else:
                BlockLines.append(codeLines[i])
            
    with open(OutputFile, 'w') as f:
        f.writelines(SaveCodes)
    print(JsonFile)
    print(TempleteFile)
     
    pass
def GetHAName(AnalysesLine,HA,Locations,Ext=[]):
    return HA["name"].replace(' ','_')


def GetConstants(AnalysesLine,HA,Locations,Ext=[]):
    strResult='' 
    for VariableName in sorted( HA["variables"].keys()):
        if HA["variables"][VariableName]["isConstant"]==True:
            for i in range(len(AnalysesLine)):
                if ('&' in AnalysesLine[i]):    
                    
                    strResult+=VariableName+'=' +HA["variables"][VariableName]["initialValue"]  
                else:
                    strResult+=AnalysesLine[i]
    return strResult


def GetVar_Conts(AnalysesLine,HA,Locations,Ext=[]):
    strResult=''     
    for VariableName in sorted( HA["variables"].keys()):
        if HA["variables"][VariableName]["isConstant"]==False:            
            ContinuourVariables.append(VariableName)
            for i in range(len(AnalysesLine)):
                if ('&' in AnalysesLine[i]):      
                    strResult+=VariableName+'=' +HA["variables"][VariableName]["initialValue"]  
                else:
                    strResult+=AnalysesLine[i]
    return strResult
def GetLocation_Init(AnalysesLine,HA,Locations,Ext=[]):
    strResult=''     
    for iLocation in Locations.keys():
        if Locations[iLocation]['boxName']:
            for i in range(len(AnalysesLine)):
                if ('&' in AnalysesLine[i]):      
                    strResult+=Locations[iLocation]['boxName']+'_FT=False'# +str(Locations[iLocation]["isInitial"]  )
                else:
                    strResult+=AnalysesLine[i]
    return strResult
def GetDicLocationName(AnalysesLine,HA,Locations,Ext=[]):
    strResult='' 
    listResult=[]
    for ik in sorted(Locations.keys()):
        if Locations[ik]['boxName']:
            for i in range(len(AnalysesLine)):
                if ('&' in AnalysesLine[i]):      
                    strResult+=str(ik)+':'+Locations[ik]['boxName']
                    if len(listResult)<len(Locations.keys())-1 and not Locations[ik+1]['boxName']=="":
                        strResult+=','
                else:
                    strResult+=AnalysesLine[i]
            listResult.append(strResult)
        strResult=''
    return ''.join(listResult)
def GetUntil(AnalysesLine,HA,Locations,Ext=[]):
    return str(Ext[0])


def GetlocationEndName(AnalysesLine,HA,Locations,Ext=[]):
    return Locations[len(Locations.keys())-1]["boxName"]
    
def GetEndlocationFuction(BlockLines,currentHA,ttol , iterations,Locations ):
    strResult='' 
    locationEndName=GetlocationEndName(BlockLines,currentHA,Locations)
    if locationEndName=="":
        return ""
    
    for i in range(len(BlockLines)):
        if '$' in BlockLines[i]:
            AnalysesLine=BlockLines[i].split("$")
            strResult=strResult+(formatCode(AnalysesLine,currentHA,Locations))
        else:
            strResult=strResult+(BlockLines[i])    
             
    
    return strResult

def GetEquations(BlockLines,HA,ttol , iterations ):
    strResult='' 
    intBlockLines=len(BlockLines)
    for boxName in sorted(HA["boxes"].keys()):
        iEditLineCount=len(HA["boxes"][boxName]["equation"])        
        for iEditLine in range(iEditLineCount): 
            strData= HA["boxes"][boxName]["equation"][iEditLine].replace('$','')
            if 1==1:
            #try:
                if '=' in strData:
                    leftEquation,rightEquation,variableName=   FormatParseLatex.formatParseLatex4Code(strData,ContinuourVariables)
                    EquationDatas={"&locationname":boxName,
                                  "&equation_Left":leftEquation,
                                  "&equation_right":rightEquation,
                                  "&variableName":variableName,
                                  "&ttol":ttol,
                                  "&iterations":iterations[0],                                  
                                  "&vtol":str(iterations[1])
                                  } 
                    if not variableName in OrderContinuourVariables:
                        OrderContinuourVariables.append(variableName)
                    for strLine in BlockLines:
                        AnalysesLine=strLine.split('$')
                        for i in range(len(AnalysesLine)):
                            if ('&' in AnalysesLine[i]):      
                                data=EquationDatas.get(AnalysesLine[i])
                                if data :   
                                    strResult+=str(data          )
                            else:
                                strResult+=AnalysesLine[i]
        
            else:
            #except Exception as e  : 
                print (str(e))
                pass           
    
    return strResult

switcher = {
        "&haname": GetHAName,
        "&constants":GetConstants,
        "&Var_Conts":GetVar_Conts,
        "&equations":GetEquations,
        '&location_Init':GetLocation_Init,
        '&DicLocationName':GetDicLocationName,
        '&locationEndName':GetlocationEndName,
        '&until':GetUntil
    } 
LineProducing=['&constants','&Var_Conts','&location_Init','&DicLocationName']
        
        
def formatCode(AnalysesLine,currentHA,Locations,Ext=[]):
    strResult=''
    for i in range(len(AnalysesLine)):
        if ('&' in AnalysesLine[i]):
            func=switcher.get(AnalysesLine[i])
            if func :
                if AnalysesLine[i] in LineProducing:
                    strResult=func(AnalysesLine,currentHA,Locations,Ext )
                    return strResult
    strResult=''   
    for i in range(len(AnalysesLine)):
        if ('&' in AnalysesLine[i]):
            func=switcher.get(AnalysesLine[i])
            if func :                
                strResult+=func(AnalysesLine,currentHA,Locations,Ext )
            else:
                func=ConvertFunction.switcherInFunction.get(AnalysesLine[i])
                if func :                
                    strResult+=func(None,OrderContinuourVariables,Locations,None)
                  
        else:
            strResult+=AnalysesLine[i]
    return strResult
def formatBlock(CurrentBlockName, BlockLines,currentHA,ttol , iterations , Locations):
    if (CurrentBlockName== "equations"):
        return GetEquations(BlockLines,currentHA,ttol , iterations )
    if (CurrentBlockName== "EndlocationFuction"):
        return GetEndlocationFuction(BlockLines,currentHA,ttol , iterations,Locations )
    elif CurrentBlockName=='locationFuction':         
        return ConvertFunction.GetLocationFunction(BlockLines,currentHA, Locations,OrderContinuourVariables)
     
    return ''
if __name__ == "__main__":
    #def convertF(JsonFile='robot.json',OutputFile='GerenateResult.py',HAName='HA1',ModelName='',ttol=10**-2, iterations=1000,TempleteFile='CodeTemplete.txt' ):

    if (len(sys.argv)>9):        
        convertF(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8],sys.argv[9]) 
    elif (len(sys.argv)>8):        
        convertF(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7],sys.argv[8]) 
    elif (len(sys.argv)>7):        
        convertF(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6],sys.argv[7]) 
    elif (len(sys.argv)>6):        
        convertF(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5],sys.argv[6] )         
    elif (len(sys.argv)>5):        
        convertF(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4],sys.argv[5] ) 
        
    elif (len(sys.argv)>4):        
        convertF(sys.argv[1],sys.argv[2],sys.argv[3],sys.argv[4] )         
    elif (len(sys.argv)>3):        
        convertF(sys.argv[1],sys.argv[2],sys.argv[3])         
    elif (len(sys.argv)>2):        
        convertF(sys.argv[1],sys.argv[2]  )         
    elif (len(sys.argv)>1):        
        convertF(sys.argv[1]  ) 
    else:
        convertF( )