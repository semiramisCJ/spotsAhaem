#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jan 18 12:08:28 2019

@author: mcastro
"""

import argparse
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inDir", help="COMPLETE path to the input directory where the *.panoctInput files are located")
parser.add_argument("-p", "--pDir", help="COMPLETE path to the panoct dir")
parser.add_argument("-s", "--strainFile", help="COMPLETE path to the genomesList_panoct.txt file")
parser.add_argument("-c", "--chromFile", help="COMPLETE path to the chromosomes.txt file")
args = parser.parse_args()



def getColumns( infileName ):
    "Parses a tab-delimited file and stores it in a list of lists. This allows to easily access the columns"
    data= open(infileName,"rU")
    lines= data.readlines()
    precolumns= [l.strip("\n") for l in lines]
    data.close()
    columns= [l.split("\t") for l in precolumns]
    return columns

def getUniqueCoreGenes(inDir, pDir, strainFile, chromFile):
    strains=getColumns(strainFile)
    strains=[c[0] for c in strains]
    
    chroms=getColumns(chromFile)
    chroms={c[0]:c[1] for c in chroms}
    
    #Load paralogs IDs
    paralogsRel=getColumns(pDir+"paralogs.txt")
    paralogsList=list(set([item for sublist in paralogsRel for item in sublist]))
    
    #Load homology groups with and without paralogs
    homologyGroups=getColumns(pDir+"matchtable.txt")
    homologyGroups={c[0]:c[1:] for c in homologyGroups}
    
    geneFamilies=open(pDir+"geneFamilies_panoct.txt", "w")
    coreOrth=open(pDir+"coreOrth_panoct.txt", "w")
    coreOrth.write("#"+"\t".join(strains)+"\n") #Write header
    
    chr_coreOrth=open(pDir+"chr_coreOrth_panoct.txt", "w")
    cplxStrains=[s+"_"+chroms[s] for s in strains]
    chr_coreOrth.write("#"+"\t".join(cplxStrains)+"\n") #Write header
    
    #Load genome coords file
    idxGene=1
    
    genomeDict={}
    for strain in strains:
        replicon=chroms[strain]
        genomeCoordsFile=inDir+strain+"_"+replicon+".panoctInput"
        bed=getColumns(genomeCoordsFile)
        bed=[c[idxGene] for c in bed if '#' not in c[0]]
        genomeDict.update({strain:bed})
    
    
    putativeUniqueGenes={}
    for ID in homologyGroups:
        allGeneFamList=[]
        chrGeneFamList=[]
        for i in range(len(strains)):
            LT=homologyGroups[ID][i]
            #strain=strains[i]
            if LT == "----------":
                LT = ''
            
            #Grow the all gene family list
            allGeneFamList.append(LT)
            strain=strains[i]
            if LT in genomeDict[strain]:
                chrGeneFamList.append(LT)
            
            if LT != '':
                #Write gene family to outfile as silix does
                geneFamilies.write(ID+"\t"+LT+"\n")
        
        #If all the LTs are in the chromosomes 
        if len(list(filter(None,chrGeneFamList))) == len(strains) and ID not in paralogsList:
            chr_coreOrth.write("\t".join(chrGeneFamList)+"\n")
        
        #If there is only one gene per gene family
        if len(list(filter(None,allGeneFamList))) == 1: 
            i=list(map(bool, allGeneFamList)).index(True)
            strain=strains[i]
            putativeUniqueGenes.update({ID:strain})
        
        #Check if this gene family is present in all strains
        #and that it has no paralogs in it
        elif len(list(filter(None,allGeneFamList))) == len(strains) and ID not in paralogsList:
            #for i in range(len(strains)):
            coreOrth.write("\t".join(allGeneFamList)+"\n")
        
        
    
    geneFamilies.close()
    coreOrth.close()
    chr_coreOrth.close()
    
    #Check if the paralogs have copies outside the same strain
    exclude=[]
    for ID in putativeUniqueGenes:
        #Get strain for this particular putative unique gene
        curr_strain=putativeUniqueGenes[ID]
        
        #If this ID is in the paralogs list
        if ID in paralogsList:
            for line in paralogsRel:
                #Look for the line it is in
                if ID in line:
                    #Iterate through all the IDs it is related to
                    for item_id in line:
                        #And check to which strain those genes are related
                        if item_id in putativeUniqueGenes and putativeUniqueGenes[item_id] != curr_strain:
                            #If the gene belongs to a different strain, then we have to exclude it
                            exclude.append(ID)
    
    uniqueGenes=open(pDir+"uniqueGenes_panoct.txt", "w")
    chr_uniqueGenes=open(pDir+"chr_uniqueGenes_panoct.txt", "w")
    note="uniqueGenesPanOCT"
    for ID in putativeUniqueGenes:
        if ID in exclude: continue
        strain=putativeUniqueGenes[ID]
        i=strains.index(strain)
        LT=homologyGroups[ID][i]
        uniqueGenes.write(strain+"\t"+LT+"\n")
        
        if LT in genomeDict[strain]:
            replicon=chroms[strain]
            chr_uniqueGenes.write(strain+"_"+replicon+"\t"+LT+"\t"+note+"\n")
    
    uniqueGenes.close()
    chr_uniqueGenes.close()
    
    return
    

getUniqueCoreGenes(args.inDir, args.pDir, args.strainFile, args.chromFile)

