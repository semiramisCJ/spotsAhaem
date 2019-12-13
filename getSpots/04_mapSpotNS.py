#!/usr/bin/env python
# coding: utf-8

import argparse
import re
parser = argparse.ArgumentParser()
parser.add_argument("-i", "--inDir", help="COMPLETE path to the input directory with chromosomes.txt and *.panoctInput files")
parser.add_argument("-r", "--ref", help="Reference genome")
parser.add_argument("-s", "--spotsDir", help="COMPLETE path to where is located the spots_LTs.txt file")
parser.add_argument("-o", "--outDir", help="COMPLETE path to the output directory; it must end with slash [/]")
args = parser.parse_args()

comillas=re.compile(r"'+")

def getColumns( infileName ):
    "Parses a tab-delimited file and stores it in a list of lists. This allows to easily access the columns"
    data= open(infileName,"rU")
    lines= data.readlines()
    precolumns= [l.strip("\n") for l in lines]
    data.close()
    columns= [l.split("\t") for l in precolumns]
    return columns

#####Load spots
def loadSpots(spotsDir):
    #####Build dictionary {strain:{spot1:[elem1, elem2..], spot2:[...], ...}
    spots=getColumns(spotsDir+"spots_LTs.txt")
    #strain	spotID	spot	size	core	accesory
    spots=[c for c in spots if '#' not in c[0]]
    spotDict_byGene={}
    for line in spots:
        strain=line[0]
        if strain not in spotDict_byGene.keys():
            spotDict_byGene.update({strain:{}})
        spotID=line[1]
        spotGenes=line[2]
        spotGenes=comillas.sub("",spotGenes)
        spotGenes=spotGenes.replace("[", "")
        spotGenes=spotGenes.replace("]", "")
        spotGenes=spotGenes.split(', ')
        for gene in spotGenes:
            spotDict_byGene[strain].update({gene:spotID})
    
    return spotDict_byGene

def writeSortedSpots(refName, sortedSpotsDict, outDir):
    ref=sortedSpotsDict[refName]
    startRef=ref[0]
    #Open outfile for sortedSpots
    sortedSpots=open(outDir+"sortedSpots.txt","w")
    sortedSpots.write("#strain\tsortedSpots\n") 
    #Write info for reference
    sortedSpots.write(refName+"\t"+" ".join(sortedSpotsDict[refName])+"\n")
    for strain in sortedSpotsDict:
        if strain == refName: continue #skip reference as we have already written that
        currList=sortedSpotsDict[strain]
        start=currList.index(startRef)
        new=currList[start:]+currList[:start]
        #Write info for current strain
        sortedSpots.write(strain+"\t"+" ".join(new)+"\n")
    
    sortedSpots.close()
    return 

#Load genome info
def loadGenomeData_writeTable(inDir, ref, spotsDir, outDir):
    #####Load list of ordered complete genomes
    #####Build dictionary {strain:[gene1, gene2, ...]}
    #####Focus on chromosomes
    mappings=getColumns(inDir+"chromosomes.txt")
    mappings={c[1]:c[0] for c in mappings} #replicon:simpleStrain
    genomeDict={}
    genomeCoords={}
    strainList=[]
    for replicon in mappings:
        simpleStrain=mappings[replicon]
        strain=simpleStrain+"_"+replicon
        #bedFile="sorted_"+strain+".bed"
        #AN4.gbff_CP031983.panoctInput
        bedFile=strain+".panoctInput"
        if strain not in genomeDict.keys():
            genomeDict.update({strain:[]})
            genomeCoords.update({strain:[]})
            strainList.append(strain)
        
        bedGenes=getColumns(inDir+bedFile)
        #start, end, product
        bedCoords={c[1]:[c[2],c[3], c[4]] for c in bedGenes if '#' not in c[0]}          
        bedGenes=[c[1] for c in bedGenes if '#' not in c[0]]
        genomeDict[strain]=bedGenes
        genomeCoords[strain]=bedCoords
    
    #Load spot/NonSpot data
    spotDict_byGene=loadSpots(spotsDir)
    
    #Name NSs
    #compIDdict={}
    sortedSpotsDict={}
    #Open outfile for sizes of spots and non spots
    sizesTable=open(outDir+"sizes.txt","w")
    sizesTable.write("strain\tcategory\tlen\n")
    for strain in spotDict_byGene:
        #Open outfile and write header
        megaTable=open(outDir+strain+"_spotNS.txt","w")
        megaTable.write("#strain\tgene\tstart\tend\tproduct\tspotID\n")
        
        #Create dictionary of sorted spots and non-spots
        sortedSpotsDict.update({strain:[]})
        
        #Create list of only spots
        onlySpots=[]
        
        #Create dict of sizes
        sizesAll={}
        
        #Create list and dictionary of compIDs
        compIDlist=[]
        compIDdict={}
        
        genomeGenes=genomeDict[strain]
        prevSpot=""
        last_prevSpot=""
        c=0
        for gene in genomeGenes:
            if gene in spotDict_byGene[strain]:
                spotID=spotDict_byGene[strain][gene]
                
                #Keep record of the order of only the spots
                if spotID not in onlySpots: onlySpots.append(spotID)
                
                prevSpot=spotID
            else:
                if prevSpot != last_prevSpot:
                    c+=1
                #spotID="NS"+str(c)+"_nextTo"+str(prevSpot)
                spotID="NS"+str(c)+"_between_"+str(prevSpot)+"_and_"
                last_prevSpot=prevSpot
            
            compID=strain+"-"+spotID
            compIDdict.update({gene:compID})
            compIDlist.append(compID)
        
        #If the chromosome starts and ends with a "NS", they must have the same name
        #So, we will re-iterate over the entire genome...
        if "NS" in compIDlist[0] and "NS" in compIDlist[-1]:
            prevSpot=""
            last_prevSpot=""
            c=0
            for gene in genomeGenes:
                if gene in spotDict_byGene[strain]:
                    spotID=spotDict_byGene[strain][gene]
                    prevSpot=spotID
                else:
                    if prevSpot != last_prevSpot:
                        c+=1
                    #spotID="NS"+str(c)+"_nextTo"+str(prevSpot)
                    spotID="NS"+str(c)+"_between_"+str(prevSpot)+"_and_"
                    last_prevSpot=prevSpot
                    
                    compID=strain+"-"+spotID
                    if compID == compIDlist[-1]: #When reaching the last NS
                        compIDdict.update({gene:compIDlist[0]}) #Correct the ID
        
        #Initial value for nextSpot set to None        
        nextSpot=None
        #Now, write table
        compIDlist=[]
        #prevSpot=onlySpots[0]
        for gene in genomeGenes:
            start=genomeCoords[strain][gene][0]
            end=genomeCoords[strain][gene][1]
            product=genomeCoords[strain][gene][2]
            compID=compIDdict[gene]
            simpleID=compID.split(strain+"-")[1]
            
            #Only for spots
            if "NS" not in simpleID and "between" not in simpleID:
                #We take the current spot as if it was the previous one (it will be if the next element is a non-spot)
                prevSpot=simpleID
                
                #We locate the index of the current spot and sum 1 unit
                i=onlySpots.index(prevSpot)+1
                
                #The 'next' spot will be the next on the list
                if i < len(onlySpots):
                    nextSpot=onlySpots[i]
                else: 
                    nextSpot=onlySpots[0]
            
            #When we find a non-spot
            else:
                #If we haven't read any Spot before
                if nextSpot == None:
                    nextSpot=onlySpots[0]
                    prevSpot=onlySpots[-1]
                    compID=compID.split("_between")[0]
                    compID=compID+"_between_"+prevSpot+"_and_"+nextSpot
                #For the next entries with the same situation as above
                elif compID.split("_between")[1] == '__and_':
                    nextSpot=onlySpots[0]
                    prevSpot=onlySpots[-1]
                    compID=compID.split("_between")[0]
                    compID=compID+"_between_"+prevSpot+"_and_"+nextSpot
                else: #Normal situations
                    compID=compID+nextSpot
            
            megaTable.write(strain+"\t"+gene+"\t"+start+"\t"+end+"\t"+product+"\t"+compID+"\n")            
            
            #Update size per strain
            if simpleID not in sizesAll:
                sizesAll.update({simpleID:0})
            sizesAll[simpleID]+=1
            
            if simpleID not in compIDlist and "between" not in simpleID:
                compIDlist.append(simpleID)
        
        megaTable.close()
        sortedSpotsDict[strain]=compIDlist
        
        #Write entries of sizes for this particular strain
        for simpleID in sizesAll:
            if "between" not in simpleID:
                category="Spot"
            else:
                category="NonSpot"
            size=str(sizesAll[simpleID])
            sizesTable.write(strain+"\t"+category+"\t"+size+"\n")
    
    
    sizesTable.close() #Close outfile after iterating through all the strains
    #Write sorted spots file
    writeSortedSpots(ref, sortedSpotsDict, outDir)
    return

loadGenomeData_writeTable(args.inDir, args.ref, args.spotsDir, args.outDir)

