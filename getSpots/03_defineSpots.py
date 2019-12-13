# -*- coding: utf-8 -*-
"""
Created on Mon Aug 27 16:58:43 2018

@author: mcastro
"""
import argparse
import re
parser = argparse.ArgumentParser()
parser.add_argument("-f", "--bbhFile", help="COMPLETE path to table with final BBH results")
parser.add_argument("-l", "--bedListFile", help="COMPLETE path to a list of sorted bed files")
parser.add_argument("-b", "--bedDir", help="COMPLETE path to sorted bed files directory")
parser.add_argument("-r", "--ref", help="Reference gbff file name without path")
parser.add_argument("-n", "--nwindow", help="Number of genes for window size to scan trough the genome")
parser.add_argument("-s", "--numStrains", help="Number of studied strains")
parser.add_argument("-sf", "--strainsListFile", help="COMPLETE path to a list of strains; it can be the fnames file used as input for the reciprocal blasts")
parser.add_argument("-o", "--outDir", help="COMPLETE path to output directory")
parser.add_argument("-i", "--silix", help="COMPLETE path to silix dictionary of gene families")
args = parser.parse_args()

def getColumns( infileName, symbol="\t" ):
    "Parses a tab-delimited file and stores it in a list of lists. This allows to easily access the columns"
    data= open(infileName,"rU")
    lines= data.readlines()
    precolumns= [l.strip("\n") for l in lines]
    data.close()
    columns= [l.split(symbol) for l in precolumns]
    return columns


def getWindow(currGenome, i, n):
    r=l=m=[]
    genomeEnd=len(currGenome)
    
    rend=2*n + i + 1
    mcoord = i+n 
    
    if rend < genomeEnd: 
        l=currGenome[ i : mcoord ]
        m=currGenome[ mcoord ]
        r=currGenome[ mcoord+1 : rend]
    
    else:
        
        if mcoord+1 < genomeEnd:
            r1B=rend - genomeEnd
            r=currGenome[mcoord+1:genomeEnd] + currGenome[0:r1B]
            
            l=currGenome[ i : mcoord ]
            m=currGenome[ mcoord ]
        
        elif mcoord+1 == genomeEnd:
            m=currGenome[mcoord]
            r=currGenome[0:n]
            l=currGenome[ i : mcoord]
        
        elif mcoord+1 == genomeEnd+1:
            mnew=mcoord - genomeEnd
            if mnew < n:
                m=currGenome[mnew]
                l=currGenome[ i : mcoord+1]
                r=currGenome[mnew+1:mnew+n+1]
        
        elif mcoord > genomeEnd:
            mnew=mcoord - genomeEnd
            m=currGenome[mnew]
            r=currGenome[mnew+1:mnew+n+1]
            l=currGenome[i:genomeEnd]+currGenome[0:mnew]
    
    return {"m":m, "r":r, "l":l}


def getEnds(start, end, geneDict):
    finals=[]
    for strain in geneDict.keys():
        finals.append(set([ geneDict[strain][start], geneDict[strain][end] ] ))
    
    return finals

def keepLargestInterval(intervals, keys=False):
    """ If keys=False, Intervals dictionary must have the following structure:
    intervals{strain: [genelist1], [genelist2], [genelist3], ... }
    Else, the structure must be 
    intervals{strain: {ID1: [genelist1]}, {ID2:[genelist2]}, {ID3:[genelist3]}, ... }"""
    intervals_proc={}
    
    for strain in intervals:
        if strain not in intervals_proc.keys():
            intervals_proc.update({strain:[]})
        if not keys:
            uniques=[list(i) for i in set(tuple(i) for i in intervals[strain])]
        else: 
            uniques=[list(i) for i in set(tuple(i) for i in intervals[strain].values())]
        
        #Sort by alphabetical order and len
        uniques=sorted(sorted(uniques, key=len, reverse=True))
        discarded=[]
        selected=[]
        notIntersected=[]
        for i in uniques:
            if i in discarded or i in intervals_proc[strain]: #Skip already analyzed elements
                continue
            for j in uniques:
                if i == j: 
                    continue
                if len(set.intersection(set(i), set(j))) == 0:
                    notIntersected.append(j)
                    continue
                chosenList=[]
                #If one is included in the other, keep only the largest one
                if sorted(list(set.intersection(set(i), set(j)))) == sorted(i):
                    chosenList=j
                    discarded.append(i)
                    selected.append(chosenList)
                elif sorted(list(set.intersection(set(i), set(j)))) == sorted(j):
                    chosenList=i
                    discarded.append(j)
                    selected.append(chosenList)
                else:
                    selected.append(i)
                    selected.append(j)
        
        #Keep unique elements from those with overlapping info
        sel_uniq=[list(i) for i in set(tuple(i) for i in selected)]
        intervals_proc[strain].extend(sel_uniq)
        
        #Add unique items with no overlapping info
        notIntersected=[list(i) for i in set(tuple(i) for i in notIntersected)]
        for i in notIntersected:
            if i not in sel_uniq and i not in discarded:
                intervals_proc[strain].append(i)
    
    return intervals_proc

def collapseOverlaps(intervals, keys=False):
    intervals_proc={}
    for strain in intervals:
        if strain not in intervals_proc.keys():
            intervals_proc.update({strain:[]})
        if not keys:
            uniques=[list(i) for i in set(tuple(i) for i in intervals[strain])]
        else: 
            uniques=[list(i) for i in set(tuple(i) for i in intervals[strain].values())]
        
        #Sort by alphabetical order and len
        uniques=sorted(sorted(uniques, key=len, reverse=True))
        
        alreadyCollapsed=[]
        notIntersected=[]
        
        for i in uniques:
            new=i
            if i in alreadyCollapsed: #Skip already analyzed elements
                continue
            for j in uniques:
                if i == j: 
                    continue
                elif len(set.intersection(set(i), set(j))) == 0:
                    notIntersected.append(j)
                    continue                        
                #If there is intersection between new and j
                else: #len(set.intersection(set(new), set(j))) > 0:
                    alreadyCollapsed.append(i)
                    alreadyCollapsed.append(j)
                    
                    
                    ################################################
                    #Analyze the starts
                    for item_n in new:
                        k=0
                        if item_n in j:
                            start=item_n #Get the first element of the intersection
                            k=j.index(start) #Get its index
                            break
                    
                    if k > 0 and set.intersection(set(j[:k]), set(new)) == set(): 
                        #If the first elements in j are not included in the original set, add them
                        new=j[:k]+new
                    
                    ################################################
                    #Analyze the ends
                    for item_n in list(reversed(new)):
                        l=0
                        if item_n in j:
                            end=item_n #Get the last element of the intersection
                            l=j.index(end) #Get its index
                            break
                    
                    if l < len(j)-1 and set.intersection(set(j[l+1:]), set(new)) == set(): 
                        #If the last elements in j are not included in the original set, add them
                        new=new+j[l+1:]
                
            ##Add new superinterval
            intervals_proc[strain].append(new)
            
        alreadyCollapsed=[list(i) for i in set(tuple(i) for i in alreadyCollapsed)]
        #Add unique items with no overlapping info
        notIntersected=[list(i) for i in set(tuple(i) for i in notIntersected)]
        for i in notIntersected:
            if i not in alreadyCollapsed and i not in intervals_proc[strain]:
                intervals_proc[strain].append(i)
    
    return intervals_proc



def getPosOrths(bbhFile, bedListFile, bedDir, ref, n, numStrains, strainsListFile, outDir):    
    n=int(n)
    numStrains=int(numStrains)
    coreBBHs=getColumns(bbhFile)
    strains=coreBBHs[0]
    
    bedFiles=getColumns(bedListFile)
    bedFiles={c[0].split('sorted_')[1].split('.bed')[0]: c[0] for c in bedFiles}
    coreByStrain={}
        
    #RefGene:line
    coreBBHs={c[0]:c[1:] for c in coreBBHs if  '#' not in c[0]} #skip header
    coreList=coreBBHs.keys()
    coreByStrain.update({ref:{c:c for c in coreList}})
    lts={c:c for c in coreList}
    
    #Build dictionary with spaces for all the other strains
    ##And all the genes present in the reference
    for strain in strains:
        strain=strain.replace('#', '')
        if strain == ref: continue #skip reference
        for refGene in coreList:
            if strain not in coreByStrain:           
                coreByStrain.update({strain:{refGene:''}})
            else: 
                coreByStrain[strain].update({refGene:''})
    
    #Get LT for core BBHs per strain 
    for refGene in coreList:
        line=coreBBHs[refGene]
        for i in range(len(strains)-1):
            LT=line[i]
            strain=strains[i+1]
            coreByStrain[strain][refGene]=LT
            lts.update({LT:refGene})
    
    
    #Get ordered core genes by each strain
    coreOrder={}
    allPos={}
    for strain in strains:
        strain=strain.replace('#', '')
        bedInfo=getColumns(bedDir+bedFiles[strain])
        order=[]
        LT=""
        for line in bedInfo:
            LT=line[3]
            start=line[1]
            end=line[2]
            allPos.update({LT:[start, end]})
            if LT in coreByStrain[strain].values():
                refGene=lts[LT]
                order.append(refGene)
        
        coreOrder.update({strain:order})
    
    
    orderedREF=coreOrder[ref]
    permissive={}
    #problematic={}
    guided_coreOrder={}
    
    #Get positional orthologs by windows of n orthologous core genes up and n downstream
    ###Use sliding windows for the reference and each of the other strains
    for i in range(len(orderedREF)):
        res=getWindow(orderedREF,i,n)
        middle=res['m']
        left=res['l']
        right=res['r']
        
        guided_coreOrder.update({i:[]})
        for strain in sorted(coreOrder.keys()):
            #Flags to tell if the upstream and downstream genes are in the same order
            l=r=False 
            if strain == ref: continue #Do not perform self-comparisons (ref vs ref)
            
            orths=coreOrder[strain]
            if middle not in orths: 
                print(middle+" not found in strain "+strain)
                continue
            k=orths.index(middle)
            
            start=k-n
            if start < 0:
                start=len(orths)+start
            
            res2=getWindow(orths, start, n)
            
            if all([p == q for p, q in zip(left, res2['l'])]): l=True
            if all([p == q for p, q in zip(right, res2['r'])]): r=True
            
            if l and r: #Strictly conserved
                guided_coreOrder[i].append(strain)
            
            #For permissive criteria: neighbourhood >= n-1
            refNeighbours=set(left+right)
            queryNeighbours=set(res2['l']+res2['r'])
            neighbourhood=len(set.intersection(refNeighbours, queryNeighbours))
            if neighbourhood >= n-1:
                if middle not in permissive.keys():
                    permissive.update({middle:{}})
                    permissive[middle].update({ref:[left+[middle]+right]})
                    permissive[middle].update({strain:[res2['l']+[res2['m']]+res2['r']] })
                else: 
                    permissive[middle].update({strain:[res2['l']+[res2['m']]+res2['r']]})
            
    
    #Get core genes with conserved position:
    p=open(outDir+"PARTIAL_conservedGenes_strict"+str(numStrains)+".txt","w") 
    #Strictly conserved in some strains, but not all
    posOrths=[]
    for i in guided_coreOrder:
        if len(guided_coreOrder[i]) < 1: continue #skip those conserved only in the reference
        res=getWindow(orderedREF,i,n)
        middle=res['m']
        left=res['l']
        right=res['r']
        
        if len(guided_coreOrder[i]) <= numStrains-1:
            p.write("The order is strictly conserved in "+str(len(guided_coreOrder[i])+1)+" strains: \n")
            p.write("Those strains are: ")
            p.write(ref+', '+", ".join(guided_coreOrder[i])+"\n")
            p.write(str(left+[middle]+right)+"\n")
        
        #Note that the reference is not included in the array
        if len(guided_coreOrder[i]) == numStrains-1: 
            posOrths.extend(left)
            posOrths.append(middle)
            posOrths.extend(right)
    
    
    #Print strictly conserved genes
    f=open(outDir+"STRICTconservedGenes_strict"+str(numStrains)+".txt","w")
    posOrths=list(set(posOrths))
    for gene in posOrths:
        #f.write(gene+"\n")
        k=orderedREF.index(gene)
        start=k-n
        if start < 0:
            start=len(orderedREF)+start
        res=getWindow(orderedREF,k,n)
        strict=res['l']+[res['m']]+res['r']
        
        f.write(str(strict)+"\n")
    
    f.close()
    #print("# of positional orthologs found: "+str(len(posOrths)))
    
    t=open(outDir+"PERMISSIVEconservedGenes_"+str(numStrains)+".txt","w")
    t.write("#middle\tstrain\tgenes\n")
    for middle in permissive.keys():
        if len(permissive[middle].keys()) == numStrains:
            for strain in permissive[middle]:
                t.write(middle+"\t"+strain+"\t"+str(permissive[middle][strain])+"\n")
    
    t.close()
    res=[permissive, coreByStrain, lts, allPos]
    return res

def getIntervals(bbhFile, bedListFile, bedDir, ref, n, numStrains,  strainsListFile, outDir, silix):
    res=getPosOrths(bbhFile, bedListFile, bedDir, ref, n, numStrains, strainsListFile, outDir)
    permissive=res[0]
    coreByStrain=res[1]
    lts=res[2] #Dictionary of LT:referenceLT
    allPos=res[3]
    numStrains=int(numStrains)
    n=int(n)
    
    #comillas=re.compile(r"'+")
    
    #Back translate permissive lists to original LTs and map them to gene families
    permissive_origLTs={}
    permissive_geneFams={}
    
    geneFamilies=getColumns(silix)
    geneFamilies={c[1]:c[0] for c in geneFamilies}
    
    #Translate them only if they have matches with ALL the other strains
    for middle in permissive.keys():
        if len(permissive[middle].keys()) == numStrains:
            permissive_origLTs.update({middle:{}})
            permissive_geneFams.update({middle:{}})
            for strain in permissive[middle].keys():
                permissive_origLTs[middle].update({strain:[]})
                permissive_geneFams[middle].update({strain:[]})
                #Map to original LT
                orderGene=[]
                orderFam=[]
                for refGene in permissive[middle][strain][0]:
                    origLT=coreByStrain[strain][refGene]
                    geneFam=geneFamilies[origLT]
                    orderGene.append(origLT)
                    orderFam.append(geneFam)
                
                permissive_origLTs[middle][strain].extend(orderGene)
                permissive_geneFams[middle][strain].extend(orderFam)
    
    #Create dictionaries and lists for interval processing
    intervals={}
        
    #Check if the ends of the interval belong to the same gene family
    for middle in permissive_origLTs.keys():
        finals=[]
        start=0
        end=len(permissive_origLTs[middle][ref])-1
        finals=getEnds(start, end, permissive_geneFams[middle])
        if set.intersection(*finals) == set.union(*finals):
            #Add the elements of the interval to dictionary only if the condition is met
            for strain in permissive_origLTs[middle]:
                #Add new strains
                if strain not in intervals.keys(): intervals.update({strain:[]})
                intervals[strain].append(permissive_origLTs[middle][strain][start:end+1]) #para abarcar el último gen
        
        #Try to reduce the interval to see if the inner genes may have equivalent extremes
        else:
            while start < n/2 and end > start:
                start=start+1
                end=end-1
                finals=getEnds(start, end, permissive_geneFams[middle])
                if set.intersection(*finals) == set.union(*finals):
                    for strain in permissive_origLTs[middle]:
                        intervals[strain].append(permissive_origLTs[middle][strain][start:end+1])
                    break
    
    #Delete redudancies: small intervals contained in bigger ones
    #Do this a few times to cover some problematic cases    
    identical=False
    r1=intervals
    r2={}
    while not identical:
        r1=keepLargestInterval(r1)
        if r1 == r2:
            identical=True
        r2=r1
    
    intervals_prenr=r2
    
    
    #Write file of intervals_prenr
    intervs=open(outDir+"intervals_prenr.txt","w")
    intervs.write("#strain\tintervalID\tinterval\tsize\n")
    for strain in intervals_prenr:
        for i in intervals_prenr[strain]:
            intervs.write(strain+"\t"+str(i)+"\t"+str(len(i))+"\n")
    intervs.close()
    
    
    #Cat subintervals into superintervals
    #Do this a few times to cover some problematic cases
    identical=False
    r1=intervals_prenr
    r2={}
    while not identical:
        r1=collapseOverlaps(r1)
        if r1 == r2:
            identical=True
        r2=r1
    
    intervals_nr=r2
    
    ####Check the role for the IDs!!!!
    
    #Create dictionary with superintervals in terms of the reference's LT
    translatedSuperIntervals={}
    for strain in intervals_nr:
        newSupIntervalID=0
        #Add strain to superintervals dictionary
        if strain not in translatedSuperIntervals.keys():
            translatedSuperIntervals.update({strain:{}})
        for superinterval in intervals_nr[strain]:
            #Translate superinterval and create new ID for it
            newSupIntervalID+=1
            tr_superinterval=[]
            for gene in superinterval:
                if gene not in lts: break
                else: val=lts[gene]
                tr_superinterval.append(val)
            
            #Add translated superinterval to dictionary but also keep the original
            translatedSuperIntervals[strain].update({newSupIntervalID:[tr_superinterval,superinterval]})
        
    #Check equivalences for each superinterval among strains
    mappings={}
    for ID1 in translatedSuperIntervals[ref]:
        qSupInt=translatedSuperIntervals[ref][ID1][0]
        mappings.update({ID1:[]})
        
        for strain in sorted(translatedSuperIntervals.keys()):
            if strain == ref: continue
            for ID2 in translatedSuperIntervals[strain]:
                currSupInt=translatedSuperIntervals[strain][ID2][0]
                interSize=len( set.intersection( set(qSupInt) , set(currSupInt) ) )
                if interSize > 0 :
                    mappings[ID1].append([strain,ID2,interSize])
                
    #Write to outfile mapped superintervals
    #And build dictionary of superintervals
    dicMapTrSupInt={}
    f=open(outDir+"superintervalMappings.txt","w")
    f.write("#RefID(size): MappedStrains(ID,size)\n")
    for ID in mappings:
        f.write(str(ID)+'('+str(len(translatedSuperIntervals[ref][ID]))+')'+": "+str(mappings[ID])+"\n")
        
        g=open(outDir+"mappedTrSuperintervals_"+str(ID)+".txt","w")
        g.write("#strain\tr_superinterval\tsize\torig_superinterval\torig_start\torig_end\n")
        dicMapTrSupInt.update({ID:{}})
        refSupInterval=translatedSuperIntervals[ref][ID][0]
        #origSupInterval=translatedSuperIntervals[ref][ID][1]
        dicMapTrSupInt[ID].update({ref:refSupInterval})
        g.write(ref+"\t"+str(refSupInterval)+"\t"+str(len(refSupInterval))+"\t")
        posFirst=allPos[refSupInterval[0]]
        posLast=allPos[refSupInterval[-1]]
        g.write(str(refSupInterval)+"\t"+str(posFirst)+"\t"+str(posLast)+"\n")
        
        for item in mappings[ID]:
            if ID not in translatedSuperIntervals[strain]:
                print("ID1 not found in translatedSuperIntervals for strain "+strain)
                continue
            
            strain=item[0]
            ID2=item[1]
            
            if ID2 not in translatedSuperIntervals[strain]:
                print("ID2 not found in translatedSuperIntervals for strain "+strain)
                continue
            
            supInterval=translatedSuperIntervals[strain][ID2][0]
            origSupInterval=translatedSuperIntervals[strain][ID][1]
            
            posFirst=allPos[origSupInterval[0]]
            posLast=allPos[origSupInterval[-1]]
            dicMapTrSupInt[ID].update({strain:supInterval})
            g.write(strain+"\t"+str(supInterval)+"\t"+str(len(supInterval))+"\t")
            g.write(str(refSupInterval)+"\t"+str(posFirst)+"\t"+str(posLast)+"\n")
                
        g.close()
    f.close()
                
    
    return [dicMapTrSupInt, coreByStrain]
    


def buildSpots(bbhFile, bedListFile, bedDir, ref, n, numStrains, strainsListFile, outDir, silix):
    comillas=re.compile(r"'+")
    res=getIntervals(bbhFile, bedListFile, bedDir, ref, n, numStrains,  strainsListFile, outDir, silix)
    
    dicMapTrSupInt=res[0]
    coreByStrain=res[1]
    
    bedFiles=getColumns(bedListFile)
    bedFiles={c[0].split('sorted_')[1].split('.bed')[0]: c[0] for c in bedFiles}
    geneFamilies=getColumns(silix)
    geneFamilies={c[1]:c[0] for c in geneFamilies}
    
    strains=getColumns(strainsListFile)
    strains=[c[0].split('.FAA')[0] for c in strains]
    strains=list(set(strains))
    
    #Get ordered core genes by each strain
    genomeOrder={}
    coords={}
    for strain in strains:
        if strain not in coords.keys():
            coords.update({strain:{}})
        bedInfo=getColumns(bedDir+bedFiles[strain])
        order=[]
        #LT=""
        for line in bedInfo:
            LT=line[3]
            coords[strain].update({LT:[line[1], line[2]]})
            order.append(LT)
        
        genomeOrder.update({strain:order})
    
    spots={}
    spots_geneFam={}
    #dicMapTrSupInt[ID]{strain:supInterval}
    for ID in dicMapTrSupInt:
        for strain in dicMapTrSupInt[ID]:
            currGenome=genomeOrder[strain]
            if strain not in spots.keys():
                spots.update({strain:[]})
                spots_geneFam.update({strain:[]})
            
            interval=dicMapTrSupInt[ID][strain]
            spot=[]
            for i in range(len(interval)-1):
                pre_currOrth= interval[i]
                pre_nextOrth= interval[i+1]
                
                #origLT=coreByStrain[strain][refGene]
                currOrth=coreByStrain[strain][pre_currOrth]
                nextOrth=coreByStrain[strain][pre_nextOrth]
                
                #if strain != ref:
                #    print(strain+"\t"+currOrth+"\t"+nextOrth)
                
                if i == 0:
                    spot.append(currOrth) #Add the first gene just in the first iteration
                
                k=currGenome.index(currOrth) #If the orths genes are contiguous
                if k+1 < len(currGenome) and currGenome[k+1] == nextOrth:
                        spot.append(nextOrth) #Just add the next gene
                elif k+1 == len(currGenome) and currGenome[0] == nextOrth:
                    spot.append(nextOrth)
                else:
                    q=currGenome.index(nextOrth)
                    missing=currGenome[k:q] #This does not includes the qth gene
                    for gene in missing:
                        if gene not in coreByStrain[strain]: #Skip other core genes
                            spot.append(gene)
                    
                    spot.append(nextOrth) #Lastly, add the next ortholog
            
            spots[strain].append([ID,spot]) #Add the spot to the final array
            #print(spot)
    
    #print(spots.values())
    
    f=open(outDir+"spots_LTs.txt","w")
    g=open(outDir+"counts_spots.txt","w")
    h=open(outDir+"coords_spots.txt","w")
    f.write("#strain\tspotID\tspot\tsize\tcore\taccesory\n")
    g.write("#strain\tspotID\tspot_size\tspot_core\tspot_accesory\tuniverse_core\tuniverse_accesory\tuniverse_size\n")
    h.write("#strain\tspotID\tstart_gene\tend_gene\tcoords_start\tcoords_end\n")
    for strain in spots:
        for item in spots[strain]:
            ID=item[0]
            spot=item[1]
            f.write(strain+"\t"+str(ID)+"\t"+str(spot)+"\t"+str(len(spot))+"\t")
            g.write(strain+"\t"+str(ID)+"\t"+str(len(spot))+"\t")
            h.write(strain+"\t"+str(ID)+"\t")
            statusDict={"core":[], "accesory":[]}
            for gene in spot:
                status=""
                if gene not in coreByStrain[strain].values():
                    status="accesory"
                else: 
                    status="core"
                
                statusDict[status].append(gene)
            
            
            f.write(str(statusDict["core"])+"\t"+str(statusDict["accesory"])+"\n")
            g.write(str(len(statusDict["core"]))+"\t"+str(len(statusDict["accesory"]))+"\t")
            g.write(str(len(coreByStrain[strain]))+"\t"+str(len(genomeOrder[strain])-len(coreByStrain[strain]))+"\t"+str(len(genomeOrder[strain]))+"\n")
            
            startCoords=str(coords[strain][spot[0]])
            startCoords=comillas.sub("",startCoords)
            endCoords=str(coords[strain][spot[-1]])
            endCoords=comillas.sub("",endCoords)
            h.write(spot[0]+"\t"+spot[-1]+"\t"+startCoords+"\t"+endCoords+"\n")
    f.close()
    g.close()
    h.close()
    return


#==============================================================================

buildSpots(args.bbhFile, args.bedListFile, args.bedDir, args.ref, 
            args.nwindow, args.numStrains, args.strainsListFile, args.outDir, 
            args.silix)
