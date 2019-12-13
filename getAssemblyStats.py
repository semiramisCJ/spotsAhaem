# -*- coding: utf-8 -*-
"""
Created on Fri Jul 31 10:39:58 2015

@author: mcastro
"""

import argparse
import os
from Bio import SeqIO
from os import walk
import numpy as np

###Argument definition
parser = argparse.ArgumentParser()
parser.add_argument("-p", "--path", help="COMPLETE path with input file(s)/folder(s). The output file will be written there. BioPython and Numpy must be installed")
parser.add_argument("-f", "--file", help="Input file name in the k-mer directories. Examples: contigs.fa, final_contigs.fasta, polished_assembly.fasta")
parser.add_argument("-m", "--min", type=int, help="Minimum contig length (# bases). Recommended: 300 as in CLC genomics workbench; NCBI accepts contigs >=200 in size")
parser.add_argument("-s", "--seq", help="Sequencing technology name. Examples: HiSeq, MiSeq, PacBio, etc. This will be used as a prefix for the output file.")
parser.add_argument("-a", "--assembler", help="Assembler name. Examples: velvet, spades, ABySS, smrtools, other. This is used as a guide for the directory structure. 'Other' supports the case when there is only one file in the main directory and there are no subdirectories, as is the case of smrtools, metassembler, unicycler")
parser.add_argument("-e", "--exclude", help="Exclude contigs with len lesser than the m parameter. Default False. If True, it will write a new fasta file with prefix minLen_")
args = parser.parse_args()

###Function definition
def getAssemblyStats(path, minContigLen, inputFileName, filterContigs):
    """For each multifasta file, get N50, number of contigs, average contig length and total length
    Gets as argument: path, min contigs length and input file name"""
    os.chdir(path)
    contigInfo={}
    totalContigs=0
    if filterContigs:
        f=open("minLen_"+str(minContigLen)+"_"+inputFileName,"w")
        for record in SeqIO.parse(inputFileName, "fasta"):
            if len(record.seq) > minContigLen:
                #contigInfo.update({record.name : len(record.seq) })
                f.write('>'+record.name+"\n"+str(record.seq)+"\n")
        f.close()
    
    for record in SeqIO.parse(inputFileName, "fasta"):
        totalContigs+=1
        if len(record.seq) > minContigLen:
            contigInfo.update({record.name : len(record.seq) })
    
    contigs=len(contigInfo)
    aveContigLen=int(np.mean(contigInfo.values()))
    
    lens=sorted(contigInfo.values(), reverse=True)
    totLen= sum(lens)
    halfLen=int(totLen/2)
    
    cumSum=0
    for l in lens:
        cumSum=cumSum+l
        if cumSum >= halfLen:
            N50=l
            break
    
    return (N50, contigs, totalContigs, aveContigLen, totLen)


def walkPaths(mypath):
    """Gets a path as argument. Returns a dictionary with a list of sub-directories and files within them"""
    org={}
    for (dirpath, dirnames, filenames) in walk(mypath):
        org.update({dirpath:[dirnames, filenames]})
    return org


def getAssemblyInfo(path, inputFileName, minContigLen, seqTechnology, assembler, exclude=False):
    """Gets as argument: path, min contigs length and input file name for getAssemblyStats. 
    Also, takes as argument sequencing technology and assembler to write output file with assembly stats"""
    
    dirs=walkPaths(path)
    assemblyInfo={}
    ################Patch for varied user input
    if exclude in ['True', 'true', 'T', 'TRUE', 'yes', 'YES', 'Y', 'y']: 
        filterContigs=True
    else: filterContigs=False
    
    for key in dirs.keys():
        if inputFileName in dirs[key][1]:
            stats=getAssemblyStats(key, minContigLen, inputFileName, filterContigs)
            kmer=key.split("/")[-1]
            
            #if (assembler == 'spades') or (assembler == 'ABySS'):
            #    kmer=kmer
            if assembler == 'velvet':
                kmer=kmer.split("_")[-1]
            elif assembler == 'smrtools':
                kmer="NA"
            else:
                kmer=kmer
            
            assemblyInfo.update({kmer: stats})
    
    os.chdir(path)
    f=open(seqTechnology+"_"+assembler+"_"+str(minContigLen)+"_"+inputFileName+".info","w")
    f.writelines("#kmerSize\tN50\tNumberOfContigsPassMinLen\tTotalContigs\tAveContigLen\tTotLen\n")
    for key in sorted(assemblyInfo.keys()):
        f.write(key+"\t")
        for i in range(len(assemblyInfo[key])):
            if i == len(assemblyInfo[key])-1:
                f.write(str(assemblyInfo[key][i])+"\n")
            else:
                f.write(str(assemblyInfo[key][i])+"\t")
    f.close()

#Call to master function with user-provided arguments
getAssemblyInfo(args.path, args.file, args.min, args.seq, args.assembler, args.exclude)
