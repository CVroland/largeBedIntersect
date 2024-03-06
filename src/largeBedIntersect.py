#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Process the intersection of two large BED files. It can be used for BEDPE files as well. Return the regions of the first file that intersect with at least one region of the second file.
Usage:
    largeBedIntersect2.py -a <bedA> -b <bedB> -o <output> [--strand] [--type first_first|second_second|first_second|second_first]] [--aType bed|bedpe] [--bType bed|bedpe] 


Auteur : Christophe Vroland
Date : 2023 02 28
"""

__authors__ = ("Christophe Vroland")
__date__ = '2023 02 28'
__email__ = ("christophe@cvroland.com")
__status__ = "Prototype"
__version__ = "0.0.1"

import bisect
import sys
import numpy as np
import pandas as pd
from pandas.io.common import get_handle
from collections import defaultdict
from typing import Literal, Generator

verbose=False


BED3_NAMES=[
    "chr",
    "start",
    "end",
]
BED3_DTYPE=defaultdict(
    lambda:"object", #everything that is not specified is not parsed/converted
    {
        0:"string",
        1:"uint32",
        2:"uint32",
    }
)
BED3_USECOLS=[
    0,
    1,
    2,
]
BED3_RENAME={
    0:"chr",
    1:"start",
    2:"end",
}

BED6_NAMES=[
    "chr",
    "start",
    "end",
    "name",
    "score",
    "strand",
]
BED6_DTYPE=defaultdict(
    lambda:"object", #everything that is not specified is not parsed/converted
    {
        0:"string",
        1:"uint32",
        2:"uint32",
        # 3:"string",  # name is not used -> do not parse it
        # 4:"float64", # score is not used -> do not parse it
        5:"string",
    }
)
BED6_USECOLS=[
    0,
    1,
    2,
    5,
]
BED6_RENAME={
    0:"chr",
    1:"start",
    2:"end",
    5:"strand",
}

# https://bedtools.readthedocs.io/en/latest/content/general-usage.html#bedpe-format
BEDPE_NAMES_FULL=[
    "chr1",
    "start1",
    "end1",
    "chr2",
    "start2",
    "end2",
    "name",
    "score",
    "strand1",
    "strand2",
]

BEDPE_DTYPE_FULL=defaultdict( 
    lambda:"string", #everything else is parsed as string
    {
        0:"string",
        1:"uint32",
        2:"uint32",
        3:"string",
        4:"uint32",
        5:"uint32",
        6:"string",
        7:"float64",
        8:"string",
        9:"string",
        10:"string",
    }
)

BEDPE_DTYPE_FIRST=defaultdict(
    lambda:"object", #everything that is not specified is not parsed/converted
    {
        0:"string",#chr1
        1:"uint32",#start1
        2:"uint32",#end1
    }
)
BEDPE_FIRST_TO_BED3_RENAME={
    0:"chr",
    1:"start",
    2:"end",
}

BEDPE_DTYPE_FIRST_STRAND=defaultdict(
    lambda:"object", #everything that is not specified is not parsed/converted
    {
        0:"string",#chr1
        1:"uint32",#start1
        2:"uint32",#end1
        8:"string",#strand1
    }
)
BEDPE_FIRST_TO_BED6_RENAME={
    0:"chr",
    1:"start",
    2:"end",
    8:"strand",
}

BEDPE_DTYPE_SECOND=defaultdict(
    lambda:"object", #everything that is not specified is not parsed/converted
    {
        3:"string",#chr2
        4:"uint32",#start2
        5:"uint32",#end2
    }
)
BEDPE_SECOND_TO_BED3_RENAME={
    3:"chr",
    4:"start",
    5:"end",
}

BEDPE_DTYPE_SECOND_STRAND=defaultdict(
    lambda:"object", #everything that is not specified is not parsed/converted
    {
        3:"string",#chr2
        4:"uint32",#start2
        5:"uint32",#end2
        9:"string",#strand2
    }
)
BEDPE_SECOND_TO_BED6_RENAME={
    3:"chr",
    4:"start",
    5:"end",
    9:"strand",
}


BEDPE_DTYPE_BOTH=defaultdict(
    lambda:"object", #everything that is not specified is not parsed/converted
    {
        0:"string",#chr1
        1:"uint32",#start1
        2:"uint32",#end1
        3:"string",#chr2
        4:"uint32",#start2
        5:"uint32",#end2
    }
)
BEDPE_DTYPE_BOTH_STRAND=defaultdict(
    lambda:"object", #everything that is not specified is not parsed/converted
    {
        0:"string",#chr1
        1:"uint32",#start1
        2:"uint32",#end1
        3:"string",#chr2
        4:"uint32",#start2
        5:"uint32",#end2
        8:"string",#strand1
        9:"string",#strand2
    }
)

RADICLSEQ_BEDPE_NAMES_FULL=[
    "chr_RNA"
    "start_RNA"
    "end_RNA"
    "chr_DNA"
    "start_DNA"
    "end_DNA"
    "qname"
    "score"
    "strand_RNA"
    "strand_DNA"
    "mapq_RNA"
    "mapq_DNA"
    "cigar_RNA"
    "cigar_DNA"
    "gene_id_RNA"
    "window_id_DNA"
    "rescue_status"
]

FileTypes=Literal["bed", "bedpe"]
BedType=Literal["BED", "BEDPE_FIRST", "BEDPE_SECOND", "BEDPE_BOTH"] 
IntersectType=Literal["first_first", "second_second", "first_second", "second_first"]
IntersectTypeWithShort=Literal["ff", "first_first", "ss", "second_second", "fs", "first_second", "sf", "second_first"]
IntersectTypeToLong={
    "ff":"first_first",
    "first_first":"first_first",
    "ss":"second_second",
    "second_second":"second_second",
    "fs":"first_second",
    "first_second":"first_second",
    "sf":"second_first",
    "second_first":"second_first",
}

# TODO : replace those functions by a "chain of responsibility" pattern
def getBedNames(bedType:BedType, stranded:bool):
    """
    Get the names of the columns of a BED file.

    Parameters
    ----------
    bedType : str
        Type of the BED file. Can be "BED" or "BEDPE".
    stranded : bool
        If the strand is considered.

    Returns
    -------
    list
        Names of the columns of the BED file.
    """
    if bedType=="BED":
        if stranded:
            return BED6_NAMES
        else:
            return BED3_NAMES
    elif bedType=="BEDPE_FIRST":
            return BEDPE_NAMES_FULL
    elif bedType=="BEDPE_SECOND":
            return BEDPE_NAMES_FULL
    else:
        raise ValueError(f"Unknown bedType {bedType}")
    
def getBedDtype(bedType:BedType, stranded:bool):
    """
    Get the dtype of the columns of a BED file.

    Parameters
    ----------
    bedType : str
        Type of the BED file. Can be "BED" or "BEDPE".
    stranded : bool
        If the strand is considered.

    Returns
    -------
    dict
        Dtypes of the columns of the BED file.
    """
    if bedType=="BED":
        if stranded:
            return BED6_DTYPE
        else:
            return BED3_DTYPE
    elif bedType=="BEDPE_FIRST":
        if stranded:
            return BEDPE_DTYPE_FIRST_STRAND
        else:
            return BEDPE_DTYPE_FIRST
    elif bedType=="BEDPE_SECOND":
        if stranded:
            return BEDPE_DTYPE_SECOND_STRAND
        else:
            return BEDPE_DTYPE_SECOND
    else:
        raise ValueError(f"Unknown bedType {bedType}")
    
def getBedUseCols(bedType:BedType, stranded:bool):
    """
    Get the usecols of the columns of a BED file.

    Parameters
    ----------
    bedType : str
        Type of the BED file. Can be "BED" or "BEDPE".
    stranded : bool
        If the strand is considered.

    Returns
    -------
    list
        Usecols of the columns of the BED file.
    """
    if bedType=="BED":
        if stranded:
            return BED6_USECOLS
        else:
            return BED3_USECOLS
    elif bedType=="BEDPE_FIRST":
        return [0,1,2]
    elif bedType=="BEDPE_SECOND":
        return [3,4,5]
    else:
        raise ValueError(f"Unknown bedType {bedType}")
    
def getBedRename(bedType:BedType, stranded:bool):
    """
    Get the rename of the columns of a BED file.

    Parameters
    ----------
    bedType : str
        Type of the BED file. Can be "BED" or "BEDPE".
    stranded : bool
        If the strand is considered.

    Returns
    -------
    dict
        Rename of the columns of the BED file.
    """
    if bedType=="BED":
        if stranded:
            return BED6_RENAME
        else:
            return BED3_RENAME
    elif bedType=="BEDPE_FIRST":
        if stranded:
            return BEDPE_FIRST_TO_BED6_RENAME
        else:
            return BEDPE_FIRST_TO_BED3_RENAME
    elif bedType=="BEDPE_SECOND":
        if stranded:
            return BEDPE_SECOND_TO_BED6_RENAME
        else:
            return BEDPE_SECOND_TO_BED3_RENAME
    else:
        raise ValueError(f"Unknown bedType {bedType}")
    




# the index clusters the intervals by chromosome and strand.

def tmpIndexInit():
    """
    Initialize the temporary index.

    Returns
    -------
    dict
        Temporary index.
    """
    return dict()

def tmpIndexAddEntry(tmpIndex:dict[str,tuple[list[int], list[int]]], chrStrand:str, start:int, end:int):
    """
    Add an entry to the temporary index.

    Parameters
    ----------
    tmpIndex : dict
        Temporary index.
    chrStrand : str
        Chromosome and strand.
    start : int
        Start position of the interval.
    end : int
        End position of the interval.
    """
    if chrStrand not in tmpIndex:
        tmpIndex[chrStrand]=([],[])
    indexStart, indexEnd=tmpIndex[chrStrand]
    idx=bisect.bisect_right(indexStart, end)
    # indexStart[idx] > end, the interval is before the current cluster and does not intersect it
    if idx==0 :
        # the interval is before all the cluster in the index : add a new cluster
        indexStart.insert(0,start)
        indexEnd.insert(0,end)
        return
    # look if the interval does not intersect the previous cluster
    preClstrIdx=idx-1
    if indexEnd[preClstrIdx] < start :
        # the interval does not intersect the previous cluster : add a new cluster
        indexStart.insert(idx, start)
        indexEnd.insert(idx, end)
        return
    else :
        # the interval intersect the previous cluster : update the end of the previous cluster
        # indexStart[idx]=min(indexStart[idx], start)
        indexEnd[preClstrIdx]=max(indexEnd[preClstrIdx], end)
        # test if it intersect the previous clusters if needed
        if indexStart[preClstrIdx] > start : # the interval start before the current cluster : it may intersect the previous clusters
            overlapsEndIdx=preClstrIdx
            overlapsBeginIdx=overlapsEndIdx
            while overlapsBeginIdx>0 and indexEnd[overlapsBeginIdx] >= start :
                overlapsBeginIdx-=1
            # update the start of the current cluster with the start of the first intersecting cluster
            indexStart[preClstrIdx]=min(start, indexStart[overlapsBeginIdx])
            # remove the intervals included in the new interval if there are any
            if (overlapsBeginIdx<overlapsEndIdx):
                indexStart[overlapsBeginIdx:overlapsEndIdx]=[]
                indexEnd[overlapsBeginIdx:overlapsEndIdx]=[]
        return

def tmpIndexPushEntry(tmpIndex:dict[str,tuple[list[int], list[int]]], chrStrand:str, start:int, end:int):
    """
    Add an entry at the end of the temporary index.
    It assumes that the start of the new entry is greater or equal to the start of the last cluster of the index.

    Parameters
    ----------
    tmpIndex : dict
        Temporary index.
    chrStrand : str
        Chromosome and strand.
    start : int
        Start position of the interval.
    end : int
        End position of the interval.
    """
    if chrStrand not in tmpIndex:
        tmpIndex[chrStrand]=([],[])
        # case where the index is empty : add a new cluster
        indexStart, indexEnd=tmpIndex[chrStrand]
        indexStart.append(start)
        indexEnd.append(end)
        return
    indexStart, indexEnd=tmpIndex[chrStrand]
    #  case if the interval does not intersect the last cluster
    if indexEnd[-1] < start :
        indexStart.append(start)
        indexEnd.append(end)
        return
    # the interval intersect the last cluster : update the end of the last cluster
    else :
        indexEnd[-1]=max(indexEnd[-1], end)
        return

def tmpIndexSingleChrStrandCreateIndex(startIndex, endIndex, previousIndex=np.empty(0, dtype=bool)):
    """
    Create the final index from the temporary index for a single chromosome and strand.

    Parameters
    ----------
    startIndex : list
        Start positions of intervals. 
    endIndex : list
        End positions of intervals.
    """
    chrStrandSize=max(endIndex[-1], len(previousIndex))
    previousIndex.resize(chrStrandSize, refcheck=False)
    # in Numpy doc : "Enlarging an array: as above, but missing entries are filled with zeros"
    for i in range(len(startIndex)):
        previousIndex[startIndex[i]:endIndex[i]]=True
    return previousIndex



def tmpIndexCreateIndex(tmpIndex:dict[str,tuple[list[int], list[int]]], previousIndex:dict[str,np.ndarray]={}):
    """
    Create the final index from the temporary index.

    Parameters
    ----------
    tmpIndex : dict
        Temporary index.
    """
    index=previousIndex
    for chrStrand, (startIndex, endIndex) in tmpIndex.items():
        index[chrStrand]=tmpIndexSingleChrStrandCreateIndex(startIndex, endIndex, previousIndex=previousIndex.get(chrStrand, np.empty(0, dtype=bool)))
    return index

def indexIntersect(index, chrStrand, start, end):
    """
    Check if an interval intersect with the index.

    Parameters
    ----------
    index : dict
        Index of intervals.
    chrStrand : str
        Chromosome and strand.
    start : int
        Start position of the interval.
    end : int
        End position of the interval.
    """
    if chrStrand not in index:
        return False
    chrStrandIndex=index[chrStrand]
    # XXX : numpy accepts indices out of bounds for slicing : end > len(chrStrandIndex) is not an error
    return chrStrandIndex[start:end].any() 

vIndexIntersect=np.vectorize(indexIntersect, excluded=[0], otypes=[np.bool_], cache=True)


# this is slower than the previous implementation
def indexBed2(bed, bedType="BED", stranded=False, chunksize=10**5):
    """
    Index a BED file.

    Parameters
    ----------
    bed : pathlike or file-like
        Path to the BED file or file-like object.
    bedType : str
        Type of the BED file.
    stranded : bool
        If the strand is considered.
    """
    if verbose : print(f"Indexing", file=sys.stderr)
    bedUseCols=getBedUseCols(bedType, stranded)
    bedDtype=getBedDtype(bedType, stranded)
    bedRename=getBedRename(bedType, stranded)
    # read the bed file by chunks
    chunkId=0
    index={}
    for chunk in pd.read_csv(
        bed, 
        sep="\t", 
        header=None, 
        dtype="object",
        usecols=bedUseCols,
        chunksize=chunksize,
        engine="c"
    ):
        chunkId+=1
        if verbose : print(f"\tchunk {chunkId}", file=sys.stderr)
        chunk=chunk.astype(bedDtype, copy=False)
        chunk.rename(columns=bedRename, inplace=True)
        if stranded:
            chunk["chrStrand"]=chunk["chr"]+chunk["strand"]
        else:
            chunk["chrStrand"]=chunk["chr"]
        # sort to reduce the cache errors
        chunk=chunk.sort_values(by=["chrStrand", "start", "end"])
        # update the index
        chrStrand=chunk["chrStrand"].to_numpy()
        start=chunk["start"].to_numpy()
        end=chunk["end"].to_numpy()
        for i in range(len(chunk)):
            chrStrandIndex=index.get(chrStrand[i], np.empty(0, np.bool_))
            newchrStrandIndexSize=max(np.max(end[i]), len(chrStrandIndex))
            if newchrStrandIndexSize > len(chrStrandIndex):
                chrStrandIndex.resize(newchrStrandIndexSize, refcheck=False)
            chrStrandIndex[start[i]:end[i]]=True
            index[chrStrand[i]]=chrStrandIndex
    return index

def indexBed(bed, bedType="BED", stranded=False, chunksize=10**5):
    """
    Index a BED file.

    Parameters
    ----------
    bed : pathlike or file-like
        Path to the BED file or file-like object.
    bedType : str
        Type of the BED file.
    stranded : bool
        If the strand is considered.
    """
    if verbose : print(f"Indexing", file=sys.stderr)
    bedUseCols=getBedUseCols(bedType, stranded)
    bedDtype=getBedDtype(bedType, stranded)
    bedRename=getBedRename(bedType, stranded)
    # read the bed file by chunks
    chunkId=0
    index={}
    for chunk in pd.read_csv(
        bed, 
        sep="\t", 
        header=None, 
        dtype="object",
        usecols=bedUseCols,
        chunksize=chunksize
    ):
        chunkId+=1
        if verbose : print(f"\tchunk {chunkId}", file=sys.stderr)
        # convert the chunk to the right dtype and rename the columns
        chunk=chunk.astype(bedDtype, copy=False)
        chunk.rename(columns=bedRename, inplace=True)
        if stranded:
            chunk["chrStrand"]=chunk["chr"]+chunk["strand"]
        else:
            chunk["chrStrand"]=chunk["chr"]
        # create a temporary index
        # sort to reduce the cache errors
        chunk=chunk.sort_values(by=["chrStrand", "start", "end"])
        chrStrand=chunk["chrStrand"].to_numpy()
        start=chunk["start"].to_numpy()
        end=chunk["end"].to_numpy()
        tmpIndex=tmpIndexInit()
        for i in range(len(chunk)):
            # because the chunk is sorted by start position, we can push the entries one by one
            tmpIndexPushEntry(tmpIndex, chrStrand[i], start[i], end[i])
        # update the index
        index=tmpIndexCreateIndex(tmpIndex, previousIndex=index)
    return index



def _intersectBed(bedA, index, bedType="BED", stranded=False, chunksize=10**5)->Generator[pd.DataFrame, None, None]:
    """
    Intersect a BED file with an index.

    Parameters
    ----------
    bedA : pathlike or file-like
        Path to the BED file or file-like object.
    index : dict
        Index of intervals.
    bedType : str
        Type of the BED file.
    stranded : bool
        If the strand is considered.
    """
    if verbose : print(f"Intersecting", file=sys.stderr)
    bedDtype=getBedDtype(bedType, stranded)
    bedRename=getBedRename(bedType, stranded)
    # read the bed file by chunks
    chunkId=0
    for chunk in pd.read_csv(
        bedA, 
        sep="\t", 
        header=None, 
        dtype="object",
        chunksize=chunksize,
        engine="c"
    ):
        if verbose : print(f"\tchunk {chunkId}", file=sys.stderr)
        chunkId+=1
        chunk=chunk.astype(bedDtype, copy=False)
        chunk.rename(columns=bedRename, inplace=True)
        maskArray=np.zeros(len(chunk), dtype=np.bool_)
        if stranded:
            chrStrandSerie=chunk["chr"]+chunk["strand"]
        else:
            chrStrandSerie=chunk["chr"]
        chrStrand=chrStrandSerie.to_numpy()
        start=chunk["start"].to_numpy()
        end=chunk["end"].to_numpy()
        maskArray=vIndexIntersect(index, chrStrand, start, end)
        yield chunk[maskArray]

def intersectBed(bedA, bedB, bedAType="BED", bedBType="BED", indexed=False, stranded=False, chunksize=10**5)->Generator[pd.DataFrame, None, None]:
    """
    Intersect two BED files.

    Parameters
    ----------
    bedA : pathlike or file-like
        Path to the first BED file or file-like object.
    bedB : pathlike or file-like
        Path to the second BED file or file-like object.
    bedAType : str
        Type of the first BED file.
    bedBType : str
        Type of the second BED file.
    stranded : bool
        If the strand is considered.
    """
    if indexed:
        index=loadIndex(bedB, bedBType)
    else :
        index=indexBed(bedB, bedBType, stranded, chunksize)
    return _intersectBed(bedA, index, bedAType, stranded, chunksize)


def _intersectTypeToBedTypes(pos:str, fileType:str)->tuple[BedType]:
    """
    Get the type of the BED files for a given intersectType.

    Parameters
    ----------
    pos : str
        Position of the BEDPE file. Can be "first" or "second". If fileType is "bed", it is ignored.
    fileType : str
        Type of the BED file. Can be "bed" or "bedpe".

    Returns
    -------
    str
        Type of the BED files.
    """
    if fileType=="bed":
        return "BED"
    if pos=="first":
        if fileType=="bedpe":
            return "BEDPE_FIRST"
        else:
            raise ValueError(f"Unknown fileType {fileType}")
    elif pos=="second":
        if fileType=="bedpe":
            return "BEDPE_SECOND"
        else:
            raise ValueError(f"Unknown fileType {fileType}")
    elif pos=="both":
        if fileType=="bedpe":
            return "BEDPE_BOTH"
        else:
            raise ValueError(f"Unknown fileType {fileType}")
    else:
        raise ValueError(f"Unknown pos {pos}")

def intersectTypeToBedTypes(intersectType:str, aFileType:str, bFileType:str)->tuple[BedType, BedType]:
    """
    Get the type of the BED files for a given intersectType.

    Parameters
    ----------
    intersectType : str
        Type of the intersection. Can be "first_first", "second_second", "first_second", "second_first".
    aFileType : str
        Type of the first BED file. Can be "bed" or "bedpe".
    bFileType : str
        Type of the second BED file. Can be "bed" or "bedpe".

    Returns
    -------
    str
        Type of the first BED file.
    str
        Type of the second BED file.
    """

    if intersectType=="first_first":
        return _intersectTypeToBedTypes("first", aFileType), _intersectTypeToBedTypes("first", bFileType)
    elif intersectType=="second_second":
        return _intersectTypeToBedTypes("second", aFileType), _intersectTypeToBedTypes("second", bFileType)
    elif intersectType=="first_second":
        return _intersectTypeToBedTypes("first", aFileType), _intersectTypeToBedTypes("second", bFileType)
    elif intersectType=="second_first":
        return _intersectTypeToBedTypes("second", aFileType), _intersectTypeToBedTypes("first", bFileType)
    else:
        raise ValueError(f"Unknown intersectType {intersectType}")

def saveIndex(output, index=None, indexFirst=None, indexSecond=None, compress=True):
    """
    Save the index to a file.

    Parameters
    ----------
    output : pathlike or file-like
        Path to the output file or file-like object.
    index : dict
        Index of intervals. Equivalent to indexFirst.
    indexFirst : dict
        Index of intervals for a single BED file or the first part of a BEDPE file.
    indexSecond : dict
        Index of intervals for the second part of a BEDPE file.

    """
    if index is not None:
        indexFirst=index
    kwds={}
    if indexFirst is not None:
        for chrStrand, idx in indexFirst.items():
            kwds[f"first_{chrStrand}"]=idx
    if indexSecond is not None:
        for chrStrand, idx in indexSecond.items():
            kwds[f"second_{chrStrand}"]=idx
    if compress:
        np.savez_compressed(output, **kwds)
    else:
        np.savez(output, **kwds)

def loadIndex(input, bedType="BED"):
    """
    Load the index from a file.

    Parameters
    ----------
    input : pathlike or file-like
        Path to the input file or file-like object.
    bType : str
        Which part of the BEDPE file to load. Can be "first" or "second".

    Returns
    -------
    dict
        Index of intervals.
    """
    if bedType in ["BED", "BEDPE_FIRST"]:
        bType="first"
    elif bedType in ["BEDPE_SECOND"]:
        bType="second"
    else:
        raise ValueError(f"Unknown bType {bType}")
    with np.load(input) as data:
        index={}
        for key in data.files:
            if key.startswith(bType+"_"):
                index[key[len(bType)+1:]]=data[key]
    return index

def indexBedAndSave(bed, output, bedType="BED", stranded=False, chunksize=10**5):
    """
    Index a BED file and save the index to a file.

    Parameters
    ----------
    bed : pathlike or file-like
        Path to the BED file or file-like object.
    output : pathlike or file-like
        Path to the output file or file-like object.
    bedType : str
        Type of the BED file.
    bType : str
        Which part of the BEDPE file to index. Can be "first", "second" or "both".
    stranded : bool
        If the strand is considered.
    """
    indexFirst=None
    indexSecond=None
    if bedType in ["BED", "BEDPE_FIRST", "BEDPE_BOTH"]:
        indexFirst=indexBed(bed, bedType, stranded, chunksize)
        print(list(indexFirst.keys()))
    if bedType in ["BEDPE_SECOND", "BEDPE_BOTH"]:
        indexSecond=indexBed(bed, bedType, stranded, chunksize)
    saveIndex(output, indexFirst=indexFirst, indexSecond=indexSecond)

def main():
    import argparse
    parser = argparse.ArgumentParser(description="Process the intersection of two large BED files. It can be used for BEDPE files as well.")
    subParsers=parser.add_subparsers(dest="command")

    # intesect command
    intersectParser=subParsers.add_parser("intersect", help="Intersect two BED files.")
    intersectParser.add_argument("-a", "--bedA", type=str, required=True, help="Path to the first BED file. If '-' or 'stdin', read from stdin.")
    intersectParser.add_argument("-b", "--bedB", type=str, required=True, help="Path to the second BED file. If '-' or 'stdin', read from stdin.")
    intersectParser.add_argument("-i", "--indexed", action="store_true", help="If the second BED file is indexed.")
    intersectParser.add_argument("-o", "--output", type=str, default="-", help="Path to the output file. If '-' or 'stdout', write to stdout. Default is '-'.")
    intersectParser.add_argument("-s", "--strand", action="store_true", help="If the strand is considered.")
    intersectParser.add_argument("--type", type=str, default="first_first", help="Which parts of the intervals to compare for BEDPE files. Can be 'first_first' (or 'ff'), 'second_second' (or 'ss'), 'first_second' (or 'fs'), 'second_first' (or 'sf'). Default is 'first_first'.")
    intersectParser.add_argument("--aType", type=str, default="bed", help="Type of the first BED file. It can be 'bed' or 'bedpe'. Default is 'bed'.")
    intersectParser.add_argument("--bType", type=str, default="bed", help="Type of the second BED file. It can be 'bed' or 'bedpe'. Default is 'bed'.")
    intersectParser.add_argument("--chunksize", type=int, default=10**5, help="Size of the chunks to read the BED files. Default is 10^5.")
    intersectParser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode.")

    # index command
    indexBedParser=subParsers.add_parser("index", help="Index a BED file. Return the index in a compressed numpy npz file.")
    indexBedParser.add_argument("bed", type=str, help="Path to the BED file. If '-' or 'stdin', read from stdin.")
    indexBedParser.add_argument("-o", "--output", type=str, default="-", help="Path to the output file. If '-' or 'stdout', write to stdout. Default is '-'.")
    indexBedParser.add_argument("-s", "--strand", action="store_true", help="If the strand is considered.")
    indexBedParser.add_argument("--bType", type=str, default="bed", help="Type of the BED file. It can be 'bed' or 'bedpe'. Default is 'bed'.")
    indexBedParser.add_argument("--type", type=str, default="first", help="Which part of the BEDPE file to index. Can be 'first', 'second' or 'both'. Default is 'first'.")
    indexBedParser.add_argument("--chunksize", type=int, default=10**5, help="Size of the chunks to read the BED files. Default is 10^5.")
    indexBedParser.add_argument("-v", "--verbose", action="store_true", help="Verbose mode.")


    args = parser.parse_args()
    global verbose
    verbose=args.verbose

    if args.command=="index":
        if args.bed in ["-", "stdin"]:
            args.bed=sys.stdin
        if args.output in ["-", "stdout"]:
            args.output=sys.stdout.buffer
        bedType=_intersectTypeToBedTypes(args.type, args.bType)
        indexBedAndSave(args.bed, args.output, bedType, args.strand, args.chunksize)

    elif args.command=="intersect":
        if args.bedA in ["-", "stdin"]:
            args.bedA=sys.stdin
        if args.bedB in ["-", "stdin"]:
            args.bedB=sys.stdin
        if args.output in ["-", "stdout"]:
            args.output=sys.stdout
        args.type=IntersectTypeToLong[args.type]
        bedAType, bedBType=intersectTypeToBedTypes(args.type, args.aType, args.bType)
        with get_handle(args.output, mode="w", compression='infer') as f:
            for chunk in intersectBed(args.bedA, args.bedB, bedAType, bedBType, args.indexed, args.strand, args.chunksize):
                chunk.to_csv(f.handle, sep="\t", index=False, mode="a", header=False)

if __name__ == "__main__":
    main()
