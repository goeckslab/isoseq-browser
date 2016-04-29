import os
import sys
import optparse
import re        # regular expressions
import string

from tt_log import logger

import numpy as np
import cPickle as pickle

import Annotations as anno
import Best        as best
import Cluster     as cl
import CigarString as cs
import pandas as pd
from sklearn.cluster import KMeans

COLORS = ['#52B3D9', '#BE90D4', '#446CB3', '#86E2D5', '#F5D76E']

MIN_REGION_SIZE = 50
FASTA_WRAP = 60                 # bases per fasta line
REGEX_NAME = re.compile ('(c\d+)')      # cluster ID in cluster name
REGEX_LEN  = re.compile ('\/(\d+)$')     # cluster length in cluster name
COMPLTAB   = string.maketrans ('ACGTacgt', 'TGCAtgca')     # for reverse-complementing reads

def getAnnotations (opt):
    omits = [] if opt.omit is None else opt.omit.split(',')            # transcripts which must not be included

    if opt.format == 'pickle':
        annotList   = anno.AnnotationList.fromPickle (opt.gtf)
    elif opt.format == 'alt':
        annotList   = anno.AnnotationList (opt.gtf, altFormat=True)
    else:     # standard format
        annotList = anno.AnnotationList (opt.gtf)

    return annotList

def getGeneFromAnnotation (opt, tranList, exonList):
    '''Add to lists of transcripts and exons: annotations for gene of interest.'''

    if opt.gtf == None:
        return tranList, exonList

    omits = [] if opt.omit is None else opt.omit.split(',')            # transcripts which must not be included

    if opt.annotations:
        annotList = opt.annotations
    else:
        if opt.format == 'pickle':
            annotList   = anno.AnnotationList.fromPickle (opt.gtf)
        elif opt.format == 'alt':
            annotList   = anno.AnnotationList (opt.gtf, altFormat=True)
        else:     # standard format
            annotList   = anno.AnnotationList (opt.gtf)

    allGenes = annotList.getGeneDict()
    if opt.gene not in allGenes:
        raise RuntimeError ('gene %s is not in the annotation file' % opt.gene)
    geneList = allGenes[opt.gene]       # a list of Annotation objects
    if len(geneList) > 1:
        logger.warning('gene %s appears %d times in annotations, first occurrence plotted' \
                           % (opt.gene, len(geneList)))
    myGene = geneList[0]

    for tran in myGene.getChildren():                       # tran is an Annotation object

        if tran.name not in omits:                          # if not in ignore list

            myTran = Transcript(tran.name, start=tran.start, end=tran.end, annot=True, ID=tran.ID)

            if hasattr(tran, 'startcodon'):
                myTran.startcodon = tran.startcodon
            if hasattr(tran, 'stopcodon'):
                myTran.stopcodon = tran.stopcodon

            for exon in tran.getChildren():                 # exon is an Annotation object
                myExon = Exon(myTran, exon.name, exon.start, exon.end, exon.strand)     # no Q score
                if hasattr (exon, 'polyAs'):
                    print exon.name
                    myExon.polyAs = exon.polyAs
                exonList.append (myExon)
                myTran.exons.append(myExon)

            tranList.append (myTran)

    return tranList, exonList

def getGeneFromMatches (opt, tranList, exonList):
    '''Add to lists of transcripts and exons: clusters which matched gene of interest.'''

    if opt.matches == None:
        return tranList, exonList

    omits = [] if opt.omit is None else opt.omit.split(',')            # clusters which must not be included
    shows = [] if opt.show is None else opt.show.split(',')            # clusters which must be included

    localList = list()                                                 # temporary list of clusters
    totClusters = 0

    for matchFile in opt.matches:                                      # --matches may have been specified more thn once

        clusterDict = cl.ClusterDict.fromPickle (matchFile)            # pickle file produced by matchAnnot.py

        for cluster in clusterDict.getClustersForGene(opt.gene):       # cluster is Cluster object
            totClusters += 1

            match = re.search (REGEX_NAME, cluster.name)

            if match is not None and match.group(1) in shows:          # if this is a force-include cluster
                localList.append ( [cluster, 'f999999p999999'] )       # fake sort key to push it to the front

            elif match is None or match.group(1) not in omits:         # shows and omits trump length filter

                full, partial = cluster.getFP()
                sortKey = 'f%06dp%06d' % (full, partial)               # single key includes full and partial counts

                matchLen = re.search(REGEX_LEN, cluster.name)          # filter by cluster length, if requested

                if matchLen is None:
                    raise RuntimeError ('no length in name: %s' % cluster.name)
                    localList.append ( [cluster, sortKey] )            # shouldn't happen -- but let it slide
                else:
                    cLen = int(matchLen.group(1))
                    if opt.minlen is None or cLen >= opt.minlen:
                        if opt.maxlen is None or cLen <= opt.maxlen:
                            localList.append ( [cluster, sortKey] )

    localList.sort(key=lambda x: x[1], reverse=True)                   # sort by full/partial counts

    if opt.nodups:                                                     # eliminate exact dups?

        tempList = list()
        uniqueClusters = set()
        totDups = 0

        for ent in localList:
            cluster = ent[0]
            key = '%9d %s | %s' % (cluster.start, cluster.cigar.prettyPrint(), cluster.cigar.MD)
            if key in uniqueClusters:
                totDups += 1
            else:
                tempList.append(ent)                                   # keep this
                uniqueClusters.add(key)                                # remember it

        localList = tempList
        logger.debug('discarded %d clusters as exact duplicates of each other' % totDups)

    if opt.howmany is not None:
        localList = localList[:opt.howmany]                            # keep the top N entries (which will include the forces)

    totFull = 0
    totPartial = 0

    for ent in localList:
        cluster = ent[0]
        myTran = Transcript(cluster.name, score=cluster.bestScore)
        myTran.chr = cluster.chr

        full, partial = cluster.getFP()
        totFull += full
        totPartial += partial

        myTran.full = full
        myTran.partial = partial


        end = 0
        start = float('inf')

        leading, trailing = cluster.cigar.softclips()

        for exonNum, exon in enumerate(cluster.cigar.exons()):         # exon is a cs.Exon object

            exonName = '%s/%d' % (myTran.name, exonNum)                # exons don't have names: make one up

            if end < exon.end:
                end = exon.end
            if start > exon.start:
                start = exon.start

            if cluster.cigar.MD is not None:                           # if MD string was supplied
                myExon = Exon(myTran, exonName, exon.start, exon.end, cluster.strand, QScore=exon.QScore(), full=full, partial=partial)
            else:
                myExon = Exon(myTran, exonName, exon.start, exon.end, cluster.strand, full=full, partial=partial)

            if exonNum == 0:
                myExon.leading = leading              # add leading softclips to first exon

            exonList.append (myExon)
            myTran.exons.append(myExon)

        myExon.trailing = trailing                    # add trailing softclips to last exon

        myTran.start = start
        myTran.end = end

        tranList.append (myTran)

        if opt.fasta is not None:
            writeFasta (opt, cluster)

    logger.debug('kept %d of %d clusters for gene %s' % (len(localList), totClusters, opt.gene))
    logger.debug('kept clusters include %d full + %d partial reads' % (totFull, totPartial))

    return tranList, exonList

def assignBlocks (opt, exonList):
    '''
    Assign exons to blocks, separated by sequence which is intronic in
    all transcripts. exonList is assumed to be sorted by ascending
    start position.
    '''

    adjust  = 0
    blockNo = 0
    exonIx  = 0
    blocks = list()

    while exonIx < len(exonList):               # can't use enumerate here, it's a double loop

        blockStart = exonList[exonIx].start     # block start = start of first exon in block
        blockEnd   = exonList[exonIx].end       # initial value, updated in the loop below
        blockStartIx = exonIx

        while exonIx < len(exonList) and exonList[exonIx].start <= blockEnd:
            myExon = exonList[exonIx]
            if myExon.end > blockEnd:
                blockEnd = myExon.end
            myExon.block = blockNo
            myExon.tran.blocks.add(blockNo)    # transcript has an exon in this block
            myExon.adjStart = myExon.start - blockStart + adjust
            exonIx += 1

        adjust += blockEnd - blockStart + 1
        blocks.append(Block(blockStart, blockEnd, adjust))
        blockNo += 1

    return blocks

def assignBlocksReverse (opt, exonList):
    '''
    Like assignblocks, but for the reverse strand, ordering blocks
    from the 5' end of the transcript. exonList is assumed to be
    sorted by decreasing exon end position.
    '''

    # I did this as a separate mirror image of assignBlocks, rather
    # than clutter the scenery with lots of forward/reverse checks.

    adjust  = 0
    blockNo = 0
    exonIx  = 0
    blocks = list()

    while exonIx < len(exonList):               # can't use enumerate here, it's a double loop

        blockStart = exonList[exonIx].end       # block start = end of last exon in block
        blockEnd   = exonList[exonIx].start     # initial value, updated in the loop below
        blockStartIx = exonIx

        while exonIx < len(exonList) and exonList[exonIx].end >= blockEnd:

            myExon = exonList[exonIx]
            if myExon.start < blockEnd:
                blockEnd = myExon.start
            myExon.block = blockNo
            myExon.tran.blocks.add(blockNo)    # transcript has an exon in this block
            myExon.adjStart = blockStart - myExon.end + adjust
            exonIx += 1

        adjust += blockStart - blockEnd + 1
        blocks.append(Block(blockStart, blockEnd, adjust))
        blockNo += 1

    return blocks

def findRegions (tranList):
    '''Find breakpoints where coverage by exons changes.'''

    # Why are we doing this? See the note in the Transcript class
    # definition below.

    breaks = list()

    for tranIx, tran in enumerate(tranList):

        for exon in tran.exons:
            breaks.append ([exon.start, 0, tranIx, tran.name, exon.name])
            breaks.append ([exon.end,   1, tranIx, tran.name, exon.name])

    breaks.sort (key=lambda x: x[0])
    curPos = breaks[0][0]
    curTranSet = set()
    region = 0

    for ix in xrange(len(breaks)):

        posit, flag, tranIx, tranName, exonName = breaks[ix]

        if posit > curPos + MIN_REGION_SIZE:             # this is a new region
            if len(curTranSet) > 0:
                for ix in curTranSet:
                    tranList[ix].regions.add(region)     # update set of regions hit by this transcript
                region += 1
            curPos = posit

        if flag == 0:                                    # exon start
####            print '%9d  start  %s' % (posit, exonName)
            curTranSet.add (tranIx)
        else:                                            # exon end
####            print '%9d  end    %s' % (posit, exonName)
            curTranSet.remove (tranIx)

    logger.debug('found %d regions' % region)

    return

def orderTranscripts (tranList):
    '''
    Order the transcripts (i,e., assign each a Y coordinate) so similar
    transcripts are close to each other.
    '''

    # The measure of similarity used here is block occupancy: The
    # distance between two transcripts is the number of blocks where
    # one transcript has exons, and the other doesn't. How many exons
    # there are, or how similar they are in length, is not looked at.

    # The ordering is done using a greedy nearest-neighbor
    # heuristic. To do it optimally turns it into a Traveling Salesman
    # problem.

    tranNames = list()
    curTran = tranList[0]            # arbitrarily start with the first transcript
    tranIx = 0

    while True:                                     # loop until break below
        if curTran.annot is False:
            tranNames.append(curTran.name)              # needed for yticks call
        else:
            tranNames.append(curTran.ID)
        curTran.tranIx = tranIx

        bestTran = best.Best(reverse=True)
        for myTran in tranList:                     # find the next closest transcript
            if myTran.tranIx is None:               # if transcript hasn't been indexed yet
                diff = len(curTran.regions.symmetric_difference(myTran.regions))
                bestTran.update(diff, myTran)

        if bestTran.which is None:                  # every transcript has its index: we're done
            break
####        else:
####            logger.debug('%2d  %s' % (bestTran.value, bestTran.which.name))

        curTran = bestTran.which
        tranIx += 1

    return tranNames

def writeFasta (opt, cluster):
    '''Write a fasta file for a cluster.'''

    # It's fairly common to want to see the sequence for interesting
    # clusters in fasta format. It's convenient to build that into
    # clusterView, since the logic for picking the most populous (or
    # otherwise interesting) clusters is already here.

    if not os.path.exists (opt.fasta):
        os.makedirs (opt.fasta)
    elif not os.path.isdir (opt.fasta):
        raise RuntimeError ('%s exists but is not a directory' % opt.fasta)

    match = re.search (REGEX_NAME, cluster.name)
    if match is None:
        raise RuntimeError ('cannot find cluster ID in %s' % cluster.name)

    if cluster.strand == '+':     # Cluster object includes bases in forward strand sense
        bases = cluster.bases
    else:
        bases = cluster.bases[::-1].translate(COMPLTAB)     # fasta file wants them in read sense

    filename = '%s/%s.fasta' % (opt.fasta, match.group(1))
    handle = open (filename, 'w')
    handle.write ('>%s\n' % cluster.name)

    for ix in xrange(0, len(bases), FASTA_WRAP):
        handle.write (bases[ix:ix+FASTA_WRAP] + '\n')

    handle.close()

def groupTran(tranList, exonList, cluster_num):
    # minVal is the minimum starting point of all the transcripts, maxVal stands for maximum
    maxVal = 0
    minVal = float('inf')
    matchTran = list()
    for tran in tranList:
        if tran.annot == False:
            if maxVal < tran.end:
                maxVal = tran.end
            if minVal > tran.start:
                minVal = tran.start
            matchTran.append(tran)

    df = pd.DataFrame(data=matchTran, columns=['tran'])
    df['min'] = minVal
    df['max'] = maxVal
    df['exons'] = df.apply(getExon, axis=1)
    df['name'] = df.apply(getName, axis=1)

    # Build a matrix contains only true and false
    #
    #   Transcript1:    -----    ----  -- -------
    #   Transcript2:  ----  ------  ------- ---
    #   Superimpose:  ---------------------------
    #   booleanTran:  FFFFFFFFFFFFFFFFFFFFFFFFFFF
    #   Overlap1:     FFTTTTTFFFFTTTTFFTTFTTTTTTT

    #     booleanTran has the same length with superimpose, which is a list of False. Change the overlap region
    #   between the booleanTran and each transcript to True.

    df['boolean'] = df.apply(toBoolean, axis=1)

    # Create a distance table that can be used in K-Means.
    #
    #                        c225/f26p50/6117  c483/f8p23/6083  c20615/f3p27/6185
    # c225/f26p50/6117           0.000000         0.029911           0.012941
    # c483/f8p23/6083            0.029911         0.000000           0.019231
    # c20615/f3p27/6185          0.012941         0.019231           0.000000

    #    The number can be interpreted as the similarity between each two transcript. 0 means they are exactly same
    #  while 1 means they have no overlap region.

    length = len(df)
    index = df['name']
    matrix = [[calcDis(df,i,j) for i in range(length)] for j in range(length)]
    distanceTable = pd.DataFrame()
    distanceTable = pd.DataFrame(matrix)
    distanceTable.columns = index
    distanceTable.index = index

    # Group transcripts, n_clusters set how mant groups should be assigned
    colorDF = pd.DataFrame()
    colorDF['name'] = df['name']
    if len(colorDF) < cluster_num:
        cluster_num = len(colorDF)
    for i in range(cluster_num):
        group = KMeans(n_clusters=i+1).fit_predict(distanceTable)
        global groupName
        groupName = 'group%s' %str(i+1)
        colorName = 'color%s' %str(i+1)
        colorDF[groupName] = group
        colorDF[colorName] = colorDF.apply(assignColor, axis=1)
    return colorDF

def getExon(row):
    startEnd = list()
    exons = row.tran.exons
    for exon in exons:
        startEnd.append((exon.start-row['min'], exon.end-row['min']))
    return startEnd

def getName(row):
    return row.tran.name

def toBoolean(row):
    booleanTran = [False for x in range(row['max']-row['min'])]
    exons = row['exons']
    for exon in exons:
        booleanTran[exon[0]:exon[1]+1] = [True for x in range(exon[1]+1-exon[0])]
    return booleanTran

def calcDis(df, i, j):
    tran1 = df.ix[i]
    tran2 = df.ix[j]
    sum1 = float(sum(tran1['boolean']))
    sum2 = float(sum(tran2['boolean']))
    overlapLength = sum([a and b for a, b in zip(tran1['boolean'], tran2['boolean'])])
    distance = (sum1 + sum2 - 2 * overlapLength) / (sum1 + sum2-overlapLength)
    return distance

def assignColor(row):
    for x in range(5):
        if row[groupName] == x:
            return COLORS[x]
            break

def changeNames(tranNames):
    newTranNames = list()
    for name in tranNames:
        if len(name) >= 30:
            splitList = name.split('|')
            newName = "|".join([splitList[0], splitList[2]])
            newTranNames.append(newName)
        else:
            newTranNames.append(name)
    return newTranNames

class Transcript (object):
    '''Just a struct actually, containing data about a transcript.'''

    def __init__ (self, name, start=None, end=None, score=None, full=None, partial=None, annot=False, ID=None, chr=None):

        self.name    = name
        self.score   = score
        self.annot   = annot            # transcript comes from annotations?
        self.tranIx  = None             # y-axis coordinate of transcript
        self.full    = full
        self.partial = partial
        self.start = start
        self.end = end
        self.ID = ID
        self.exons   = list()           # Exon objects for this transcript
        self.blocks  = set()            # blocks where this transcript has exon(s)
        self.regions = set()            # regions where this transcript has exon(s)
        self.chr = chr
        # What's the difference between a block and a region? Every
        # exon boundary defines a new region. A new block occurs only
        # when exon coverage transitions from 0 to >0. The example
        # below comprises 5 regions, but only one block.

        #    ==============
        #           ==============
        #                      ==================
        #    |      |     |    | |              |

        # I originally used block occupancy to group similar
        # transcripts in orderTranscripts. But that turned out not to
        # work very well, so I invented regions as a finer-grain
        # version of that idea.

class Exon (object):
    '''Struct containing data about an exon.'''

    def __init__ (self, tran, name, start, end, strand, QScore=None, full=None, partial=None):

        self.tran     = tran            # Transcript object containing this exon
        self.name     = name
        self.start    = start
        self.end      = end
        self.strand   = strand
        self.QScore   = QScore
        self.full     = full
        self.partial  = partial
        self.block    = None            # block number where this exon resides
        self.adjStart = None            # start of exon in phony x-axis coordinates
        self.leading  = 0               # number of leading softclipped bases
        self.trailing = 0               # number of trailing softclipped bases

class Block (object):
    '''Struct for plot block.'''

    # A plot block is a vertical span representing a range contiguous
    # bases. Plot blocks are separated by vertical lines representing
    # regions of the reference, of unspecified length, which contain
    # no exons.

    # one implication of that scheme is that the x axis of the plot is
    # meaningless: it represents neither genomic nor RNA sequence range.

    def __init__ (self, start, end, boundary):

        self.start    = start          # actual genomic start coord
        self.end      = end            # actual genomic end coord
        self.boundary = boundary       # right-hand boundary x-coord in phony space
        self.annot    = False          # block contains annotation exons?
