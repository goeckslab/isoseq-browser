import os, sys
sys.path.append(os.getcwd())

import getGene
from getGene import *
from bokeh.plotting import Figure
from bokeh.models import ColumnDataSource, HoverTool, HBox, VBoxForm, PanTool, WheelZoomTool, BoxZoomTool, ResetTool, ResizeTool, PreviewSaveTool, Range1d, LinearAxis
from bokeh.io import curdoc
from bokeh.models.widgets import Slider, Select, TextInput, Button, VBox, PreText, DataTable, TableColumn, AutocompleteInput
from sklearn.cluster import KMeans
import pandas as pd
from collections import Counter

PLOT_WIDTH = 1200
TITLE_FONT_SIZE = "25pt"
REFERENCE_COLOR = "#22313F"
MATCH_COLOR = "#52B3D9"

def selectGene(opt):
    tranList = list()                                      # list of Transcript objects
    exonList = list()                                      # list of Exon objects
    global isAnnot, isMatch
    if opt.gtf is not None:
        try:
            getGeneFromAnnotation (opt, tranList, exonList)    # lists will be changed
            isAnnot = True
        except (RuntimeError, IOError):
            Console.text = 'Console:\ngene %s is not in the \nannotation file' % opt.gene
            isAnnot = False
    if opt.matches is not None:
        try:
            Console.text = 'Console:\nReading pickle file...'
            getGeneFromMatches (opt, tranList, exonList)       # lists will be changed
            isMatch = True
        except (RuntimeError, IOError):
            Console.text = 'Console:\ngene %s is not in the \nmatch file' % opt.gene
            isMatch = False
    if len(exonList) == 0:
        Console.text = 'Console:\nno exons found for gene %s \nin annotation or match files' % opt.gene
        raise RuntimeError ('Console:\nno exons found for gene %s \nin annotation or match files' % opt.gene)
        return
    return tranList, exonList

def getLength(row):
    return len(row['Clusters'])

def getChromosome(tranList):
    for tran in tranList:
        if tran.annot is False:
            chromosome = tran.chr
            break
    return chromosome

def updateMatch():
    global matchList
    try:
        if matchList == Matches.value.strip().split(','):
            return False
        else:
            matchList = Matches.value.strip().split(',')
            return True
    except NameError:
        matchList = Matches.value.strip().split(',')
        return True

def updateGene(attrname, old, new):
    source.data=dict(x=[], y=[], color=[], line_alpha=[],
                                 QScore=[], start=[], end=[])
    blockSource.data=dict(top=[], bottom=[], left=[], right=[], exon=[],
                        start=[], end=[], chromosome=[], xs=[], ys=[])
    p.title = ""

    allGenes = Counter()
    if updateMatch():
        matchList == Matches.value.strip().split(',')
        for matchFile in matchList:
            geneDict = dict()
            try:
                clusterDict = cl.ClusterDict.fromPickle(matchFile)            # pickle file produced by matchAnnot.py
            except IOError:
                Console.text = 'Console:\nNo such file: %s\n change gene to restart' %matchFile
                break
            myDict = clusterDict.getGeneDict()
            for key, val in myDict.iteritems():
                geneDict.setdefault(key, len(val))
            geneDict = Counter(geneDict)
            allGenes = allGenes + geneDict
        Auto.completions = list(allGenes.keys())
        geneSource.data = dict(
            Gene = allGenes.keys(),
            Cluster = allGenes.values(),
        )

    global Annotations, gtf, opt
    Console.text = 'Console:\nReading annotation file...'
    try:
        Annotations, gtf
        if gtf != GTF.value.strip():
            try:
                Annotations = getAnnotations(getParams(GTF.value.strip(), None, None))
            except IOError:
                Console.text = 'Console:\nNo such file: %s' %str(gtf)
                Annotations = None
            gtf = GTF.value.strip()
    except NameError:
        gtf = GTF.value.strip()
        print gtf
        try:
            Annotations = getAnnotations(getParams(gtf, None, None))
        except IOError:
            Console.text = 'Console:\nNo such file: %s' %str(gtf)
            Annotations = None

    opt = getParams(GTF.value.strip(), matchList,
            Gene.value.strip(), forMat=Format.value.strip(), annotations=Annotations)
    opt.gene = Gene.value.strip()

    global tranList
    tranList, exonList = selectGene(opt)
    chromosome = getChromosome(tranList)

    forwardStrand = '-' if opt.flip else '+'
    if exonList[0].strand == forwardStrand:
        exonList.sort(key=lambda x: x.start)               # sort the list by start position
        blocks = assignBlocks (opt, exonList)              # assign each exon to a block
    else:
        exonList.sort(key=lambda x: x.end, reverse=True)   # sort the list by decreasing end position
        blocks = assignBlocksReverse (opt, exonList)       # assign each exon to a block -- backwards

    findRegions (tranList)                       # determine regions occupied by each transcript
    tranNames = orderTranscripts (tranList)
    tranNames = changeNames(tranNames)
    global length, df, colorDF,boundaryDF, height
    length = len(tranNames)
    Console.text = 'Console:\nCreating plot...'

    maxVal = 0
    minVal = float('inf')
    for tran in tranList:
        if tran.annot == False:
            if maxVal < tran.end:
                maxVal = tran.end
            if minVal > tran.start:
                minVal = tran.start

    p.plot_height = 40*(length+4)
    p.title = "Transcript of %s" % opt.gene
    p.y_range.factors = tranNames[::-1]

    Console.text = 'grouping...'
    if Group.value == "on" and isMatch is True:
        colorDF = groupTran(tranList, exonList, 5)
    else:
        colorDF = None

    getExonData(exonList, colorDF)

    source.data = dict(
        xs=df['xs'],
        ys=df['ys'],
        color=df['colors'],
        line_alpha=df['alpha'],
        x=df['circlex'],
        y=df['circley'],
        QScore=df['QScore'],
        start=df['start'],
        end=df['end'],
        width = df['width'],
    )
    getBoundaryData(blocks)

    xs = list(boundaryDF['xs'])
    ys= list(boundaryDF['ys'])
    xs.insert(0, (0, 0))
    ys.insert(0, (0, length+1))

    blockNum = len(boundaryDF)
    right = list(boundaryDF['boundary'])
    right.insert(0, 0)
    del right[-1]
    blockSource.data = dict(
         top = [(length+1) for x in range(blockNum)],
         bottom = [0 for x in range(blockNum)],
         right = boundaryDF['boundary'],
         left = right,
         start = boundaryDF['start'],
         end = boundaryDF['end'],
         exon = [x+1 for x in range(blockNum)],
         chromosome = [chromosome for x in range(blockNum)],
         xs=xs,
         ys=ys,
    )
    if isAnnot == False:
        Console.text = 'Console:\nSuccess! Annotation\n file is missing.'
    elif isMatch == False:
        Console.text = 'Console:\nSuccess! Match file\n is missing.'
    else:
        Console.text = 'Console:\nSuccess!'

def updateFP(attrname, old, new):
    global df
    df['alpha'] = df.apply(greaterFP, axis=1)
    source.data = dict(
        xs=df['xs'],
        ys=df['ys'],
        color=df['colors'],
        line_alpha=df['alpha'],
        x=df['circlex'],
        y=df['circley'],
        QScore=df['QScore'],
        start=df['start'],
        end=df['end'],
        width=df['width'],
    )

def updateGroup(attrname, old, new):
    colors = list()
    global colorDF, df
    for index, row in df.iterrows():
        color = getColor(row['tran'], colorDF)
        colors.append(color)
    df['colors'] = colors
    source.data = dict(
        xs=df['xs'],
        ys=df['ys'],
        color=df['colors'],
        line_alpha=df['alpha'],
        x=df['circlex'],
        y=df['circley'],
        QScore=df['QScore'],
        start=df['start'],
        end=df['end'],
        width=df['width'],
    )

def updateWidth(attrname, old, new):
    df['width'] = Width.value
    p.plot_height = Width.value*2*(length+4)
    source.data = dict(
        xs=df['xs'],
        ys=df['ys'],
        color=df['colors'],
        line_alpha=df['alpha'],
        x=df['circlex'],
        y=df['circley'],
        QScore=df['QScore'],
        start=df['start'],
        end=df['end'],
        width=df['width'],
    )


def getExonData(exonList, colorDF):
    global df
    df = pd.DataFrame()
    columns = ['name', 'xs', 'ys', 'colors', 'circlex', 'circley', 'QScore',
                'start', 'end', 'tran', 'full', 'partial', 'annot']
    df = pd.DataFrame(columns=columns)
    for myExon in exonList:
        exonSize = myExon.end - myExon.start + 1
        adjStart = myExon.adjStart
        if colorDF is not None:
            color = getColor(myExon.tran.name, colorDF)
        else:
            if myExon.tran.annot:
                color = REFERENCE_COLOR
            else:
                color = MATCH_COLOR
        xs = [adjStart, adjStart+exonSize]
        ys = [length-(myExon.tran.tranIx), length-(myExon.tran.tranIx)]
        circlex = (adjStart + adjStart+exonSize)/2
        circley = length-(myExon.tran.tranIx)
        data = pd.Series([myExon.name, xs, ys, color, circlex, circley,
                        myExon.QScore, myExon.start, myExon.end,
                        myExon.tran.name, myExon.full, myExon.partial,
                        myExon.tran.annot], index=[columns])

        df = df.append(data,ignore_index=True)
    df['alpha'] = 1
    df['width'] = 20

def getColor(exonName, colorDF):
    if exonName not in list(colorDF.name):
        color = '#22313F'
    else:
        row = colorDF.loc[colorDF['name'] == exonName]
        colorName = 'color%s' %str(Cluster.value)
        try:
            color = row[colorName]
        except ValueError:
            color = row['color1']
        except KeyError:
            color = row['color1']
    return color

def getBoundaryData(blocks):
    boundaryDF = pd.DataFrame()
    boundary = list()
    start = list()
    end = list()
    for bound in blocks:
        boundary.append(bound.boundary)
        start.append(bound.start)
        end.append(bound.end)
    global boundaryDF
    boundaryDF['boundary'] = boundary
    boundaryDF['xs'] = zip(boundary, boundary)
    boundaryDF['ys'] = [(0, length+1) for x in range(len(boundaryDF))]
    boundaryDF['start'] = start
    boundaryDF['end'] = end
    boundaryDF = boundaryDF.sort('boundary')

def greaterFP(row):
    alphaVal = Alpha.value
    if row['annot'] is False:
        if row['full'] < Full.value or row['partial'] < Partial.value:
            return alphaVal
    else:
        return 1.0

def saveFasta(attrname, old, new):
    opt.fasta = Save.value.strip()
    tranList = list()
    exonList = list()
    getGeneFromMatches (opt, tranList, exonList)
    opt.fasta = None

def createPlot(df, boundaryDF):
    p = Figure(plot_height=900, plot_width=PLOT_WIDTH, title="", y_range=[],
                title_text_font_size=TITLE_FONT_SIZE)
    p.xgrid.grid_line_color = None
    p.ygrid.grid_line_color = None

    quad = p.quad(top="top", bottom="bottom", left="left", right="right", source=blockSource,
        fill_color="grey", hover_fill_color="firebrick",
        fill_alpha=0.05, hover_alpha=0.3,
        line_color=None, hover_line_color="white")
    p.multi_line(xs="xs", ys="ys", source=blockSource, color="black",
                line_width=2, line_alpha=0.4, line_dash="dotted")
    p.multi_line(xs="xs", ys="ys", source=source, color="color",
                line_width="width", line_alpha='line_alpha')

    p.add_tools(HoverTool(tooltips=[("chromosome", "@chromosome"),("exon", "@exon"),
                ("start", "@start"), ("end", "@end")], renderers=[quad]))
    return p

class getParams(object):
    def __init__(self, gtf, matches, gene, forMat="standard", omit=None,
                 show=None, howmany=None, nodups=None, minlen=None, maxlen=None, output="exon.png",
                 flip=None, yscale=1.0, details=None, fasta=None, title=None, notes=None, full=None,
                 partial=None, highsupport=None, annotations=None):
        self.gtf = gtf
        self.matches = matches
        self.gene = gene
        self.format = forMat                        #format is a keyword
        self.omit = omit
        self.show = show
        self.howmany = howmany
        self.nodups = nodups
        self.minlen = minlen
        self.maxlen = maxlen
        self.output = output
        self.flip = flip
        self.yscale = yscale
        self.details = details
        self.fasta = fasta
        self.title = title
        self.notes = notes
        self.highsupport = highsupport
        self.full = full
        self.partial = partial
        self.annotations = annotations

GTF = TextInput(title="Enter the name of annotation file", value="gencode.vM8.annotation.gtf")
Format = TextInput(title="Enter the format of annotation file, standard is gtf", value="standard")
Matches = TextInput(title="Enter the name of pickle files from MatchAnnot, e.g. of multiple files: match1.pickle,match2.pickle", value="matches.pickle")
Gene = TextInput(title="Select gene to visualize")
Alpha = Slider(title="Alpha value of exons", value=1.0, start=0, end=1.0, step=0.1)
Full = Slider(title="Full support threshold", value=0, start=0, end=30, step=1.0)
Partial = Slider(title="Partial support threshold", value=0, start=0, end=50, step=1.0)
Group = Select(title="Group isoform or not", value="on", options=["on", "off"])
Cluster = Slider(title="The number of groups", value=3, start=1, end=5, step=1.0)
Width = Slider(title="The width of transcripts", value=20, start=5, end=30, step=1)
Save = TextInput(title="Enter the folder name to data in Fasta", value=None)
Auto = AutocompleteInput(title="AutoComplete", completions=[])

blockSource = ColumnDataSource(data=dict(top=[], bottom=[], left=[], right=[], exon=[],
                                start=[], end=[], chromosome=[], xs=[], ys=[]))
source = ColumnDataSource(data=dict(x=[], y=[], color=[], line_alpha=[], width=[],
                            QScore=[], start=[], end=[]))
geneSource = ColumnDataSource(data=dict(Gene=[], Cluster=[]))

df = pd.DataFrame()
boundaryDF = pd.DataFrame()
colorDF = pd.DataFrame()
outDF = pd.DataFrame()


Console = PreText(text='Console:\nStart visualize by entering \nannotations, pickle file and gene.\nPress Enter to submit.\n',
                    width=250, height=100)
p = createPlot(df, boundaryDF)

Gene.on_change('value', updateGene)
Full.on_change('value', updateFP)
Partial.on_change('value', updateFP)
Alpha.on_change('value', updateFP)
Cluster.on_change('value', updateGroup)
Save.on_change('value', saveFasta)
Width.on_change('value', updateWidth)

dataColumns = [
        TableColumn(field="Gene", title="Gene"),
        TableColumn(field="Cluster", title="Cluster"),
    ]
data_table = DataTable(source=geneSource, columns=dataColumns, width=200, height=1200)

paramSource = ColumnDataSource(data=dict(
        Parameter=['annotation', 'format', 'matches', 'gene', 'full', 'alpha', 'partial', 'group', '# of groups', 'fasta'],
        Description=['annotations file, in format specified by --format, Reload page to update' ,
                    'format of annotation file: standard, alt, pickle (default=standard: gtf)',
                    'pickle file from matchAnnot.py, if there are multiple files, seperate them with comma no space, Reload page to update',
                    'gene to plot (required)',
                    'add alpha channel to transcripts without low full/partial support, change full/partial to update',
                    'full support threshold, work together with alpha',
                    'partial support threshold, work together with alpha',
                    'group the trascript on similarity',
                    'assign transcripts into how many groups',
                    'out put the .fasta files to a folder, enter the folder name'],
))

paramColumns = [
        TableColumn(field="Parameter", title="Parameter"),
        TableColumn(field="Description", title="Description"),
    ]
param_table = DataTable(source=paramSource, columns=paramColumns, width=1200, height=700)
main = VBoxForm(p, param_table)

controls = [Console, GTF, Format, Matches, Auto, Gene, Width, Alpha, Full, Partial, Group, Cluster, Save]
inputs = HBox(VBoxForm(*controls), width=250)
curdoc().add_root(HBox(inputs, main, data_table, width=1800))
