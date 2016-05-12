import os, sys
sys.path.append(os.getcwd())

import getGene
from getGene import *
from bokeh.plotting import Figure
from bokeh.models import ColumnDataSource, HoverTool, HBox, VBoxForm, PanTool, WheelZoomTool, BoxZoomTool, ResetTool, ResizeTool, PreviewSaveTool, Range1d, LinearAxis
from bokeh.io import curdoc
from bokeh.models.widgets import Slider, Select, TextInput, Button, VBox, PreText, DataTable, TableColumn
from sklearn.cluster import KMeans
import pandas as pd
from collections import Counter

PLOT_WIDTH = 1200                                           # the width of plott
TITLE_FONT_SIZE = "25pt"                                    # font size of plot title
REFERENCE_COLOR = "#22313F"                                 # the color of Reference transcripts
MATCH_COLOR = "#52B3D9"                                     # the color of matched isoforms

# Select isoforms of a particular gene
def selectGene(opt):
    tranList = list()                                      # list of Transcript objects
    exonList = list()                                      # list of Exon objects
    global isAnnot, isMatch                                # whether annotation transcripts or matched isoforms are in the inputs

    if opt.gtf is not None:
        try:
            getGeneFromAnnotation (opt, tranList, exonList)                 # find transcripts, exons of the gene in the reference database
            isAnnot = True
        except RuntimeError:                                                # if gene is not in the annotation file
            Console.text = 'Console:\ngene %s is not in the \nannotation file' % opt.gene
            isAnnot = False
        except IOError:                                                     # if annotation file is not found
            Console.text = 'Console:\nannotations file %s is not found' % opt.gtf
            isAnnot = False
    if opt.matches is not None:
        try:
            Console.text = 'Console:\nReading pickle file...'
            getGeneFromMatches (opt, tranList, exonList)                    # find transcripts, exons of the gene in the matched isoforms
            isMatch = True
        except RuntimeError:                                                # if gene is not in the matched file
            Console.text = 'Console:\ngene %s is not in the \nmatch file' % opt.gene
            isMatch = False
        except IOError:                                                     # if mached file is not found
            Console.text = 'Console:\nannotations file %s is not found' % opt.matches
            isMatch = False

    if len(exonList) == 0:                                                  # if there is no exons in the annotation or match files
        Console.text = 'Console:\nno exons found for gene %s \nin annotation or match files' % opt.gene
        raise RuntimeError ('Console:\nno exons found for gene %s \nin annotation or match files' % opt.gene)
        return
    return tranList, exonList

def getChromosome(tranList):                                # find out which chromosome does the gene locate on
    chromosome = None
    for tran in tranList:
        if tran.annot is False:
            chromosome = tran.chr
            break
    return chromosome

def updateMatch():                                      # update visuzlization  if the match files changes, return True if update is needed
    checkUpdate = Matches.value.strip().split(',')
    checkUpdate = ['../' + x for x in checkUpdate]
    try:
        if matchList == checkUpdate:
            return False
        else:
            return True
    except NameError:
        return True

def updateGene(attrname, old, new):                     # update visualization if the gene changes
    source.data=dict(x=[], y=[], color=[], line_alpha=[],                       # the data source for plotting transcripts
                                 QScore=[], start=[], end=[])
    blockSource.data=dict(top=[], bottom=[], left=[], right=[], exon=[],        # the data source for plotting block boundary
                        start=[], end=[], chromosome=[], xs=[], ys=[])
    p.title = "Transcript of %s" % Gene.value.strip()                           # update the title of plot

    # create and display a list of all the genes
    global matchList
    if updateMatch():                                               # if matches files are updated
        matchList = Matches.value.strip().split(',')
        allGenes = Counter()                                        # create a counter hastable(dictionary) object
        for matchFile in matchList:                                 # there can be multiple match files as inputs
            print '----------------------------------------------'
            print matchFile
            try:
                clusterDict = cl.ClusterDict.fromPickle(matchFile)            # pickle file produced by matchAnnot.py
                print clusterDict
            except IOError:
                Console.text = 'Console:\nNo such file: %s\n change gene to restart' %matchFile
                break
            geneDict = dict()                                               # how many isoforms for each gene
            myDict = clusterDict.getGeneDict()
            for key, val in myDict.iteritems():
                geneDict.setdefault(key, len(val))
            geneDict = Counter(geneDict)
            allGenes = allGenes + geneDict                                  # add information of new match file
        geneSource.data = dict(                                             # update the table in visualization
            Gene = allGenes.keys(),
            Isoforms = allGenes.values(),
        )

    # load the reference transcripts
    global Annotations, gtf, opt
    Console.text = 'Console:\nReading annotation file...'
    try:                                                                    # if it's not the first time make the plot
        Annotations, gtf
        if gtf != GTF.value.strip():                                        # and if the annotation file changes
            try:
                Annotations = getAnnotations(getParams(GTF.value.strip(),   # get the lists of transcripts in annot file
                                None, None))                                # hold it in RAM for quicker respond when changing the gene
            except IOError:                                               # annotation file not found
                Console.text = 'Console:\nNo such file: %s' %str(gtf)
                Annotations = None
            gtf = GTF.value.strip()
    except NameError:                                   # if it's the first time make the plot, do the same
        gtf = GTF.value.strip()
        try:
            Annotations = getAnnotations(getParams(gtf, None, None))
        except IOError:
            Console.text = 'Console:\nNo such file: %s' %str(gtf)
            Annotations = None

    opt = getParams('../' + GTF.value.strip(), matchList,               # get the options from inputs
            Gene.value.strip(), forMat=Format.value.strip(), annotations=Annotations)
    opt.gene = Gene.value.strip()

    global tranList
    tranList, exonList = selectGene(opt)                    # select transcripts by gene
    chromosome = getChromosome(tranList)                    # find out which chromosome does the gene locate

    forwardStrand = '-' if opt.flip else '+'
    if exonList[0].strand == forwardStrand:
        exonList.sort(key=lambda x: x.start)               # sort the list by start position
        blocks = assignBlocks (opt, exonList)              # assign each exon to a block
    else:
        exonList.sort(key=lambda x: x.end, reverse=True)   # sort the list by decreasing end position
        blocks = assignBlocksReverse (opt, exonList)       # assign each exon to a block -- backwards

    findRegions (tranList)                       # determine regions occupied by each transcript
    tranNames = orderTranscripts (tranList)                 # get the names of transcripts, placed them in the right order
    tranNames = changeNames(tranNames)                      # if the length of name is too long, reduce it

    global length, df, colorDF,boundaryDF
    length = len(tranNames)
    Console.text = 'Console:\nCreating plot...'

    p.plot_height = Width.value*2*(length+4)                   # set the height of plot according to the length of transcripts

    p.y_range.factors = tranNames[::-1]             # set the y axis tick to the transcripts names

    Console.text = 'grouping...'
    if Group.value == "on" and isMatch is True:
        colorDF = groupTran(tranList, exonList, 5)          # group the transcripts by similarities
    else:
        colorDF = None

    getExonData(exonList, colorDF)                          # get the data of each isoform that can be directly used to plot, put it into a global dataframe

    source.data = dict(                                     # update the data for plot
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
    getBoundaryData(blocks)                                 # get the data of each block that can be directly used to plot, put it into a global dataframe

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

# update the number of groups
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

# update the width of the each isoform
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
    df['width'] = Width.value

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

# find out the position of boundaries
def getBoundaryData(blocks):
    boundaryDF = pd.DataFrame()
    boundary = list()
    start = list()
    end = list()
    for bound in blocks:
        boundary.append(bound.boundary)
        start.append(bound.start)
        end.append(bound.end)
    # add data used for plotting to pandas DataFrame
    global boundaryDF
    boundaryDF['boundary'] = boundary
    boundaryDF['xs'] = zip(boundary, boundary)
    boundaryDF['ys'] = [(0, length+1) for x in range(len(boundaryDF))]
    boundaryDF['start'] = start
    boundaryDF['end'] = end
    boundaryDF = boundaryDF.sort('boundary')

# find out which isoform is below full/partial threshold, which isoform is not
def greaterFP(row):
    alphaVal = Alpha.value
    if row['annot'] is False:
        if row['full'] < Full.value or row['partial'] < Partial.value:
            return alphaVal
    else:
        return 1.0

# save the transcripts to .fasta file, the function is copied from MatchAnnot
def saveFasta(attrname, old, new):
    opt.fasta = Save.value.strip()
    tranList = list()
    exonList = list()
    getGeneFromMatches (opt, tranList, exonList)
    opt.fasta = None

# create the visualization plot
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

# a class containing all the input parameters
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

# create all kinds of widgets
GTF = TextInput(title="Enter the name of annotation file", value="gencode.vM9.annotation.gtf")
Format = TextInput(title="Enter the format of annotation file, standard is gtf", value="standard")
Matches = TextInput(title="Enter the name of pickle files from MatchAnnot, e.g. of multiple files: a.pickle,b.pickle", value="matches.pickle")
Gene = TextInput(title="Select gene to visualize")
Alpha = Slider(title="Alpha value of exons", value=1.0, start=0, end=1.0, step=0.1)
Full = Slider(title="Full support threshold", value=0, start=0, end=30, step=1.0)
Partial = Slider(title="Partial support threshold", value=0, start=0, end=50, step=1.0)
Group = Select(title="Group isoform or not", value="on", options=["on", "off"])
Cluster = Slider(title="The number of groups", value=3, start=1, end=5, step=1.0)
Width = Slider(title="The width of transcripts", value=20, start=5, end=30, step=1)
Save = TextInput(title="Enter the folder name to data in Fasta", value=None)

# the data used for plotting isoforms, boundaries and gene
blockSource = ColumnDataSource(data=dict(top=[], bottom=[], left=[], right=[], exon=[],
                                start=[], end=[], chromosome=[], xs=[], ys=[]))
source = ColumnDataSource(data=dict(x=[], y=[], color=[], line_alpha=[], width=[],
                            QScore=[], start=[], end=[]))
geneSource = ColumnDataSource(data=dict(Gene=[], Isoforms=[]))

# the console box
Console = PreText(text='Console:\nStart visualize by entering \nannotations, pickle file and gene.\nPress Enter to submit.\n',
                    width=250, height=100)
# the visualization plot
p = createPlot(pd.DataFrame(), pd.DataFrame())

# a table of with all the genes in the match files, and how many isoforms in each gene
geneColumns = [
        TableColumn(field="Gene", title="Gene"),
        TableColumn(field="Isoforms", title="Isoforms"),
    ]
geneCountTable = DataTable(source=geneSource, columns=geneColumns, width=200, height=1200)

# a description of all the parameters
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
paramTable = DataTable(source=paramSource, columns=paramColumns, width=1200, height=700)

# make changes to the plot when widgets are updated
Gene.on_change('value', updateGene)
Full.on_change('value', updateFP)
Partial.on_change('value', updateFP)
Alpha.on_change('value', updateFP)
Cluster.on_change('value', updateGroup)
Save.on_change('value', saveFasta)
Width.on_change('value', updateWidth)

# the position of plot and widgets
controls = [Console, GTF, Format, Matches, Gene, Width, Alpha, Full, Partial, Group, Cluster, Save]
main = VBoxForm(p, paramTable)
inputs = HBox(VBoxForm(*controls), width=250)
curdoc().add_root(HBox(inputs, main, geneCountTable, width=1800))
