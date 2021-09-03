#
# Part of this  code is taken from https://groups.google.com/forum/#!topic/etetoolkit/cYTHHsL21KY
# written by Jaime Huerta-Cepas

"""Ecrit un  arbre au format XML .
Usage:
  genere_xml.py <treeFile> <speciesDicoFile> <syntenyLeftFile> <syntenyRightFile>

Positional arguments:
  treeFile              family name / phylogenic trees in newick format
  speciesDicoFile       sequence name / species
  synLeftFile           sequence name / list of families associated to the genes on the left
  synRightFile          sequence name / list of families associated to the genes on the right

"""
import sys
import re
import random
import os
from Bio import Phylo
from io import StringIO
# from cStringIO import StringIO
import time
from ete3 import Phyloxml, phyloxml

from lxml import etree
from xml.etree import ElementTree
from xml.dom import minidom

from lxml.etree import XMLParser, parse

from docopt import docopt
import sqlite3
import zlib
import base64

args = docopt(__doc__)

tree = args["<treeFile>"]
speciesDico = args["<speciesDicoFile>"]
synLeft = args["<syntenyLeftFile>"]
synRight = args["<syntenyRightFile>"]
arguments = docopt(__doc__, version='1.0.0rc2')

sys.setrecursionlimit(15000)
print(sys.getrecursionlimit())

def loadDico(fileDico):
    file = open(fileDico,"r")
    speciesDico = {}
    for line in file:
        tline = re.split('\|',line)
        cds=tline[0].strip(' ')
        spec=tline[-1].strip(' \n')
        speciesDico[cds]=spec
    file.close()
    return speciesDico

def loadSynte(fileSynte):
    file = open(fileSynte,"r")
    synteDico = {}
    for line in file:
        tline = re.split(' ',line)
        cds=tline[0].strip('\n')
        i = 1
        synte = ""
        while (i < len(tline)) :
            synte = synte + " "+tline[i]
            #print "debug"+tline[i]
            if 'seqdefdico' in globals():
                tcds = re.split(":",tline[i])
                if (len(tcds) > 1):
                    #print tcds[1]
                    if tcds[1] in seqdefdico:
                        #print seqdefdico[tcds[1]]
                        seqdef=seqdefdico[tcds[1]]

                        synte=synte +"|"+seqdef.replace(" ","_")
                        print("SYNTE="+synte)
            i = i + 1
        synteDico[cds]=synte.strip(' \n')
    file.close()
    return synteDico

def createPhyloXML(fam,newick):
    handle = StringIO(newick)
    trees = Phylo.read(handle, 'newick')
    rd = str(random.randint(0,1000))
    Phylo.write([trees], 'tmpfile-'+rd+'.xml', 'phyloxml')
    file = open('tmpfile-'+rd+'.xml', 'r')
    text = file.read()
    file.close()
    os.remove('tmpfile-'+rd+'.xml')
    p = XMLParser(huge_tree=True)
    text = text.replace("phy:", "")
    text = re.sub("b'([^']*)'", "\\1", text)
    text = re.sub('branch_length_attr="[^"]+"', "", text)
    header = "<phyloxml>"

    text = re.sub('<phyloxml[^>]+>', header, text)
    text = text.replace('Phyloxml', 'phyloxml')
    tree = etree.fromstring(text,parser=p)
    # ajout du nom d'arbre
    treename = etree.Element("name")
    treename.text = fam
    ins = tree.find('phylogeny')
    ins.append(treename)

    clade = tree.xpath("/phyloxml/phylogeny/clade")
    subtree = tree.xpath("/phyloxml")
    nbfeuille = 0
    famspecies = {}

    for element in clade[0].iter('clade'):
        enom=element.find('name')
        if (enom is not None) :
            nbfeuille = nbfeuille + 1
            cds = enom.text
            sp = dico.get(cds)
            if (not  sp):
                print ("undefined species for "+ cds)
                sp = "undefined"
            famspecies[sp] = 1
            synLeft=synteLeftDico.get(cds)
            if (not synLeft):
                synLeft = ""
            synRight=synteRightDico.get(cds)
            if (not synRight):
                synRight = ""

            evrec = etree.Element("eventsRec")
            leaf = etree.Element("leaf")
            if 'crossdico' in globals():
                crossref = etree.Element("crossref")
                if cds in crossdico[0]:
                    response = crossdico[0][cds]
                    tabbuf = response.split("|")
                    print (tabbuf)
                    for buf in tabbuf:
                        print ("process "+buf)
                        tabcross = buf.split(":")
                        if len(tabcross) > 1 :
                            crossref.set(tabcross[0], tabcross[1])
            leaf.set('speciesLocation', sp)
            if 'seqdefdico' in globals():
                if cds in seqdefdico:
                    leaf.set('defintiion', seqdefdico[cds])
            leaf.set('syntenyLeft', synLeft)
            leaf.set('syntenyRight', synRight)
            if 'crossdico' in globals():
                leaf.append(crossref)
            evrec.append(leaf)
            element.append(evrec)
    print ("Number of leaves : ")
    print (nbfeuille)
    nbspecies = len(famspecies)
    print ("Number of species : ")
    print (nbspecies)

    treesize =  etree.Element("size")
    treesize.set('leaves',str(nbfeuille))
    treesize.set('species',str(nbspecies))
    e=subtree[0].find('phylogeny')
    e.append(treesize)
    text =  minidom.parseString(ElementTree.tostring(subtree[0])).toprettyxml()
    # remove blank lines
    cleantext = "\n".join([ll.rstrip() for ll in text.splitlines() if ll.strip()])
    return cleantext

print ("Loading species...")
dico = loadDico(speciesDico)
print ("OK")

print ("Loading synteny left... ")
synteLeftDico =  loadSynte(synLeft)
print ("OK")

print ("Loading synteny right... ")
synteRightDico =  loadSynte(synRight)
print ("OK")

#Creates empty phyloxml document
#project = Phyloxml()   a decommenter si on veut un fichier xml unique

# Loads newick tree
treefile = open(tree,"r")
xmlouputfile = open(tree+".xml","w")


for line in treefile:
    tline = re.split(' ',line)
    newick=tline[1]
    fam=tline[0]
    phyloxmltree = createPhyloXML(fam,newick)
    xmlouputfile.write(phyloxmltree)
    print ("Famille "+fam+" OK")
