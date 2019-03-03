#!/usr/bin/env python3

import cairo
import argparse
from re import finditer
from Bio import SeqIO

# Will contain a list of all motifs
motifs = []

# Used to generate a list of all regular expressions
regex = []

# Contains a dictionary of regular expressions for each nucleotide
nucleotide_regex = {"t":"[Tt]", "a":"[Aa]", "c":"[Cc]", "g":"[Gg]", "u":"[Tt]", "w":"[AaTt]", "s":"[CcGg]", 
                    "m":"[AaCc]", "k":"[GgTt]", "r":"[AaGg]", "y":"[CcTt]", "b":"[CcGgTt]", "d":"[AaGgTt]", 
                    "h":"[AaCcTt]", "v":"[AaCcGg]", "n":"[AaCcTtGg]"}

# Will contain all sequences to be identified with their name
seq_dict = {}

# Contains RGB values
# Red, Orange, Yellow, Green, Cyan, Blue, Violet, Purple-Red, Light Gray, Dark Gray       
RGBs = [[0.8,0,0],[0.8,0.4,0],[0.8,0,0.8],[0,0.8,0],[0,0.8,0.8],
        [0,0,0.8],[0.6,0,0.6],[0.6,0,0.3],[0.25,0.25,0.25],[0.5,0.5,0.5]]


def get_arguments():
	'''Labels motifs on a given exon and associated introns.'''
	parser = argparse.ArgumentParser(description = "Labels motifs on a given exon and associated introns.")
	parser.add_argument("-m", "--motifs", help = "Motif txt file. Max of 10 motifs must be separated with a space or newline.",\
		required=False, type=str, default = "motifs.txt")
	parser.add_argument("-f", "--fasta", help= "Fasta file with one exon per entry. Exon must be in uppercase, intron in lowercase.",\
		required=False, type=str, default = "INSR.fasta")
	parser.add_argument("-s", "--scaling", help = "Sets scaling factor. Image can handle 780 nucleotides when scaling is one. Default scaling is 3.",\
		required=False, type=str, default = 3)
	parser.add_argument("-o", "--output", help= "Name of outputted SVG file. Default is 'output.svg'",\
		required=False, type=str, default = "output.svg")
	return parser.parse_args()

args = get_arguments()

# Argparse functions will eventually go here
motifs_file = args.motifs
fasta_file = args.fasta
title = args.output
scaling = args.scaling

def getMotifs(motifs_file):
    '''This function generates a list of all desired motifs from a text file.'''
    with open(motifs_file, "r") as fh:
        for line in fh:
            line = line.strip("\n").split(" ")
            for item in line:
                if item != "":
                    motifs.append(item.lower())
            
def getRegex(motif):
    '''Returns a regular expressions (regex) for motif marking.'''
    expression = ""
    for nucleotide in range(len(motif)):
        expression += nucleotide_regex[motif[nucleotide]]
    regex.append(expression)
    
def parseFasta(fasta_file):
    '''Adds Fasta file to dictionary. Query is the key and sequence is the value.'''
    fasta_sequences = SeqIO.parse(open(fasta_file),'fasta')
    for fasta in fasta_sequences:
        description, sequence = fasta.description, str(fasta.seq)
        seq_dict[description] = sequence

def getTrons(sequence):
    '''Returns the positions of the end of the first intron, end of the exon, and end of the second intron.'''
    exon_nucleotides = "[ACTG]+"
    for match in finditer(exon_nucleotides, sequence):
        intron1 = match.start()
        exon = match.end()
        intron2 = len(sequence)
    return intron1, exon, intron2
    
def getMotifPositions(query, sequence):
    '''Returns position of the motif to be drawn'''
    start = []
    length = []
    for match in finditer(query, sequence):
        start.append(match.start())
        length.append(match.end() - match.start())
    return start, length

def getNumMotifsPerTron(intron1, exon, intron2, query, sequence):
    '''Returns the number of matches per intron/exon'''
    intron1Matches = 0
    exonMatches = 0
    intron2Matches = 0
    for match in finditer(query, sequence[0:intron1]):
        intron1Matches += 1
    for match in finditer(query, sequence[intron1:exon]):
        exonMatches += 1
    for match in finditer(query, sequence[exon:intron2]):
        intron2Matches += 1
    return intron1Matches, exonMatches, intron2Matches

def drawMotifs(regex, seq_dict, title):
    '''Draws the introns, exon, and sequence motifs.'''
    width = 800
    height = (len(seq_dict.keys()) // 5 + 1) * 800
    surface = cairo.SVGSurface(title, width, height)
    ctx = cairo.Context(surface) 
    ctx.select_font_face("Courier", cairo.FONT_SLANT_NORMAL, cairo.FONT_WEIGHT_NORMAL)
    ctx.set_line_width(1)
    loop = -1
    for key in seq_dict.keys():
        loop += 1
        ctx.set_source_rgb(0,0,0)
        ctx.set_font_size(13)
        ctx.move_to(20, 30 + (loop * 180))
        ctx.show_text(key)
        intron1, exon, intron2 = getTrons(seq_dict[key])
        ctx.move_to(20, 80 + (loop * 180))
        ctx.line_to(20 + (intron2 // scaling), 80 + (loop * 180))  #(x,y)
        ctx.stroke()
        ctx.rectangle(20 + (intron1 // scaling), 70 + (loop * 180), 
                     ((exon - intron1) // scaling), 20)   #(x,y,width,height)
        ctx.fill()
        motif_color = -1
        for query in regex:
            motif_color += 1
            start, length = getMotifPositions(query, seq_dict[key])
            firstPass = True
            for pos in range(len(start)):
                ctx.set_source_rgb(RGBs[motif_color][0], RGBs[motif_color][1], RGBs[motif_color][2])
                ctx.rectangle(20 + (start[pos] // scaling), 40 + (loop * 180), 
                              (length[pos] // scaling), 20)
                ctx.fill()
                if firstPass == True: 
                    ctx.set_font_size(10)
                    ctx.move_to(20, 100 + (loop * 180) + (motif_color * 10))
                    I1Matches, EMatches, I2Matches = getNumMotifsPerTron(intron1, exon, intron2, query, seq_dict[key])
                    string1 = motifs[motif_color] + " - First Intron: " + str(I1Matches) + ", Exon: " + str(EMatches) 
                    string2 = string1 + ", Second Intron: " + str(I2Matches)
                    ctx.show_text(string2)
                    firstPass = False
            
def main():
    '''Draws 2 introns and 1 exon from FASTA file, and marks its motifs from a text file.'''
    getMotifs(motifs_file)
    for motif in motifs:
        getRegex(motif)
    parseFasta(fasta_file)
    drawMotifs(regex, seq_dict, title)
    
main()
