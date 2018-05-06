#!/usr/bin/python
## Only works for python2, not Python3 - unsure why????

import csv
import gffutils
from Bio.Seq import Seq
from Bio.Alphabet import generic_dna

HEP_HighMeth = open('/Users/rudyfabunan/Documents/Final_DataSets/MethylationAnalysis/Gene/GeneBody_Analysis/Genes_HEP_CG_HighlyMethylatedGenes.txt','rU')
# LEP_HighMeth = open('/Users/rudyfabunan/Documents/Final_DataSets/MethylationAnalysis/Gene/GeneBody_Analysis/Genes_LEP_CG_HighlyMethylatedGenes.txt','rU')
# HEP_LowMeth = open('/Users/rudyfabunan/Documents/Final_DataSets/MethylationAnalysis/Gene/GeneBody_Analysis/Genes_HEP_CG_LowMethylatedGenes.txt','rU')
# LEP_LowMeth = open('/Users/rudyfabunan/Documents/Final_DataSets/MethylationAnalysis/Gene/GeneBody_Analysis/Genes_LEP_CG_LowMethylatedGenes.txt','rU')

HEP_HM_r = HEP_HighMeth.readlines()
# LEP_HM_r = LEP_HighMeth.readlines()
# HEP_LM_r = HEP_LowMeth.readlines()
# LEP_LM_r = LEP_LowMeth.readlines()

HEP_fa = 'HEP_recode2.fa'
# LEP_fa = 'LEP_recode2.fa'

out_HEP = open('HEP_out.txt','w')
# out_LEP = open('LEP_out.txt','w')
oHEP = csv.writer(out_HEP, delimiter = '\t')
# oLEP = csv.writer(out_LEP, delimiter = '\t')
# out_HEP_l = open('HEPl_out.txt','w')
# out_LEP_l = open('LEPl_out.txt','w')
# oHEPl = csv.writer(out_HEP_l, delimiter = '\t')
# oLEPl = csv.writer(out_LEP_l, delimiter = '\t')

## Create the database:
## Will create another file in your directory titled 'Araport_genes'
## Once this is made, you do not have to recreate it

# db = gffutils.create_db('Araport11_genes.20151202.gff3','Araport_genes')


## Once the db is created, you can call it like so
db = gffutils.FeatureDB('Araport_genes')

## Create empty list for each populaton
## Empty list will contain the gene name, counts of CG in the exons, and the length of the concatenated exons
exons_HEP = []
exons_LEP = []


## Function to extract features according to GFF from the fasta file
## This function takes in three arguments:
## type - the fature needed (ex. exons, etc)
## pop_fasta - the fasta file
## list - the empty list to append
def get_features(type,pop_fasta,list):
    for t in db.features_of_type('mRNA', order_by='start'): ## The GFF is heirachal and one of the top levels in my GFF is mRNA
        Exons_only = '' ## empty string which will be used for exon sequences
        for i in db.children(t, featuretype=type, order_by='start'): ## Exon features are a lower level within mRNA and this is to go through all the exons per gene
            seq = i.sequence(pop_fasta, use_strand=True) ## Extracting the sequence feature from the fasta file
            Exons_only += seq ## adding on to the string 'Exons_only'
        Exons_only = Seq(Exons_only, generic_dna) ## making the string into a Seq class which has string methods. This case we are using generic DNA, but if I wanted to get CHG or CHH sites I would use ambiguous_dna instead of generic_dna
        list.append([t.id,Exons_only.count('CG'),len(Exons_only)]) ## appending the Seq class string to the list so I can iterate through it later

## Function to make a table and output into a csv file
## This function takes in three arguments
## input - this is file with genes I want to look for in the list
## output - This is the output file
## list - The list to iterate through to find the genes from the input file

def maketable(input,output,list):
    for gene in input:
        gene_name = gene.strip('\n') + '.1' ## The gff contains multiple splice forms for each gene. I only want the first one so I need to add the '.1' at the end of the gene name
        for i in list:
            if gene_name == i[0]:
                output.writerow(i)

## Calling the functions above
get_features('exon',HEP_fa,exons_HEP)
# get_features('exon',LEP_fa,exons_LEP)

maketable(HEP_HM_r,oHEP,exons_HEP)
# maketable(LEP_HM_r,oLEP,exons_LEP)
# maketable(HEP_LM_r,oHEPl,exons_HEP)
# maketable(LEP_LM_r,oLEPl,exons_LEP)
#
## closing the files
out_HEP.close()
# out_LEP.close()
# out_HEP_l.close()
# out_LEP_l.close()
