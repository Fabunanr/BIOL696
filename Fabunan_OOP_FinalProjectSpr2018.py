## RENAME this file YourLastName_OOP_FinalProjectSpr2018.py

##Assignment: Add to the constructor and methods of a base class and superclasses
##            which inherit the base class properties.

## Begin with the base Seq class and the DNA class we created in lecture below.

### Seq Class
#
#  Constructor:
#  (1) Use the string functions upper and string to clean up self.sequence.
#  (2) Add a variable self.kmers to the constructor and make it equal to an empty list.

#  Methods:
#  (1) Add a method called make_kmers that makes kmers of a given length from self.sequence
#      appends these to self.kmers. Default kmer paramter=3.
#  (2) Add a method called fasta that returns a fasta formatted string like this:
#      >species gene
#      AGATTGATAGATAGATAT

class Seq:

    def __init__(self,sequence,gene,species):
        self.sequence=sequence.strip(' ').upper()
        self.gene=gene
        self.species=species
        self.kmers=[]

    def __str__(self):
        return self.sequence

    def print_record(self):
        print(self.species + " " + self.gene + ": " + self.sequence)

    def make_kmers(self, k=3):
        n= len(self.sequence)
        for i in range(0, n-k+1):
           self.kmers.append(self.sequence[i:i+k])
        return self.kmers

    def fasta(self):
        return '>' + self.species + ' ' + self.gene + '\n' + self.sequence

import re

## DNA Class: INHERITS Seq class

 # Constructor:
 # Use re.sub to change any non nucleotide characters in self.sequence into an 'N'.
 #     re.sub('[^ATGCU]','N',sequence) will change any character that is not a
 #     capital A, T, G, C or U into an N. (Seq already uppercases and strips.)
 #
 # Methods:
 # Add a method called print_info that is like print_record, but adds geneid and an
 #     empty space to the beginning of the string.


class DNA(Seq):

    def __init__(self,sequence,gene,species,geneid,**kwargs):
        super().__init__(sequence,gene,species)
        self.sequence=re.sub('[^ATGCU]','N',self.sequence)
        self.geneid=geneid
        # self.sequence=re.sub('[^ATGCU]','N',sequence)

    def analysis(self):
        gc=len(re.findall('G',self.sequence) + re.findall('C',self.sequence))
        return gc

    def print_info(self):
       print(self.geneid + ' ' + self.species + " " + self.gene + ": " + self.sequence)

## RNA Class:  INHERITS DNA class

 # Construtor:
 # Use the super() function (see DNA Class example).
 # (1) Automatically change all Ts to Us in self.sequence.
 # (2) Add self.codons equals to an empty list
 #
 # Methods:
 # (1) Add make_codons which breaks the self.sequence into 3 letter codons
 #     and appends these codons to self.codons unless they are less than 3 letters long.
 # (2) Add translate which uses the Global Variable standard_code below to
 #     translate the codons in self.codons and returns a protein sequence.

standard_code = {
     "UUU": "F", "UUC": "F", "UUA": "L", "UUG": "L", "UCU": "S",
     "UCC": "S", "UCA": "S", "UCG": "S", "UAU": "Y", "UAC": "Y",
     "UAA": "*", "UAG": "*", "UGA": "*", "UGU": "C", "UGC": "C",
     "UGG": "W", "CUU": "L", "CUC": "L", "CUA": "L", "CUG": "L",
     "CCU": "P", "CCC": "P", "CCA": "P", "CCG": "P", "CAU": "H",
     "CAC": "H", "CAA": "Q", "CAG": "Q", "CGU": "R", "CGC": "R",
     "CGA": "R", "CGG": "R", "AUU": "I", "AUC": "I", "AUA": "I",
     "AUG": "M", "ACU": "T", "ACC": "T", "ACA": "T", "ACG": "T",
     "AAU": "N", "AAC": "N", "AAA": "K", "AAG": "K", "AGU": "S",
     "AGC": "S", "AGA": "R", "AGG": "R", "GUU": "V", "GUC": "V",
     "GUA": "V", "GUG": "V", "GCU": "A", "GCC": "A", "GCA": "A",
     "GCG": "A", "GAU": "D", "GAC": "D", "GAA": "E", "GAG": "E",
     "GGU": "G", "GGC": "G", "GGA": "G", "GGG": "G"}

class RNA(DNA):

    def __init__(self,sequence,gene,species,geneid,**kwargs):
        super().__init__(sequence,gene,species,geneid)
        self.sequence = re.sub('T','U',self.sequence)
        self.codons = []

    def make_codons(self):
        x = re.findall('.{1,3}',self.sequence)
        for i in x:
            if len(i) == 3:
                self.codons.append(i)
            else:
                pass

    def translate(self):
        seq = ''
        for i in self.codons:
            try:
                AA = standard_code[i]
                seq += AA
            except KeyError:
                AA = 'X'
                seq += AA
        print(seq)

# ## Protein Class: INHERITS Seq class

 # Construtor:
 # Use the super() function (see DNA Class example).
 # (1) Use re.sub to change any non LETTER characters in self.sequence into an 'X'.
 # (2) self.aa_counts has already been added - it is a dictionary with all 21
 #     amino acids as keys and values set to 0
 #
 # Methods:
 # Add tabulate_amino_acids, which counts the amino acids in self.sequence
 #     and puts these in the dictionary self.counts

class Protein(Seq):
    def __init__(self,sequence,gene,species,geneid,**kwargs):
        super().__init__(sequence,gene,species)
        self.sequence =  re.sub('[^A-Z]','X',self.sequence.rstrip())
        self.aa_counts={'A':0,'C':0,'D':0,'E':0,'F':0,'G':0,'H':0,'I':0,'K':0,'L':0,'M':0,'P':0,'Q':0,'R':0,'S':0,'T':0,'U':0,'V':0,'W':0,'X':0,'Y':0}

    def tabulate_amino_acids(self):
        for aa in self.sequence:
            if aa in self.aa_counts:
                self.aa_counts[aa] +=1
            else:
                self.aa_counts[aa] = 1
        return(self.aa_counts)
