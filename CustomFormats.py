# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 20:31:93 2017

@author: Sunghee Woo
"""


def initForwardCodon():
    return {"TTT": "F", "TTC": "F", "TTA": "L", "TTG": "L",
            "TCT": "S", "TCC": "S", "TCA": "S", "TCG": "S",
            "TAT": "Y", "TAC": "Y", "TAA": "X", "TAG": "X",
            "TGT": "C", "TGC": "C", "TGA": "X", "TGG": "W",
            "CTT": "L", "CTC": "L", "CTA": "L", "CTG": "L",
            "CCT": "P", "CCC": "P", "CCA": "P", "CCG": "P",
            "CAT": "H", "CAC": "H", "CAA": "Q", "CAG": "Q",
            "CGT": "R", "CGC": "R", "CGA": "R", "CGG": "R",
            "ATT": "I", "ATC": "I", "ATA": "I", "ATG": "M",
            "ACT": "T", "ACC": "T", "ACA": "T", "ACG": "T",
            "AAT": "N", "AAC": "N", "AAA": "K", "AAG": "K",
            "AGT": "S", "AGC": "S", "AGA": "R", "AGG": "R",
            "GTT": "V", "GTC": "V", "GTA": "V", "GTG": "V",
            "GCT": "A", "GCC": "A", "GCA": "A", "GCG": "A",
            "GAT": "D", "GAC": "D", "GAA": "E", "GAG": "E",
            "GGT": "G", "GGC": "G", "GGA": "G", "GGG": "G",
            }


def initReverseCodon():
    return {"AAA": "F", "AAG": "F", "AAT": "L", "AAC": "L",
            "AGA": "S", "AGG": "S", "AGT": "S", "AGC": "S",
            "ATA": "Y", "ATG": "Y", "ATT": "X", "ATC": "X",
            "ACA": "C", "ACG": "C", "ACT": "X", "ACC": "W",
            "GAA": "L", "GAG": "L", "GAT": "L", "GAC": "L",
            "GGA": "P", "GGG": "P", "GGT": "P", "GGC": "P",
            "GTA": "H", "GTG": "H", "GTT": "Q", "GTC": "Q",
            "GCA": "R", "GCG": "R", "GCT": "R", "GCC": "R",
            "TAA": "I", "TAG": "I", "TAT": "I", "TAC": "M",
            "TGA": "T", "TGG": "T", "TGT": "T", "TGC": "T",
            "TTA": "N", "TTG": "N", "TTT": "K", "TTC": "K",
            "TCA": "S", "TCG": "S", "TCT": "R", "TCC": "R",
            "CAA": "V", "CAG": "V", "CAT": "V", "CAC": "V",
            "CGA": "A", "CGG": "A", "CGT": "A", "CGC": "A",
            "CTA": "D", "CTG": "D", "CTT": "E", "CTC": "E",
            "CCA": "G", "CCG": "G", "CCT": "G", "CCC": "G",
            }


def initChr():
    return ['chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7',
            'chr8', 'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14',
            'chr15', 'chr16', 'chr17', 'chr18', 'chr19', 'chr20',
            'chr21', 'chr22', 'chr23', 'chrX', 'chrY']


def initMutationData():
    return {'chr1': {}, 'chr2': {}, 'chr3': {}, 'chr4': {}, 'chr5': {},
            'chr6': {}, 'chr7': {}, 'chr8': {}, 'chr9': {}, 'chr10': {},
            'chr11': {}, 'chr12': {}, 'chr13': {}, 'chr14': {}, 'chr15': {},
            'chr16': {}, 'chr17': {}, 'chr18': {}, 'chr19': {}, 'chr20': {},
            'chr21': {}, 'chr22': {}, 'chr23': {}, 'chrX': {}, 'chrY': {}}


def formatChrNameForHuman(reference_name_str):
    if reference_name_str == "1":
        return "chr1"
    elif reference_name_str == "2":
        return "chr2"
    elif reference_name_str == "3":
        return "chr3"
    elif reference_name_str == "4":
        return "chr4"
    elif reference_name_str == "5":
        return "chr5"
    elif reference_name_str == "6":
        return "chr6"
    elif reference_name_str == "7":
        return "chr7"
    elif reference_name_str == "8":
        return "chr8"
    elif reference_name_str == "9":
        return "chr9"
    elif reference_name_str == "10":
        return "chr10"
    elif reference_name_str == "11":
        return "chr11"
    elif reference_name_str == "12":
        return "chr12"
    elif reference_name_str == "13":
        return "chr13"
    elif reference_name_str == "14":
        return "chr14"
    elif reference_name_str == "15":
        return "chr15"
    elif reference_name_str == "16":
        return "chr16"
    elif reference_name_str == "17":
        return "chr17"
    elif reference_name_str == "18":
        return "chr18"
    elif reference_name_str == "19":
        return "chr19"
    elif reference_name_str == "20":
        return "chr20"
    elif reference_name_str == "21":
        return "chr21"
    elif reference_name_str == "22":
        return "chr22"
    elif reference_name_str == "23":
        return "chr23"
    elif reference_name_str == "X":
        return "chrX"
    elif reference_name_str == "Y":
        return "chrY"
    elif reference_name_str.startswith("chr"):
        return reference_name_str
    else:
        return -1


def initDnaFaFile(chr_str):
    return {'chr1': 'Homo_sapiens.GRCh38.dna.chromosome.1.fa', 
    'chr2': 'Homo_sapiens.GRCh38.dna.chromosome.2.fa', 
    'chr3': 'Homo_sapiens.GRCh38.dna.chromosome.3.fa', 
    'chr4': 'Homo_sapiens.GRCh38.dna.chromosome.4.fa', 
    'chr5': 'Homo_sapiens.GRCh38.dna.chromosome.5.fa', 
    'chr6': 'Homo_sapiens.GRCh38.dna.chromosome.6.fa', 
    'chr7': 'Homo_sapiens.GRCh38.dna.chromosome.7.fa', 
    'chr8': 'Homo_sapiens.GRCh38.dna.chromosome.8.fa', 
    'chr9': 'Homo_sapiens.GRCh38.dna.chromosome.9.fa', 
    'chr10': 'Homo_sapiens.GRCh38.dna.chromosome.10.fa', 
    'chr11': 'Homo_sapiens.GRCh38.dna.chromosome.11.fa', 
    'chr12': 'Homo_sapiens.GRCh38.dna.chromosome.12.fa', 
    'chr13': 'Homo_sapiens.GRCh38.dna.chromosome.13.fa', 
    'chr14': 'Homo_sapiens.GRCh38.dna.chromosome.14.fa', 
    'chr15': 'Homo_sapiens.GRCh38.dna.chromosome.15.fa', 
    'chr16': 'Homo_sapiens.GRCh38.dna.chromosome.16.fa', 
    'chr17': 'Homo_sapiens.GRCh38.dna.chromosome.17.fa', 
    'chr18': 'Homo_sapiens.GRCh38.dna.chromosome.18.fa', 
    'chr19': 'Homo_sapiens.GRCh38.dna.chromosome.19.fa', 
    'chr20': 'Homo_sapiens.GRCh38.dna.chromosome.20.fa', 
    'chrX': 'Homo_sapiens.GRCh38.dna.chromosome.X.fa', 
    'chrY': 'Homo_sapiens.GRCh38.dna.chromosome.Y.fa'}


def initGenePFile(chr_str):
    return {'chr1': 'genes_mu_chr1_test.p', 'chr2': 'genes_mu_chr2.p', 'chr3': 'genes_mu_chr3.p', 'chr4': 'genes_mu_chr4.p', 'chr5': 'genes_mu_chr5.p', 'chr6': 'genes_mu_chr6.p', 'chr7': 'genes_mu_chr7.p', 'chr8': 'genes_mu_chr8.p', 'chr9': 'genes_mu_chr9.p', 'chr10': 'genes_mu_chr10.p', 'chr11': 'genes_mu_chr11.p', 'chr12': 'genes_mu_chr12.p', 'chr13': 'genes_mu_chr13.p', 'chr14': 'genes_mu_chr14.p', 'chr15': 'genes_mu_chr15.p', 'chr16': 'genes_mu_chr16.p', 'chr17': 'genes_mu_chr17.p', 'chr18': 'genes_mu_chr18.p', 'chr19': 'genes_mu_chr19.p', 'chr20': 'genes_mu_chr20.p', 'chrX': 'genes_mu_chrX.p', 'chrY': 'genes_mu_chrY.p'}


def initGFFFile(chr_str):
    return {'chr1': 'Homo_sapiens.GRCh38.88.chr1.gff3', 'chr2': 'Homo_sapiens.GRCh38.91.chromosome.2.chr.gff3', 'chr3': 'Homo_sapiens.GRCh38.91.chromosome.3.chr.gff3', 'chr4': 'Homo_sapiens.GRCh38.91.chromosome.4.chr.gff3', 'chr5': 'Homo_sapiens.GRCh38.91.chromosome.5.chr.gff3', 'chr6': 'Homo_sapiens.GRCh38.91.chromosome.6.chr.gff3', 'chr7': 'Homo_sapiens.GRCh38.91.chromosome.7.chr.gff3', 'chr8': 'Homo_sapiens.GRCh38.91.chromosome.8.chr.gff3', 'chr9': 'Homo_sapiens.GRCh38.91.chromosome.9.chr.gff3', 'chr10': 'Homo_sapiens.GRCh38.91.chromosome.10.chr.gff3', 'chr11': 'Homo_sapiens.GRCh38.91.chromosome.11.chr.gff3', 'chr12': 'Homo_sapiens.GRCh38.91.chromosome.12.chr.gff3', 'chr13': 'Homo_sapiens.GRCh38.91.chromosome.13.chr.gff3', 'chr14': 'Homo_sapiens.GRCh38.91.chromosome.14.chr.gff3', 'chr15': 'Homo_sapiens.GRCh38.91.chromosome.15.chr.gff3', 'chr16': 'Homo_sapiens.GRCh38.91.chromosome.16.chr.gff3', 'chr17': 'Homo_sapiens.GRCh38.91.chromosome.17.chr.gff3', 'chr18': 'Homo_sapiens.GRCh38.91.chromosome.18.chr.gff3', 'chr19': 'Homo_sapiens.GRCh38.91.chromosome.19.chr.gff3', 'chr20': 'Homo_sapiens.GRCh38.91.chromosome.20.chr.gff3', 'chrX': 'Homo_sapiens.GRCh38.91.chromosome.X.chr.gff3', 'chrY': 'Homo_sapiens.GRCh38.91.chromosome.Y.chr.gff3'}


def writeListToFile(self, outFile, in_list):
    write_list = '\t'.join(in_list) + '\n'
    outFile.write(write_list)
