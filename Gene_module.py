# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 20:04:58 2017

@author: Sunghee Woo
"""
import os
# import sys
# import gc
# import jsonpickle
import cPickle as pickle
import Graph_module  # , trie, tree, parameter, tree_dna
import Parameters as param
import CustomFormats as formats
from os.path import dirname  # , abspath
script_dir = dirname(__file__)
# sys.setrecursionlimit(30000)


class ConstructGraphPerChr(object):

    def __init__(self, chr_str):
        print ('Processing:', chr_str)
        self.GGraph = Graph_module.GGraph()
        self.GGraph.addNode(-1, 's', 0, 'source')
        self.chr = ''
        self.seq = ''
        self.chr_str = chr_str
        self.dna_fa_file_dic = formats.initDnaFaFile(chr_str)
        self.gene_p_file_dic = formats.initGenePFile(chr_str)
        self.gff_file_dic = formats.initGFFFile(chr_str)
        self.fasta = os.path.join(
            script_dir, 'dna', self.dna_fa_file_dic[chr_str])
        self.p_gene = os.path.join(
            script_dir, 'pickle', self.gene_p_file_dic[chr_str])
        self.gff = os.path.join(script_dir, 'gff', self.gff_file_dic[chr_str])
        self.vcf = os.path.join(script_dir, 'vcf', param.vcf_file_name())
        self.debug_trans_seq = ''
        self.Process()

    def Process(self):
        print (' -Reading reference DNA')
        self.parseDNA()
        print (' -Creating transcript graph from GFF')
        self.parseGFF()
        print (' -Creating variant graph fvrom VCF')
        self.parseVCF()
        print (' -Merging graph')
        self.GGraph.mergeGraph()
        self.debugPrint()
        # print (' -Pickling')
        # pickle.dump(self.GGraph, open(self.p_gene, "wb"))

    def parseDNA(self):
        # , open("F:\\Ensembl\\dna_chr1.txt",'w') as o_fasta:
        with open(self.fasta, 'r') as in_fasta:
            header = in_fasta.readline()
            self.chr = header.strip().split(' ')[0]
            seqlist = in_fasta.read().splitlines()
            self.seq = ''.join(seqlist)
            #pickle.dump(self.seq, open("F:\\Ensembl\\dna_chr1.p", "wb"))
            # o_fasta.write(self.seq)
#            for line in in_fasta:
#                if line.startswith('>'):
#                    self.chr = line[1:].split(' ')[0]
#                else:
#                    self.seq += line.strip()

#    def createDNA(self):
#        self.tree = tree_dna.SuffixTree(self.prmtr)
#        self.tree.makeST(self.tree, self.seq)

    def parseGFF(self):
        in_transcript = False
        trans_id_cds = -1
        tmp_cnt = 0
#        if self.trie_flag:
#            self.tree = trie.SuffixTrie(self.prmtr)
#        else:
#            self.tree = tree.SuffixTree(self.prmtr)
        trans_seq = ''
        list_coor_tree = []
        with open(self.gff, 'r') as in_gff:
            if self.debug_trans_seq != '':
                debug_seq = open(self.debug_trans_seq, 'w')
            for line in in_gff:
                if line.startswith('#'):
                    continue
                data = line.split('\t')
                if data[0] != self.chr_str and 'chr'+data[0] != self.chr_str:
                    continue
                if data[2] == 'mRNA':
                    #chr_id = data[0]
                    tags = data[8].split(';')
                    for i in tags:
                        #if i.startswith('Parent=gene:'):
                        #    gene_id = i.replace('Parent=gene:','')
                        if i.startswith('ID=transcript:'):
                            trans_id = i.replace('ID=transcript:', '')
                    in_transcript = True
                    node_from = -1
                    tmp_cnt += 1
                    if list_coor_tree != []:
                        #self.tree.makeST(self.tree, trans_seq, list_coor_tree)
                        if self.debug_trans_seq != '':
                            #debug_seq.write(trans_seq+'\n'+','.join(str(x) for x in list_coor_tree)+'\n')
                            debug_seq.write(trans_id_cds+'\n'+trans_seq+'\n')
                        list_coor_tree = []
                        trans_seq = ''
                    continue
                elif data[2] == 'exon':
                    continue
                elif data[2].find('UTR') > -1:
                    continue
                elif data[2] == 'CDS' and in_transcript == True:
                    tags = data[8].split(';')
                    for i in tags:
                        #                        if i.startswith('Parent=gene:'):
                        #                            gene_id = i.replace('Parent=gene:','')
                        if i.startswith('Parent=transcript:'):
                            trans_id_cds = i.replace('Parent=transcript:', '')
                            break
                    if trans_id != trans_id_cds:
                        print ('gff3 parsing is going wrong',
                               trans_id, trans_id_cds)
                    start = int(data[3])-1
                    end = int(data[4])
                    seq = self.getSeq(start, end)
                    for i in range(0, len(seq)):
                        node_idx = start + i
                        # trans_id_cds)
                        self.GGraph.addNode(node_idx, seq[i], 0, 'REF')
                        self.GGraph.addEdge(node_from, node_idx, 0)
                        node_from = node_idx
                        list_coor_tree.append(node_idx)
                    trans_seq += seq
                else:
                    in_transcript = False
                    continue
#        pickle.dump(self.GGraph, open("F:\\Ensembl\\genes_chr1.p", "wb"))
#        pickle.dump(self.tree, open(self.p_tree, "wb"))
        if self.debug_trans_seq != '':
            debug_seq.close()
#        pickle.dump(self.GGraph, open(self.p_gene, "wb"))
#        frozen = jsonpickle.encode(self.GGraph)
        print ('num genes:', tmp_cnt)
        print ('G len ref:', self.GGraph.graph_node_len)

    def parseVCF(self):
        with open(self.vcf, 'r') as in_vcf:
            for line in in_vcf:
                if line.startswith('#'):
                    continue
                data = line.strip().split('\t')
                if data[0] != self.chr_str and 'chr'+data[0] != self.chr_str:
                    continue
#                chr_id = data[0]
                start = int(data[1])-1
                ref = data[3]
                mut_list = data[4].split(',')
                for mut in mut_list:
                    if len(ref) == 1:  # SNP / insertion
                        if len(mut) == 1:  # SNP
                            self.GGraph.addSNP(start, ref, mut, 0)
                        else:  # insertion
                            self.GGraph.addINS(start, ref, mut, 0)
                    elif len(ref) > 1:  # deletion / substitution
                        if len(mut) == len(ref):  # substitution
                            self.GGraph.addSUB(start, ref, mut, 0)
                        else:  # deletion
                            self.GGraph.addDEL(start, ref, mut, 0)
                    else:
                        print ('unexpected format in VCF:', line)

#        pickle.dump(self.GGraph, open(self.p_gene, "wb"))
        print ('G len:', self.GGraph.graph_node_len)
        print ('SNP:', self.GGraph.SNP)
        print ('SUB:', self.GGraph.SUB)
        print ('INS:', self.GGraph.INS)
        print ('DEL:', self.GGraph.DEL)

    def getSeq(self, start, end):
        return self.seq[start:end]

    def debugPrint(self):
        with open('debug_log.txt', 'w') as dFile:
            for idx in self.GGraph.G:
                node = self.GGraph.G[idx]
                write_str = ''
                write_str += str(idx) + '\t'
                write_str += str(node.base) + '\t'
                write_str += str(node.start) + '\t'
                write_str += str(node.end) + '\t'
                write_str += str(node.type) + '\t'
                write_str += ','.join("::".join((str(k), str(v)))
                                      for k, v in node.from_edges.items()) + '\t'
                write_str += ','.join("::".join((str(k), str(v)))
                                      for k, v in node.to_edges.items()) + '\n\n'
                dFile.write(write_str)
            for node in self.GGraph.M:
                write_str = ''
                write_str += str(node.base) + '\t'
                write_str += str(node.start) + '\t'
                write_str += str(node.end) + '\t'
                write_str += str(node.type) + '\t'
                write_str += ','.join("::".join((str(k), str(v)))
                                      for k, v in node.from_edges.items()) + '\t'
                write_str += ','.join("::".join((str(k), str(v)))
                                      for k, v in node.to_edges.items()) + '\n\n'
                dFile.write(write_str)


class unpickleGraph(object):

    def __init__(self, chr_str):
        self.chr_str = chr_str
        self.dna_fa_file_dic = formats.initDnaFaFile(chr_str)
        self.gene_p_file_dic = formats.initGenePFile(chr_str)
        self.gff_file_dic = formats.initGFFFile(chr_str)
        self.fasta = os.path.join(
            script_dir, 'dna', self.dna_fa_file_dic[chr_str])
        self.p_gene = os.path.join(
            script_dir, 'pickle', self.gene_p_file_dic[chr_str])
        self.gff = os.path.join(script_dir, 'gff', self.gff_file_dic[chr_str])
        self.vcf = os.path.join(script_dir, 'vcf', param.vcf_file_name())
        #        self.prmtr = parameter.Parameter()
        print ('loading gene pickle')
        self.GGraph = pickle.load(open(self.p_gene, "rb"))
#        print ('loading tree pickle')
#        self.tree = pickle.load( open("F:\\Ensembl\\trie_chr1.p", "rb" ) )
        print ('pickle loaded')


class testConstructGraphPerChr(object):

    def __init__(self):
        #        self.prmtr = parameter.Parameter()
        self.GGraph = Graph_module.GGraph()
        self.GGraph.addNode(-1, 's', 0, 'source')
        self.Process()

    def Process(self):
        self.parseGFF()

    def parseGFF(self):
        cnt = 0
#        if self.trie_flag:
#            self.tree = trie.SuffixTrie(self.prmtr)
#        else:
#            self.tree = tree.SuffixTree(self.prmtr)
        trans_seq = ''
        list_coor_tree = []
        with open(self.gff, 'r') as in_gff:
            for line in in_gff:
                if line.startswith('A') or line.startswith('T') or line.startswith('C') or line.startswith('G'):
                    trans_seq = line.strip()
                else:
                    list_coor_tree = [int(x) for x in line.strip().split(',')]
                    self.tree.makeST(self.tree, trans_seq,
                                     list_coor_tree[:len(trans_seq)])
                    cnt += 1
                    if cnt % 50 == 0:
                        print (cnt, 'trans indexed')

if __name__ == '__main__':
    parse_obj = ConstructGraphPerChr('chr1')
