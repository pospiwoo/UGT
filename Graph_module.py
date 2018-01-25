# -*- coding: utf-8 -*-
'''
Created on Mon Jul 17 20:45:30 2017

@author: Sunghee Woo
'''
from collections import defaultdict


class GraphNode:

    def __init__(self, base_str, coverage, node_type):
        self.base = base_str
        self.cov = coverage
        self.from_edges = defaultdict(int)
        self.to_edges = defaultdict(int)
        self.visited = 0
        self.type = node_type


class GraphNodeIdx:

    def __init__(self, idx, base_str, coverage, node_type):
        self.start = idx
        self.end = idx + 1
        self.base = base_str
        self.cov = coverage
        self.from_edges = defaultdict(int)
        self.to_edges = defaultdict(int)
        self.visited = 0
        self.type = node_type


class GGraph:

    def __init__(self):
        self.G = {}
        self.M = []
        self.SNP = 0
        self.SUB = 0
        self.INS = 0
        self.DEL = 0
        self.graph_node_len = 0
        self.idx_to_delete = []

    def addNode(self, idx, base_str, coverage, node_type):
        try:
            idx = int(idx)
        except ValueError:
            return False
        if base_str == 'N':
            print ('this base is masked', idx, base_str)
        if idx in self.G:
            # print 'Node already exist, adding coverage to existing node',
            # idx, base_str, coverage
            self.G[idx].cov += coverage
        else:
            self.G[idx] = GraphNodeIdx(idx, base_str, coverage, node_type)
            self.graph_node_len += 1
        return True

    def addEdge(self, fro_idx, to_idx, cov):
        try:
            fro_idx = int(fro_idx)
            to_idx = int(to_idx)
        except ValueError:
            print ('Node indexes should be INT', fro_idx, to_idx)
            return False
        if fro_idx in self.G and to_idx in self.G:
            # Assign edge from previous node
            self.G[fro_idx].to_edges[to_idx] += cov
            # Assign edge to next node
            self.G[to_idx].from_edges[fro_idx] += cov
        else:
            print ('Nodes not exist', fro_idx, to_idx)
            return False
        return True

    def addMutation(self, base_str, coverage, node_type, fro_idx, to_idx, ref_coor, cov):
        node_inst = GraphNodeIdx(-ref_coor, base_str, coverage, node_type)
        try:
            fro_idx = int(fro_idx)
            to_idx = int(to_idx)
        except ValueError:
            print ('Node indexes should be INT', fro_idx, to_idx)
            return False
        node_inst.to_edges[to_idx] += cov  # Assign edge from previous node
        node_inst.from_edges[fro_idx] += cov  # Assign edge to next node
        # Assign REF node edges
        node_inst_from = self.G[fro_idx]
        node_inst_to = self.G[to_idx]
        virtual_coor = len(self.M)
        self.M.append(node_inst)
        if node_type == 'SNP':
            self.SNP += 1
        elif node_type == 'SUB':
            self.SUB += 1
        elif node_type == 'INS':
            self.INS += 1
        elif node_type == 'DEL':
            self.DEL += 1
        self.graph_node_len += 1
        node_inst_from.to_edges[-virtual_coor] += cov  # Assign edge from previous node
        node_inst_to.from_edges[-virtual_coor] += cov  # Assign edge to next node
        node_inst.start = -virtual_coor
        node_inst.end = -virtual_coor + len(base_str)
        return True

    def removeGraphNode(self, node_inst):
        node_inst.type += '_removed'
        for idx, cov in node_inst.from_edges.iteritems():
            from_node = self.G[idx]
            from_node.to_edges.pop(node_inst.idx, None)
        for idx, cov in node_inst.to_edges.iteritems():
            to_node = self.G[idx]
            to_node.from_edges.pop(node_inst.idx, None)

    def merge2Nodes(self, node_idx_1, node_idx_2):
        node_inst_1 = self.G[node_idx_1]
        node_inst_2 = self.G[node_idx_2]
        node_inst_1.base += node_inst_2.base
        node_inst_1.end = node_inst_2.start + 1
        node_inst_1.to_edges = node_inst_2.to_edges
        # node_inst_2.to_edges = {} # Leave this for tracking
        node_inst_2.from_edges = {}
        node_inst_2.base = 'N'
        node_inst_2.type = 'X'
        self.idx_to_delete.append(node_idx_2)

    def mergeGraph(self):
        for idx in self.G:
            node = self.G[idx]
            if node.type != 'REF' or len(node.to_edges) != 1:
                continue
            idx_next = next(iter(node.to_edges.keys()))
            next_node = self.G[idx_next]
            prev_idx = idx
            while 1:
                if next_node.type != 'REF' or len(next_node.from_edges) != 1 \
                        or idx_next - prev_idx != 1:
                    break
                self.merge2Nodes(idx, idx_next)
                # print len(next_node.to_edges), next_node.to_edges
                if len(next_node.to_edges) != 1:
                    break
                prev_idx = idx_next
                idx_next = next(iter(next_node.to_edges.keys()))
                next_node = self.G[idx_next]
        # Delete merged leftover nodes
        for idx in self.idx_to_delete:
            self.G.pop(idx, None)
        self.idx_to_delete = []

    def removeEdge(self, fro_idx, to_idx):
        try:
            fro_idx = int(fro_idx)
            to_idx = int(to_idx)
        except ValueError:
            print ('Node indexes should be INT', fro_idx, to_idx)
            return False
        if fro_idx in self.G and to_idx in self.G:
            # Remove edge from previous node
            self.G[fro_idx].to_edges.pop(to_idx, None)
            # Remove edge to next node
            self.G[to_idx].to_edges.pop(fro_idx, None)
        else:
            print ('removeEdge(): Nodes not exist', fro_idx, to_idx)
            return False
        return True

    def lookupMutationGraphNode(self, original_seq_len, node_from, node_to, snp_char, cov, node_type):
        if self.graph_node_len == original_seq_len:
            # Graph does not have any mutation node so we add
            new_ind = self.graph_node_len
            self.addNode(new_ind, snp_char, cov, node_type)
            self.addEdge(node_from, new_ind, cov)
            self.addEdge(new_ind, node_to, cov)
            return True
        elif self.graph_node_len < original_seq_len:
            print ('graph structure was not valid upon initiation')
            return False
        for idx, node in self.G.iteritems():
            if idx < original_seq_len:
                continue
            if node_from in node.from_edges and node_to in node.to_edges and node.base == snp_char:
                node.cov += cov
                self.addEdge(node_from, idx, cov)
                self.addEdge(idx, node_to, cov)
                return True
        # We don't have identical existing mutation node so we add
        new_ind = self.graph_node_len
        self.addNode(new_ind, snp_char, cov, node_type)
        self.addEdge(node_from, new_ind, cov)
        self.addEdge(new_ind, node_to, cov)
        return True

    def refNodeExists(self, ref_coor, ref_base, prev_ind, next_ind):
        if prev_ind not in self.G or next_ind not in self.G:
            return False
        for i in range(0, len(ref_base)):
            if ref_coor+i not in self.G:
                return False
            if self.G[ref_coor+i].base != ref_base[i]:
                #                print ('VCF coordinate does not match with reference',\
                #                ref_coor, self.G[ref_coor+i].base, ref_base)
                return False
        return True

    def addSNP(self, ref_coor, ref_base, mut_base, cov):
        prev_ind = ref_coor - 1
        next_ind = ref_coor + 1
        if not self.refNodeExists(ref_coor, ref_base, prev_ind, next_ind):
            return False
        self.addMutation(mut_base, cov, 'SNP', prev_ind,
                         next_ind, ref_coor, cov)
        return True

    def addINS(self, ref_coor, ref_base, mut_base, cov):
        prev_ind = ref_coor
        next_ind = ref_coor + 1
        if not self.refNodeExists(ref_coor, ref_base, prev_ind, next_ind):
            return False
        self.addMutation(mut_base[1:], cov, 'INS',
                         prev_ind, next_ind, ref_coor, cov)
        return True

    def addSUB(self, ref_coor, ref_base, mut_base, cov):
        prev_ind = ref_coor - 1
        next_ind = ref_coor + len(mut_base)
        if not self.refNodeExists(ref_coor, ref_base, prev_ind, next_ind):
            return False
        self.addMutation(mut_base, cov, 'SUB', prev_ind,
                         next_ind, ref_coor, cov)
        return True

    def addDEL(self, ref_coor, ref_base, mut_base, cov):
        prev_ind = ref_coor
        next_ind = ref_coor + len(ref_base)
        if not self.refNodeExists(ref_coor, ref_base, prev_ind, next_ind):
            return False
        self.addMutation('', cov, 'DEL', prev_ind, next_ind, ref_coor, cov)
        return True
