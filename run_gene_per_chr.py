#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 16:21:14 2018

@author: swoo
"""
import os
import sys
import pickle
import gc
import Parameters as param
import CustomFormats as formats
import Gene_module as gene


if sys.argv[1] == '0':
    parse_obj = gene.ConstructGraphPerChr('chr1')
    print (' -Reading reference DNA')
    parse_obj.parseDNA()
    print (' -Creating transcript graph from GFF')
    parse_obj.parseGFF()
elif sys.argv[1] == '1':
    parse_obj = gene.unpickleGraph('chr1')
    for idx in parse_obj.GGraph.G:
        node = parse_obj.GGraph.G[idx]
        print idx, node.base, \
            node.start, node.end, \
            node.type
    # print (' -Creating variant graph fvrom VCF')
    # parse_obj.parseVCF()
