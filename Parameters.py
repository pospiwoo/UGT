# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 20:11:43 2017

@author: Sunghee Woo
"""
import os
from os.path import dirname  # , abspath
script_dir = dirname(__file__)

# Parameters
# SPL_generator_module.py
empty_spl_ = os.path.join(script_dir, 'spl', 'empty.spl')
spl_out_file_name_ = os.path.join(script_dir, 'spl', 'vcf_somatic_test_2.spl')
vcf_spl_bam_dir_ = os.path.join(script_dir, 'vcf')
# ACGT_module.py
ms2db_file_name_ = os.path.join(
    script_dir, 'ms2db', 'vcf_somatic_test_2.ms2db')
fasta_out_file_ = os.path.join(script_dir, 'fa', 'vcf_somatic_test_gae_2.fa')
fasta_out_file_ = os.path.join(script_dir, 'fa', 'dbSNP_clinvar.fa')
max_part_AA_len_ = 30
# gene_module.py
vcf_file_name_ = 'common_all_20170403_0.vcf'
# Location_module.py
ref_fa_file_name_ = os.path.join(
    script_dir, 'ref_fa', 'uniprot_canonical.fasta')
ms_result_file_ = os.path.join(
    script_dir, 'ms_result', '20180126_ALDER_G2957_470_560_directDIA_20180125F1_Peptide.xls')
location_file_ = os.path.join(
    script_dir, 'loc', '20180126_ALDER_G2957_470_560_directDIA_20180125F1_Peptide.txt')
# Event_module.py
gff_file_name_ = os.path.join(script_dir, 'gff', 'Homo_sapiens.GRCh38.88.gff3')
event_filename_ = os.path.join(
    script_dir, 'event', location_file_.replace('.txt', '_event.txt'))


def empty_spl():
    return empty_spl_


def ref_fa_file_name():
    return ref_fa_file_name_


def spl_out_file_name():
    return spl_out_file_name_


def vcf_spl_bam_dir():
    return vcf_spl_bam_dir_

# ACGT_module.py


def ms2db_file_name():
    return ms2db_file_name_


def fasta_out_file():
    return fasta_out_file_


def max_part_AA_len():
    return max_part_AA_len_

# Gene_module.py


def vcf_file_name():
    return vcf_file_name_

# Location_module.py


def location_file():
    return location_file_


def ms_result_file():
    return ms_result_file_

# Location_module.py


def gff_file_name():
    return gff_file_name_


def event_filename():
    return event_filename_


# debug_output_on_ = False # False # True
#debug_output_file_name_ = os.path.join(path.out_dir, 'debug_output.txt')
 # normal
# self.parent_dir = dirname(dirname(dirname(abspath(__file__)))) # py2exe

#        self.out_dir = os.path.join(self.parent_dir, 'Output')
#        if not os.path.exists(self.out_dir):
#            os.makedirs(self.out_dir)
#
#        self.input_dir = os.path.join(self.parent_dir, 'Input')
#        self.encoding_dir = os.path.join(self.parent_dir, 'encoding')

#    def makeNewOutputDir(self, new_output_dir):
#        self.out_dir = os.path.join(self.out_dir, new_output_dir)
#        if not os.path.exists(self.out_dir):
#            os.makedirs(self.out_dir)
