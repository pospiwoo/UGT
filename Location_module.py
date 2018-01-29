# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 21:54:28 2017

@author: Sunghee Woo
"""
import os
import sys
import time
import re
from datetime import timedelta
import Parameters as param
import CustomFormats as formats


class GenomicLocation:

    def __init__(self, master_seq, ref_fa_file, fasta_database_name, location_file_name, ms_result_file_name):
        self.pep_seq_col = 0
        self.ref_fa_file = ref_fa_file
        self.fasta_database_name = fasta_database_name
        self.location_file_name = location_file_name
        self.ms_result_file_name = ms_result_file_name
        self.pepdic = {}
        self.cnt_peptides = 0
        self.master_seq = master_seq
        self.initDBFileDIR()
        print 'Processing determinNovelty....'
        self.determinNovelty()

    def initDBFileDIR(self):
        if os.path.isdir(self.fasta_database_name):
            self.db_list = fasta_database_name.split(',')
        else:
            self.db_list = [self.fasta_database_name]

    def determinNovelty(self):
        cnt_all = 0
        cnt_novel = 0
        with open(self.ms_result_file_name, 'r') as inResultFile:
            for line in inResultFile:
                if line == "":
                    continue
                if line.startswith("#") or line.startswith('PEP'):
                    # caption = line
                    continue
                data = line.split('\t')
                if len(data) < 3:
                    continue
                # TLHNLVIQYASQGR', 'NaN', 'sequence_992465', '157936.015625
                cleaned_seq = formats.CleanPeptideString(
                    data[self.pep_seq_col])
                # gene_id = data[1]
                # index_id = data[2]
                # quant_val = float(data[3])
                cnt_all += 1
                if formats.replaceIdenMassAA(cleaned_seq[1:-2]) in self.master_seq:
                    continue
                # print cleaned_seq
                cnt_novel += 1
                if cleaned_seq in self.pepdic:
                    self.pepdic[cleaned_seq][0].append(line)
                    self.pepdic[cleaned_seq][
                        3] = self.pepdic[cleaned_seq][3] + 1
                    if float(data[15]) < self.pepdic[cleaned_seq][4]:
                        self.pepdic[cleaned_seq][4] = float(data[15])
                else:
                    # [[spectra file information],
                    # [location information],
                    # [extra needed],
                    # peptides count,
                    # fdr,
                    # location count]
                    self.pepdic[cleaned_seq] = [
                        [], {}, [], 1, float(0), 0]
                    self.pepdic[cleaned_seq][0].append(line)
                    self.cnt_peptides += 1
            print '# of peptides', cnt_all
            print '# novel peptides', cnt_novel
            # for i in self.pepdic:
            #     print i, self.pepdic[i]

    def binarySearch(self, array, key, imin, imax):
        if (imax < imin):
            return imax
        else:
            imid = (imin+imax)/2
            if array[imid] > key:
                return self.binarySearch(array, key, imin, imid-1)
            elif array[imid] < key:
                return self.binarySearch(array, key, imid+1, imax)
            else:
                return imid

    def get_location(self, start_seq, end_seq, start, end, length, strand, chrNum):
        i = -1
        j = -1
        while start_seq >= 0:
            i += 1
            start_seq -= length[i]
        while end_seq > 0:
            j += 1
            end_seq -= length[j]
        if strand == 0:
            tmp_start = start[i] - start_seq
            tmp_end = start[j] - end_seq
        else:
            tmp_start = end[i] + start_seq
            tmp_end = end[j] + end_seq
        # print ' start,end,i,j ',tmp_start,tmp_end,i,j,
        location = [[], [], []]
        if strand == 0:
            location.append(tmp_end)
            location.append(tmp_start)
            if j-i <= 0:
                location[0].append(tmp_end)
                location[1].append(tmp_start)
                location[2].append(tmp_start-tmp_end)
            else:
                location[0].append(start[i])
                location[1].append(tmp_start)
                location[2].append(tmp_start-start[i])
                for k in range(j-i-1):
                    location[0].append(start[k+i+1])
                    location[1].append(end[k+i+1])
                    location[2].append(end[k+i+1]-start[k+i+1])
                location[0].append(tmp_end)
                location[1].append(end[j])
                location[2].append(end[j]-tmp_end)
        else:
            location.append(tmp_start)
            location.append(tmp_end)
            if j-i <= 0:
                location[0].append(tmp_start)
                location[1].append(tmp_end)
                location[2].append(tmp_end-tmp_start)
            else:
                location[0].append(tmp_start)
                location[1].append(end[i])
                location[2].append(end[i]-tmp_start)
                for k in range(j-i-1):
                    location[0].append(start[k+i+1])
                    location[1].append(end[k+i+1])
                    location[2].append(end[k+i+1]-start[k+i+1])
                location[0].append(start[j])
                location[1].append(tmp_end)
                location[2].append(tmp_end-start[j])
        location.append(strand)
        location.append(chrNum)
        return location

    def SpliceDBRead(self):
        for db in self.db_list:
            print 'Reading database: ', db
            fasta_database = open(db, 'r')
            dummy, fileExtension = os.path.splitext(fasta_database_name)

            if fileExtension.strip() == '.txt':
                file_case = 0
            elif fileExtension.strip() == '.fa':
                file_case = 1
            elif fileExtension.strip() == '.fasta':
                file_case = 1

            if file_case == 1:
                fasta = fasta_database.readlines()
                fasta_info = []
                index_info = [0]
                fasta_seq = ''
                #sequence_info = []

                count = 0
                for i in range(len(fasta)):
                    if count % 2 == 0:
                        fasta_info.append(fasta[i])
                    else:
                        # sequence_info.append(fasta[i].strip())
                        fasta_seq += fasta[i].strip() + 'X'
                        index_info.append(len(fasta_seq))
                    count += 1

                del fasta
            else:
                fasta_info = []
                index_info = [0]
                fasta_seq = ''
                for file in fasta_database:
                    input = open(file.strip(), 'r')
                    fasta = input.readlines()
                    count = 0
                    for i in range(len(fasta)):
                        if count % 2 == 0:
                            fasta_info.append(fasta[i])
                        else:
                            # sequence_info.append(fasta[i].strip())
                            fasta_seq += fasta[i].strip() + 'X'
                            index_info.append(len(fasta_seq))
                        count += 1

                    del fasta
            #>Splice@chr1@176@15444@0;20082126-20082213;20073655-20073753;20072948-20073092;20072024-20072144;20070131-20070219;  example of splice_info

            pep_list = self.pepdic.keys()
            for pep in pep_list:
                pep = formats.CleanPeptideString(pep)
                location_index = [m.start()
                                  for m in re.finditer(pep, fasta_seq)]
                pep_location = []
                for location in location_index:
                    index = self.binarySearch(
                        index_info, location, 0, len(index_info)-1)
                    splice_info = fasta_info[index].split('@')
                    #full_seq = sequence_info[index]
                    #start_seq = full_seq.find(pep)
                    start_seq = location - index_info[index]
                    end_seq = start_seq + len(pep)
                    start_seq = start_seq * 3
                    end_seq = end_seq * 3
                    if splice_info[0].find('Splice') > -1 or splice_info[0].find('Varient') > -1:
                        chrNum = splice_info[1]
                        splice_in = splice_info[4].split(';')
                        # print start_seq,end_seq, sequence, full_seq,
                        start = []
                        end = []
                        length = []
                        for i, splice in enumerate(splice_in):
                            if i == 0:
                                strand = int(splice)
                            elif i == len(splice_in)-1:
                                continue
                            else:
                                if splice.find('/') > -1:
                                    splice = splice.split('/')
                                else:
                                    splice = splice.split('-')
                                start.append(int(splice[0]))
                                end.append(int(splice[1]))
                                length.append(int(splice[1])-int(splice[0]))
                    else:
                        strand = int(splice_info[3])
                        start = [int(splice_info[1])]
                        end = [int(splice_info[2])]
                        length = [end[0]-start[0]]
                        chrNum = splice_info[0][1:]

                    pep_location.append(self.get_location(
                        start_seq, end_seq, start, end, length, strand, chrNum))

                prev = 0
                i = 0
                while i < len(pep_location):  # delete same location seq
                    if prev == pep_location[i][3]:
                        del pep_location[i]
                    else:
                        prev = pep_location[i][3]
                        i += 1

                # save pep location in pepdic with checking the occurence in
                # previous pepdic
                for loc in pep_location:
                    chrNum = loc[6]  # [-1]
                    if not chrNum in self.pepdic[pep][1]:
                        self.pepdic[pep][1][chrNum] = [loc]
                    else:
                        list_in_pepdic = self.pepdic[pep][1].get(chrNum)
                        is_in_list = 0
                        for i in list_in_pepdic:
                            if i[3] == loc[3]:
                                is_in_list = 1
                        if is_in_list == 0:
                            self.pepdic[pep][1][chrNum].append(loc)
        del fasta_info
        del index_info
        del fasta_seq

    def LocationCounts(self):
        # location count fill in
        for pep in self.pepdic:
            count = 0
            for chr in self.pepdic[pep][1]:
                count += len(self.pepdic[pep][1][chr])
            if len(self.pepdic[pep]) == 5:
                self.pepdic[pep].append(count)
            else:
                self.pepdic[pep][5] = count

        # copy list
        coordi_list = {}
        for pep in self.pepdic:
            for chrNum in self.pepdic[pep][1]:
                if not chrNum in coordi_list.keys():
                    coordi_list[chrNum] = {}
                for list in self.pepdic[pep][1][chrNum]:
                    coordi_list[chrNum][list[3]] = pep
        #coordi_list = {'chr1':{100000:'PEPTIDE',1000100:'PEPDIDE'}}
        # list_in_chr = [1000000,10000100,...] coordinates of peptide in each chromosome
        # grouping_inf = [0,1,3,6,10,...] index of beginning groups

        with open(self.location_file_name, 'w') as outFile:
            outFile.write(
                '#chr\tstart-end\tPEP\tspec_count\tlocation_count\tFDR\tSprob\n')

            for chrNum in coordi_list:
                list_in_chr = coordi_list[chrNum].keys()
                list_in_chr.sort()
                grouping_info = []
                gstart = -3000
                for index in range(len(list_in_chr)):
                    if gstart + 5000 < list_in_chr[index]:
                        grouping_info.append(index)
                    gstart = list_in_chr[index]

                Sprob = []
                for group in range(len(grouping_info)-1):
                    prob = 1
                    for index in range(grouping_info[group], grouping_info[group+1]):
                        prob *= (1 - (1-self.pepdic[coordi_list[chrNum][list_in_chr[index]]][
                                 4]) / self.pepdic[coordi_list[chrNum][list_in_chr[index]]][5])
                    prob = 1 - prob
                    Sprob.append(prob)
                if len(grouping_info) == 0:
                    Sprob.append(0)
                else:
                    prob = 1
                    for index in range(grouping_info[len(grouping_info)-1], len(coordi_list[chrNum])):
                        prob *= (1 - (1-self.pepdic[coordi_list[chrNum][list_in_chr[index]]][
                                 4]) / self.pepdic[coordi_list[chrNum][list_in_chr[index]]][5])
                    prob = 1 - prob
                    Sprob.append(prob)

                for coord in list_in_chr:
                    Sp_index = self.binarySearch(grouping_info, list_in_chr.index(
                        coord), 0, len(grouping_info)-1)
                    pep = coordi_list[chrNum][coord]
                    outFile.write(chrNum+'\t')
                    for i in self.pepdic[pep][1].get(chrNum):
                        if i[3] != coord:
                            continue
                        for j in range(len(i[0])):
                            outFile.write(str(i[0][j])+'/'+str(i[1][j])+';')
                        outFile.write(
                            '\t'+pep+'\t'+str(self.pepdic[pep][3])+'\t'+str(self.pepdic[pep][5])+'\t')
                        outFile.write(str(self.pepdic[pep][4])+'\t')
                        # outFile.write(str(Sprob[Sp_index])+'\n')
                        outFile.write(
                            str(Sprob[Sp_index])+'\t'+str(i[-2])+'\n')

    def executionTime(func):  # timer decorator
        def wrapper(*args, **kwargs):
            begin = time.time()
            func(*args, **kwargs)
            finish = time.time()
            time_in_sec = finish - begin
            print 'Excecution time:', str(timedelta(seconds=time_in_sec))
        return wrapper

    @executionTime
    def process(self):
        print 'Processing Locations....'
        self.SpliceDBRead()
        self.LocationCounts()


if __name__ == "__main__":
    if len(sys.argv) == 4:
        ref_fa_file_name = sys.argv[1]
        fasta_database_name = sys.argv[2]
        location_file_name = sys.argv[3]
        ms_result_file_name = sys.argv[4]
    elif len(sys.argv) == 0:
        ref_fa_file_name = ''
        fasta_database_name = '/home/s3cha/data/SpliceDB/IG_VU_DB/fasta_file/UNCID_1582766.974b1fa9-8ee7-4f02-b0ff-221fc98abe5f.sorted_genome_alignments20150521-950-hevqjo_IGH_AND_UNMAPPED_splicegraph_1.fa'
        location_file_name = '/home/s3cha/data/SpliceDB/IG_VU_DB/temp_location.txt'
        ms_result_file_name = '/home/s3cha/data/SpliceDB/IG_VU_DB/VU_IG_DB_result_SPEC01.p'
    else:
        ref_fa_file_name = param.ref_fa_file_name()
        fasta_database_name = param.fasta_out_file()
        location_file_name = param.location_file()
        ms_result_file_name = param.ms_result_file()

    master_seq = formats.parseRefProtFASTA(ref_fa_file_name)
    loc_obj = GenomicLocation(master_seq,
                              ref_fa_file_name,
                              fasta_database_name,
                              location_file_name,
                              ms_result_file_name)
    loc_obj.process()
    del loc_obj
