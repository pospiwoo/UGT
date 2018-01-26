# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 21:54:28 2017

@author: Sunghee Woo
"""
import os
import sys
import time
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
        # self.parseRefProtFASTA()
        self.determinNovelty()
        # db_list = fasta_database_name.split(',')

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
                print cleaned_seq
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
            print cnt_all, cnt_novel
            for i in self.pepdic:
                print i, self.pepdic[i]

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
        for db in db_list:
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

            pep_list = pepdic.keys()
            for pep in pep_list:
                pep = CleanPeptideString(pep)
                location_index = [m.start()
                                  for m in re.finditer(pep, fasta_seq)]
                pep_location = []
                for location in location_index:
                    index = binary_search(
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

                    pep_location.append(get_location(
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
                    if not pepdic[pep][1].has_key(chrNum):
                        pepdic[pep][1][chrNum] = [loc]
                    else:
                        list_in_pepdic = pepdic[pep][1].get(chrNum)
                        is_in_list = 0
                        for i in list_in_pepdic:
                            if i[3] == loc[3]:
                                is_in_list = 1
                        if is_in_list == 0:
                            pepdic[pep][1][chrNum].append(loc)
        del fasta_info
        del index_info
        del fasta_seq

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
        with open(ms2db_file_name, 'r') as inFile:
            for line in inFile:
                if line.find('<Database') > -1:
                    continue
                if line.find('<Gene Name=') > -1:
                    tmp = line.split('"')
                    self.geneName = tmp[1]
                    self.ExonCount = tmp[3]
                    self.ForwardFlag = int(tmp[7])
                    self.exon_length = []
                    self.exon_array = []
                    self.next_exon = []
                    self.previous_exon = []
                    self.coordi_start = []
                    self.coordi_end = []
                    self.source_exon = []
                    self.incoming_indicator = []
                    for i in range(int(self.ExonCount)):
                        self.exon_length.append(0)
                        self.exon_array.append('')
                        self.next_exon.append({})
                        self.previous_exon.append({})
                        self.incoming_indicator.append(1)
                        self.coordi_start.append(0)
                        self.coordi_end.append(0)
                    continue
                if line.find('<Exon Index') > -1:
                    tmp = line.split('"')
                    self.current_exon = int(tmp[1])
                    self.coordi_start[self.current_exon] = int(tmp[3])
                    self.coordi_end[self.current_exon] = int(tmp[5])
                    continue
                if line.find('<ExonSequence') > -1:
                    tmp = line.split('"')
                    self.exon_length[self.current_exon] = int(tmp[1])
                    tmp1 = tmp[2].split('>')
                    tmp2 = tmp1[1].split('<')
                    self.exon_array[self.current_exon] = tmp2[0]
                    if self.coordi_end[self.current_exon] < 0:
                        self.coordi_start[self.current_exon] = \
                            self.coordi_end[self.current_exon] - \
                            self.exon_length[self.current_exon]
                    continue
                if line.find('<ExtendsExon Index=') > -1:
                    tmp = line.split('"')
                    self.previous_exon[self.current_exon][int(tmp[1])] = 0
                    continue
                if line.find('<LinkFrom Index=') > -1:
                    tmp = line.split('"')
                    self.previous_exon[self.current_exon][int(tmp[1])] = 1
                    continue
                if line.find('<LinkFromInsertion Index=') > -1:
                    tmp = line.split('"')
                    self.previous_exon[self.current_exon][int(tmp[1])] = 2
                    continue
                if line.find('<LinkFromDeletion Index=') > -1:
                    tmp = line.split('"')
                    self.previous_exon[self.current_exon][int(tmp[1])] = 3
                    continue
                if line.find('<LinkFromMutation Index=') > -1:
                    tmp = line.split('"')
                    self.previous_exon[self.current_exon][int(tmp[1])] = 4
                    continue
                if line.find('</Gene>') == -1:
                    continue
                if line.find('</Gene>') > -1:
                    error = 0
                    for i in range(len(self.previous_exon)):
                        tmp = self.previous_exon[i].keys()
                        for k in tmp:
                            if k >= len(self.previous_exon):
                                error = 1
                                print 'Unfeasible link in: ' + str(self.current_exon)
                                break
                        if error == 1:
                            break
                        if self.previous_exon[i] == {}:
                            self.source_exon.append(i)
                        else:
                            for j in self.previous_exon[i].keys():
                                self.next_exon[j][i] = self.previous_exon[i][j]
                    if error == 1:
                        continue
                    for depth in xrange(int(self.ExonCount)):
                        if self.incoming_indicator[depth] != 1:
                            continue
                        if self.ForwardFlag == 0:
                            coordi_x = self.coordi_end[depth]
                        else:
                            coordi_x = self.coordi_start[depth]
                        self.check_walk = [0]
                        self.MPS(depth, [], self.exon_array, coordi_x)
                    continue

    def writeAA(self, output_array, path, coordi_x, coordi_y):
        if 'X' in output_array or 'N' in output_array:
            return
        AA_array = self.getAA(output_array)
        tmp_s = []
        tmp_e = []
        check_deletion = []
        for i in range(len(path)):
            if path[i] != path[-1]:
                if self.next_exon[path[i]][path[i+1]] == 3:
                    check_deletion.append(path[i])

        for i in path:
            tmp_s.append(self.coordi_start[i])
            tmp_e.append(self.coordi_end[i])
        if self.ForwardFlag == 0:
            tmp_e[0] = coordi_x
            tmp_s[len(tmp_s)-1] = coordi_y
            if len(tmp_s) == 1:
                tmp_s[0] = coordi_y
                tmp_e[0] = coordi_x
        else:
            tmp_s[0] = coordi_x
            tmp_e[len(tmp_e)-1] = coordi_y
        if self.ForwardFlag == 0:
            i = 0
            while i < len(tmp_s)-1:
                while tmp_s[i] == tmp_e[i+1]:
                    tmp_s[i] = tmp_s[i+1]
                    del tmp_s[i+1]
                    del tmp_e[i+1]
                    if i == len(tmp_s)-1:
                        break
                i += 1
        else:
            i = 0
            while i < len(tmp_s)-1:
                while tmp_e[i] == tmp_s[i+1]:
                    tmp_e[i] = tmp_e[i+1]
                    del tmp_s[i+1]
                    del tmp_e[i+1]
                    if i == len(tmp_s)-1:
                        break
                i += 1
        code_len = 0
        s_len = tmp_e[0] - tmp_s[0]
        e_len = tmp_e[len(tmp_s)-1] - tmp_s[len(tmp_s)-1]
        for j in range(len(tmp_s)):
            code_len += (tmp_e[j] - tmp_s[j])
        for j in range(3):
            self.outFile.write('>Splice@'+self.geneName+'@' +
                               str(self.counter_idx)+'@'+str(self.ForwardFlag)+';')
            self.counter_idx += 1
            for i in range(len(tmp_s)):
                if self.ForwardFlag == 0:
                    if i == 0:
                        if s_len == 1:
                            if j == 0:
                                self.outFile.write(
                                    str(tmp_s[i])+'/'+str(tmp_e[i]-j)+';')
                            else:
                                continue
                        elif s_len == 2:
                            if j == 0 or 1:
                                self.outFile.write(
                                    str(tmp_s[i])+'/'+str(tmp_e[i]-j)+';')
                            else:
                                continue
                        else:
                            self.outFile.write(
                                str(tmp_s[i])+'/'+str(tmp_e[i]-j)+';')
                    elif i == len(tmp_s)-2 and e_len == 1 and (code_len-j) % 3 == 2:
                        self.outFile.write(
                            str(tmp_s[i]+1)+'/'+str(tmp_e[i])+';')
                        break
                    elif i == len(tmp_s)-1:
                        self.outFile.write(
                            str(tmp_s[i]+(code_len-j) % 3)+'/'+str(tmp_e[i])+';')
                    else:
                        if s_len == 1 and j == 2 and i == 1:
                            self.outFile.write(
                                str(tmp_s[i])+'/'+str(tmp_e[i]-1)+';')
                        else:
                            self.outFile.write(
                                str(tmp_s[i])+'/'+str(tmp_e[i])+';')
                else:
                    if i == 0:
                        if s_len == 1:
                            if j == 0:
                                self.outFile.write(
                                    str(tmp_s[i])+'/'+str(tmp_e[i]-j)+';')
                            else:
                                continue
                        elif s_len == 2:
                            if j == 0 or 1:
                                self.outFile.write(
                                    str(tmp_s[i]+j)+'/'+str(tmp_e[i])+';')
                            else:
                                continue
                        else:
                            self.outFile.write(
                                str(tmp_s[i]+j)+'/'+str(tmp_e[i])+';')
                    elif i == len(tmp_s)-2 and e_len == 1 and (code_len-j) % 3 == 2:
                        self.outFile.write(
                            str(tmp_s[i])+'/'+str(tmp_e[i]-1)+';')
                        break
                    elif i == len(tmp_s)-1:
                        self.outFile.write(
                            str(tmp_s[i])+'/'+str(tmp_e[i]-(code_len-j) % 3)+';')
                    else:
                        if s_len == 1 and j == 2 and i == 1:
                            self.outFile.write(
                                str(tmp_s[i]+1)+'/'+str(tmp_e[i])+';')
                        else:
                            self.outFile.write(
                                str(tmp_s[i])+'/'+str(tmp_e[i])+';')

            for i in check_deletion:
                if self.ForwardFlag == 0:
                    self.outFile.write(
                        '<deletion>' + str(self.coordi_start[i])+' ')
                else:
                    self.outFile.write(
                        '<deletion>' + str(self.coordi_end[i])+' ')
            self.outFile.write('\n')
            self.outFile.write(AA_array[j] + '\n')

            if (os.path.getsize(self.out_file_name) > 100*1024*1024):
                self.file_counter_idx += 1
                self.out_file_name = self.file_n+"_" + \
                    str(self.file_counter_idx)+self.ext_n
                self.outFile.close()
                self.outFile = open(self.out_file_name, 'w')

    def checkLoop(self, node_index, path):
        if node_index in path:
            print 'loop exist in: ' + self.geneName + ' at path ',
            for i in path:
                print str(i)+' ',
            print node_index
            return True
        return False

    def MPS2(self, node_index, L_tmp, path, tmp_array, coordi_x):
        self.check_walk[-1] += self.exon_length[node_index]
        if self.checkLoop(node_index, path):
            return path
        output_array = ''
        tmp_path = []
        if L_tmp - self.exon_length[node_index] < 0 or node_index > self.depth_cutoff:
            path.append(node_index)
            for i in range(len(path)):
                if i != len(path)-1:
                    output_array += tmp_array[path[i]]
                    if output_array != '':
                        tmp_path.append(path[i])
                else:
                    output_array = output_array + tmp_array[path[i]][:L_tmp]
                    tmp_path.append(path[i])
            if self.ForwardFlag == 0:
                coordi_y = self.coordi_end[node_index] - L_tmp
            else:
                coordi_y = self.coordi_start[node_index] + L_tmp
            self.writeAA(output_array, tmp_path, coordi_x, coordi_y)
            path.pop()
            return
        elif len(self.next_exon[node_index]) == 0 or node_index > self.depth_cutoff:
            path.append(node_index)
            for i in range(len(path)):
                output_array += tmp_array[path[i]]
                if output_array != '':
                    tmp_path.append(path[i])
            if self.ForwardFlag == 0:
                coordi_y = self.coordi_start[node_index]
            else:
                coordi_y = self.coordi_end[node_index]
            self.writeAA(output_array, tmp_path, coordi_x, coordi_y)
            path.pop()
            return
        else:
            L_tmp = L_tmp - self.exon_length[node_index]
            tmp = self.next_exon[node_index].keys()
            tmp.sort()
            for k, j in enumerate(tmp):
                if k == 0:
                    if self.check_walk[-1] < self.min_exon_len and self.coordi_end[j] < 0:
                        continue
                    self.check_walk.append(self.check_walk[-1])
                    if self.next_exon[node_index][j] != 0:
                        self.check_walk[-1] = 0
                    path.append(node_index)
                    self.MPS2(j, L_tmp, path, tmp_array, coordi_x)
                    self.check_walk.pop()
                    path.pop()
                else:
                    self.check_walk.append(self.check_walk[-1])
                    if self.next_exon[node_index][j] != 0 and self.check_walk[-1] < self.min_exon_len:
                        self.check_walk.pop()
                        continue
                    elif self.next_exon[node_index][j] != 0:
                        self.check_walk[-1] = 0
                    path.append(node_index)
                    tmp_sum = self.nu_length
                    tmp_array = self.exon_array[:]
                    for i in range(len(path)):
                        tmp_sum = tmp_sum - \
                            self.exon_length[path[len(path)-i-1]]
                        if tmp_sum <= 0 and tmp_sum + self.exon_length[path[len(path)-i-1]] > 0:
                            tmp_array[
                                path[len(path)-i-1]] = tmp_array[path[len(path)-i-1]][-tmp_sum:]
                            if self.ForwardFlag == 0:
                                coordi_x = self.coordi_start[
                                    path[len(path)-i-1]] + (tmp_sum + self.exon_length[path[len(path)-i-1]])
                            else:
                                coordi_x = self.coordi_end[
                                    path[len(path)-i-1]] - (tmp_sum + self.exon_length[path[len(path)-i-1]])
                        elif tmp_sum <= 0:
                            tmp_array[path[len(path)-i-1]] = ''
                    self.MPS2(j, L_tmp, path, tmp_array, coordi_x)
                    self.check_walk.pop()
                    path.pop()

    def MPS(self, node_index, path, temp_array, coordi_x):
        self.check_walk[-1] += self.exon_length[node_index]
        output_array = ''
        temp_path = []
        if self.incoming_indicator[node_index] == 2:
            self.MPS2(node_index, self.nu_length, path, temp_array, coordi_x)
        elif len(self.next_exon[node_index]) == 0 or node_index > self.depth_cutoff:
            self.incoming_indicator[node_index] = 2
            path.append(node_index)
            for i in range(len(path)):
                output_array += temp_array[path[i]]
                if output_array != '':
                    temp_path.append(path[i])
            if self.ForwardFlag == 0:
                coordi_y = self.coordi_start[node_index]
            else:
                coordi_y = self.coordi_end[node_index]
            self.writeAA(output_array, temp_path, coordi_x, coordi_y)
            path.pop()
            return
        else:
            self.incoming_indicator[node_index] = 2
            temp = self.next_exon[node_index].keys()
            temp.sort()
            for k, j in enumerate(temp):
                if k == 0:
                    if self.check_walk[-1] < self.min_exon_len and self.coordi_end[j] < 0:
                        continue
                    self.check_walk.append(self.check_walk[-1])
                    if self.next_exon[node_index][j] != 0:
                        self.check_walk[-1] = 0
                    path.append(node_index)
                    self.MPS(j, path, temp_array, coordi_x)
                    self.check_walk.pop()
                    path.pop()
                else:
                    self.check_walk.append(self.check_walk[-1])
                    if self.next_exon[node_index][j] != 0 and self.check_walk[-1] < self.min_exon_len:
                        # print geneName,path,node_index,self.check_walk
                        self.check_walk.pop()
                        continue
                    elif self.next_exon[node_index][j] != 0:
                        self.check_walk[-1] = 0
                    path.append(node_index)
                    tmp_sum = self.nu_length
                    temp_array = self.exon_array[:]
                    for i in range(len(path)):
                        tmp_sum = tmp_sum - \
                            self.exon_length[path[len(path)-i-1]]
                        if tmp_sum <= 0 and tmp_sum + self.exon_length[path[len(path)-i-1]] > 0:
                            temp_array[
                                path[len(path)-i-1]] = temp_array[path[len(path)-i-1]][-tmp_sum:]
                            if self.ForwardFlag == 0:
                                coordi_x = self.coordi_start[path[len(path)-i-1]] \
                                    + tmp_sum \
                                    + self.exon_length[path[len(path)-i-1]]
                            else:
                                coordi_x = self.coordi_end[path[len(path)-i-1]] \
                                    - tmp_sum \
                                    - self.exon_length[path[len(path)-i-1]]
                        elif tmp_sum <= 0:
                            temp_array[path[len(path)-i-1]] = ''
                    self.MPS(j, path, temp_array, coordi_x)
                    self.check_walk.pop()
                    path.pop()


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
                              ref_fa_file_name, fasta_database_name, location_file_name, ms_result_file_name)
    # loc_obj.process()
    del loc_obj
