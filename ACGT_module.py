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
sys.setrecursionlimit(3000)


class FASTAConvert:

    def __init__(self, ms2db_file_name, out_file_prefix, input_length):
        self.aa_length = 29
        self.nu_length = self.aa_length * 3
        self.min_exon_len = int(input_length)
        self.depth_cutoff = 2000
        self.file_count = 1
        self.counter_idx = 0
        self.ms2db_file_name = ms2db_file_name
        (self.file_n, self.ext_n) = os.path.splitext(out_file_prefix)
        self.out_file_name = self.file_n + "_" + \
            str(self.file_count) + self.ext_n
        self.outFile = open(self.out_file_name, 'w')
        self.ForwardCode = formats.initForwardCodon()
        self.ReverseCode = formats.initReverseCodon()

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

    def getAA(self, output_array):
        output_seq = []
        for j in range(3):
            AA_seq = ''
            if self.ForwardFlag == 0:
                for i in range(len(output_array[j:])/3):
                    AA_seq += self.ForwardCode[output_array[3*i+j:3*(i+1)+j]]
            else:
                for i in range(len(output_array[j:])/3):
                    AA_seq += self.ForwardCode[output_array[3*i+j:3*(i+1)+j]]
            output_seq.append(AA_seq)
        return output_seq

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
    if len(sys.argv) > 2:
        ms2db_file_name = sys.argv[1]
        out_file = sys.argv[2]
        input_length = int(sys.argv[3])
    else:
        ms2db_file_name = param.ms2db_file_name()
        out_file = param.fasta_out_file()
        input_length = param.max_part_AA_len()

    fa_convert_obj = FASTAConvert(ms2db_file_name, out_file, input_length)
    fa_convert_obj.process()
    del fa_convert_obj
