# -*- coding: utf-8 -*-
"""
Created on Fri Jun 23 21:13:85 2017

@author: Sunghee Woo
"""
import os
import sys
import time
import re
from datetime import timedelta
import Parameters as param
import CustomFormats as formats


class SPLConvert:

    def __init__(self, GFFinput_file_name, input_file_folder, output_file_name):
        self.min_DP = 2
        self.Species = 'Human'  # 'CElegans', 'Human'
        self.check_sample = [0, 1000]  # 0: check off, 1: check on
        self.GFFinput_file_name = GFFinput_file_name
        self.input_file_folder = input_file_folder
        self.output_file_name = output_file_name
        self.file_mapping_info = {}
        self.chr_ref = formats.initChr()
        self.splice_data = formats.initMutationData()
        self.in_data = formats.initMutationData()
        self.del_data = formats.initMutationData()
        self.mu_data = formats.initMutationData()
        self.initFileNames()

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
        self.processFiles()
        self.writeSPL()

    def writeSPL(self):
        with open(output_file_name, 'w') as oFile:
            self.writeFileMappingHeader(oFile)
            self.writeSingleVariant(oFile, self.splice_data, 'Splice')
            self.writeSingleVariant(oFile, self.in_data, 'Insertion')
            self.writeSingleVariant(oFile, self.del_data, 'Deletion')
            self.writeSingleVariant(oFile, self.mu_data, 'Mutation')

    def writeFileMappingHeader(self, oFile):
        oFile.write('#')
        for i in self.file_mapping_info:
            oFile.write(i+':'+str(self.file_mapping_info[i]))
            if i != self.file_mapping_info.keys()[-1]:
                oFile.write(',')
        oFile.write('\n')

    def writeSingleVariant(self, oFile, variant_data, variant_type):
        tmp_chr = variant_data.keys()
        # read depth might not always same for each vcf file format. this is for TCGA samples.
#       temp = line[7].split(';')
#       DP = 1
#       for checkdp in temp:
#           if checkdp.startswith('DP'):
#               DP = int(checkdp.split('=')[1])
        tmp_chr.sort()
        oFile.write('#'+variant_type+'\n')
        for chromosome in tmp_chr:
            data = variant_data[chromosome]
            for key, value in data.items():
                for list in value:
                    #temp_sum = sum(list[1].values())
                    # if temp_sum < min_DP:
                    #   continue
                    oFile.write(chromosome+'\t'+str(key) +
                                '\t'+str(list[0])+'\t')
                    for i in list[1]:
                        oFile.write(str(i)+':'+str(list[1][i]))
                        if i != list[1].keys()[len(list[1])-1]:
                            oFile.write(',')
                    oFile.write('\t'+list[2]+'\n')

    def initFileNames(self):
        self.in_filelist = []
        if os.path.isdir(self.input_file_folder):
            for path, subdirs, files in os.walk(self.input_file_folder):
                for name in files:
                    self.in_filelist.append(os.path.join(path, name))
        else:
            self.in_filelist.append(self.input_file_folder)

    def processFiles(self):
        for afile in self.in_filelist:
            afile = afile.strip()
            afile_name = os.path.split(afile)[1]
            # print 'Read file: '+file
            if self.file_mapping_info.has_key(afile):
                print afile, 'has been read previously'
                continue
            SAM_file_idx = len(self.file_mapping_info)
            print 'Read file: ' + afile_name, SAM_file_idx
            dummy, fileExtension = os.path.splitext(afile)
            if fileExtension.strip() == '.sam':
                file_case_intxt = 1
            elif fileExtension.strip() == '.gff':
                file_case_intxt = 2
            elif fileExtension.strip() == '.bam':
                file_case_intxt = 3
            elif fileExtension.strip() == '.tmp' or \
                    fileExtension.strip() == '.spl':
                file_case_intxt = 4
            elif fileExtension.strip() == '.vcf':
                file_case_intxt = 5
            elif fileExtension.strip() == '.maf':
                file_case_intxt = 6
            else:
                print 'Input file name is not readable: ' + fileExtension
                continue
            if file_case_intxt != 4:
                self.file_mapping_info[afile_name] = SAM_file_idx
            if file_case_intxt != 3:
                SourceFile = open(afile, 'r')
            if file_case_intxt == 5:
                self.addVCF(SourceFile, SAM_file_idx)
#            if file_case_intxt == 1:
#                        addSAM(SourceFile,ID)
#            elif file_case_intxt == 2:
#                        addGFF(SourceFile,ID)
#            elif file_case_intxt == 3:
#                        addBAM(ID)
#            elif file_case_intxt == 4:
#                        addTMP(SourceFile,ID,file_mapping_info)
#            elif file_case_intxt == 6:
#                        addMAF(SourceFile,ID)

    def baseconvert(self, n, base):
        """convert positive decimal integer n to equivalent in another base (2-36)"""
        digits = "0123456789abcdefghijklmnopqrstuvwxyz"
        try:
            n = int(n)
            base = int(base)
        except:
            return ""
        if n < 0 or base < 2 or base > 36:
            return ""

        s = ""
        while 1:
            r = n % base
            s = digits[r] + s
            n = n / base
            if n == 0:
                break
        return s

    def checkParentID(self, part1, part2):
        if part1[0] == -1 or part2[0] == -1:
            return False
        tmp1 = part1[-1].split(";")[0].split("=")[1]
        tmp2 = part2[-1].split(";")[0].split("=")[1]
        if tmp1 == tmp2:
            return True
        else:
            return False

    def decomp(self, CIGAR, CIGAR_begin, sign, strand):
        del_location = []
        seq_location = []
        pattern = re.compile('S|N|M|D|I|[0-9]+')
        CIGAR_info = [(m.group()) for m in pattern.finditer(CIGAR)]
        CIGAR_interval = [int(n)
                          for i, n in enumerate(CIGAR_info) if i % 2 == 0]
        index_D = [i for i, n in enumerate(CIGAR_info) if n == sign]
        for i in index_D:
            del_interval = int(CIGAR_info[i-1])
            del_start = CIGAR_begin
            for j in range((i-1)/2):
                if CIGAR_info[2*j+1] != 'I' and CIGAR_info[2*j+1] != 'S':
                    del_start += CIGAR_interval[j]
            #del_start = CIGAR_begin + sum(CIGAR_interval[0:(i-1)/2])
            del_end = del_start + del_interval
            del_location.append((del_start, del_end))
            seq_finder = 0
            for j in range((i-1)/2):
                if CIGAR_info[2*j+1] != 'D' and CIGAR_info[2*j+1] != 'N':
                    seq_finder += CIGAR_interval[j]
            seq_location.append(
                (seq_finder, seq_finder+CIGAR_interval[(i-1)/2]))
        return [del_location, seq_location]

    def addVCF(self, SourceFile, SAM_file_idx):
        for line in SourceFile:
            if line.startswith('#'):
                continue
            line = line.split('\t')
            if str(line[0].strip()).startswith('chr'):
                chromosome = line[0].strip()
            else:
                chromosome = 'chr'+str(line[0].strip())
            if chromosome not in self.chr_ref:
                # print 'Chromosome: ', chromosome,' is not readable'
                continue
            # read depth might not always same for each vcf file format. this is for TCGA samples.
#       temp = line[7].split(';')
#       DP = 1
#       for checkdp in temp:
#           if checkdp.startswith('DP'):
#               DP = int(checkdp.split('=')[1])
            if 'DP=' in line[7]:
                tmp = line[7].split(';')
                tmp_index = 0
                for index in xrange(len(tmp)):
                    if tmp[index].startswith('DP='):
                        tmp_index = index
    #            print line[7], tmp_index
                DP = int(line[7].split(';')[tmp_index].split('=')[1])
            else:
                DP = 1
                continue
            ##

            strand = '.'
            coordi = int(line[1])-1
            REF = line[3]
            ALT_g = line[4].split(',')
            for ALT in ALT_g:
                if len(REF) == 1 and len(ALT) == 1:
                    data = self.mu_data.get(chromosome)
                    if coordi in data:
                        check = 0
                        for i, check_end in enumerate(data[coordi]):
                            if check_end[0] == ALT:
                                data[coordi][i][1][SAM_file_idx] = DP
                                check = 1
                        if check == 0:
                            data[coordi].append(
                                [ALT, {SAM_file_idx: DP}, strand])
                    else:
                        data[coordi] = [[ALT, {SAM_file_idx: DP}, strand]]
                elif len(REF) > len(ALT):
                    data = self.del_data.get(chromosome)
                    coordi = coordi + 1
                    coordi_end = coordi+len(REF)-len(ALT)
                    if coordi in data:
                        check = 0
                        for i, check_end in enumerate(data[coordi]):
                            if check_end[0] == coordi_end:
                                data[coordi][i][1][SAM_file_idx] = DP
                                check = 1
                        if check == 0:
                            data[coordi].append(
                                [coordi_end, {SAM_file_idx: DP}, strand])
                    else:
                        data[coordi] = [
                            [coordi_end, {SAM_file_idx: DP}, strand]]
                elif len(REF) < len(ALT):
                    data = self.in_data.get(chromosome)
                    coordi = coordi + len(REF)
                    if coordi in data:
                        for i, check_end in enumerate(data[coordi]):
                            if check_end[0] == ALT[len(REF):]:
                                data[coordi][i][1][SAM_file_idx] = DP
                                check = 1
                            if check == 0:
                                data[coordi].append(
                                    [ALT[len(REF):], {SAM_file_idx: DP}, strand])
                    else:
                        data[coordi] = [
                            [ALT[len(REF):], {SAM_file_idx: DP}, strand]]
        return

    def addBAM(self, SAM_file_idx):
        id_num = 1
        for line in sys.stdin:  # SourceFile:
            if line == "":
                continue
            if line.startswith(">"):
                continue
            if line.startswith("@"):
                continue
            if line.startswith("track"):
                continue
            data = line.strip()
            #data = data.split('\t')
            data = data.split()
            # if not data[2].startswith('chr'):
            #        continue

            QNAME = data[0]
            FLAG = int(data[1])
            RNAME = data[2]
            POS = int(data[3])-1
#            MAPQ = int(data[4])
            CIGAR = data[5]
#            RNEXT = data[6]
#            PNEXT = int(data[7])
            TLEN = int(data[8])
#            SEQ = data[9]
#            QUAL = data[10]
#            p = re.compile('[0-9]+M')
            ###############################
            # if RNAME != chr:
            #        continue
            ###############################
            if self.check_sample[0] == 1:
                if 'N' in CIGAR:
                    id_num += 1
                    if id_num > self.check_sample[1]:
                        break

            tmp_flag = self.baseconvert(FLAG, 2)
            if tmp_flag == "0" or tmp_flag == 0:
                strand = "+"
            elif len(tmp_flag) < 5:
                strand = "+"
            elif int(tmp_flag[-5]) == 0:
                strand = "+"
            elif int(tmp_flag[-5]) == 1:
                strand = "-"
            else:
                print QNAME, CIGAR, POS, TLEN, "can't parse strand from: ", FLAG
                continue
            #chr_name = formatChrNameForCElegans(RNAME)
            if self.Species == 'Human':
                chr_name = formats.formatChrNameForHuman(RNAME)
                #chr_name = RNAME
#            elif self.Species == 'CElegans':
#                    chr_name = self.formatChrNameForCElegans(RNAME)
            if chr_name not in self.chr_ref:
                continue
            #chr_name = formats.formatChrNameForHuman(RNAME)
            if chr_name == -1:
                continue
            #############################################
            CIGAR_begin = POS
            #############################################
            if 'S' in CIGAR:
                continue
            if 'N' in CIGAR:
                sign = 'N'
                [data_location, seq_location] = self.decomp(
                    CIGAR, CIGAR_begin, sign, strand)
                data = self.splice_data.get(chr_name)
                for coord_info in data_location:
                    if coord_info[0] in data:
                        check = 0
                        for i, check_end in enumerate(data[coord_info[0]]):
                            if check_end[0] == coord_info[1] and strand == check_end[2]:
                                if SAM_file_idx in check_end[1]:
                                    data[coord_info[0]][i][
                                        1][SAM_file_idx] += 1
                                else:
                                    data[coord_info[0]][i][1][SAM_file_idx] = 1
                                #data[coord_info[0]][i][1] += 1
                                check = 1
                        if check == 0:
                            data[coord_info[0]].append(
                                [coord_info[1], {SAM_file_idx: 1}, strand])
                    else:
                        data[coord_info[0]] = [
                            [coord_info[1], {SAM_file_idx: 1}, strand]]
                self.splice_data[chr_name] = data

    def addTMP(self, SourceFile, file_name_info):
        ignore_list = []
        table_file_name = {}
        input_file_name_info = {}
        file_temp = SourceFile.readline()
        file_temp = file_temp[1:]
        if file_temp != '':
            file_temp = file_temp.strip()
            file_temp = file_temp.split(',')
            for i in file_temp:
                i = i.split(':')
                input_file_name_info[i[0]] = int(i[1])
        print input_file_name_info
        for new_file_name in input_file_name_info:
            if new_file_name in file_name_info:
                print new_file_name + ' already in file'
                ignore_list.append(input_file_name_info[new_file_name])
            else:
                file_name_number = len(file_name_info)
                file_name_info[new_file_name] = file_name_number
                table_file_name[input_file_name_info[
                    new_file_name]] = file_name_number
        # table_file_name >> convert from input file name info to new added file name info
        # ignore_list >> contain information of overlapped file name info

        id_num = 0
        for line in SourceFile:
            if self.check_sample[0] == 1:
                id_num += 1
                if id_num > self.check_sample[1]:
                    break
            line = line.strip()
            if re.match(line, '#Splice'):
                case = 0
                continue
            elif re.match(line, '#Deletion'):
                case = 1
                continue
            elif re.match(line, '#Insertion'):
                case = 2
                continue
            elif re.match(line, '#Mutation'):
                case = 3
                continue
            line = line.split('\t')
            chromosome = line[0]
            first = int(line[1])
            if case != 2 and case != 3:
                second = int(line[2])
            else:
                second = line[2]
            files = line[3]
            strand = line[4]
            files = files.split(',')
            file_info = {}
            for i in files:
                i = i.split(':')
                file_info[int(i[0])] = int(i[1])
            if case == 0:
                tmp_data = self.splice_data[chromosome]
            elif case == 1:
                tmp_data = self.del_data[chromosome]
            elif case == 2:
                tmp_data = self.in_data[chromosome]
            elif case == 3:
                tmp_data = self.mu_data[chromosome]
            new_file_info = {}
            for i in file_info:
                if i not in ignore_list:
                    new_file_info[table_file_name[i]] = file_info[i]
            if first in tmp_data:
                check = 0
                for j, check_end in enumerate(tmp_data[first]):
                    if check_end[0] == second and check_end[2] == strand:
                        for i in new_file_info:
                            # print
                            # ':j:',j,':new_file_info:',new_file_info,':i:',i,':tmp_data[first][j]',tmp_data[first][j]
                            tmp_data[first][j][1][i] = new_file_info[i]
                        check = 1
                if check == 0:
                    tmp_data[first].append([second, new_file_info, strand])
            else:
                tmp_data[first] = [[second, new_file_info, strand]]


if __name__ == "__main__":
    if len(sys.argv) == 4:
        GFFinput_file_name = sys.argv[1]
        input_file_folder = sys.argv[2]
        output_file_name = sys.argv[3]
    elif len(sys.argv) == 3:
        GFFinput_file_name = param.empty_spl()
        input_file_folder = sys.argv[1]
        output_file_name = sys.argv[2]
    else:
        GFFinput_file_name = param.empty_spl()
        input_file_folder = param.vcf_spl_bam_dir()
        output_file_name = param.spl_out_file_name()

    spl_obj = SPLConvert(GFFinput_file_name,
                         input_file_folder, output_file_name)
    spl_obj.process()
    del spl_obj

    ##### deprecated #####
    '''
            if 'D' in CIGAR:
                    sign = 'D'
                    [del_location,dummy] = decomp(CIGAR,CIGAR_begin,sign,strand)
                    del_temp = del_data.get(chr_name)
                    for coord_info in del_location:
                            if del_temp.has_key(coord_info[0]):
                                    check = 0
                                    for i,check_end in enumerate(del_temp[coord_info[0]]):
                                            if check_end[0] == coord_info[1] and strand == check_end[2]:
                                                    if check_end[1].has_key(SAM_file_idx):
                                                            del_temp[coord_info[0]][i][1][SAM_file_idx] += 1
                                                    else:
                                                            del_temp[coord_info[0]][i][1][SAM_file_idx] = 1
                                                    check = 1
                                    if check == 0:
                                            del_temp[coord_info[0]].append([coord_info[1],{SAM_file_idx:1},strand])
                            else:
                                    del_temp[coord_info[0]] = [[coord_info[1],{SAM_file_idx:1},strand]]
                    del_data[chr_name] = del_temp
                    #for i in del_location:
                    # i[0] is starting position of del
                    # i[1] is ending position of del
                    # make a database for a del and compare the del position if they are overlap

            if 'I' in CIGAR:
                    #in_string = 'Bla Bla,' # calculate insertion string from location information and input strng
                    sign = 'I'
                    [in_location,seq_location] = decomp(CIGAR,CIGAR_begin,sign,strand)
                    in_temp = in_data.get(chr_name)
                    #print in_temp,in_location,CIGAR,chr_name
                    for i,coord_info in enumerate(in_location):
                            in_string = SEQ[seq_location[i][0]:seq_location[i][1]]
                            #print in_temp
                            if in_temp.has_key(coord_info[0]):
                                    check = 0
                                    for i,check_end in enumerate(in_temp[coord_info[0]]):
                                            #print in_temp[coord_infom[0]], check_end, coord_infom, CIGAR
                                            if check_end[0] == in_string and strand == check_end[2]:
                                                    if check_end[1].has_key(SAM_file_idx):
                                                            in_temp[coord_info[0]][i][1][SAM_file_idx] = in_temp[coord_info[0]][i][1].get(SAM_file_idx)+1
                                                    else:
                                                            in_temp[coord_info[0]][i][1][SAM_file_idx] = 1
                                                    #in_temp[coord_infom[0]][i][1] += 1
                                                    check = 1
                                    if check == 0:
                                            in_temp[coord_info[0]].append([in_string,{SAM_file_idx:1},strand])
                            else:
                                    in_temp[coord_info[0]] = [[in_string,{SAM_file_idx:1},strand]]

                    in_data[chr_name] = in_temp

            #mutation detection and recording
                if 'M' in CIGAR:
                        sign = 'M'
                        [M_location, seq_location] = decomp(CIGAR,CIGAR_begin,sign,strand)
                        mu_temp = mu_data.get(chr_name)
                    for i,coord_info in enumerate(M_location):
                                M_string = SEQ[seq_location[i][0]:seq_location[i][1]]
                                DNA_string = DNA_data[M_location[i][0]:M_location[i][1]]
                                mu_info = check_mutation(M_string,DNA_string,M_location[i][0],line)
                                if mu_info == []:
                                        continue
                                for j in mu_info:
                                        if mu_temp.has_key(j[0]):
                                                check = 0
                                                for k,check_end in enumerate(mu_temp[j[0]]):
                                                        if check_end[0] == j[1] and strand == check_end[2]:
                                                                if check_end[1].has_key(SAM_file_idx):
                                                                        mu_temp[j[0]][k][1][SAM_file_idx] += 1
                                                                else:
                                                                        mu_temp[j[0]][k][1][SAM_file_idx] = 1
                                                                #in_temp[coord_infom[0]][i][1] += 1
                                                                check = 1
                                                if check == 0:
                                                        mu_temp[j[0]].append([j[1],{SAM_file_idx:1},strand])
                                        else:
                                                mu_temp[j[0]] = [[j[1],{SAM_file_idx:1},strand]]

                    mu_data[chr_name] = mu_temp
    '''
