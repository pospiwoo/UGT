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


class PGEvent:

    def __init__(self, master_seq, ref_fa_file, fasta_database_name, location_file_name, ms_result_file_name):
        self.ForwardCode = formats.initForwardCodon()
        self.ReverseCode = formats.initReverseCodon()
        self.novelgene
        self.fusiongene
        self.alternativesplice
        self.novelsplice
        self.insertion
        self.mutation
        self.deletion
        self.translatedutr
        self.exonboundary
        self.novelexon
        self.geneboundary
        self.frameshift
        self.reversestrand
        self.pseudo
        self.IG
        self.indelcount
        self.na

    def CompareLocation(self, FS, FE, SS, SE):
        case = 0
        if FE < SS or FS > SE:
            case = 10
        elif FS > SS:
            if FE > SE:
                case = 1
            elif FE < SE:
                case = 2
            elif FE == SE:
                case = 3
        elif FS < SS:
            if FE > SE:
                case = 4
            elif FE < SE:
                case = 5
            elif FE == SE:
                case = 6
        elif FS == SS:
            if FE > SE:
                case = 7
            elif FE < SE:
                case = 8
            elif FE == SE:
                case = 9

        return case

    def CheckStopCodon(self, dna_string, strand):
        if strand == 0:
            dna_string = dna_string[::-1]
        for i in range(len(dna_string)/3):
            if strand == 1:
                if self.ForwardCode.get(dna_string[i*3:(i+1)*3]) == 'X':
                    return False
            else:
                if self.ReverseCode.get(dna_string[i*3:(i+1)*3]) == 'X':
                    return False
        return True

    def compare_start(self, item1, item2):
        if int(item1[0]) < int(item2[0]):
            return -1
        elif int(item1[0]) > int(item2[0]):
            return 1
        else:
            return 0

    def compare_pseudo(self, item1, item2):
        if item1[0][0] < item2[0][0]:
            return -1
        elif item1[0][0] > item2[0][0]:
            return 1
        else:
            return 0

    def sort_location(self, item1, item2):
        if item1[-1] == 1:
            position1 = item1[0][0]
        else:
            position1 = item1[0][-1]
        if item2[-1] == 1:
            position2 = item2[0][0]
        else:
            position2 = item2[0][-1]

        if position1 < position2:
            return -1
        elif position1 > position2:
            return 1
        else:
            return 0

    def sort_by_chromosome(self, item1, item2):
        position1 = item1[1]
        position2 = item2[1]

        if position1 < position2:
            return -1
        elif position1 > position2:
            return 1
        else:
            return 0

    def WriteLocation(self, location):
        beg = location[0][0]
        end = location[0][1]
        string = ''
        for i in range(len(beg)):
            if beg[i] < 0:
                if beg[i] < -1000:
                    string += ('M' + str(end[i] - beg[i]))
                else:
                    string += ('I' + str(end[i] - beg[i]))
            else:
                string += (str(beg[i]) + '-' + str(end[i]))
            if i != len(beg) - 1:
                string += (';')
        return string

    def CheckRelation(self, start, end, gframe, ref_seq):
        exon_case = 10
        for cds in ref_seq[4]:
            if end < cds[3]:
                continue
            elif start > cds[4]:
                continue
            exon_case = self.CompareLocation(start, end, cds[3], cds[4])
            if cds[7] == '.':
                cds_gframe = 0
            elif ref_seq[3][6] == '+':
                cds_gframe = (cds[3] % 3 + int(cds[7])) % 3
            else:
                cds_gframe = (cds[3] %
                              3 + (cds[4] - cds[3] - int(cds[7])) % 3) % 3
            if cds_gframe != gframe:
                exon_case = exon_case - 10
        return exon_case

    def Bigger(self, x, y):
        if x > y:
            return x
        else:
            return y

    def Recruit(self, current_location, unknown_nonsplice, known, unknown_splice, unknown_indel, indicator, strand):
        group = [[current_location], []]
        NearConstant = 1000
        start = current_location[0][0][0]
        end = current_location[0][1][-1]
        start_gframe = start % 3
        end_gframe = (current_location[0][
                      0][-1] - (sum(current_location[0][1][:-1]) - sum(current_location[0][0][:-1]) % 3)) % 3

        # if start < 0 or end < 0:
        #    print 'negative start or end error',
        #    print group,start,end
        #    return group

        if start < 0:
            start = current_location[0][0][1]
            start_gframe = (
                start - (current_location[0][1][0] - current_location[0][0][0])) % 3
        elif end < 0:
            end = current_location[0][1][-2]
            end_gframe = (current_location[0][
                          0][-1] - (sum(current_location[0][1][:-2]) - sum(current_location[0][0][:-2]) % 3)) % 3
        # grouping the splice event ( if the two different peptide share it's
        # junction, combine them together )
        if indicator == 1 and unknown_splice != None:
            index = 0
            start_c = start
            end_c = end
            while index < len(unknown_splice):
                candidate = unknown_splice[index]
                if candidate[5] != strand:
                    index += 1
                    continue
                if start_c > candidate[0][1][-1] + NearConstant:
                    index += 1
                    continue
                if end_c < candidate[0][0][0] - NearConstant:
                    break
                isOverlap = True
                check = False
                for i in range(len(current_location[0][0])-1 if len(current_location[0][0]) < len(candidate[0][0]) else len(candidate[0][0])-1):
                    if current_location[0][1][i] != candidate[0][1][i] or current_location[0][0][i+1] != candidate[0][0][i+1]:
                        isOverlap = False
                    else:
                        check = True
                if isOverlap and check and current_location != candidate and start_gframe == candidate[0][0][0] % 3:
                    group[0].append(candidate)
                    unknown_splice.remove(candidate)
                    index -= 1
                else:
                    index += 1

        elif indicator == 2 and unknown_indel != None:
            index = 0
            start_c = start
            end_c = end
            while index < len(unknown_indel):
                candidate = unknown_indel[index]
                if candidate[5] != strand:
                    index += 1
                    continue
                if start_c > candidate[0][1][-1] + NearConstant:
                    index += 1
                    continue
                if end_c < candidate[0][0][0] - NearConstant:
                    break
                isOverlap = True
                check = False
                for i in range(len(current_location[0][0])-1 if len(current_location[0][0]) < len(candidate[0][0]) else len(candidate[0][0])-1):
                    if current_location[0][1][i] != candidate[0][1][i] or current_location[0][0][i+1] != candidate[0][0][i+1]:
                        isOverlap = False
                    else:
                        check = True
                if isOverlap and check and current_location != candidate and start_gframe == candidate[0][0][0] % 3:
                    group[0].append(candidate)
                    unknown_indel.remove(candidate)
                    index -= 1
                else:
                    index += 1

        # grouping the event
        # find relative unknown sequence
        index = 0
        start_c = start
        end_c = end
        if unknown_nonsplice != None:
            while index < len(unknown_nonsplice):
                candidate = unknown_nonsplice[index]
                if candidate[5] != strand:
                    index += 1
                    continue
                candidate_end_gframe = (candidate[0][
                                        0][-1] - (sum(candidate[0][1][:-1]) - sum(candidate[0][0][:-1]) % 3)) % 3
                if start_c > candidate[0][1][-1] + NearConstant:
                    index += 1
                    continue
                if end_c < candidate[0][0][0] - NearConstant:
                    break
                if candidate[0][1][-1] < current_location[0][1][0] and start_gframe == candidate[0][0][0] % 3:
                    group[0].append(candidate)
                    unknown_nonsplice.remove(candidate)
                    start_c = candidate[0][0][0]
                    index -= 1
                elif candidate[0][0][0] > current_location[0][0][-1] and end_gframe == candidate_end_gframe:
                    group[0].append(candidate)
                    unknown_nonsplice.remove(candidate)
                    end_c = candidate[0][1][-1]
                else:
                    index += 1

        index = 0
        start_c = start
        end_c = end
        while index < len(known):
            candidate = known[index]
            candidate_end_gframe = (candidate[0][
                                    0][-1] - (sum(candidate[0][1][:-1]) - sum(candidate[0][0][:-1]) % 3)) % 3
            if start_c > candidate[0][1][-1] + NearConstant:
                index += 1
                continue
            if end_c < candidate[0][0][0] - NearConstant:
                break
            if candidate[0][1][-1] < current_location[0][1][0] and start_gframe == candidate[0][0][0] % 3:
                group[1].append(candidate)
                # known.remove(candidate)
                # start_c = candidate[0][0][0]
            elif candidate[0][0][0] > current_location[0][0][-1] and end_gframe == candidate_end_gframe:
                group[1].append(candidate)
                # known.remove(candidate)
                # end_c = candidate[0][1][-1]
            index += 1
        return group

    def FindNonSpliceEvent(self, current_location, ref_seq, event):
        start = current_location[0][0][0]
        end = current_location[0][1][-1]
        gframe = start % 3
        if ref_seq[3][6] == '+':
            ref_strand = 1
        else:
            ref_strand = 0
        if ref_strand != int(current_location[5]):
            event[-1] = self.Bigger(event[-1], 1)
            return event  # reverse
        case = self.CompareLocation(current_location[0][0][0], current_location[
            0][1][-1], ref_seq[0], ref_seq[1])
        exon_case = self.CheckRelation(start, end, gframe, ref_seq)

        if case == 1 or case == 5:
            event[6] = 1
        if len(ref_seq[6]) > 0:
            for utr in ref_seq[6]:
                # if utr[3] < start < utr[4] or utr[3] < end < utr[4]:
                if utr[3] <= start and end <= utr[4]:
                    event[7] = 1
        if len(ref_seq[5]) > 0:
            for utr in ref_seq[5]:
                if utr[3] <= start and end <= utr[4]:
                    # if utr[3] < start < utr[4] or utr[3] < end < utr[4]:
                    event[7] = 1
        if exon_case in [0, 10]:
            event[9] = 1
        elif exon_case in [1, 4, 5, 6, 7]:
            utr = ref_seq[:]
            utr[4] = ref_seq[6]
            if self.CheckRelation(start, end, gframe, utr) not in [0, 10]:
                event[7] = 1
            utr[4] = ref_seq[5]
            if self.CheckRelation(start, end, gframe, utr) not in [0, 10]:
                event[7] = 1
            event[8] = 1
        elif exon_case < 0:
            event[10] = 2
        else:
            # print 'exon_case error', exon_case, current_location, ref_seq,
            # event
            event[10] = 0
        # print event, case,exon_case, ref_seq,current_location
        return event

    def PepString(self, peptide, location):
        beg = location[0][0]
        end = location[0][1]
        index = 0
        pivot = 2
        if location[-1] == 1:
            new_string = peptide[0:2]
            for i in range(len(beg) - 1):
                index += (end[i] - beg[i])
                if index < 3 or end[i] < 0:
                    continue
                if peptide[0] != '-' and i == 0:
                    index = index - 3
                new_string = new_string + peptide[pivot:index / 3 + 2] + ':'
                pivot = index / 3 + 2
            new_string += peptide[pivot:]
        else:
            new_string = peptide[-2:]
            for i in range(len(beg) - 1):
                index += (end[i] - beg[i])
                if index < 3 or end[i] < 0:
                    continue
                if peptide[-1] != '-' and i == 0:
                    index = index - 3
                new_string = ':' + \
                    peptide[-(index / 3 + 2):-pivot] + new_string
                pivot = index / 3 + 2
            new_string = peptide[:-pivot] + new_string

        return new_string

    def CheckPseudo(self, event_group, chr):
        ret_val = [False, '', '']
        if pseudo_gff_dic.get(chr) == None:
            return ret_val
        pseudo_ref = self.pseudo_gff_dic.get(chr)
        for pseudo in pseudo_ref:
            if pseudo[1][-1] < event_group[0][0][0][0][0]:
                continue
            elif pseudo[0][0] > event_group[0][0][0][1][-1]:
                break
            for index in range(len(pseudo[0])):
                for splice_index in range(len(event_group[0][0][0][0])):
                    # print
                    # index,pseudo,pseudo[0][index],pseudo[1][index],event_group[0][0][0][0][splice_index],event_group[0][0][0][1][splice_index]
                    pseudo_event = self.CompareLocation(event_group[0][0][0][0][splice_index], event_group[
                        0][0][0][1][splice_index], pseudo[0][index], pseudo[1][index])
                    if pseudo_event != 10:
                        ret_val[0] = True
                        ret_val[1] = pseudo[2]
                        ret_val[2] = pseudo[3]

        return ret_val

    def ChooseEvent(self, event, event_group, chr, IG_check):
        if IG_check:
            sevent = 'IG gene'
            self.IG += 1
        # and event[-1] != 2:
        elif sum(event[:]) == 0 or (sum(event[:]) - event[1]) == 0:
            sevent = 'novel gene'
            pseudogene = self.CheckPseudo(event_group, chr)
            if pseudogene[0]:
                sevent = pseudogene[2]
                sevent = 'pseudogene'
                self.pseudo += 1
            else:
                self.novelgene += 1
        # event[0] in [9, 10]:# and event[1] == 0 and event[2] == 0:
        elif event[1] == 3:
            sevent = 'fusion gene'
            self.fusiongene += 1
        elif event[7] != 0 and event[8] != 0:
            sevent = 'translated UTR'
            self.translatedutr += 1
        elif event[8] != 0:
            sevent = 'exon boundary'
            self.exonboundary += 1
        elif event[10] != 0:
            sevent = 'frame shift'
            self.frameshift += 1
        elif event[0] == 5:
            sevent = 'alternative splice'
            self.alternativesplice += 1

        elif event[7] != 0:
            sevent = 'translated UTR'
            self.translatedutr += 1
        elif (event[0] > 0)and event[5] == 0:
            sevent = 'novel splice'
            self.novelsplice += 1
        elif event[3] != 0:
            sevent = 'insertion'
            self.indelcount += 1
            self.insertion += 1
        elif event[4] != 0:
            sevent = 'mutation'
            self.indelcount += 1
            self.mutation += 1
        elif event[5] != 0:
            sevent = 'deletion'
            self.indelcount += 1
            self.deletion += 1

        elif event[9] != 0:
            sevent = 'novel exon'
            pseudogene = self.CheckPseudo(event_group, chr)
            if pseudogene[0]:
                sevent = pseudogene[2]
                sevent = 'pseudogene'
                self.pseudo += 1
            else:
                self.novelexon += 1
        elif event[6] != 0:
            sevent = 'gene boundary'
            self.geneboundary += 1
        # elif event[10] != 0:
        #    sevent = 'frame shift'
        #    frameshift += 1
        elif event[11] != 0:
            sevent = 'reverse strand'
            pseudogene = self.CheckPseudo(event_group, chr)
            if pseudogene[0]:
                sevent = pseudogene[2]
                sevent = 'pseudogene'
                self.pseudo += 1
            else:
                self.reversestrand += 1
        else:
            sevent = 'NA'
            self.na += 1
        return sevent

    def GbrowserLink(self, tmpGbrowserCoor):
        tmpGbrowserLink = "http://genome.ucsc.edu/cgi-bin/hgTracks?db="
        tmpGbrowserLink += db_str
        tmpGbrowserLink += "&position="
        tmpGbrowserLink += tmpGbrowserCoor
        tmpGbrowserLink += "&hgt.customText=http://bix.ucsd.edu/tmp/"
        tmpGbrowserLink += usa_gffile_str
        return tmpGbrowserLink

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
    if len(sys.argv) == 5:
        novel_location_filename = sys.argv[1]
        known_location_filename = sys.argv[2]
        gff_filename = sys.argv[3]
        event_filename = sys.argv[4]
    elif len(sys.argv) == 0:
        novel_location_filename = 'Location_Normal_colon_SV_merged_somatic.txt'
        known_location_filename = 'Location_Combine123_RefSeq_recalculated_Correction.txt'
        gff_filename = 'Homo_sapiens.GRCh37.88_formatted.gff'
        event_filename = sys.argv[4]
    else:
        novel_location_filename = param.location_file()
        known_location_filename = param.location_file()
        gff_filename = param.gff_file_name()
        event_filename = param.event_filename()

    event_obj = PGEvent(novel_location_filename,
                        known_location_filename,
                        gff_filename,
                        event_filename)
    event_obj.process()
    del event_obj
