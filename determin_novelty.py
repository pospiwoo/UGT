import sys

# ref_fa_file = sys.argv[1]
# MS_result_file = sys.argv[2]
ref_fa_file = 'ref_fa/uniprot_canonical.fasta'
ref_fa_file = 'ref_fa/uniprot_can_iso.fasta'
MS_result_file = 'result/20180126_ALDER_G2957_470_560_directDIA_20180125F1_Peptide.xls'

master_seq = ''


def replaceIdenMassAA(seq_str):
    return seq_str.replace('I', 'L')
    # return seq_str.replace('L', 'I')

with open(ref_fa_file, 'r') as inFaFile:
    for line in inFaFile:
        if line.startswith('>'):
            cur_header = line[1:]
            master_seq += 'X'
            continue
        master_seq += replaceIdenMassAA(line.strip())

# print len(master_seq), master_seq[:1000]

cnt_all = 0
cnt_novel = 0
with open(MS_result_file, 'r') as inResultFile:
    for line in inResultFile:
        data = line.split('\t')
        seq = data[0]
        cnt_all += 1
        if replaceIdenMassAA(seq[1:-2]) not in master_seq:
            print seq
            cnt_novel += 1
        # if cnt_all % 1000 == 1:
        #     print cnt_all, cnt_novel
        # gene_id = data[1]
        # index_id = data[2]
        # quant_val = float(data[3])
# TLHNLVIQYASQGR', 'NaN', 'sequence_992465', '157936.015625
print cnt_all, cnt_novel
