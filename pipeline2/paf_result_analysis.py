import csv
import json
import sys

from pafpy import PafFile


# Input(argv1): R1 paf File
# Input(argv2): R2 paf File
# Output(argv3): R1 Read Counts
# Output(argv4): R2 Read Counts
# Output(argv5): sequence_FPKM csv file
# Output(argv6): cazyfamily_FPKM csv file


def get_brief_records(path, cazyfamilies, sequences):
    reads = {}
    with open(path) as fileobj:
        paf = PafFile(fileobj)
        for record in paf:
            sequence = record.tname.split('|')[0]
            if sequence not in sequences['name_list']:
                seq_info = {}
                seq_info['seq_read_count'] = 0
                seq_info['seq_length'] = record.tlen
                sequences['name_list'][sequence] = seq_info
                sequences['seq_amount'] += 1

            multi_cazyfamilies = record.tname.split('|')[1:]
            if len(multi_cazyfamilies) > 1:
                del multi_cazyfamilies[-1]

            for cazyfamily in multi_cazyfamilies:
                if cazyfamily in cazyfamilies['name_list']:
                    seq_array = [sequence, record.tlen]
                    if seq_array not in cazyfamilies['name_list'][cazyfamily]['cazyfamily_sequence_list']:
                        cazyfamilies['name_list'][cazyfamily]['cazyfamily_sequence_list'].append(seq_array)
                        cazyfamilies['name_list'][cazyfamily]['seq_amount'] += 1
                else:
                    cazy_info = {}
                    cazy_info['cazyfamily_read_count'] = 0
                    seq_array = [sequence, record.tlen]
                    cazy_info['cazyfamily_sequence_list'] = [seq_array]
                    cazy_info['seq_amount'] = 1
                    cazyfamilies['name_list'][cazyfamily] = cazy_info
                    cazyfamilies['family_amount'] += 1

            read_id = record.qname.split('/')[0]
            if read_id in reads.keys():
                reads[read_id]['seq_id'].append(record.tname.split('|')[0])
                reads[read_id]['cazy_family'] += multi_cazyfamilies
            else:
                read = {}
                read['cazy_family'] = multi_cazyfamilies
                read['seq_id'] = [record.tname.split('|')[0]]
                reads[read_id] = read

    output_name = path.split('.paf')[0] + '.json'
    with open(output_name, 'w') as f:
        json.dump(reads, f)
    f.close()

    return reads


def cazy_family_counter(cazyfamilies, reads_R1, reads_R2, seq_id_intersec, seq_id_list_R1_subtraction,
                        seq_id_list_R2_subtraction):
    for read_id in seq_id_intersec:
        read_in_R1 = reads_R1[read_id]
        for read1_cazyfamily in read_in_R1['cazy_family']:
            cazyfamilies['name_list'][read1_cazyfamily]['cazyfamily_read_count'] += 1
        read_in_R2 = reads_R2[read_id]
        for read2_cazyfamily in read_in_R2['cazy_family']:
            if read2_cazyfamily not in read_in_R1['cazy_family']:
                cazyfamilies['name_list'][read2_cazyfamily]['cazyfamily_read_count'] += 1

    for read_id in seq_id_list_R1_subtraction:
        read_in_R1 = reads_R1[read_id]
        for read1_cazyfamily in read_in_R1['cazy_family']:
            cazyfamilies['name_list'][read1_cazyfamily]['cazyfamily_read_count'] += 1

    for read_id in seq_id_list_R2_subtraction:
        read_in_R2 = reads_R2[read_id]
        for read2_cazyfamily in read_in_R2['cazy_family']:
            cazyfamilies['name_list'][read2_cazyfamily]['cazyfamily_read_count'] += 1

    for cazyfamily_name in cazyfamilies['name_list']:
        cazyfamily = cazyfamilies['name_list'][cazyfamily_name]
        seq_amount = len(cazyfamily['cazyfamily_sequence_list'])
        sum_seq_length = 0
        for seq in cazyfamily['cazyfamily_sequence_list']:
            sum_seq_length += seq[1]
        average_length = sum_seq_length / seq_amount
        cazyfamily['average_seq_length'] = average_length


def sequence_counter(sequences, reads_R1, reads_R2, seq_id_intersec, seq_id_list_R1_subtraction,
                     seq_id_list_R2_subtraction):
    for read_id in seq_id_intersec:
        read_in_R1 = reads_R1[read_id]
        read1_sequences = read_in_R1['seq_id']
        for read1_sequence in read1_sequences:
            sequences['name_list'][read1_sequence]['seq_read_count'] += 1

        read_in_R2 = reads_R2[read_id]
        read2_sequences = read_in_R2['seq_id']
        for read2_sequence in read2_sequences:
            if read2_sequence not in read1_sequence:
                sequences['name_list'][read2_sequence]['seq_read_count'] += 1

    for read_id in seq_id_list_R1_subtraction:
        read_in_R1 = reads_R1[read_id]
        read1_sequences = read_in_R1['seq_id']
        for read1_sequence in read1_sequences:
            sequences['name_list'][read1_sequence]['seq_read_count'] += 1

    for read_id in seq_id_list_R2_subtraction:
        read_in_R2 = reads_R2[read_id]
        read2_sequences = read_in_R2['seq_id']
        for read2_sequence in read2_sequences:
            sequences['name_list'][read2_sequence]['seq_read_count'] += 1


def sequence_FPKM(sequences, amount_all_reads):
    for sequence_name in sequences['name_list'].keys():
        sequence = sequences['name_list'][sequence_name]
        convert_all_reads = amount_all_reads / pow(10, 6)
        convert_length = sequence['seq_length'] / 1000
        FPKM = sequence['seq_read_count'] / (convert_all_reads * convert_length)
        sequence['FPKM'] = FPKM

    # with open(sys.argv[3], 'w') as f:
    #     json.dump(sequences, f)
    # f.close()
    data_file = open(sys.argv[5], 'w')
    csv_writer = csv.writer(data_file)
    count = 0

    seq_ids = sequences['name_list']
    for seq_id in seq_ids.keys():
        if count == 0:
            # Writing headers of CSV file
            header = ['seq_id']
            header += seq_ids[seq_id].keys()
            csv_writer.writerow(header)
            count += 1
        row = [seq_id] + list(seq_ids[seq_id].values())
        csv_writer.writerow(row)
    data_file.close()


def cazyfamily_FPKM(cazyfamilies, amount_all_reads):
    for cazyfamily_name in cazyfamilies['name_list'].keys():
        cazyfamily = cazyfamilies['name_list'][cazyfamily_name]
        convert_all_reads = amount_all_reads / pow(10, 6)
        convert_length = cazyfamily["average_seq_length"] / 1000
        FPKM = cazyfamily["cazyfamily_read_count"] / (convert_length * convert_all_reads)
        cazyfamily['FRKM'] = FPKM

    # with open(sys.argv[4], 'w') as f:
    #     json.dump(cazyfamilies, f)
    # f.close()
    data_file = open(sys.argv[6], 'w')
    csv_writer = csv.writer(data_file)
    count = 0

    cazyfamily_names = cazyfamilies['name_list']
    for cazyfamily_name in cazyfamily_names.keys():
        if count == 0:
            # Writing headers of CSV file
            header = ['cazyfamily_name']
            cazy_title = list(cazyfamily_names[cazyfamily_name].keys())
            del cazy_title[1]
            header += cazy_title
            csv_writer.writerow(header)
            count += 1
        cazy_info = list(cazyfamily_names[cazyfamily_name].values())
        del cazy_info[1]
        row = [cazyfamily_name] + cazy_info
        csv_writer.writerow(row)
    data_file.close()


def main():
    cazyfamilies = {}
    cazyfamilies['name_list'] = {}
    cazyfamilies['family_amount'] = 0

    sequences = {}
    sequences['name_list'] = {}
    sequences['seq_amount'] = 0

    path_R1 = sys.argv[1]
    reads_R1 = get_brief_records(path_R1, cazyfamilies, sequences)

    path_R2 = sys.argv[2]
    reads_R2 = get_brief_records(path_R2, cazyfamilies, sequences)

    seq_id_list_R1 = reads_R1.keys()
    seq_id_list_R2 = reads_R2.keys()
    seq_id_intersec = list(set(seq_id_list_R1).intersection(set(seq_id_list_R2)))
    seq_id_list_R1_subtraction = list(set(seq_id_list_R1).difference(set(seq_id_list_R2)))
    seq_id_list_R2_subtraction = list(set(seq_id_list_R2).difference(set(seq_id_list_R1)))

    amount_all_reads = int(sys.argv[3])
    cazy_family_counter(cazyfamilies, reads_R1, reads_R2, seq_id_intersec, seq_id_list_R1_subtraction,
                        seq_id_list_R2_subtraction)
    sequence_counter(sequences, reads_R1, reads_R2, seq_id_intersec, seq_id_list_R1_subtraction,
                     seq_id_list_R2_subtraction)
    sequence_FPKM(sequences, amount_all_reads)
    cazyfamily_FPKM(cazyfamilies, amount_all_reads)


if __name__ == '__main__':
    main()
