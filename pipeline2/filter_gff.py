import sys
import time


# Input(argv1): prodigal.gff file from dbCAN
# Input(argv2): overview.txt file from dbCAN
# Output(argv3): a filtered gff file

def gff_filter(gff_file, overview_file, output_file):
    with open(gff_file, 'r') as f:
        data = f.readlines()
    f.close()

    with open(overview_file, 'r') as f:
        raw_refs = f.readlines()
    f.close()

    ref_id = []
    for raw_ref in raw_refs[1:]:
        id = raw_ref.split()[0]
        ref_id.append(id)

    result = []
    for row in data:
        if row.startswith('#'):
            continue
        else:
            brief = row.split('\n')[0]
            columns = brief.split()
            des = columns[-1]
            id = des.split(';partial')[0]
            id_suffix = id.split('_')[1]
            new_id = columns[0] + '_' + id_suffix
            if new_id in ref_id:
                new_brief = brief.replace(columns[0], new_id, 1) + '\n'
                result.append(new_brief)

    with open(output_file, 'w') as output:
        output.writelines(result)
    output.close()


if __name__ == '__main__':
    gff_filter(gff_file=sys.argv[1], overview_file=sys.argv[2], output_file=sys.argv[3])
