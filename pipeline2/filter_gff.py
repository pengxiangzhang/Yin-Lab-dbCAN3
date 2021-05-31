import sys


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

    ref_id = {}
    for raw_ref in raw_refs[1:]:
        id = raw_ref.split()[0]
        HMMER = raw_ref.split()[1]
        Hotpep = raw_ref.split()[2]
        DIAMOND = raw_ref.split()[3]
        ref_id[id] = {'HMMER': HMMER, 'Hotpep': Hotpep, 'DIAMOND': DIAMOND}

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

            if new_id in ref_id.keys():
                # print(new_id)
                reference = ref_id[new_id]
                output_id = 'ID=' + new_id + '|HMMER=' + reference['HMMER'] + '|Hotpep=' + reference[
                    'Hotpep'] + '|DIAMOND=' + reference['DIAMOND']
                new_brief = brief.replace(id, output_id, 1) + '\n'
                result.append(new_brief)

    with open(output_file, 'w') as file:
        file.writelines(result)
    file.close()


def check_result(new_id):
    with open('uniInput', 'r') as f:
        raw_checks = f.readlines()
    f.close()

    checks = []
    for raw_check in raw_checks:
        if raw_check.startswith('>'):
            raw_id = raw_check.split(' # ')[0]
            id = raw_id.split('>')[1]
            checks.append(id)

    # if new_id in checks:
    #     print('checked:' + new_id)
    # else:
    #     print(new_id)


if __name__ == '__main__':
    gff_filter(gff_file=sys.argv[1], overview_file=sys.argv[2], output_file=sys.argv[3])
