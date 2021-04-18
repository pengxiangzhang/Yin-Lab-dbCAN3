import csv
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import sys


# Input(argv1): a paf csv
# Input(argv2): a paf csv


def FPKM_statistic_options(target_file, option):
    with open(target_file, 'r') as f:
        data = f.readlines()
    f.close()

    cazyfamily_dic = {}
    cazyFamily_count = {}
    for row in data[2:]:
        if row.startswith(option):
            row = row.replace('\n', '')
            row_infos = row.split(',')
            split_tag = option + "="
            cazyFamily_info = row_infos[0].split(split_tag)[1]
            cazyFamilys = cazyFamily_info.split('+')
            if cazyFamilys[0] == '-':
                continue
            for cazyFamily in cazyFamilys:
                if cazyFamily in cazyfamily_dic.keys():
                    cazyFamily = cazyFamily.split('(')[0]
                    cazyFamily_count[cazyFamily] += 1
                    cazyfamily_dic[cazyFamily] = cazyfamily_dic[cazyFamily] + float(row_infos[-1])
                else:
                    cazyFamily = cazyFamily.split('(')[0]
                    cazyfamily_dic[cazyFamily] = float(row_infos[-1])
                    cazyFamily_count[cazyFamily] = 1

    cazyFamily_pairs = []
    for cazyFamily in cazyfamily_dic.keys():
        pair = [cazyFamily, cazyfamily_dic[cazyFamily] / cazyFamily_count[cazyFamily]]
        cazyFamily_pairs.append(pair)
    df = pd.DataFrame(cazyFamily_pairs)
    df.to_csv('FPKM_statistic_options.csv', index=False, header=['cazyfamily_name', 'FPKM'])


def FPKM_statistic(target_file):
    with open(target_file, 'r') as f:
        data = f.readlines()
    f.close()

    FPKM_dic = {}
    cazyFamily_count = {}
    for row in data[2:]:
        if 'FPKM' in row:
            trans_info = row.split(';')[1]
            if 'HMMER=-' not in trans_info:
                cazyFamily_info = trans_info.split('HMMER=')[1].split('(')[0]
                FPKM_value = float(row.split('FPKM \"')[1].split('\";')[0])
                if cazyFamily_info in FPKM_dic:
                    FPKM_dic[cazyFamily_info] += FPKM_value
                    cazyFamily_count[cazyFamily_info] += 1
                else:
                    FPKM_dic[cazyFamily_info] = FPKM_value
                    cazyFamily_count[cazyFamily_info] = 1

    gff_statistic = []
    for cazyFamily in FPKM_dic.keys():
        processed_FPKM = FPKM_dic[cazyFamily] / cazyFamily_count[cazyFamily]
        gff_statistic.append([cazyFamily, processed_FPKM])
    df = pd.DataFrame(gff_statistic)
    df.to_csv('FPKM_statistic.csv', index=False, header=['cazyFamily', 'FPKM'])


def compare_FPKM():
    with open('FPKM_statistic_options.csv', 'r') as f1:
        refes = pd.read_csv(f1)
    f1.close()
    with open(sys.argv[1], 'r') as f2:
        mines = pd.read_csv(f2)
    f2.close()

    common_cazyFamilys = list(set(refes['cazyfamily_name']).intersection(set(mines['cazyfamily_name'])))
    diffs = []
    y = []
    try:
        for common_cazyFamily in common_cazyFamilys:
            refe_FPKM = refes.loc[refes['cazyfamily_name'] == common_cazyFamily]
            mine_FPKM = mines.loc[mines['cazyfamily_name'] == common_cazyFamily]
            dif = float(mine_FPKM['FPKM']) - float(refe_FPKM['FPKM'])
            diffs.append([common_cazyFamily, dif])
            y.append(dif)
        result = pd.DataFrame(diffs)
        result.to_csv('FPKM_diff_o.csv', index=False, header=['cazyfamily_name', 'FPKM_diff'])

        x = np.arange(0, len(diffs), 1)
        plt.plot(x, y, 'r--')
        plt.show()
        # print(list(result[1]))
    except:
        print('error from compare')


if __name__ == '__main__':
    target_file = sys.argv[2]
    FPKM_statistic_options(target_file, 'Hotpep')
    compare_FPKM()
