'''
DESCRIPTION
This script takes the data output by autoccivr2.py and counts the gene pairs with biotype combinations including coding_coding, coding_pseudogene, pseudogene_coding and psuedogene_psuedogene.

PROCEDURE
1. Run this script after running autoccivr2.py in your project directory.
2. The tables of gene pairs with each biotype combination are stored in 'data/[dataset name]/biotype_combination'.
3. When all analyses have been completed, the results are summarised for each species and mode and stored in 'result'. The name of summary files end with '_pg'.
'''

import pandas as pd
import os
import re
from autoccivr2 import RANGE_DICT

# dirpath直下のsubdir pathを一括取得してリストとして返す関数
def get_subdir_path(dirpath):
    objects = os.listdir(dirpath)
    subdir = [os.path.join(dirpath, f) for f in objects if os.path.isdir(os.path.join(dirpath, f))]
    return subdir

# data直下のdirpathを取得
sp_dirlist = get_subdir_path('data')
print(sp_dirlist)

# resultディレクトリをdata直下に作成
result_dir = 'result'
os.makedirs(result_dir, exist_ok=True)

def biotypes_count(mode):

    # 各species dirに対して
    for i in sp_dirlist:

        print('----------')
        # 結果出力先の準備
        # result直下にspecies名のディレクトリを作成
        data_name = os.path.splitext(os.path.basename(i))[0]
        print('Analysing '+data_name)

        # i直下に biotype_combination というフォルダを作る
        biocombi_dir = os.path.join(i,'biotype_combination')
        os.makedirs(os.path.join(biocombi_dir),exist_ok=True)
        mode_dir = os.path.join(biocombi_dir,mode)
        os.makedirs(os.path.join(mode_dir),exist_ok=True)
    
        # speciesディレクトリ直下に4つのディレクトリ（coding_coding, coding_pseudogene, pseudogene_coding, pseudogene_pseudogene）を作成
        cc_dir = os.path.join(mode_dir,'coding_coding')
        cp_dir = os.path.join(mode_dir,'coding_pseudogene')
        pc_dir = os.path.join(mode_dir,'pseudogene_coding')
        pp_dir = os.path.join(mode_dir,'pseudogene_pseudogene')
        os.makedirs(cc_dir,exist_ok=True)
        os.makedirs(cp_dir,exist_ok=True)
        os.makedirs(pc_dir,exist_ok=True)
        os.makedirs(pp_dir,exist_ok=True)

        # species summury を作成する
        sp_summary = pd.DataFrame(columns=['coding_coding','coding_pseudogene','pseudogene_coding','pseudogene_pseudogene'])

        # 解析用データを取得
        # mode2直下のdirpathを取得
        output_dirlist = get_subdir_path(os.path.join(i,mode))
        
        # 各階層dirに対して
        for j in output_dirlist:
            # 階層名取得、名前変更 -> table_mode2_(xxxx,xxxx)
            table_filename = os.path.basename(j).replace('ccivr_output_','table_') + '.csv'
            range_name = re.findall('(?<=\().+?(?=\))',table_filename)[0]
            range_idx = [k for k, v in RANGE_DICT.items() if v == range_name][0]
            # Table.csv のpathを取得し、df化
            table_path = os.path.join(j,'Table.csv')
            df = pd.read_csv(table_path, header=0)
            # gene biotype でソート
            df_cc = df[(df['gene_biotype']=='protein_coding') & (df['_gene_biotype']=='protein_coding')].reset_index()
            df_cp = df[(df['gene_biotype']=='protein_coding') & (df['_gene_biotype'].str.contains('pseudogene'))].reset_index()
            df_pc = df[(df['gene_biotype'].str.contains('pseudogene')) & (df['_gene_biotype']=='protein_coding')].reset_index()
            df_pp = df[(df['gene_biotype'].str.contains('pseudogene')) & (df['_gene_biotype'].str.contains('pseudogene'))].reset_index()
            # カウント
            sp_summary.at[range_idx,'coding_coding'] = df_cc['id'].nunique()
            sp_summary.at[range_idx,'coding_pseudogene'] = df_cp['id'].nunique()
            sp_summary.at[range_idx,'pseudogene_coding'] = df_pc['id'].nunique()
            sp_summary.at[range_idx,'pseudogene_pseudogene'] = df_pp['id'].nunique()
            # csv出力
            df_cc.to_csv(os.path.join(cc_dir,table_filename),index=False)
            df_cp.to_csv(os.path.join(cp_dir,table_filename),index=False)
            df_pc.to_csv(os.path.join(pc_dir,table_filename),index=False)
            df_pp.to_csv(os.path.join(pp_dir,table_filename),index=False)

            print('COUNT FINISH '+j)

        sp_summary.sort_index().to_csv(os.path.join(result_dir,'summary_'+data_name+'_'+mode+'_pg.csv'))

def main():
    biotypes_count('mode2')
    biotypes_count('mode3')

main()