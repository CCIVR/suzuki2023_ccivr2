'''
DESCRIPTION
This script automatically performs both mode2 and mode3 of CCIVR2 analyses, with different ranges, on multiple datasets.
Each dataset is analysed repeatedly for each class of the analysis range (i.e. AS_TSS or AS_TES position range) defined in "RANGE_DICT".

REQUIREMENT
*Python 3.8+
*pandas

PROCEDURE
1. If you want to change the range classes for your analysis, edit the value of "RANGE_DICT" in this script (otherwise it is unnecessary).
2. Create a directory named "data" under your project directory and place the datasets in the "data" directory. For details on the format of the dataset, see the README of CCIVR2.
3. Run this script in your project directory.
4. Under the "data" directory, a new directory is created with the same name as each dataset, and the results of each analysis will be placed in it.
5. The percentages of the number of genes extracted per class in each mode are compiled in a summary file for each species and stored in the directory named "result".
'''

import subprocess
import shutil
import glob
import os
import pandas as pd

RANGE_DICT = {
    1:'-5000,-4501',
    2:'-4500,-4001',
    3:'-4000,-3501',
    4:'-3500,-3001',
    5:'-3000,-2501',
    6:'-2500,-2001',
    7:'-2000,-1501',
    8:'-1500,-1001',
    9:'-1000,-501',
    10:'-500,-1',
    11:'0,500',
    12:'501,1000',
    13:'1001,1500',
    14:'1501,2000',
    15:'2001,2500',
    16:'2501,3000',
    17:'3001,3500',
    18:'3501,4000',
    19:'4001,4500',
    20:'4501,5000'
}

def culc_value(file,output,mode,range):
    stdin = mode + '\n' + range
    subprocess.run(['ccivr2',file,'-o',output],input=stdin,text=True)

def one_species(file_path):

    # species data 名のディレクトリを作成する
    data_name = os.path.splitext(os.path.basename(file_path))[0]
    dir_path = os.path.join('data', data_name)
    os.makedirs(dir_path, exist_ok=True)

    # species summury を作成する
    sp_summary = pd.DataFrame(columns=['%mode2','%mode3'])

    # modeごとディレクトリ作成
    mode2_dir = os.path.join(dir_path,'mode2')
    mode3_dir = os.path.join(dir_path,'mode3')
    os.makedirs(mode2_dir, exist_ok=True)
    os.makedirs(mode3_dir, exist_ok=True)


    # 各階級ごとにmode2でccivr2をまわす
    for c in RANGE_DICT.keys():
        # 計算
        culc_value(file=file_path,output=mode2_dir,mode='2',range=RANGE_DICT[c])
        # Summary出力パスを取得
        sum_path2 = os.path.join(mode2_dir,'ccivr_output_mode2_('+RANGE_DICT[c]+')','Summary.csv')
        # ％値を取得してspecies summary に書き込み
        if os.path.exists(sum_path2):
            summary2 = pd.read_csv(sum_path2,header=None,index_col=0)
            rate2 = str(summary2.at['extracted genes',2])
            sp_summary.at[c,'%mode2'] = rate2
        else:
            sp_summary.at[c,'%mode2'] = 0
        print('-----------------------')

    # 同様にmode3でccivr2をまわす
    for c in RANGE_DICT.keys():

        culc_value(file=file_path,output=mode3_dir,mode='3',range=RANGE_DICT[c])
        sum_path3 = os.path.join(mode3_dir,'ccivr_output_mode3_('+RANGE_DICT[c]+')','Summary.csv')

        if os.path.exists(sum_path3):
            summary3 = pd.read_csv(sum_path3,header=None,index_col=0)
            rate3 = (summary3.at['extracted genes',2])
            sp_summary.at[c,'%mode3'] = rate3
        else:
            sp_summary.at[c,'%mode3'] = 0
        print('-----------------------')

    # species summary をcsv出力
    result_dir = 'result'
    os.makedirs(result_dir, exist_ok=True)
    sp_summary.to_csv(os.path.join(result_dir,'summary_'+data_name+'.csv'))

    # 解析済みファイルを各speciesのディレクトリへ移動
    new_path = shutil.move(file_path,dir_path)

    print('FINISH '+data_name)


def main():

    # 全speciesのファイルを取得
    path_list = glob.glob(os.path.join('data','*.csv'))
    print(path_list)
    for i in path_list:
        print('***********************')
        print('Reading '+i)
        
        one_species(i)

main()