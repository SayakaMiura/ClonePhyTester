import sys
import os
import shutil
import Functions
import glob
import pandas as pd
import ast

ID2Res = {}
for f in sys.argv[2:]:
    l = open(f).readline().strip().split('\t')
    ID2Res[l[0]+'\t'+l[1]] = ast.literal_eval(l[2])

Input = sys.argv[1]
SeqLs, Rpath, MLTEDpath, Summary, InputParse = Functions.parse_input(Input)

#####MLTED, TreeVec, RF, and error rate of ordering mutations (scores) reported in ref [1] in README.txt####
xl = pd.ExcelFile('Summary.xlsx')

####make output excel with score of "YournMethod"####
#Out: output file name
Out=Input[:-6]+'_summary.xlsx'
Functions.make_output_excel(Out,ID2Res,xl,Rpath)

#make plot
os.system(Rpath+' Plots.r')
shutil.copy2('plot.jpg',Out[:-5]+'.jpg')
os.remove('plot.jpg')
os.remove('All.xlsx')