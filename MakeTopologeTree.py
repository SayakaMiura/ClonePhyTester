import Functions
import sys
import os

Nwk=sys.argv[1]
Rscript=sys.argv[2]
RunR='library(phytools)\nlibrary(phangorn)\nItree=read.newick(\''+Nwk.replace('\\','\\\\')+'\')\nItree1=di2multi(Itree)\nwrite.tree(Itree1,file=\''+Nwk.replace('\\','\\\\')[:-4]+'multi.nwk'+'\')\n'
Functions.GetOut('RunR.r',RunR)		
os.system('\"'+Rscript+'\" RunR.r')	