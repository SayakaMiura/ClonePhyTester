import sys
import os
import shutil
import Functions
import glob
import pandas as pd

####Parse input file####
#SeqLs: list lines in "Data" file given in your input file. See SequenceList.txt for the format of "Data."
#Rpath: path to Rscript
#MLTEDpath: path to ./main (MLTED software)
#Summary: output errors in the input file
#InputParse: 'y' or 'n', in which 'y' indictes that there is errors in the input file and the clonephytester program terminate here.
Input=sys.argv[1] #*.input
SeqLs, Rpath, MLTEDpath, Summary, InputParse = Functions.parse_input(Input)

####File of true clone sequences####
#Sim2TrueFile: a dictionary to assign a file of true clone sequences for each dataset.
#Sim2TrueFile: a key is "Sim" in your "Data" file. See SequenceList.txt for the format of "Data."
#Sim2TrueFile: "ID" is "Data" in you "Data" file.
Sim2TrueFile={'SimA':'ID'}

if InputParse!='y': print (Summary)
else:
            ####Scores for inferred clone phylogeny is stored in ID2Res####
            #Scores are computed for each given inference in the "Data" file.
            ID2Res={}
            for In0 in SeqLs:
               ####Parse a file of inferred clone sequences####
               #Sim: "Sim" in your "Data" file. See SequenceList.txt for the format of "Data."
               #ID: "Data" in your "Data" file.
               #TruMeg: A file (.meg format) with true clone sequences for a dataset.
               #InfMeg: A file (.meg format) with inferred clone sequences for a dataset.
               #Summary: if there are errors in files of true and inferred clone sequences, the errors are shown.
               #InfParse: y' or 'n', in which 'y' indictes that there is errors in the input files and scores are not computed for this dataset.
               Sim,ID,TruMeg,InfMeg,Summary,InfParse = Functions.parse_inferred_file(In0,Summary,Sim2TrueFile)
               if InfParse!='y': print (Summary)
               if InfParse=='y':
                      ####Scores for this dataset are stored in ResAll####
                      ResAll={}						  
                      ####Copy mega files####
                      #InfMegTMP: Copied file of inferred clone sequences. Redundant sequences are removed, if any.
                      #TruMegTMP: Copied file of true clone sequences.
                      InfMegTMP=InfMeg.split('/')[-1]
                      TruMegTMP=TruMeg.split('/')[-1]
                      shutil.copy2(InfMeg, InfMegTMP)
                      shutil.copy2(TruMeg, TruMegTMP)
                      Functions.RmRedunSeq(InfMegTMP)
                      InfMegTMP=InfMegTMP[:-4]+'_NoRedun.meg'
                      
                      ####Annotate inferred clones####
                      #inferred clone sequences are paired with the most similar true clone sequences.
                      #This clone annotation is output in the same directory of the mega file of inferred clone sequences.
                      #PairList: list the line of clone annotation file
                      #Inf2TruClo: dictionary of inferred clones matched with true clones
                      #TrueCloneLs: list of true clone
                      #Tru2Inf2GE: the number of SNV assignment errors (GE) for each inferred clone
                      #TrueClone_Ls: list of true clones.
                      PairList,Inf2TruClo,TrueCloneLs,Tru2Inf2GE,TrueClone_Ls=Functions.CloneAnno(TruMegTMP,InfMegTMP)
                      shutil.copy2(InfMegTMP[:-4]+'_TrueCloneAnnotation.txt',InfMeg[:-4]+'_TrueCloneAnnotation.txt')
                      
                      ####infer true and inferred clone phylogeny####
                      #use MEGACC to infer bifurcating phylogeny.
                      #use phangorn pachage in R to infer multifurcating phylogeny from bifurcating phylogeny.
                      #TrueMultiTree: nwk format true clone phylogeny.
                      #InfMultiTree: nwk format inferred clne phylogeny. The number of clones was made to be the same as that in its true clone phylogeny.
                      TrueMultiTree, InfMultiTree=Functions.Build_true_and_inferred_phylogeny(TrueClone_Ls,PairList,TruMegTMP,InfMegTMP,Rpath)
                      
                      ####Compute RF####
                      #RF: RF distance calcuated by using PhyloNet_3.6.1.jar
                      if TrueMultiTree!='NA' and InfMultiTree!='NA':
                          RF=Functions.Compute_RF(TrueMultiTree, InfMultiTree)
                          ResAll['RF']=RF
             
                      ####Compute mutation order error####
                      #Sequential_Error,Parallel_Error,Concurrent_Error: see ref [1] for th detail of thse error rates
                      Sequential_Error,Parallel_Error,Concurrent_Error = Functions.Compute_mut_order_error(InfMegTMP,TruMegTMP)
                      ResAll['AncError']=Sequential_Error
                      ResAll['SibError']=Parallel_Error
                      ResAll['CluError']=Concurrent_Error

                      ####Compute MLTED####
                      #MLTED: MLTED score computed by using the 'main' software.
                      MLTED=Functions.Compute_MLTED(TrueMultiTree,InfMegTMP,TruMegTMP,Rpath,MLTEDpath)
                      ResAll['MLTED']=MLTED

                      ####Compute TreeVec####
                      #TreeVec: TreeVec distance by using treespace package in R.
                      TreeVec=Functions.Compute_TreeVec(TrueMultiTree,InfMegTMP,Inf2TruClo,TrueCloneLs,Tru2Inf2GE,Rpath)
                      ResAll['TreeVec']=TreeVec

                      ####Store all scores####
                      ID2Res[Sim+'\t'+ID]=ResAll

                      ####Delete unnecessary files####
                      #os.remove('test.out')
                      #os.remove('test.nwk')
                      #os.remove('test1.nwk')
                      os.remove(InfMegTMP[:-4]+'_TrueCloneAnnotation.txt')
                      os.remove(InfMegTMP)
                      os.remove(TruMegTMP)

#####MLTED, TreeVec, RF, and error rate of ordering mutations (scores) reported in ref [1] in README.txt####
xl = pd.ExcelFile('Example/Tips_for_modification_of_clonephytester/Example_A-user_generated_data/Summary_example_A.xlsx')

####make output excel with score of "YournMethod"####
#Out: output file name
Out=Input[:-6]+'_summary.xlsx'
Functions.make_output_excel(Out,ID2Res,xl,Rpath)
os.remove('All.xlsx')
#make plot
#os.system(Rpath+' Plots.r')
#shutil.copy2('plot.jpg',Out[:-5]+'.jpg')
#os.remove('plot.jpg')
