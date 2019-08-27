import sys
import os
import shutil
import Functions
import glob
import pandas as pd
Input=sys.argv[1] #*.input
Sim2TrueFile={'G7':'G7\\True\\G7-ID.meg','G12':'G12\\True\\G12-ID.meg','P10':'P10\\True\\ID_True.meg','MA':'MA\\True\\ID_True_NoRedun.meg'}
xl = pd.ExcelFile('Summary.xlsx')#MA50xMAbest_summary.xlsx')#PhyloWGSSum3.xlsx')#
Out=Input[:-6]+'_summary.xlsx'
Input=open(Input,'r').readlines()
SeqLs=''
Rpath=''
MLTEDpath=''
Summary=''
for i in Input:
    i=i.strip().split('\t')
    if i[0]=='Data': SeqLs=i[1]
    elif i[0]=='R': Rpath='\"'+i[1]+'\"'
    elif i[0]=='MLTED': MLTEDpath=i[1]
    else: 
       print 'incorrect information ID',i[0]
       Summary+='incorrect information ID: '+i[0]+'\n'	   
Go='y'
if SeqLs=='':
    print 'Data file is not assigned. Data was not analyzed.'
    Summary+='Data file is not assigned. Data was not analyzed.'	
    Go='n'	
if Rpath=='':
    print 'R path is not assigned. Data was not analyzed.'
    Summary+='R path is not assigned. Data was not analyzed.'
    Go='n'
if MLTEDpath=='':
    print 'MLTED path is not assigned. Data was not analyzed.'
    Summary+='MLTED path is not assigned. Data was not analyzed.'	
    Go='n'
print Rpath
#open('A','r').readlines()	
Sheet2In={}
#Sheet2In['Clone_count']=['CloneC'] ##############
#Sheet2In['Ancestral_clone_count']=['AncNotFound'] #############
#Sheet2In['clone_genotype_error']=['GE'] ##################
Sheet2In['MLTED']=['MLTED']#,'MutPosiCountTrue','MutPosiCountInf']	
Sheet2In['RF']=['RF']
Sheet2In['TreeVec']=['TreeVec']
Sheet2In['Sequential_MutPair_count']=['AncError']#['AncIc','AncError']
Sheet2In['parallel_MutPair_count']=['SibError']#['SibIc','SibError']
Sheet2In['concurrent_MutPair_count']=['CluError']#['CluIc','CluError']
SheetLs=Sheet2In.keys()
def CleanFile(Target):	
    FiLs=glob.glob(Target)
    for Fi in FiLs:
        os.remove(Fi)	
if Go=='y':
    if os.path.exists(SeqLs)!=True:
        print 'Data file does not exists. ',SeqLs,' Data was not analyzed.'
        Summary+='Data file does not exists. '+SeqLs+' Data was not analyzed.'
    else:	
        ID2Res={}	
        SeqLs=open(SeqLs,'r').readlines()
        Head=SeqLs[0].strip()
        if	Head!='Sim\tData\tInferred Seq':
            print 'Data file should contain \"Sim\tData\tInferred Seq\" column. Data was not analyzed.'
            Summary+='Data file should contain \"Sim\tData\tInferred Seq\" column. Data was not analyzed.'
        else:
            SeqLs=SeqLs[1:]		
            for In0 in SeqLs:
               In=In0.strip().split('\t')
               if len(In)<3: 
                  print 'This line in the Data file is incorrect.',In0
                  Summary+='This line in the Data file is incorrect. '+In0+' This line was skipped.\n'
               else:				  
                  Sim=In[0]
                  ID=In[1]
                  TruMeg=Sim2TrueFile[Sim].replace('ID',ID)				  
                  InfMeg=In[2]	
		  
                  if os.path.exists(InfMeg)!=True:
                      print 'Inferred sequence file does not exists.',InfMeg
                      Summary+='Inferred sequence file does not exists. '+InfMeg+'\n'	
                  if os.path.exists(TruMeg)!=True:
                      print 'True sequence file does not exists.',TruMeg
                      Summary+='True sequence file does not exists. '+TruMeg+'\n'	
                  if os.path.exists(InfMeg)==True and os.path.exists(TruMeg)==True:
                      ResAll={}						  
                      print InfMeg,TruMeg
                      InfMegTMP=InfMeg.split('\\')[-1]
                      TruMegTMP=TruMeg.split('\\')[-1]					  
                      shutil.copy2(InfMeg, InfMegTMP)	
                      shutil.copy2(TruMeg, TruMegTMP)						  
                      os.system('python RmRedunSeq.py '+InfMegTMP)	 #*_NoRedun.meg	
                      os.remove(InfMegTMP)					  
                      InfMegTMP=InfMegTMP[:-4]+'_NoRedun.meg'					  
                      print 'compute GE'
                      PairList, TruCloCount, InfCloCount, AveGE, MedGE=Functions.PairClone(TruMegTMP,InfMegTMP)	
                      TmutNum=Functions.CountMutPosi(TruMegTMP)
                      ImutNum=Functions.CountMutPosi(InfMegTMP)	
                      ResAll['MutPosiCountTrue']=TmutNum
                      ResAll['MutPosiCountInf']=ImutNum					  
                      Functions.GetOut_from_list(PairList,'True clone\tInf clone\tGE\tAnnotation\n',InfMegTMP[:-4]+'_TrueCloneAnnotation.txt')  					  
                      print InfCloCount, AveGE	
                      ResAll['CloneC']=InfCloCount
                      ResAll['GE']=AveGE					  
                      print 'compute RF'
                      TrueFound_Ls, CorrectInf_Ls  = Functions.compute_E(PairList,99999999)			
          
                      RmNum, Compare2TreeIn, TrueMultiTree, InfMultiTree=Functions.MakeCompare2TreeIn(TrueFound_Ls,PairList,TruMegTMP,InfMegTMP,9999999,999999,'AA',Rpath)
                      print TrueMultiTree, InfMultiTree
                      if TrueMultiTree!='NA' and InfMultiTree!='NA':
                          RFin=Functions.DoRF_repeatPrune(TrueMultiTree, InfMultiTree,0)
                          RFin=RFin.split('\t')
                          print RFin					  
                          RF=1.0*(int(RFin[1])+int(RFin[2]))/(int(RFin[3])+int(RFin[4]))					  
                          print RF
                          ResAll['RF']=RF
                      CleanFile('Inf_ancestral_states_*.txt')
                      CleanFile('True_ancestral_states_*.txt')	
                      os.remove('test.nwk')	
                      os.remove('test1.nwk')
                      os.remove('True.meg')
                      os.remove('Inf.meg')					  
                      print 'Compute mutation order error'
                      os.system('python MutOrder.py '+InfMegTMP) #*_MutOrder.txt
                      os.system('python MutOrder.py '+TruMegTMP)	
                      TCla2Posi=Functions.Classify(TruMegTMP[:-4]+'_MutOrder.txt')
                      ICla2Posi=Functions.Classify(InfMegTMP[:-4]+'_MutOrder.txt')
         
                      AncDecNotFound, AncDecIncorrect = Functions.Count1(TCla2Posi, ICla2Posi, 'AncDec')	
                      DecAncNotFound, DecAncIncorrect = Functions.Count1(TCla2Posi, ICla2Posi, 'DecAnc')			  
                      SibNotFound, SibIncorrect = Functions.Count1(TCla2Posi, ICla2Posi, 'Sib')
                      CluNotFound, CluIncorrect = Functions.Count1(TCla2Posi, ICla2Posi, 'Clu')
                      AncTc=len(TCla2Posi['AncDec'])+len(TCla2Posi['DecAnc'])
                      AncIc=len(ICla2Posi['AncDec'])+len(ICla2Posi['DecAnc'])
                      AncFP=AncDecIncorrect+DecAncIncorrect
                      AncFN=AncDecNotFound+DecAncNotFound
                      ResAll['AncIc']=AncIc
                      if AncIc==0: 	ResAll['AncError']=	'Infinite'			  
                      else: ResAll['AncError']=((1.0*AncFP/AncIc)+(1.0*AncFN/AncTc))/2					  
                      SibTc=len(TCla2Posi['Sib'])
                      SibIc=len(ICla2Posi['Sib'])
                      SibFP=SibIncorrect
                      SibFN=SibNotFound	
                      ResAll['SibIc']=SibIc
                      if SibIc==0: ResAll['SibError']='Infinite'				  
                      else: ResAll['SibError']=((1.0*SibFP/SibIc)+(1.0*SibFN/SibTc))/2						  
                      CluTc=len(TCla2Posi['Clu'])
                      CluIc=len(ICla2Posi['Clu'])
                      CluFP=CluIncorrect
                      CluFN=CluNotFound						  
                      ResAll['CluIc']=CluIc
                      if CluIc==0: ResAll['CluError']='Infinite'					  
                      else: ResAll['CluError']=((1.0*CluFP/CluIc)+(1.0*CluFN/CluTc))/2	
                      print ResAll
                      os.remove(TruMegTMP[:-4]+'_MutOrder.txt')
                      os.remove(InfMegTMP[:-4]+'_MutOrder.txt')					  
                      print 'compute MLTED'
                      Itree=Functions.MakeMLTEin(InfMegTMP,'NA',InfMegTMP[:-4]+'_MLTEin.txt',Rpath)	
                      Ttree=Functions.MakeMLTEin(TruMegTMP,'NA',TruMegTMP[:-4]+'_MLTEin.txt',Rpath)	
                      if os.path.exists(InfMegTMP[:-4]+'_MLTEin.txt')==True and os.path.exists(TruMegTMP[:-4]+'_MLTEin.txt')==True: 
                             os.system(MLTEDpath+' '+TruMegTMP[:-4]+'_MLTEin.txt'+ ' ' +InfMegTMP[:-4]+'_MLTEin.txt'+' >'+'red.txt')
                             Res=open('red.txt','r').readlines()
                             for i in Res:
                                   
                                      if i.find('Normalized Similarity = ')!=-1: ResAll['MLTED']=1.0-float(i.split(' ')[-1])
                             os.remove('red.txt')	
                             os.remove(InfMegTMP[:-4]+'_MLTEin.txt')	
                             os.remove(TruMegTMP[:-4]+'_MLTEin.txt')								 
                      print ResAll	
                      print 'annotate ancestral clone'
                      Inf2TruClo,TrueCloneLs,Tru2Inf2GE=Functions.CloneAnno(InfMegTMP[:-4]+'_TrueCloneAnnotation.txt')	

                      MissTrueCloLs=Functions.FindMissingTruClo(InfMegTMP,Inf2TruClo,TrueCloneLs,Tru2Inf2GE)
                      AncCloTruLs=Functions.FindAnc(TruMegTMP)
   
                      NotFoundLs=[]
                      for Anc in AncCloTruLs:
                              if MissTrueCloLs.count(Anc)!=0: NotFoundLs.append(Anc)
				  
                      ResAll['AncNotFound']=len(NotFoundLs)
                      print 'Compute TreeVec'                      
                      InfTree, InfCloneLs=Functions.DupRenameClone_MPtree(InfMegTMP,Inf2TruClo,TrueCloneLs,Tru2Inf2GE,Rpath)	
                
                      print TrueCloneLs,InfCloneLs	
                      print Ttree,InfTree					  
      				  
                      Dist=Functions.Dotreespace(Ttree,InfTree,TrueCloneLs,InfCloneLs,Rpath)
                      ResAll['TreeVec']=Dist
                      print  ResAll 
                      os.remove('RunR.r')	
                      os.remove('RunTreeSpace.r')	
                      os.remove('test.out')	
                      os.remove('test.nwk')
                      os.remove('test1.nwk')
                      shutil.copy2(InfMegTMP[:-4]+'_TrueCloneAnnotation.txt',InfMeg[:-4]+'_TrueCloneAnnotation.txt')						  
                      os.remove(InfMegTMP[:-4]+'_TrueCloneAnnotation.txt')		
                      os.remove(InfMegTMP)	
                      os.remove(TruMegTMP)					  
                      ID2Res[Sim+'\t'+ID]=ResAll                      					  
                     
print 'make output excel'	
				  
writer = pd.ExcelWriter(Out, engine='xlsxwriter') 
for Sheet in Sheet2In:
    InLs=Sheet2In[Sheet]
    df = xl.parse(Sheet) 
    ColLs= list(df)
	
    Len=len(df['Sim'].tolist())
    AllDic={}
    for Col in ColLs:
        AllDic[Col]=df[Col].tolist()	
    InLen=len(InLs)
    InC=0
    while InC<InLen:
       In=InLs[InC]
       ColLs.append('YourMethod')#In)	   
       C=0
       ValLs=[]	   
       while C<Len:
          ID=df['Sim'].tolist()[C]+'\t'+df['Data'].tolist()[C]
          if ID2Res.has_key(ID)==True: 
                 if ID2Res[ID].has_key(In)==True: ValLs.append(ID2Res[ID][In])
                 else: 	ValLs.append('NA')			 
          else: ValLs.append('NA')
          C+=1
     #  AllDic[In]=ValLs
       AllDic['YourMethod']=ValLs	   
       InC+=1
    df = pd.DataFrame(AllDic)
    df = df[ColLs]	
    df.to_excel(writer, sheet_name=Sheet)
writer.save()		   
CleanFile('*_summary.txt')
CleanFile('*_NoRedun_ancestral_states.txt')
CleanFile('*_MLTEin.txt')

shutil.copy2(Out,'All.xlsx')
os.system(Rpath+' Plots.r')
os.remove('All.xlsx')
shutil.copy2('plot.jpg',Out[:-5]+'.jpg')
os.remove('plot.jpg')	  
	