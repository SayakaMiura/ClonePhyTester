import os
import shutil
import glob
from shutil import copy
import numpy
from ete3 import Tree
from Bio import Phylo
from itertools import combinations
from io import StringIO
from Bio.Phylo.Consensus import _BitString
import pandas as pd

###Example of adding your own score computation###
def Compute_clone_count_error(TruMegFile,InfMegFile):
    TcloLs, Tcloseq = ReadMegSeq(TruMegFile)
    IcloLs, Icloseq = ReadMegSeq(InfMegFile)
    CCerror=len(IcloLs)-len(TcloLs)
    return CCerror
##################################################

def parse_input(Input):
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
            print ('incorrect information ID',i[0])
            Summary+='incorrect information ID: '+i[0]+'\n'
    Go='y'
    if SeqLs=='':
        print ('Data file is not assigned. Data was not analyzed.')
        Summary+='Data file is not assigned. Data was not analyzed.'
        Go='n'
    if os.path.exists(SeqLs)!=True:
        print ('Data file does not exists. ',SeqLs,' Data was not analyzed.')
        Summary+='Data file does not exists. '+SeqLs+' Data was not analyzed.'
        Go='n'
    if os.path.exists(SeqLs)==True:
        SeqLs=open(SeqLs,'r').readlines()
        Head=SeqLs[0].strip()
        if    Head!='Sim\tData\tInferred Seq':
            print ('Data file should contain \"Sim\tData\tInferred Seq\" column. Data was not analyzed.')
            Summary+='Data file should contain \"Sim\tData\tInferred Seq\" column. Data was not analyzed.'
            Go='n'
    if Rpath=='':
        print ('R path is not assigned. Data was not analyzed.')
        Summary+='R path is not assigned. Data was not analyzed.'
        Go='n'
    if MLTEDpath=='':
        print ('MLTED path is not assigned. Data was not analyzed.')
        Summary+='MLTED path is not assigned. Data was not analyzed.'
        Go='n'
    return SeqLs[1:], Rpath, MLTEDpath, Summary, Go

def parse_inferred_file(In0,Summary,Sim2TrueFile):
    In=In0.strip().split('\t')
    Go='y'
    if len(In)<3:
        print ('This line in the Data file is incorrect.',In0)
        Summary+='This line in the Data file is incorrect. '+In0+' This line was skipped.\n'
        Go='n'
        Sim=''
        ID=''
        TruMeg=''
        InfMeg=''
    else:
        Sim=In[0]
        ID=In[1]
        TruMeg=Sim2TrueFile[Sim].replace('ID',ID)
        InfMeg=In[2]
        if os.path.exists(InfMeg)!=True:
                print ('Inferred sequence file does not exists.',InfMeg)
                Summary+='Inferred sequence file does not exists. '+InfMeg+'\n'
                Go='n'
        if os.path.exists(TruMeg)!=True:
                print ('True sequence file does not exists.',TruMeg)
                Summary+='True sequence file does not exists. '+TruMeg+'\n'
                Go='n'
    return Sim,ID,TruMeg,InfMeg,Summary,Go

def Clean(Target):
	filelist = glob.glob(Target)
	for f in filelist:
		os.remove(f)
		
def RmRedunSeq(Meg):
    Out2=Meg[:-4]+'_NoRedun.meg'
    NameOrder, Name2Seq=ReadMegSeq(Meg)
    out2='#MEGA\n!Title SNVs;\n!Format datatype=dna;\n\n'
    c=0
    RmSeq=[]
    Name2IdenLs={}
    IdenTar=[]
    SeqNum=len(NameOrder)
    Len=len(Name2Seq[NameOrder[0]])
    while c<SeqNum:
        Ref=NameOrder[c]
        RefSeq=Name2Seq[Ref]
        Name2IdenLs[Ref]=[Ref]
        Tc=0
        while Tc<SeqNum:
            Tar=NameOrder[Tc]
            TarSeq=Name2Seq[Tar]
            DifC=0
            Pc=0
            while Pc<Len:
                TarNuc=TarSeq[Pc]
                RefNuc=RefSeq[Pc]
                if TarNuc!=RefNuc: DifC+=1
                Pc+=1
            if DifC==0:
                    RmSeq.append(Tar)
                    Name2IdenLs[Ref].append(Tar)
            Tc+=1
        c+=1
    Done=[]
    for Name in Name2Seq:
        Code=Name in Done
        if Code!=True:
         IdenLs=Name2IdenLs[Name]
         Seq=Name2Seq[Name]
         MutC=0
         c=0
         while c<Len:
             if Seq[c]=='T': MutC+=1
             c+=1
         if MutC!=0:
           out2+=Name+'\n'+Name2Seq[Name]+'\n'
           Done+=IdenLs
    GetOut(Out2,out2)

def ReadMegSeq(Meg): 
  Meg=open(Meg,'r').readlines()
  Read='s'
  NameOrder=[]
  Name2Seq={}
  for i in Meg:
    if i[0]=='#' and i.strip()!='#MEGA' and i.strip()!='#mega' :
        Read='n'
        Name=i.strip()
        NameOrder.append(Name)
        Name2Seq[Name]=''
    elif Read=='n': Name2Seq[Name]+=i.strip()
  return NameOrder, Name2Seq

def UpMeg(Name2Seq,NameLs,Out):
    out='#MEGA\n!Title SNVs;\n!Format datatype=dna;\n\n'
    for Name in NameLs:
        if Name[0]!='#': Name='#'+Name
        out+=Name+'\n'+Name2Seq[Name]+'\n'
    GetOut(Out,out)

def all_parents(tree):
    parents = {}
    for clade in tree.find_clades(order='level'):
        for child in clade:
            parents[child] = clade
    return parents

def Nwk2NodeMap(OriTree):
    tree = Tree(OriTree)
    edge = 0
    for node in tree.traverse():
       if not node.is_leaf():
          node.name = "NODE_%d" %edge
          edge += 1
    InTree= tree.write(format=8)	
    tree1 = Phylo.read(StringIO(InTree), 'newick')	
    Child2Parent = all_parents(tree1)
    Anc2Dec={}
    TipLs=[]
    for Chi in Child2Parent:
        if (Child2Parent[Chi].name in Anc2Dec)!=True: Anc2Dec[Child2Parent[Chi].name]=[]
        if Chi.name!='hg19': Anc2Dec[Child2Parent[Chi].name].append(Chi.name)
        if Chi.name!=None and Chi.name.find('NODE_')==-1 and Chi.name!='hg19': TipLs.append(Chi.name)
    return  Anc2Dec,TipLs

def ListClade(Anc2Dec,TipLs):
    Clade2Clone={}
    for Anc in Anc2Dec:
        if Anc==None: Clade2Clone[None]=TipLs
        else:
         All=Anc2Dec[Anc]
         Clone=[]
         NodeLs=[]	 
         while All!=[]:	 
          for Item in All:
             if Item.find('NODE_')!=-1: NodeLs.append(Item)
             else: Clone.append(Item)
          All=[]
          for Node in NodeLs:
             All+=Anc2Dec[Node]	  
          NodeLs=[]
         Clade2Clone[Anc]=Clone
    return Clade2Clone

def MakeTopologeTree(Nwk,Rscript,intag):
    RunR='library(phytools)\nlibrary(phangorn)\nItree=read.newick(\''+Nwk+'\')\nItree1=di2multi(Itree)\nwrite.tree(Itree1,file=\''+Nwk[:-4]+'multi.nwk'+'\')\n'
    rfile = intag+"_runr.r"
    GetOut(rfile,RunR)
    os.system('\"'+Rscript+'\" '+rfile)
    os.remove(rfile)

def Meg2MP(Meg,RootTaxa,Rpath,intag):
    MegID=Meg[:-4]
    CloLs,Clo2Seq=ReadMegSeq(Meg)
    if CloLs.count(RootTaxa)==0: Clo2Seq[RootTaxa]='A'*(len(Clo2Seq[CloLs[0]]))
    out='#MEGA\n!Title SNVs;\n!Format datatype=dna;\n\n'
    for Clo in Clo2Seq:
        out+=Clo+'\n'+Clo2Seq[Clo]+'\n'
    GetOut(MegID+'_WithOut.meg',out)		
    os.system('megacc -a infer_MP_nucleotide.mao -d '+MegID+'_WithOut.meg'+' -o '+MegID+'.nwk')
    os.remove(MegID+'_summary.txt')
    os.remove(MegID+'_WithOut.meg')
    if os.path.exists(MegID+'.nwk')==True:
        MakeTopologeTree(MegID+'.nwk',Rpath,intag)
        os.remove(MegID+'.nwk')
        Clean(MegID+'_changes_list*.txt')
        Clean(MegID+'_ancestral_states*.txt')
    if os.path.exists(MegID+'multi.nwk')==True:  			 
       Nwkstr=open(MegID+'multi.nwk','r').readlines()[0]
       os.remove(MegID+'multi.nwk')
       Tree=RootTree(Nwkstr, RootTaxa[1:],intag)
    else: Tree='NA'
    return Tree
	
def RootTree(OriNwk, RootTaxa, intag):	
        treename = intag+"_test.nwk"
        treename1 = intag+"_test1.nwk"
        OutF=open(treename,'w')
        OutF.write(OriNwk)
        OutF.close()		
        trees = list(Phylo.parse(treename, 'newick'))
        for tree in trees:
           tree = tree.root_with_outgroup({'name': RootTaxa})
        Phylo.write(trees, treename1, "newick")	
        Tree_res=open(treename1,'r').readlines()[0]
        os.remove(treename)
        os.remove(treename1)
        return Tree_res

def Build_true_and_inferred_phylogeny(TrueFound_Ls,PairList,TrueMeg,InfMeg,Rpath,intag):
   InfCloLs, InfClo2Seq =ReadMegSeq(InfMeg)
   CorrectSeq_rename={'#Normal':'A'*len(InfClo2Seq[InfCloLs[0]])}
   for Pair in PairList:
            Pair=Pair.split('\t')
            if Pair[-1].strip()=='Initial' and TrueFound_Ls.count(Pair[0])!=0:			
             InfSeq=InfClo2Seq[Pair[1]]
             CorrectSeq_rename[Pair[0]]	=InfSeq
   TruCloLs,TruClo2Seq=ReadMegSeq(TrueMeg)
   TruClo2Seq['#Normal']='A'*len(TruClo2Seq[TruCloLs[0]])
   TrueFound_Ls.append('#Normal')
   truefile=intag+"_True.meg"
   inffile = intag+"_Inf.meg"
   UpMeg(TruClo2Seq,TrueFound_Ls,truefile)
   UpMeg(CorrectSeq_rename,TrueFound_Ls,inffile)
   TrueMultiTree=Meg2MP(truefile,'#Normal',Rpath,intag)
   InfMultiTree=Meg2MP(inffile,'#Normal',Rpath,intag)
   os.remove(truefile)
   os.remove(inffile)
   if TrueMultiTree=='' or InfMultiTree=='': return 'NA','NA'
   else: return TrueMultiTree, InfMultiTree

def PairClone(TreuMeg,InfMeg):
    TruCloLs, TruClo2Seq = ReadMegSeq(TreuMeg)
    InfCloLs, InfClo2Seq = ReadMegSeq(InfMeg)
    DoneInf=[]
    PairLs=[]	
    TruCloCount=0
    TruCloLs1=[]
    for TrueClo in TruCloLs:
        TruSeq=TruClo2Seq[TrueClo]	
        TruCloCount+=1		
        if TrueClo!='#Normal' and TrueClo!='#hg19':
            BestInf=''
            BestGE=999999999999
            for InfClo in InfCloLs:
                 if InfClo!='#Normal' and InfClo!='#hg19':
                       InfSeq=InfClo2Seq[InfClo]				 
                       GE=CountDifNum(TruSeq,InfSeq)
                       if GE<BestGE:
                            BestGE=GE
                            BestInf=InfClo
            DoneInf.append(BestInf)
            TruCloLs1.append(TrueClo)
            PairLs.append(TrueClo+'\t'+BestInf+'\t'+str(BestGE)+'\tInitial')
    for InfClo in InfCloLs:
        if 	InfClo!='#Normal' and InfClo!='#hg19' and DoneInf.count(InfClo)==0:
            InfSeq=InfClo2Seq[InfClo]
            BestTru=''
            BestGE=999999999999	
            for TrueClo in TruCloLs:			
                if TrueClo!='#Normal' and TrueClo!='#hg19':
                      TruSeq=TruClo2Seq[TrueClo]
                      GE=CountDifNum(TruSeq,InfSeq)
                      if GE<BestGE:
                            BestGE=GE
                            BestTru=TrueClo
            PairLs.append(BestTru+'\t'+InfClo+'\t'+str(BestGE)+'\tExtraInf')
    return 	PairLs, TruCloCount, TruCloLs1

def CountDifNum(Seq0,Seq1):
            Len=len(Seq0)		
            Dif=0
            c=0
            while c<Len:
                if Seq0[c]!=Seq1[c]: Dif+=1
                c+=1
            return Dif
					  
def Classify(File):
    In=open(File,'r').readlines()
    Cla2Posi={'DecAnc':[],'AncDec':[],'Sib':[],'Clu':[]}	
    for i in In:
       i=i.strip().split('\t')
       if (i[1] in Cla2Posi)!=True: Cla2Posi[i[1]]=[]
       Cla2Posi[i[1]].append(i[0])	   
    return Cla2Posi

def Count1(Tdic, Idic, Cate):
      if (Cate in Tdic)!=True: Tposi=[]
      else: Tposi=Tdic[Cate]
      if (Cate in Idic)!=True: Iposi=[]
      else: Iposi=Idic[Cate]
      NotFouncd=0
      for Posi in Tposi:
         if Iposi.count(Posi)==0: NotFouncd+=1
      Incorrect=0
      for Posi in Iposi:
         if Tposi.count(Posi)==0: Incorrect+=1	 
      return NotFouncd, Incorrect

def GetCommon(CLs,C2S):
    All=list(C2S.keys())
    Len=len(C2S[All[0]])
    c=0
    In=''	
    while c<Len:
        Good='y'
        for Clo in CLs:
             if C2S['#'+Clo][c]!='T' : 
               Good='n'
        for Clo in All:
             if C2S[Clo][c]=='T' and CLs.count(Clo[1:])==0: Good='n'        			 
        if Good=='y':
            In+=str(c)+ ','
        c+=1
    if In!='': In=In[:-1]
    return In

def MakeMLTEin(Meg,OriTree,OutMLTEDfile,Rpath,intag):
  if OriTree=='NA':
       OriTree=Meg2MP(Meg,'#hg19',Rpath,intag)
  if OriTree.strip()=='NA': return  OriTree
  else:  
    Anc2Dec,TipLs = Nwk2NodeMap(OriTree)
    Clade2Clone = ListClade(Anc2Dec,TipLs)	
    CloLs,Clo2Seq=ReadMegSeq(Meg)
    Len=len(Clo2Seq[list(Clo2Seq.keys())[0]])
    Anc2DecIn=''
    CladeIn={}
    AssignMut=[]
    if len(Anc2Dec[None])==1: 
       NoRoot='y'
       RootClade=Anc2Dec[None][0]
    else: 
       NoRoot=''	
       RootClade='Root'
    RmTip=['hg19']	
    for Tip in TipLs:
       MutIn=GetCommon([Tip],Clo2Seq)
       if MutIn=='': RmTip.append(Tip)	   
       else: 
           CladeIn[Tip]=Tip+'='+MutIn+'\n'	
           AssignMut+=MutIn.split(',')		   
    for Anc in Anc2Dec:
       DecLs=Anc2Dec[Anc]
       if Anc==None : Ain=RootClade+':'	   
       else:Ain=Anc+':'
       Added='n'	   
       for Dec in DecLs:
          if RmTip.count(Dec)==0:	   
             Ain+=Dec+','
             Added='y'			 
       if Added=='y':Ain=Ain[:-1]
       else: pass
       if Anc==None and NoRoot=='y': pass	   
       elif Added=='n': pass	   
       else: Anc2DecIn+=Ain+'\n'
       CloneLs=Clade2Clone[Anc]
       MutIn=GetCommon(CloneLs,Clo2Seq)
       if Anc==None and NoRoot!='y': pass	   
       elif Anc==None: 
           CladeIn[RootClade]=RootClade+'='+MutIn+'\n'	
           AssignMut+=MutIn.split(',')			  
       else: 
           CladeIn[Anc]=Anc+'='+MutIn+'\n'
           AssignMut+=MutIn.split(',')	
    UnAssign=''
    c=0
    while c<Len:
        if AssignMut.count(str(c))==0: UnAssign+=str(c)+',' 	
        c+=1
    if UnAssign!='': UnAssign=UnAssign[:-1]
    CladeIn[RootClade]=	CladeIn[RootClade][:-1]+','+UnAssign+'\n'
    CladeInStr=''
    for Clade in CladeIn:
        CladeInStr+=CladeIn[Clade]	
    GetOut(OutMLTEDfile,CladeInStr+Anc2DecIn)
    return  OriTree

def CloneAnno(TruMegTMP,InfMegTMP):
            CloAnno=InfMegTMP[:-4]+'_TrueCloneAnnotation.txt'
            PairList, TruCloCount, TrueClone_Ls =PairClone(TruMegTMP,InfMegTMP)
            GetOut_from_list(PairList,'True clone\tInf clone\tGE\tAnnotation\n',CloAnno)
            InfClo2TrueClo2GE={}
            TrueClo2InfClo2GE={}			
            CloAnno=open(CloAnno,'r').readlines()[1:]
            AllTru=[]			
            for Line in CloAnno:
              Line=Line.strip().split('\t')
              Inf=Line[1][1:]			  
              Tru=Line[0][1:]
              AllTru.append(Tru)			  
              GE=float(Line[2])			  
              if (Inf in InfClo2TrueClo2GE)!=True: InfClo2TrueClo2GE[Inf]={}
              InfClo2TrueClo2GE[Inf][Tru]=GE
              if (Tru in TrueClo2InfClo2GE)!=True: TrueClo2InfClo2GE[Tru]={}
              TrueClo2InfClo2GE[Tru][Inf]=GE
            InfClo2TrueClo={}
            for InfClo in InfClo2TrueClo2GE:
                 TrueClo2GE=InfClo2TrueClo2GE[InfClo]
                 TruCloLs=list(TrueClo2GE.keys())
                 if len(TruCloLs)==1: InfClo2TrueClo[InfClo]=TruCloLs[0]
                 else:
                    Best=99999999999999
                    Bclo=''
                    for Tclo in TrueClo2GE:
                         GE=TrueClo2GE[Tclo]
                         if GE<Best: 
                               Best=GE
                               Bclo=Tclo
                    InfClo2TrueClo[InfClo]=Bclo
            return PairList,InfClo2TrueClo,AllTru,TrueClo2InfClo2GE,TrueClone_Ls

def DupRenameClone_MPtree(InfFile,InfClo2TrueClo,AllTru,TrueClo2InfClo2GE,Rpath,intag):
            CloLs, Clo2Seq=ReadMegSeq(InfFile)
            AllTru=list(set(AllTru))
            Tru2Count={}	
            InfCloLs=[]
            outMeg='#MEGA\n!Title SNVs;\n!Format datatype=dna;\n\n'
            for Clo in Clo2Seq:
               if Clo!='#hg19':			
                  TruClo=InfClo2TrueClo[Clo[1:]]
                  if (TruClo in Tru2Count)!=True: Tru2Count[TruClo]=1
                  else:
                      Tru2Count[TruClo]+=1					  
                      TruClo+='Dup'+str(Tru2Count[TruClo]-1)
                  outMeg+='#'+TruClo+'\n'+	Clo2Seq[Clo]+'\n'
                  InfCloLs.append(TruClo)
            MissAdd=''
            for Tc in AllTru:
                if (Tc in Tru2Count)!=True:
                      InfDic=TrueClo2InfClo2GE[Tc]
                      InfLs=list(InfDic.keys())
                      Bclo=''					  
                      if len(InfLs)==1: Bclo=InfLs[0]
                      else:
                         Best=99999999999999
                         for Iclo in InfDic:
                              GE=InfDic[Iclo]
                              if GE<Best: 
                                    Best=GE
                                    Bclo=Tclo					  
                      outMeg+='#'+Tc+'\n'+Clo2Seq['#'+Bclo]+'\n'			
                      InfCloLs.append(Tc)
            GetOut(InfFile[:-4]+'_tmp.meg',outMeg)
            InfTree=Meg2MP(InfFile[:-4]+'_tmp.meg','#hg19',Rpath,intag)
            os.remove(InfFile[:-4]+'_tmp.meg')		
            return InfTree.strip(),InfCloLs

def PruneTree(KeepCloSet, Nwk_str):
    from ete3 import Tree
    if Nwk_str=='NA': return 'NA'
    else:
      t=Tree(Nwk_str)
      t.prune(KeepCloSet)
      t1=t.write(format=1)
      return t1

def DoRF_fromtree(TrueTree,InfTree,nametag):
     PhyloNet_dir=''
     PhyloNet='PhyloNet_3.6.1.jar'
     Current_dir=os.getcwd()
     RF=[]	
     Res=''	 
     PhyloNETin='#NEXUS\n\n'
     PhyloNETin+='BEGIN NETWORKS;\n\n' 
     PhyloNETin+='Network tree1 = '+TrueTree
     PhyloNETin+='Network tree2 = '+InfTree
     PhyloNETin+='END;\n\n'
     PhyloNETin+='BEGIN PHYLONET;\n\n'
     PhyloNETin+='SymmetricDifference tree1 tree2;\n\n'
     PhyloNETin+='END;\n'
     phynex = Current_dir+"/"+nametag+"_PhyloNetIn.nex"
     rfout = Current_dir+'/'+nametag+"_RF.out"
     GetOut(phynex,PhyloNETin)
     os.system('java -jar '+PhyloNet +' '+ phynex+'>'+rfout)
     os.chdir(Current_dir)
     if os.path.exists(rfout)==True:
       RF=open(rfout,'r').readlines()[2:]
     for i in RF:
         Res+='\t'+i.split(':')[-1].strip()
     import time		 
     time.sleep(1)
     if os.path.exists(rfout)==True: os.remove(rfout)
     if os.path.exists(phynex)==True: os.remove(phynex)
     return Res

def Compute_RF(TrueTree,InfTree,nametag):
    if TrueTree!='NA':
     Ttree=Phylo.read(StringIO(TrueTree), "newick")
     TrueFound_Ls0=Ttree.get_terminals()
    # TrueFound_Ls=[]
     PrCloSet=[]
     for i in TrueFound_Ls0:
          if i.name!='Normal':
            PrCloSet.append(i.name.replace('#',''))
    if TrueTree=='NA' or InfTree=='NA': return '\tno tree'
    else:
        BestOut=''
        TruPr=PruneTree(PrCloSet, TrueTree)
        InfPr=PruneTree(PrCloSet, InfTree)
        Res=DoRF_fromtree(TruPr,InfPr,nametag)
        if Res!='':
                    BestOut='\t'+Res.split('\t')[1]+'\t'+Res.split('\t')[2]+'\t'+Res.split('\t')[3]+'\t'+Res.split('\t')[4]
        if BestOut=='': return 'NA'
        else:
         RFin=BestOut.split('\t')
         RF=1.0*(int(RFin[1])+int(RFin[2]))/(int(RFin[3])+int(RFin[4]))
         return RF

def Compute_mut_order_error(InfMegTMP,TruMegTMP):
    MutOrder(InfMegTMP)
    MutOrder(TruMegTMP)

    TCla2Posi=Classify(TruMegTMP[:-4]+'_MutOrder.txt')
    ICla2Posi=Classify(InfMegTMP[:-4]+'_MutOrder.txt')
    AncDecNotFound, AncDecIncorrect = Count1(TCla2Posi, ICla2Posi, 'AncDec')
    DecAncNotFound, DecAncIncorrect = Count1(TCla2Posi, ICla2Posi, 'DecAnc')
    SibNotFound, SibIncorrect = Count1(TCla2Posi, ICla2Posi, 'Sib')
    CluNotFound, CluIncorrect = Count1(TCla2Posi, ICla2Posi, 'Clu')
    AncTc=len(TCla2Posi['AncDec'])+len(TCla2Posi['DecAnc'])
    AncIc=len(ICla2Posi['AncDec'])+len(ICla2Posi['DecAnc'])
    AncFP=AncDecIncorrect+DecAncIncorrect
    AncFN=AncDecNotFound+DecAncNotFound
    if AncIc==0: AncError=    'Infinite'
    else: AncError=((1.0*AncFP/AncIc)+(1.0*AncFN/AncTc))/2
    SibTc=len(TCla2Posi['Sib'])
    SibIc=len(ICla2Posi['Sib'])
    SibFP=SibIncorrect
    SibFN=SibNotFound
    if SibIc==0: SibError='Infinite'
    else: SibError=((1.0*SibFP/SibIc)+(1.0*SibFN/SibTc))/2
    CluTc=len(TCla2Posi['Clu'])
    CluIc=len(ICla2Posi['Clu'])
    CluFP=CluIncorrect
    CluFN=CluNotFound
    if CluIc==0: CluError='Infinite'
    else: CluError=((1.0*CluFP/CluIc)+(1.0*CluFN/CluTc))/2
    os.remove(TruMegTMP[:-4]+'_MutOrder.txt')
    os.remove(InfMegTMP[:-4]+'_MutOrder.txt')
    return AncError,SibError,CluError

def Compute_MLTED(Ttree,InfMegTMP,TruMegTMP,Rpath,MLTEDpath,intag):
    if Ttree=='NA': return 'NA'
    else:
        InfMegFn = InfMegTMP[:-4]+"_MLTEin.txt"
        TruMegFn = TruMegTMP[:-4]+"_MLTEin.txt"
        redFn = InfMegTMP[:-4]+"_red.txt"
        Itree=MakeMLTEin(InfMegTMP,'NA',InfMegFn,Rpath,intag)
        Ttree=Ttree.replace('Normal:','hg19:')
        Ttree_redun=MakeMLTEin(TruMegTMP,Ttree,TruMegFn,Rpath,intag)
        MLTED='NA'
        if Itree=='NA': return 'NA'
        if os.path.exists(InfMegFn)==True and os.path.exists(TruMegFn)==True:
            os.system(MLTEDpath+' '+TruMegFn+ ' ' +InfMegFn+' >'+redFn)
            Res=open(redFn,'r').readlines()
            for i in Res:
                if i.find('Normalized Similarity = ')!=-1: MLTED=1.0-float(i.split(' ')[-1])
            os.remove(redFn)
            os.remove(InfMegFn)
            os.remove(TruMegFn)
            return MLTED

def Compute_TreeVec(Ttree,InfMegTMP,Inf2TruClo,TrueCloneLs,Tru2Inf2GE,Rpath,intag):
    if Ttree=='NA': return 'NA'
    else:
      Ttree=Ttree.replace('Normal:','hg19:')
      InfTree, InfCloneLs=DupRenameClone_MPtree(InfMegTMP,Inf2TruClo,TrueCloneLs,Tru2Inf2GE,Rpath,intag)
      if InfTree=='NA': return 'NA'
      else:
        TreeVec=Dotreespace(Ttree,InfTree,TrueCloneLs,InfCloneLs,Rpath,intag)
        # os.remove('RunTreeSpace.r')
        return TreeVec

def Dotreespace(TTre,ITre,TCloLs,ICloLs,Rpath,intag):
    CatLs=''
    TipLs=''
    for Clo in TCloLs:
       CatLs+='\''+Clo+'\','
       TipLs+='\''+Clo+'\','
    for Clo in ICloLs:
       if Clo.find('Dup')!=-1:
           Cat=Clo.split('Dup')[0]
           CatLs+='\''+Cat+'\','
           TipLs+='\''+Clo+'\','           		   
    CatLs=CatLs[:-1]
    TipLs=TipLs[:-1]
    CurDir=os.getcwd()	
    testout = intag+"_test.out"
    treespaceIn='library(treespace)\n\n'
    treespaceIn+='Tree1 <- read.tree(text=\"'+TTre+'\")\n'
    treespaceIn+='Tree2 <- read.tree(text=\"'+ITre.strip()+'\")\n'
    treespaceIn+='trees <- list(Tree1,Tree2)\n'
    treespaceIn+='df <- data.frame (category = c('+CatLs+'),'+'tiplabels = c('+TipLs+'))\n'
    treespaceIn+='dists <- relatedTreeDist(trees,df)\n'
    #treespaceIn+='write.table(dists[1],\''+CurDir+'/test.out'+'\')\n'
    treespaceIn+='write.table(dists[1],\''+testout+'\')\n'
    runtreespacef = intag+"_RunTreeSpace.r"
    GetOut(runtreespacef,treespaceIn)
    os.system(Rpath+' '+runtreespacef)
    if os.path.exists(testout)==True:
          Dist=open(testout,'r').readlines()[1].strip().split(' ')[-1]
          os.remove(testout)

    else: Dist='NA'
    os.remove(runtreespacef)
    return Dist	

def MutOrder(Meg):
    Out=Meg[:-4]+'_MutOrder.txt'
    CloLs, Clo2Seq=ReadMegSeq(Meg)
    Len1=len(Clo2Seq[CloLs[0]])
    c0=0
    Len0=Len1-1
    out=''
    while c0<Len0:
        c1=c0+1
        while c1<Len1:
            Ls1Uni=[]
            Ls0Uni=[]
            Com=[]
            ComUn=[]
            for Clo in Clo2Seq:
              if Clo!='#hg19':
                Nuc0=Clo2Seq[Clo][c0]
                Nuc1=Clo2Seq[Clo][c1]
                if Nuc0=='T' and Nuc1!='T': Ls0Uni.append(Clo)
                elif Nuc0!='T' and Nuc1=='T': Ls1Uni.append(Clo)
                elif Nuc0=='T' and Nuc1=='T': Com.append(Clo)
                elif Nuc0!='T' and Nuc1!='T': ComUn.append(Clo)
                else:
                    print ('??',c0,c1,Clo)
                    open('S','r').readlines()
            if Ls1Uni==[] and Ls0Uni==[] and Com!=[] : out+=str(c0)+'-'+str(c1)+'\tClu\n'
            elif Ls1Uni!=[] and Ls0Uni==[] and Com!=[] : out+=str(c0)+'-'+str(c1)+'\tDecAnc\n'
            elif Ls1Uni==[] and Ls0Uni!=[] and Com!=[] : out+=str(c0)+'-'+str(c1)+'\tAncDec\n'
            elif Ls1Uni!=[] and Ls0Uni!=[] : out+=str(c0)+'-'+str(c1)+'\tSib\n'
            elif Ls1Uni==[] and Ls0Uni==[] and Com==[] and ComUn!=[]: out+=str(c0)+'-'+str(c1)+'\tUnassign\n'
            elif (Ls1Uni!=[] or Ls0Uni!=[]) and Com==[] : out+=str(c0)+'-'+str(c1)+'\tUnassign\n'
            else:
                out+=str(c0)+'-'+str(c1)+'\tUnassign\n'
                print ('???',c0,c1,Ls1Uni, Ls0Uni, Com, ComUn)
                open('S','r').readlines()
            c1+=1
        c0+=1
    GetOut(Out,out)

def GetOut(File,In):
    OutF=open(File,'w')
    OutF.write(In)
    OutF.close()

def GetOut_from_list(List,Head,OutF):
     out=Head
     for i in List:
        out+=i+'\n'
     GetOut(OutF,out)

def make_output_excel(Out,ID2Res,xl,Rpath):
 SheetLs=xl.sheet_names
 #print (SheetLs)
 writer = pd.ExcelWriter(Out, engine='xlsxwriter')
 Sheet2In={}
 Sheet2In['MLTED']=['MLTED']
 Sheet2In['RF']=['RF']
 Sheet2In['TreeVec']=['TreeVec']
 Sheet2In['Sequential_MutPair_count']=['AncError']
 Sheet2In['parallel_MutPair_count']=['SibError']
 Sheet2In['concurrent_MutPair_count']=['CluError']
 for Sheet in SheetLs:
     Pre=Sheet in Sheet2In
     if Pre!=True:
        Sheet2In[Sheet]=[Sheet]
 for Sheet in Sheet2In:
  if SheetLs.count(Sheet)!=0:
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
                if (ID in ID2Res)==True:
                    if (In in ID2Res[ID])==True: ValLs.append(ID2Res[ID][In])
                    else:     ValLs.append('NA')
                else: ValLs.append('NA')
                C+=1
        AllDic['YourMethod']=ValLs
        InC+=1
    df = pd.DataFrame(AllDic)
    df = df[ColLs]
    df.to_excel(writer, sheet_name=Sheet)
 writer.save()
 Clean('*_summary.txt')
 Clean('*_NoRedun_ancestral_states.txt')
 Clean('*_MLTEin.txt')
 Clean('*.nwk')
 Clean('*.meg')
 Clean('Inf_changes_list_*.txt')
 Clean('True_changes_list_*.txt')
 shutil.copy2(Out,'All.xlsx')

