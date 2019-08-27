import os
import shutil
import glob
from shutil import copy
import numpy
from ete2 import Tree
from Bio import Phylo
from itertools import combinations
from StringIO import StringIO
from Bio.Phylo.Consensus import _BitString

def Clean(Target):
	filelist = glob.glob(Target)
	for f in filelist:
		os.remove(f)
		

def FindRoot(Dec2Anc,OutGroup):
     Root=''
     Dec=Dec2Anc.keys()[0]
     while Dec2Anc.has_key(Dec)==True:	 
         Root=Dec2Anc[Dec]
         Dec=Root		 
     return Root			  
     	 

def CountIden(CDic):
    C=0
    for i in CDic:
         if CDic[i]=='y': C+=1
    return C		 
	

def ReadMegSeq(Meg): 
  Meg=open(Meg,'r').readlines()
  Read='s'
  out2=''
  NameOrder=[]
  Name2Seq={}
  for i in Meg:
    if i[0]=='#' and i.strip()!='#MEGA' and i.strip()!='#mega' :
        Read='n'
        Name=i.strip()
        NameOrder.append(Name)
        Name2Seq[Name]=''
    elif Read=='n': Name2Seq[Name]+=i.strip()
    elif Read=='s': out2+=i
  return NameOrder, Name2Seq, out2
  
def GetMutPos(Seq):
    TMP=[]
    Len=len(Seq)
    c=0
    while c<Len:
        if Seq[c]=='T': TMP.append(c)
        c+=1
    return TMP 

def UpMeg(Name2Seq,NameLs,AA,Out):
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
def Nwk2NodeMap(OriTree): #Nwk
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
    Dec2Anc={}
    TipLs=[]
    for Chi in Child2Parent:
       # print Chi.name, Child2Parent[Chi].name
        if Anc2Dec.has_key(Child2Parent[Chi].name)!=True: Anc2Dec[Child2Parent[Chi].name]=[]
        if Chi.name!='hg19': Anc2Dec[Child2Parent[Chi].name].append(Chi.name)	
        Dec2Anc[Chi.name]=Child2Parent[Chi].name
        if Chi.name!=None and Chi.name.find('NODE_')==-1 and Chi.name!='hg19': TipLs.append(Chi.name)	
    	
    	
    return  Anc2Dec,Dec2Anc,TipLs
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
          print All	  
         Clade2Clone[Anc]=Clone
    
    	 
    return Clade2Clone			
def Meg2MP(Meg,RootTaxa,Rpath):
    MegID=Meg[:-4]
    CloLs,Clo2Seq,out=ReadMegSeq(Meg)
    if CloLs.count('#hg19')==0: Clo2Seq['#hg19']='A'*(len(Clo2Seq[CloLs[0]]))
    for Clo in Clo2Seq:
        out+=Clo+'\n'+Clo2Seq[Clo]+'\n'
    GetOut(MegID+'_WithOut.meg',out)		
    os.system('megacc -a infer_MP_nucleotide.mao -d '+MegID+'_WithOut.meg'+' -o '+MegID+'.nwk') 
  
    if os.path.exists(MegID+'.nwk')==True:      
        os.system('python MakeTopologeTree.py ' +MegID+'.nwk '+Rpath)
    if os.path.exists(MegID+'multi.nwk')==True:  			 
       Nwkstr=open(MegID+'multi.nwk','r').readlines()[0]
              
       Tree=RootTree(Nwkstr, RootTaxa[1:])
       os.remove(MegID+'.nwk')	 
       os.remove(MegID+'multi.nwk')	
       os.remove(MegID+'_summary.txt')	
       AncFile=glob.glob(MegID+'_ancestral_states_*.txt')
       for File in AncFile:
           os.remove(File) 		   
    else: Tree='NA\n'	   
    os.remove(MegID+'_WithOut.meg')	

	   
    return Tree 		
def do_MP(InMeg,Rpath):
    Tree=''
    if os.path.exists(InMeg[:-4]+'.nwk')!=True:
      CloLs, Clo2Seq, AA=ReadMegSeq(InMeg)
      Clo2Seq['#Normal']='A'*len(Clo2Seq[CloLs[0]]) 
      UpMeg(Clo2Seq,Clo2Seq.keys(),'AA',InMeg[:-4]+'_withNarmal.meg')
      os.system('megacc -a infer_MP_nucleotide.mao -d '+InMeg[:-4]+'_withNarmal.meg'+' -o '+InMeg[:-4]+'.nwk')

    if os.path.exists(InMeg[:-4]+'.nwk')==True and os.path.exists(InMeg[:-4]+'multi.nwk')!=True:# or InfRenameMeg.find('CloneFinder')!=-1:		
             os.system('python MakeTopologeTree.py ' +InMeg[:-4]+'.nwk'+' '+Rpath)
    if os.path.exists(InMeg[:-4]+'multi.nwk')==True:			 
             Tree_str=open(InMeg[:-4]+'multi.nwk','r').readlines()[0]
             Tree=RootTree(Tree_str, 'Normal')
    if os.path.exists(InMeg[:-4]+'.nwk')==True: os.remove(InMeg[:-4]+'.nwk')
    if os.path.exists(InMeg[:-4]+'_withNarmal.meg')==True:  os.remove(InMeg[:-4]+'_withNarmal.meg')	
    if os.path.exists(InMeg[:-4]+'multi.nwk')==True:  os.remove(InMeg[:-4]+'multi.nwk')	
    return Tree
	
def RootTree(OriNwk, RootTaxa):	
        OutF=open('test.nwk','w')
        OutF.write(OriNwk)
        OutF.close()		
        trees = list(Phylo.parse('test.nwk', 'newick'))
        for tree in trees:
           tree = tree.root_with_outgroup({'name': RootTaxa})
        Phylo.write(trees, 'test1.nwk', "newick")	

        return open('test1.nwk','r').readlines()[0]	


def FindBestTrueClone(PairList):
    Inf2Best={}
    Inf2TruLs=ListTrueClone(PairList)
    for Inf in Inf2TruLs:	
     BestGE=9999999999
     BestClo=''	
     for i in PairList:
        i=i.split('\t')
        Inf0=i[1]
        Tru=i[0]
        GE=float(i[2])		
        if Inf==Inf0 and GE<BestGE: 
             BestGE=GE
             BestClo=Tru			  
     Inf2Best[Inf]=[BestClo]
    return Inf2Best	
	
def GE_Prune_AllInf(TrueFound_Ls,PairList,TrueMeg,InfMeg,GEcut):
         AveGE_AllInf_PruneTrue='NA'	
         Inf2TruLs={}  
         BestGEls=[]		 
         for Pair in PairList:
            Pair=Pair.split('\t')
            if TrueFound_Ls.count(Pair[0])!=0:			

               Inf=	Pair[1]
               Tru=	Pair[0]
               GE=int(Pair[2])
               if Inf2TruLs.has_key(Inf)!=True: Inf2TruLs[Inf]=[GE,Tru]
               if Inf2TruLs[Inf][0]>GE:
                    Inf2TruLs[Inf]=[GE,Tru]

         InfCloLs, InfClo2Seq, AA=ReadMegSeq(InfMeg)
	 
         for InfB in Inf2TruLs:

             BestGEls.append(Inf2TruLs[InfB][0])			 
	 
         AveGE_AllInf_PruneTrue=1.0*sum(BestGEls)/len(BestGEls)/len(InfClo2Seq[InfCloLs[0]])
         AveGE_AllInf_PruneTrue=str(AveGE_AllInf_PruneTrue)
 
         return AveGE_AllInf_PruneTrue		
def MakeCompare2TreeIn(TrueFound_Ls,PairList,TrueMeg,InfMeg,CutGE,pruneNum,ID,Rpath):
   print TrueMeg,InfMeg

   InfCloLs, InfClo2Seq, AA=ReadMegSeq(InfMeg)
   CorrectSeq_rename={'#Normal':'A'*len(InfClo2Seq[InfCloLs[0]])}
   InfInc=[]
   for Pair in PairList:
            Pair=Pair.split('\t')
            if Pair[-1].strip()=='Initial' and TrueFound_Ls.count(Pair[0])!=0:			
             InfSeq=InfClo2Seq[Pair[1]]
             CorrectSeq_rename[Pair[0]]	=InfSeq
             InfInc.append(Pair[1])
   InfInc=list(set(InfInc))
   InfRm=len(InfCloLs)-len(InfInc)
   print 'num inf removed',InfRm	

   Go='y'	 
   if Go=='y':#else:	

    TruCloLs,TruClo2Seq,AA=ReadMegSeq(TrueMeg)
    TruClo2Seq['#Normal']=	'A'*len(TruClo2Seq[TruCloLs[0]])
    TrueFound_Ls.append('#Normal')		 
      
    PruneC=0
    UpMeg(TruClo2Seq,TrueFound_Ls,'AA','True.meg')
    UpMeg(CorrectSeq_rename,TrueFound_Ls,'AA','Inf.meg')
    TrueMultiTree=do_MP('True.meg',Rpath) 
    InfMultiTree=do_MP('Inf.meg',Rpath)
    if TrueMultiTree=='' or InfMultiTree=='': return str(InfRm),'NA','NA','NA'
    else: 
      In='java -classpath PhyloCore.jar:parallelcolt-0.9.4.jar treecomparison.Compare2Trees '
      In+='-s \"'+TrueMultiTree.strip()+'\" '
      In+='\"'+InfMultiTree.strip()+'\" >Clone2CompareOut.txt'
 
      return str(InfRm),In, TrueMultiTree, InfMultiTree	  


def GetSharePosi1(Ori2Seq,ShareNuc):
  for i in Ori2Seq:
    Name=i
  SharePosi=[]
  Len=len(Ori2Seq[Name])
  c=0
  while c<Len:
    AllMut='y'
    for Ori in Ori2Seq:
      if Ori!='#hg19' and Ori!='#Normal':
       Nuc=Ori2Seq[Ori][c]
       if Nuc!=ShareNuc: AllMut='n'
    if AllMut=='y': 
        SharePosi.append(c)
    c+=1
  return SharePosi
def CountMutPosi(MegFile):
    Cls,Cseq,AA=ReadMegSeq(MegFile)
    Cseq1={}
    for Seq in Cseq:
       Cseq1[Seq]=Cseq[Seq].replace('?','A')	
    AposiLs=GetSharePosi1(Cseq1,'A')
    MutNum=len(Cseq[Seq])-len(AposiLs)
    return MutNum	
def PairClone(TreuMeg,InfMeg): #[TrueClone\tInfClone\tGE, ...]	
    TruCloLs, TruClo2Seq, AA = ReadMegSeq(TreuMeg)
    InfCloLs, InfClo2Seq, AA = ReadMegSeq(InfMeg)
    DoneInf=[]
    PairLs=[]	
    TruCloCount=0
    InfCloCount=0
    BestGELs=[]	
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
            PairLs.append(TrueClo+'\t'+BestInf+'\t'+str(BestGE)+'\tInitial')
            BestGELs.append(BestGE)			
    for InfClo in InfCloLs:
        if 	InfClo!='#Normal' and InfClo!='#hg19' : InfCloCount+=1	
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
            BestGELs.append(BestGE)	
    AveGE=1.0*sum(BestGELs)/len(BestGELs)
    PerAveGE=AveGE/len(TruSeq)	
    MedGE=numpy.median(BestGELs)/len(TruSeq)	

    return 	PairLs, TruCloCount, InfCloCount, PerAveGE,MedGE#,PruneGE 
def compute_E(PairList,CutGE):
    TrueDone=[]
    InfDone=[]	
    CorrectInfC=0
    TrueFoundC=0
    for Pair in PairList:
        Pair=Pair.split('\t')
        Tru=Pair[0]
        Inf=Pair[1]
        GE=int(Pair[2])
        if GE<=CutGE:
            TrueDone.append(Tru)
            InfDone.append(Inf)
    TrueDone=list(set(TrueDone))
    InfDone=list(set(InfDone))
    return TrueDone,InfDone	
def CountDifNum(Seq0,Seq1):
            Len=len(Seq0)		
            Dif=0
            c=0
            while c<Len:
                if Seq0[c]!=Seq1[c]: Dif+=1
                c+=1
            return Dif

def ReadNodeMap_fromNwk(Nwk): #str
    tree = Tree(Nwk)#('(((A,B,E),C),D);')
        
    edge = 0
    for node in tree.traverse():
       if not node.is_leaf():
          node.name = "NODE%d" %edge  
          edge += 1
    	  
    	     
    Dec2Anc={}	  
    for node in tree.traverse():
       if not node.is_leaf():

          for Chi in node.children:	  
              Dec2Anc[Chi.name]=node.name
    print Dec2Anc
	
    Node2AncClo={}	  
    for node in tree.traverse():
        if node.is_leaf()==True:	

            if node.dist<=0.0001:
                Anc=Dec2Anc[node.name]
                Node2AncClo[Anc]=node.name	
    print Node2AncClo
	
    Dec2AncSum={}
    for Dec in Dec2Anc:
        Anc=Dec2Anc[Dec]	
        if Node2AncClo.has_key(Dec)==True: 
            Dec=	Node2AncClo[Dec]	

        if Node2AncClo.has_key(Anc)==True:
            Anc=Node2AncClo[Anc]	
					
        if Anc!=Dec and Anc!='Normal' and Dec!='Normal':			

           Dec2AncSum[Dec]=Anc
		   
    print Dec2AncSum

    return Dec2Anc,Node2AncClo,Dec2AncSum	


	
def GetRoot(Dec2Anc,Node):			  
    Code=Node in Dec2Anc				
    while Code==True: 
                    Node=Dec2Anc[Node]
                    Code=Node in Dec2Anc				
    return Node	
	
def GetPair(Ind,A2D,D2A,Index2Clo):		 
             Decs=A2D[D2A[Ind]]
             Pair=''
             for Dec in Decs: 
                if Dec!=Ind: Pair=Index2Clo[Dec]
             return Pair
			 
def CountBad(File):
    File=open(File,'r').readlines()[1:]
    Bad=0
    Tot=0	
    for i in File:
       if i.find('Good')==-1: Bad+=1
       Tot+=1	   
    return Bad,Tot	  
	
def GetHead(Head):
    Head=Head.strip().split('\t')
    Len=len(Head)
    c=0
    Name2Col={}
    NameOrder=[]	
    while c<Len:
        Name2Col[Head[c]]=c
        NameOrder.append(Head[c])		
        c+=1
    return NameOrder,Name2Col 	

def GetHeadObsFreqTaHead(Head): 
  Head=Head.strip().split('\t')
  SampOrder=[]
  Name2Col={}
  c=0
  Len=len(Head)
  while c<Len:
    i=Head[c]
    if i.find(':')!=-1:
        Samp=i.split(':')[0]
        Code=Samp in SampOrder
        if Code!=True:
             SampOrder.append(Samp)
        Name2Col[i]=c
    c+=1
  return SampOrder,Name2Col
 
def ListColStr(File):
  File=open(File,'r').readlines()
  NameOrder,Name2Col=GetHead(File[0])
  File=File[1:]
  Tu2Freq={}
  for Tu in NameOrder:
    Tu2Freq[Tu]=[]
  for i in File:
    i=i.strip().split('\t')
    for Tu in Name2Col:
        Tu2Freq[Tu].append(i[Name2Col[Tu]])
  return Tu2Freq	  
  
					  
def Classify(File):
    In=open(File,'r').readlines()
    Posi2Cla={}
    Cla2Posi={'DecAnc':[],'AncDec':[],'Sib':[],'Clu':[]}	
    for i in In:
       i=i.strip().split('\t')
       Posi2Cla[i[0]]=i[1]
       if Cla2Posi.has_key(i[1])!=True: Cla2Posi[i[1]]=[]
       Cla2Posi[i[1]].append(i[0])	   
    return Cla2Posi
def Count1(Tdic, Idic, Cate):
      if Tdic.has_key(Cate)!=True: Tposi=[]
      else: Tposi=Tdic[Cate]
      if Idic.has_key(Cate)!=True: Iposi=[]	  
      else: Iposi=Idic[Cate]
      NotFouncd=0
      for Posi in Tposi:
         if Iposi.count(Posi)==0: NotFouncd+=1
      Incorrect=0
      for Posi in Iposi:
         if Tposi.count(Posi)==0: Incorrect+=1	 
      return NotFouncd, Incorrect	
def GetCommon(CLs,C2S):
    print 'test',CLs
    All=C2S.keys()
	
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
def MakeMLTEin(Meg,OriTree,OutMLTEDfile,Rpath): #OriTree='(out,(G,(A,H)),((B,C),(D,E,F)));'
  if OriTree=='NA':
       OriTree=Meg2MP(Meg,'#hg19',Rpath)
  if OriTree.strip()=='NA': return  OriTree
  else:  
    Anc2Dec,Dec2Anc,TipLs = Nwk2NodeMap(OriTree)
    Clade2Clone = ListClade(Anc2Dec,TipLs)	
    CloLs,Clo2Seq,AA=ReadMegSeq(Meg)
    Len=len(Clo2Seq[Clo2Seq.keys()[0]])	
    Anc2DecIn=''
    CladeIn={}
    AssignMut=[]	
    print 'root',len(Anc2Dec[None]),Anc2Dec[None]
	
    if len(Anc2Dec[None])==1: 
       NoRoot='y'
       RootClade=Anc2Dec[None][0]
    else: 
       NoRoot=''	
       RootClade='Root'	
    print RootClade
	
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
       else:
          print Anc,DecLs,RmTip	   
   
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
    print UnAssign	
    CladeIn[RootClade]=	CladeIn[RootClade][:-1]+','+UnAssign+'\n'
    CladeInStr=''
    for Clade in CladeIn:
        CladeInStr+=CladeIn[Clade]	
    GetOut(OutMLTEDfile,CladeInStr+Anc2DecIn) 
    return  OriTree 
def CloneAnno(CloAnno):

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
              if InfClo2TrueClo2GE.has_key(Inf)!=True: InfClo2TrueClo2GE[Inf]={}
              InfClo2TrueClo2GE[Inf][Tru]=GE
              if TrueClo2InfClo2GE.has_key(Tru)!=True: TrueClo2InfClo2GE[Tru]={}
              TrueClo2InfClo2GE[Tru][Inf]=GE			  

            InfClo2TrueClo={}
            for InfClo in InfClo2TrueClo2GE:
                 TrueClo2GE=InfClo2TrueClo2GE[InfClo]
                 TruCloLs=TrueClo2GE.keys()
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
            return InfClo2TrueClo,AllTru,TrueClo2InfClo2GE
def DupRenameClone_MPtree(InfFile,InfClo2TrueClo,AllTru,TrueClo2InfClo2GE,Rpath):			
            CloLs, Clo2Seq, outMeg=ReadMegSeq(InfFile)
            AllTru=list(set(AllTru))
            Tru2Count={}	
            InfCloLs=[]			
            for Clo in Clo2Seq:
               if Clo!='#hg19':			
                  TruClo=InfClo2TrueClo[Clo[1:]]
                  if Tru2Count.has_key(TruClo)!=True: Tru2Count[TruClo]=1
                  else:
                      Tru2Count[TruClo]+=1					  
                      TruClo+='Dup'+str(Tru2Count[TruClo]-1)
				  
                  outMeg+='#'+TruClo+'\n'+	Clo2Seq[Clo]+'\n'
                  InfCloLs.append(TruClo)
            MissAdd=''
            for Tc in AllTru:
                if Tru2Count.has_key(Tc)!=True: 
                      InfDic=TrueClo2InfClo2GE[Tc]
                      InfLs=InfDic.keys()
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
            InfTree=Meg2MP(InfFile[:-4]+'_tmp.meg','#hg19',Rpath)
            os.remove(InfFile[:-4]+'_tmp.meg')			
            return InfTree.strip(),InfCloLs
def CheckOver(NewCloSet,PreBestSetLs):
    Good='n'
    for Pre in PreBestSetLs:
       	 
         All='y'
         for NClo in NewCloSet:
             if Pre.count(NClo)==0: All='n'
         if All=='y': Good='y'
    return Good					
def PruneTree(KeepCloSet, Nwk_str):
    from ete2 import Tree
    t=Tree(Nwk_str)
    t.prune(KeepCloSet)
    t1=t.write(format=1)	

    return t1
	  
def GetPruneComb(input0,pruneNum):
    Tot=len(input0)
    input=[]
    for i in input0:
        input.append(i.replace('#',''))	
    PruC2Comb={0:[input]}	
    output = sum([map(list, combinations(input, i)) for i in range(len(input) + 1)], [])

    PrunC=1
    while PrunC<=pruneNum:
        NumEl=Tot-PrunC
        Ls=[]
        for i in output:
            if len(i)==NumEl: Ls.append(i)		
        PruC2Comb[PrunC]=Ls
        PrunC+=1		
    return PruC2Comb	
def DoRF_fromtree(TrueTree,InfTree):

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
     GetOut(Current_dir+'PhyloNetIn.nex',PhyloNETin)	


     os.system('java -jar '+PhyloNet +' '+ Current_dir+'PhyloNetIn.nex>'+Current_dir+'RF.out')
     os.chdir(Current_dir)
     if os.path.exists(Current_dir+'RF.out')==True:	 
       RF=open(Current_dir+'RF.out','r').readlines()[2:]
	 
     for i in RF:
         Res+='\t'+i.split(':')[-1].strip()
     import time		 
     time.sleep(1)		 
   		 
     if os.path.exists(Current_dir+'RF.out')==True: os.remove(Current_dir+'RF.out')
     if os.path.exists(Current_dir+'PhyloNetIn.nex')==True: os.remove(Current_dir+'PhyloNetIn.nex')	 
     return Res		
def DoRF_repeatPrune(TrueTree,InfTree,pruneNum):
    from Bio import Phylo
    from cStringIO import StringIO	
    Ttree=Phylo.read(StringIO(TrueTree), "newick")
    TrueFound_Ls0=Ttree.get_terminals()
    TrueFound_Ls=[]
    for i in TrueFound_Ls0:
        if i.name!='Normal': TrueFound_Ls.append('#'+i.name)	
    print TrueFound_Ls

    PruneC2CombLs=GetPruneComb(TrueFound_Ls,pruneNum)	

    TrueFound_Ls.append('#Normal')		 
      
    PruneC=0

    RFoutIn=''  
    PreBestSet=PruneC2CombLs[0] 
    if TrueTree=='' or InfTree=='': return '\tno tree'	
    else:	
     while PruneC<=pruneNum:
        PruneCombLs=PruneC2CombLs[PruneC]

        BestRF=9999
        BestOut=''
        BestSet=[]		
        for PrCloSet in PruneCombLs:
           Good=CheckOver(PrCloSet,PreBestSet)
		   
           if Good=='y':	
	
            TruPr=PruneTree(PrCloSet, TrueTree) 
            InfPr=PruneTree(PrCloSet, InfTree) 
		
            Res=DoRF_fromtree(TruPr,InfPr)

            if Res!='':
                RF=float(Res.split('\t')[7])

                if RF<BestRF: 
                    BestRF=RF
                    BestOut='\t'+Res.split('\t')[1]+'\t'+Res.split('\t')[2]+'\t'+Res.split('\t')[3]+'\t'+Res.split('\t')[4]
                    BestSet=[PrCloSet]
                elif RF==BestRF: BestSet.append(PrCloSet)
        PreBestSet=BestSet				
        if BestOut=='': BestOut='\tNA'*4
        RFoutIn+=BestOut
        PruneC+=1
	 
     return RFoutIn			
def Dotreespace(TTre,ITre,TCloLs,ICloLs,Rpath):
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
    treespaceIn='library(treespace)\n\n'
    treespaceIn+='Tree1 <- read.tree(text=\"'+TTre+'\")\n'
    treespaceIn+='Tree2 <- read.tree(text=\"'+ITre.strip()+'\")\n'
    treespaceIn+='trees <- list(Tree1,Tree2)\n'
    treespaceIn+='df <- data.frame (category = c('+CatLs+'),'+'tiplabels = c('+TipLs+'))\n'
    treespaceIn+='dists <- relatedTreeDist(trees,df)\n'
    treespaceIn+='write.table(dists[1],\''+CurDir.replace('\\','\\\\')+'\\\\test.out'+'\')\n'	
    GetOut('RunTreeSpace.r',treespaceIn)
    os.system(Rpath+' RunTreeSpace.r')
    Dist=open('test.out','r').readlines()[1].strip().split(' ')[-1]	
	
    print 'dist',Dist,'\n'	
    return Dist	
def FindMissingTruClo(InfFile,InfClo2TrueClo,AllTru,TrueClo2InfClo2GE):			
            CloLs, Clo2Seq, outMeg=ReadMegSeq(InfFile)
            AllTru=list(set(AllTru))
            Tru2Count={}	
            InfCloLs=[]			
            for Clo in Clo2Seq:
               if Clo!='#hg19':			
                  TruClo=InfClo2TrueClo[Clo[1:]]
                  if Tru2Count.has_key(TruClo)!=True: Tru2Count[TruClo]=1
                  else:
                      Tru2Count[TruClo]+=1					  
                      TruClo+='Dup'+str(Tru2Count[TruClo]-1)
				  
                  outMeg+='#'+TruClo+'\n'+	Clo2Seq[Clo]+'\n'
                  InfCloLs.append(TruClo)
            MissAdd=''
            MissTruCloLs=[]			
            for Tc in AllTru:
                if Tru2Count.has_key(Tc)!=True:
                      MissTruCloLs.append(Tc)				

	
            return MissTruCloLs
def FindAnc(Meg):
     CloLs,Clo2Seq,AA=ReadMegSeq(Meg)
     AncLs=[]
     for Clo1 in CloLs:
       if Clo1!='#hg19':	 
         Seq1=Clo2Seq[Clo1]
         Anc='n'
         Len=len(Seq1)		 
         for Clo2 in CloLs:
             if Clo1!=Clo2:
                  Seq2=Clo2Seq[Clo2]
                  c=0
                  Uni='n'
                  while c<Len:
                     if Seq1[c]=='T' and Seq2[c]=='A': Uni='y'
                     c+=1
                  if Uni=='n': Anc='y'
         if Anc=='y': AncLs.append(Clo1[1:])				  
     return AncLs	
def GetOut(File,In):
    OutF=open(File,'w')
    OutF.write(In)
    OutF.close()	
def GetOut_from_list(List,Head,OutF):
     out=Head
     for i in List:
        out+=i+'\n'
     GetOut(OutF,out)		
	 