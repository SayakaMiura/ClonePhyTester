import Functions
import sys

Meg=sys.argv[1]
Out=Meg[:-4]+'_MutOrder.txt'

CloLs, Clo2Seq, AA=Functions.ReadMegSeq(Meg)
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
                print '??',c0,c1,Clo
                open('S','r').readlines()
        if 	Ls1Uni==[] and Ls0Uni==[] and Com!=[] : out+=str(c0)+'-'+str(c1)+'\tClu\n'
        elif 	Ls1Uni!=[] and Ls0Uni==[] and Com!=[] : out+=str(c0)+'-'+str(c1)+'\tDecAnc\n'		
        elif 	Ls1Uni==[] and Ls0Uni!=[] and Com!=[] : out+=str(c0)+'-'+str(c1)+'\tAncDec\n'			
        elif 	Ls1Uni!=[] and Ls0Uni!=[] : out+=str(c0)+'-'+str(c1)+'\tSib\n'
        elif 	Ls1Uni==[] and Ls0Uni==[] and Com==[] and ComUn!=[]: out+=str(c0)+'-'+str(c1)+'\tUnassign\n'
        elif 	(Ls1Uni!=[] or Ls0Uni!=[]) and Com==[] : out+=str(c0)+'-'+str(c1)+'\tUnassign\n'		
        else:
                out+=str(c0)+'-'+str(c1)+'\tUnassign\n'		
                print '???',c0,c1,Ls1Uni, Ls0Uni, Com, ComUn
                open('S','r').readlines()
        c1+=1
    c0+=1

Functions.GetOut(Out,out)