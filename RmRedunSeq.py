import sys
import Functions
Meg=sys.argv[1]
Out2=Meg[:-4]+'_NoRedun.meg'
NameOrder, Name2Seq, out2=Functions.ReadMegSeq(Meg)

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

    out2+=Name+'\n'+Name2Seq[Name]+'\n'		
    Done+=IdenLs

Out2F=open(Out2,'w')
Out2F.write(out2)
Out2F.close()
