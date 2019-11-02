import Functions


Input = config['filename']
SeqLs, Rpath, MLTEDpath, Summary, InputParse = Functions.parse_input(Input)
excel_out = Input[:Input.rfind('.')]+"_summary.xlsx"
if InputParse != 'y':
    print (Summary)
    exit(1)

Sim2TrueFile={'G7':'G7/True/G7-ID.meg','G12':'G12/True/G12-ID.meg','P10':'P10/True/ID_True.meg','MA':'MA/True/ID_True_NoRedun.meg','G7cna':'G7/True/G7-ID.meg','MA50':'MA/True/ID_True_NoRedun.meg','TGlinear':'TGlinear/True/IDCells_True.meg','TGstep':'TGstep/True/IDCells_True.meg','TGconst':'TGconst/True/IDCells_True.meg'}

file_list = [l.strip().split('\t')[2] for l in SeqLs]
mid_list = [l[:-4] for l in file_list]
param_files = [l+"_start_params.txt" for l in mid_list]
#seq_idx = range(1,len(SeqLs))

#print (open(param_files[0]).readlines())
#print (os.path.getsize(param_files[0]))

for i in range(len(mid_list)):
    param_fn = param_files[i]
    param_f = open(param_fn,'w')
    param_f.write(SeqLs[i].strip()+'\n')
    param_f.close()


rule all:
    input: 
        excel_out

rule runseq:
    input:
        data="{filename}_start_params.txt"
    output:
        "{filename}_result.txt"
    params:
        conf=Input
    shell:
        "python3 clonephytester_snakemake.py {params.conf} {input.data} {output}"

rule test:
    input: 
        "{filename}.meg"
    output:
        "{filename}_test.txt"
    shell:
        "cat {input} > {output}"

rule summarize:
    input:
        expand("{filename}_result.txt",filename=mid_list)
    output:
        excel_out
    params:
        conf=Input
    shell:
        "python3 collect_results.py {params.conf} {input}"