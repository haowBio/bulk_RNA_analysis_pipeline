import os
import subprocess
import pandas as pd

def extra_ge_tr_count(prepde_prog,gtf_list,out_dir):
    '''
    extra all samples gene and transcript read count using prepDE.py3
    '''
    ge_out = os.path.join(out_dir,"gene_count_matrix.csv")
    tr_out = os.path.join(out_dir,"transcript_count_matrix.csv")
    extra_cmd=f"{prepde_prog} -i {gtf_list} -g {ge_out} -t {tr_out}"
    subprocess.run(extra_cmd,shell=True,check=True)


def extra_transcript_tpm(gtf,out_file):
    '''
    extra transcript TPM expression for each sample
    '''
    transcript_tpm = {}
    with open(gtf,"r") as file:
        for line in file:
            if line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if parts[2] == 'transcript':
                attributes = parts[8].split(";")
                for attr in attributes:
                    attr = attr.strip()
                    if attr.startswith('transcript_id'):
                        transcript_id = attr.split(' ')[1].strip('"')
                    elif attr.startswith('TPM'):
                        tpm = float(attr.split(' ')[1].strip('"'))
                if transcript_id and tpm:
                    transcript_tpm[transcript_id] = tpm
    with open(out_file,"w") as out:
        out.write("Transcript_id\tTPM\n")
        for transcript_id, tpm in transcript_tpm.items():
            out.write(transcript_id+"\t"+str(tpm)+"\n")






if __name__ == '__main__':
    prepde_prog = '/home/haowu/software/stringtie/prepDE.py3'
    extra_ge_tr_count(prepde_prog,gtf_list,out_dir)
    extra_transcript_tpm(gtf,out_file)
