import os
import subprocess

def stringtie_assembly(stringtie_prog, in_bam, in_gtf, out_gtf, thread):
    '''
    StringTie version=2.1.8 to assembly
    '''
    assembly_cmd = f"{stringtie_prog}  {in_bam} -G {in_gtf} -o {out_gtf} -p {thread}"
    subprocess.run(assembly_cmd, shell=True, check=True)

def stringtie_merge(stringtie_prog,guide_gtf,output_dir,samples, merged_gtf):
    '''
    merge all gtf
    '''
    gtf_files = [os.path.join(output_dir,sp, f"{os.path.basename(sp)}_stringtie_raw.gtf") for sp in samples]
    merge_cmd = f"{stringtie_prog} --merge -G {guide_gtf} -o {merged_gtf} {' '.join(gtf_files)}"
    subprocess.run(merge_cmd, shell=True, check=True)

def stringtie_quant(stringtie_prog, in_bam, in_gtf,out_gtf,gene_table,thread):
    '''
    quantification of transcripts
    '''
    quant_cmd = f"{stringtie_prog}  {in_bam} -A {gene_table} -e -G {in_gtf} -o {out_gtf} -p {thread}"
    subprocess.run(quant_cmd, shell=True, check=True)

if __name__ == '__main__':
    stringtie_prog = "/home/haowu/software/stringtie/stringtie"
    in_gtf = "/home/haowu/reference/Gencode/Human/Gencode_v38/gencode.v38.primary_assembly.annotation.gtf"
    thread = 8
    stringtie_assembly(stringtie_prog, in_bam, in_gtf, out_gtf, thread)
    stringtie_merge(stringtie_prog,output_dir,input_dirs, merged_gtf)