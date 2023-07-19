import argparse
import os
import sys
import datetime
import glob
import pandas as pd
from Trim import trim_reads
from Align import star_align
from Assembly import stringtie_assembly,stringtie_merge,stringtie_quant
from Quantity import extra_ge_tr_count,extra_transcript_tpm
_MAN = \
'''
##########################################
Analysis illumina RNA seq Data
##########################################
'''

def illumina_RNA_pipeline(args):
    input_dirs = os.listdir(args.input_dir)  # 获取输入目录下的所有文件夹
    for input_dir in input_dirs:
        # 构建每个样本的输入和输出路径
        sample_input_dir = os.path.join(args.input_dir, input_dir)
        sample_output_dir = os.path.join(args.output_dir, input_dir)
    
        # trim reads
        print(f"### remove low quantity reads and adapters using trimmomatic {input_dir}",
              datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        trim_reads(trim_prog=args.trimmomatic_path, raw_data_dir=sample_input_dir, adapter=args.adapter_path,
                   out_dir=sample_output_dir, thread=args.threads)
        
        # align reads to genome
        print(f"### align reads to genome using STAR {input_dir}",
              datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        star_align(star_prog=args.star_path, index=args.index, gtf=args.gtf, infq_dir=sample_output_dir, thread=args.threads)

        # transcript assembly
        print(f"### transcript assembly using StingTie {input_dir}",
              datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        in_bam = os.path.join(sample_output_dir,
                              f"{os.path.basename(sample_output_dir)}_alignmentAligned.sortedByCoord.out.bam")
        out_gtf = os.path.join(sample_output_dir,
                               f"{os.path.basename(sample_output_dir)}_stringtie_raw.gtf")
        stringtie_assembly(stringtie_prog=args.stringtie_path, in_bam=in_bam, in_gtf=args.gtf, out_gtf=out_gtf,
                           thread=args.threads)
    
     # merge gtf
    print("### merge gtf",
              datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
    stringtie_merge(stringtie_prog=args.stringtie_path,guide_gtf=args.gtf,output_dir=args.output_dir,samples=input_dirs,merged_gtf=args.merged_gtf)

    for input_dir in input_dirs:
        sample_output_dir = os.path.join(args.output_dir, input_dir)
        print(f"### transcript quantification using StingTie {input_dir}",
              datetime.datetime.now().strftime("%Y-%m-%d %H:%M:%S"))
        in_bam = os.path.join(sample_output_dir,
                              f"{os.path.basename(sample_output_dir)}_alignmentAligned.sortedByCoord.out.bam")
        out_gtf = os.path.join(sample_output_dir,
                               f"{os.path.basename(sample_output_dir)}_stringtie_final.gtf")
        gene_table = os.path.join(sample_output_dir,
                               f"{os.path.basename(sample_output_dir)}_gene_abund.tab")
        stringtie_quant(stringtie_prog=args.stringtie_path, gene_table= gene_table, in_bam=in_bam, in_gtf=args.merged_gtf, out_gtf=out_gtf,
                           thread=args.threads)
        
    # prepare gtf_list file for running prepDE
    with open(os.path.join(args.output_dir,"all_gtf"),"w") as gtf_file:
        for input_dir in input_dirs:
            sample_output_dir = os.path.join(args.output_dir, input_dir)
            out_gtf = os.path.join(sample_output_dir,
                               f"{os.path.basename(sample_output_dir)}_stringtie_final.gtf")
            gtf_file.write(input_dir+'\t'+out_gtf+'\n') 


    # extra gene and transcript count matrix
    print(f"### extra gene and transcript count matrix using prepDE",
          datetime.datetime.now().strftime("%Y-%m-%d %H:%M%S"))
    gtf_list = os.path.join(args.output_dir,"all_gtf")
    extra_ge_tr_count(prepde_prog=args.prepDE,gtf_list=gtf_list,out_dir=args.output_dir)

    # extra transcript TPM for each sample
    for input_dir in input_dirs:
        sample_output_dir = os.path.join(args.output_dir, input_dir)
        in_gtf = os.path.join(sample_output_dir,
                            f"{os.path.basename(sample_output_dir)}_stringtie_final.gtf")
        out_file = os.path.join(sample_output_dir,
                            f"{os.path.basename(sample_output_dir)}_transcript_TPM.txt")
        extra_transcript_tpm(gtf= in_gtf,out_file=out_file)

    # merge all samples transcript TPM matrix
    sample_files = glob.glob(args.output_dir + '/' + '**/*_transcript_TPM.txt') #**表示递归匹配任意目录。当使用 ** 时，它将匹配目标目录及其所有子目录中的文件。
    dataframes = []
    for file in sample_files:
        sample = os.path.basename(file).split("_")[0]
        df = pd.read_csv(file,sep="\t")
        df.rename(columns={'TPM':sample},inplace=True)
        dataframes.append(df)
    merge_df = None
    for df in dataframes:
        if merge_df is None:
            merge_df = df
        else:
            merge_df = pd.merge(merge_df,df,on='Transcript_id',how='outer')
    merge_df = merge_df.fillna(0)
    out_file = os.path.join(args.output_dir,"transcript_TPM_matrix.csv")
    merge_df.to_csv(out_file,index=False)
    
     # merge all samples gene TPM matrix
    sample_files = glob.glob(args.output_dir + '/' + '**/*_gene_abund.tab') #**表示递归匹配任意目录。当使用 ** 时，它将匹配目标目录及其所有子目录中的文件。
    dataframes = []
    for file in sample_files:
        sample = os.path.basename(file).split("_")[0]
        df = pd.read_csv(file,sep="\t")
        df = df[["Gene ID","Gene Name","TPM"]]
        df.rename(columns={'TPM':sample},inplace=True)
        dataframes.append(df)
    merge_df = None
    for df in dataframes:
        if merge_df is None:
            merge_df = df
        else:
            merge_df = pd.merge(merge_df,df,on=['Gene ID',"Gene Name"],how='outer')
    merge_df = merge_df.fillna(0)
    out_file = os.path.join(args.output_dir,"gene_TPM_matrix.csv")
    merge_df.to_csv(out_file,index=False)
    





if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=_MAN)
    parser.add_argument("-t", "--trimmomatic_path", help="path of trimmomatic bin file",
                        default="/home/public/softwares/Trimmomatic-0.39/trimmomatic-0.39.jar")
    parser.add_argument("-s", "--star_path", help="path of STAR bin file", default="/usr/bin/STAR")
    parser.add_argument("--stringtie_path", help="path of StringTie bin file",
                        default="/home/haowu/software/stringtie/stringtie")
    parser.add_argument("-i", "--input_dir", help="directory of all input samples ", required=True)
    parser.add_argument("--index", help="index of STAR")
    parser.add_argument("-g", "--gtf", help="path of gtf")
    parser.add_argument("-a", "--adapter_path", help="path of adapter file",
                        default="/home/public/softwares/Trimmomatic-0.39/adapters/PE_adapters.fa")
    parser.add_argument("-o", "--output_dir", help="directory of out files", required=True)
    parser.add_argument("-m", "--merged_gtf", help="path of stringtie merged gtf", required=True)
    parser.add_argument("-@", "--threads", type=int, help="number of threads to use", default=4)
    parser.add_argument("-p", "--prepDE", help="path of prepDE.py3 file", 
                        default="/home/haowu/software/stringtie/prepDE.py3")

    args = parser.parse_args()

    illumina_RNA_pipeline(args)
