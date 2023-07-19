import subprocess
import os

def star_align(star_prog, index, gtf, infq_dir, thread):
    '''
    STAR version=2.7.9a align to genome
    '''

    in_fqs = [
        os.path.join(infq_dir, f"{os.path.basename(infq_dir)}_filtered_R1.fastq.gz"),
        os.path.join(infq_dir, f"{os.path.basename(infq_dir)}_filtered_R2.fastq.gz")
        ]
    prefix = os.path.join(infq_dir, f"{os.path.basename(infq_dir)}_alignment")

    align_cmd = f"{star_prog} --twopassMode Basic --twopass1readsN -1 --runThreadN {thread} --genomeDir {index} --readFilesIn {' '.join(in_fqs)} --readFilesCommand zcat --sjdbGTFfile {gtf} --outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverReadLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000 --outFileNamePrefix {prefix} --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM 79711714937"
    subprocess.run(align_cmd, shell=True, check=True)

if __name__ == '__main__':
    star_prog = "/usr/bin/STAR"
    index = "/home/haowu/reference/Gencode/Human/Gencode_v38/index/STAR"
    gtf = "/home/haowu/reference/Gencode/Human/Gencode_v38/gencode.v38.primary_assembly.annotation.gtf"
    infq_dir = "/path/to/input"
    thread = 8
    star_align(star_prog, index, gtf, infq_dir, thread)
