import subprocess
import os

def trim_reads(trim_prog, raw_data_dir, adapter, out_dir,thread):
    '''
    Use trimmomatic-0.39.jar to remove low quantity reads and adapters
    '''

    in_fqs = [
        os.path.join(raw_data_dir, f"{os.path.basename(raw_data_dir)}_1.fq.gz"),
        os.path.join(raw_data_dir, f"{os.path.basename(raw_data_dir)}_2.fq.gz")
        ]
    out_fqs = [
        os.path.join(out_dir, f"{os.path.basename(raw_data_dir)}_filtered_R1.fastq.gz"),
        os.path.join(out_dir, f"{os.path.basename(raw_data_dir)}_unpaired_R1.fastq.gz"),
        os.path.join(out_dir, f"{os.path.basename(raw_data_dir)}_filtered_R2.fastq.gz"),
        os.path.join(out_dir, f"{os.path.basename(raw_data_dir)}_unpaired_R2.fastq.gz")
    ]
    trim_log = os.path.join(out_dir, "read_surviving_stat.txt")

    os.makedirs(out_dir, exist_ok=True)

    trim_cmd = f"java -jar {trim_prog} PE {' '.join(in_fqs)} {' '.join(out_fqs)} ILLUMINACLIP:{adapter}:2:30:10 LEADING:3 TRAILING:3 SLIDINGWINDOW:4:15 MINLEN:36 -threads {thread} 2> {trim_log}"
    subprocess.run(trim_cmd, shell=True, check=True)

if __name__ == '__main__':
    trim_prog = "/home/public/softwares/Trimmomatic-0.39/trimmomatic-0.39.jar"
    #raw_data_dir = "/path/to/raw_data"
    adapter = "/home/public/softwares/Trimmomatic-0.39/adapters/PE_adapters.fa"
    #out_dir = "/path/to/output"
    thread = 8
    trim_reads(trim_prog,raw_data_dir,adapter,out_dir,thread)
