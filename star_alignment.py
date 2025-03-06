'''Aligning Hough samples
author: Carsten Knutsen
Date: March 13 2023
conda environment: "alignment"
STAR v2.10.a
'''

import subprocess as sp
import os
import shutil
import pandas as pd

def map_sample(rootdir, outputdir):
    samples = sorted(set([x.split('.')[0].split('_')[0] for x in os.listdir(rootdir)]))
    for sample in samples:
        print(sample)
        fq1 = f'{rootdir}/{sample}_R1_001.fastq.gz'
        fq2 = f'{rootdir}/{sample}_R2_001.fastq.gz'
        tmp_dir = f'{outputdir}/{sample}/STARtmp'
        os.makedirs(f'{outputdir}/{sample}', exist_ok=True)
        if os.path.exists(tmp_dir):
            shutil.rmtree(tmp_dir)

        try:
            if not os.path.exists(f'{outputdir}/{sample}/ReadsPerGene.out.tab'):
                call = [os.getenv('STAR', 'STAR'),
                        '--runThreadN', '6'
                        '--runMode', 'alignReads',
                        '--genomeDir', '/media/carsten/hdd/genomes/mouse/GRCm39_index',
                        '--genomeLoad', 'LoadAndKeep',
                        '--outTmpDir', f'{outputdir}/{sample}/STARtmp',
                        '--readFilesIn', fq1, fq2,
                        '--readFilesCommand', 'zcat',
                        '--outFilterType', 'BySJout',
                        '--outFilterMultimapNmax', '20', #default
                        '--alignSJoverhangMin', '8', #default
                        '--alignSJDBoverhangMin', '1', #default
                        '--outFilterMismatchNmax 999', #default
                        '--outFilterMismatchNoverLmax 0.04', #default
                        '--alignIntronMin', '20', #default
                        '--alignIntronMax', '1000000', #default
                        '--alignMatesGapMax', '1000000', #default
                        '--outSAMstrandField', 'intronMotif',
                        '--outFileNamePrefix', f'{outputdir}/{sample}/',
                        '--outSAMtype', 'BAM', 'Unsorted',
                        '--outSAMattributes', 'NH', 'HI', 'AS', 'NM',
                        '--outFilterMatchNminOverLread', '0.4',
                        '--outFilterScoreMinOverLread', '0.4',
                        '--outReadsUnmapped', 'Fastx',
                        '--quantMode', 'GeneCounts'
                ]
                sp.run(call, check=True)
            else:
                continue
        except:
            raise OSError(f'STAR mapping failed on {dir}')
            continue

    call2 = ['STAR',
             '--runMode', 'alignReads',
             '--genomeDir', '/media/carsten/hdd/genomes/mouse/GRCm39_index',
             '--genomeLoad', 'Remove']
    sp.run(call2, check=True)
def generate_csv(directory, fn = 'counts.csv'):
    ##Need to change the pandas csv stitching for larger runs, works for now
    for index, dir in enumerate(os.listdir(output_dir)):
        if index == 0:
            read = f'{output_dir}/{dir}/ReadsPerGene.out.tab'

            table = pd.read_csv(read,
                                sep='\t',
                                header=None,
                                index_col=0)
            counts_df = pd.DataFrame(table[1])
            counts_df.rename(columns={1: f'{dir}'}, inplace=True)
        else:
            break

    for index, dir in enumerate(os.listdir(output_dir)):
        if index == 0:
            continue
        else:
            read = f'{output_dir}/{dir}/ReadsPerGene.out.tab'
            table = pd.read_csv(read,
                                sep='\t',
                                header=None,
                                index_col=0)
            column = pd.DataFrame(table[1])
            column.rename(columns={1: f'{dir}'}, inplace=True)
            counts_df = pd.concat([counts_df, column], axis=1, join='inner')
    counts_df.to_csv(os.path.join(directory, fn))
if __name__ =='__main__':
    input_dir = f'data/raw/Juvenile_endothelial_ bulk_seq_2022'
    output_dir = f'data/aligned_sequencing'
    os.makedirs(output_dir, exist_ok=True)
    map_sample(input_dir, output_dir)
    generate_csv(output_dir)



