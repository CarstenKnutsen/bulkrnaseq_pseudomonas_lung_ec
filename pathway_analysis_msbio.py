'''Aligning Hough samples
author: Carsten Knutsen
Date: March 19 2023
conda environment: "msbio"
MSBio v3.5.230101
'''

import subprocess as sp
import os
import shutil
import pandas as pd
import numpy as np
import time
import json


deg_fn = 'data/output_data/deg_tests.xlsx'
output_fol = 'data/output_data/metascape_pathway'

def deg_to_metascape(deg_fn,
                         output_fol,
                         msbio_fol = '/home/carsten/MSBio/msbio_v3.5.20230101/',
                         species='10090', # Takes taxonomy id 10090 is mouse
                         fc_name='logFC',
                         fc_threshold=1,
                         pvalue_name='FDR',
                         pvalue_threshold=0.05):
    '''Takes a .xlsx file with DEGs and runs metascape on up and downregulated genes from each sheet'''

    ## tmp file paths for msbio
    print('writing tmp files')
    input_fol = f'{msbio_fol}data/tmp/input/'
    tmp_output_fol = f'{msbio_fol}data/tmp/output/'
    job_file = f'{input_fol}batch.job'
    ### write csv files for gene lists
    for path in [input_fol, tmp_output_fol]:
        if os.path.exists(path):
            shutil.rmtree(path)
        os.makedirs(path, exist_ok=True)
    degs = pd.read_excel(deg_fn, sheet_name=None, index_col=0)
    for key in degs.keys():
        df = degs[key]
        df = df.loc[(abs(df[fc_name]) > fc_threshold) &
                    (df[pvalue_name] < pvalue_threshold)]
        df.loc[df[fc_name] > 0].index.to_series(name='Gene').to_csv(f'{input_fol}{key}_up.csv', index=False) ## write the upregulated genes tmp files
        df.loc[df[fc_name] < 0].index.to_series(name='Gene').to_csv(f'{input_fol}{key}_down.csv', index=False) ## write the downregulated genes tmp files
    ### make job file
    json_dict = {}
    for fn in os.listdir(input_fol):
        tmp_output_fol2 = f'{tmp_output_fol}{fn.split(".")[0]}'
        os.makedirs(tmp_output_fol2, exist_ok=True)
        tmp_dict = {}
        tmp_dict['input'] = f'/{"/".join(f"{input_fol}{fn}".split("/")[-4:])}'
        tmp_dict['output'] =f'/{"/".join(f"{tmp_output_fol2}".split("/")[-4:])}'
        # tmp_dict['single'] = 'true'.strip("'")
        # tmp_dict['source_tax_id'] = species
        json_dict[fn.split('.')[0]] = tmp_dict
    with open(job_file, "w") as outfile:
        outfile.write(
            '\n'.join(json.dumps(json_dict[k]) for k in sorted(json_dict.keys()))
        )
    ## run MSBIO
    print('running MSBio')
    # commands
    up = f'bin/up.sh'
    job_run = f'bin/ms.sh /data/tmp/input/batch.job -u -S {species}'
    down = f'bin/down.sh'
    unlock = f'sudo chown -R $USER ~/MSBio'
    # run commands
    sp.run(up, cwd = msbio_fol, shell=True)
    sp.run(job_run, cwd = msbio_fol, shell=True)
    sp.run(down, cwd = msbio_fol, shell=True)
    sp.run(unlock, cwd = msbio_fol, shell=True)
    # copy to final location
    shutil.rmtree(output_fol)
    shutil.copytree(tmp_output_fol, output_fol)
    for path in [input_fol, tmp_output_fol]:
        if os.path.exists(path):
            shutil.rmtree(path)

if __name__ =='__main__':
    deg_to_metascape(deg_fn, output_fol)



