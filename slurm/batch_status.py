import os
from pathlib import Path
import re
import numpy as np
import json
import h5py
from itertools import groupby

def all_equal(iterable):
    g = groupby(iterable)
    return next(g, True) and not next(g, False)

def get_h5_status(filename, batch):
    if os.path.isfile(filename):
        expected_dset_paths = [
            'common/finished_all',
            'fLBIT/state_real/measurements',
            'fLBIT/state_real/status',
            'fLBIT/state_real/mem_usage',
            'fLBIT/state_real/number_probabilities',
        ]

        try:
            with h5py.File(filename, 'r') as h5file:
                expected_dsets = [h5file.get(path) for path in expected_dset_paths]
                missing_dsets = [path for dset,path in zip(expected_dsets,expected_dset_paths) if dset is None]
                if len(missing_dsets) > 0:
                    return f"FAILED|missing datasets:{missing_dsets}"
                len_of_dsets = [len(expected_dsets[1]), len(expected_dsets[2]), len(expected_dsets[3]), np.shape(expected_dsets[4])[-1]]
                has_equal_iters = all_equal(len_of_dsets)
                if not has_equal_iters:
                    return f"FAILED|(unequal iters:{len_of_dsets})"

                has_finished_all   = expected_dsets[0][()]
                has_expected_iters = len(expected_dsets[1]) == batch['time_steps']

                if has_expected_iters and has_finished_all:
                    return "FINISHED"
                if has_expected_iters and not has_finished_all:
                    return "TIMEOUT"
                if not has_expected_iters and has_finished_all:
                    return "TIMEOUT"
                return "FAILED"
        except:
            return "FAILED"
    else:
        return "MISSING"


def write_batch_status(batch, output_base):
    config_file = f'{batch["config_file"]}'
    output_path = f'{output_base}/{batch["projectname"]}/{batch["output_path"]}'
    status_file = f'status/{Path(config_file).stem}.status'
    Path(status_file).parent.mkdir(parents=True, exist_ok=True)
    with open(status_file, 'w') as sf:
        for offset, extent in zip(batch['seed_offset'], batch['seed_extent']):
            extent_size = len(batch['seed_extent'])
            offset_size = len(batch['seed_offset'])
            if offset_size != extent_size:
                raise ValueError(
                    f"offset:{offset_size} and extent:{extent_size} are not equal lengths")

            ## Start checking the files
            for seed in range(offset, offset+extent):
                filename = f'{output_path}/{batch["output_stem"]}_{seed}.h5'
                h5status = get_h5_status(filename=filename,batch=batch)
                if not "FINISHED" in h5status:
                    print(f'{config_file}: {seed}|{h5status}')
                sf.write(f'{seed}|{h5status}\n')

    return 0


output_base = '/mnt/WDB-AN1500/mbl_transition'
config_dirs = [
    'config-L[12,16,20]',
    'config-L[24]',
    'config-L[28]',
    'config-L[32]',
    ]

for config_dir in config_dirs:
    for batch in sorted(Path(config_dir).rglob('*.json')):
        with open(batch, 'r') as fp:
            write_batch_status(json.load(fp), output_base)
