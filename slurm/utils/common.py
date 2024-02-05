import argparse
import json
import linecache
import platform
import re
from copy import copy
from functools import cache
from itertools import groupby
from multiprocessing import Pool
from pathlib import Path
import numpy as np
from numba import njit, prange
from rclone_python import rclone




@cache
def get_config_linenumber(config, key):
    with open(config, "r") as file:
        for idx, line in enumerate(file):
            if line.startswith(key):
                return idx + 1


def get_config_value(config, key):
    lineno = get_config_linenumber(config, key)
    if lineno is None:
        raise LookupError(f"Could not find [{key=}] in  [{config=}]")
    return linecache.getline(config, lineno=lineno).split()[2]


def all_equal(iterable):
    g = groupby(iterable)
    return next(g, True) and not next(g, False)


@njit(cache=True)
def arrays_equal(a, b):
    if a.shape != b.shape:
        return False
    for i in prange(a.size):
        if a[i] != b[i]:
            return False
    return True


def get_entry_indices(dset, expected_num, kind=256):
    return np.where(dset['event'] == kind)[0][:expected_num]




def write_batch_status_single_thread(batch):
    parser = argparse.ArgumentParser(description='SLURM batch generator')
    parser.add_argument('--update-status', action='store_true', help='Update status files', default=False)
    args = parser.parse_args()
    output_prfx = f'{batch["output_prfx"]}'
    config_file = f'{batch["config_file"]}'
    status_file = f'{batch["status_dir"]}/{Path(config_file).stem}.status'
    seed_status = copy(batch.get('seed_status'))  # Old one
    batch['seed_status'] = []
    if platform.node() != "neumann" and args.update_status:
        raise AssertionError("--update-status is only valid on neumann")

    if platform.node() == "neumann" and args.update_status:
        Path(status_file).parent.mkdir(parents=True, exist_ok=True)
        output_path = f'{output_prfx}/{batch["projectname"]}/{batch["output_path"]}'
        print(f"Updating status: {status_file}")
        with open(status_file, 'w') as sf:
            status_count = 0
            for sidx, (offset, extent) in enumerate(zip(batch['seed_offset'], batch['seed_extent'])):
                extent_size = len(batch['seed_extent'])
                offset_size = len(batch['seed_offset'])
                if offset_size != extent_size:
                    raise ValueError(
                        f"offset:{offset_size} and extent:{extent_size} are not equal lengths")

                ## Start checking the h5 files
                is_finished = True
                for seed in range(offset, offset + extent):
                    h5status = None
                    if seed_status is not None:
                        if seed_status[sidx] == "FINISHED":
                            h5status = seed_status[sidx]  # Short circuit
                    if h5status is None:
                        filename = f'{output_path}/{batch["output_stem"]}_{seed}.h5'
                        h5status = get_h5_status(filename=filename, batch=batch)
                    if not "FINISHED" in h5status:
                        print(f'{config_file}: {seed}|{h5status}')
                        is_finished = False
                    sf.write(f'{seed}|{h5status}\n')
                if is_finished:
                    batch['seed_status'].append("FINISHED")
                else:
                    batch['seed_status'].append("PENDING")
        # We can now check if there are stray simulations
        strays_file = f'strays/{Path(config_file).stem}.strays'
        Path(strays_file).parent.mkdir(parents=True, exist_ok=True)
        with open(status_file, 'r') as sf, open(strays_file, 'w') as st:
            for h5file in sorted(Path(output_path).rglob('*.h5')):
                if 'save' in h5file.name:
                    continue
                seed_found = int(re.split('_|\.', h5file.name)[1])
                seed_stray = True
                for offset, extent in zip(batch['seed_offset'], batch['seed_extent']):
                    if offset <= seed_found < offset + extent:
                        seed_stray = False
                if seed_stray is True:
                    print(f"Stray seed found for {batch['config_file']}: {seed_found}")
                    st.write(f'{h5file}\n')
    else:
        if platform.node() == "neumann":
            status_file = "{}/{}/{}".format(batch['output_prfx'], batch['projectname'], status_file)
        else:
            rclone.copyto(f'neumann:{batch["output_prfx"]}/{batch["projectname"]}/{status_file}',
                          f'{status_file}',
                          args=['-L', '--update'])
        if not os.path.isfile(status_file):
            raise FileNotFoundError(f"{status_file}")

        status_count = 0
        for offset, extent in zip(batch['seed_offset'], batch['seed_extent']):
            extent_size = len(batch['seed_extent'])
            offset_size = len(batch['seed_offset'])
            if offset_size != extent_size:
                raise ValueError(
                    f"offset:{offset_size} and extent:{extent_size} are not equal lengths")
            is_finished = True
            for idx, seed in enumerate(range(offset, offset + extent)):
                sfline = linecache.getline(status_file, idx + status_count + 1).rstrip()
                # print(idx, seed, sfline)
                sfseed, sfstatus = sfline.split('|', maxsplit=1)
                if seed != int(sfseed):
                    raise ValueError(f'seed mismatch [{seed=}] != [{sfseed=}]')
                if sfstatus != "FINISHED":
                    is_finished = False
                    break
                # print(linecache.getline(status_file, idx+status_count))
            status_count += extent
            if is_finished:
                batch['seed_status'].append('FINISHED')
            else:
                batch['seed_status'].append('PENDING')

    return batch



def update_batch_status_single_thread(config_paths):
    for batch_filename in sorted(Path(config_paths['config_dir']).rglob('*.json')):
        with open(batch_filename, 'r') as fp:
            batchjson = write_batch_status_single_thread(json.load(fp))

        with open(batch_filename, 'w') as fp:
            json.dump(batchjson, fp, sort_keys=True, indent=4)

