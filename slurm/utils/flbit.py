from multiprocessing import Pool
import numpy as np
import os
import platform
import re
import json
import h5py
import argparse
from pathlib import Path
from itertools import groupby
import linecache
from copy import copy
from rclone_python import rclone
import traceback
import linecache
from functools import cache
from numba import njit,prange

def get_unique_config_string(d: dict, dl: dict, delim: str):
    str_L = str(d['model::model_size'])
    str_J = "[{}±{}_{}±{}_{}±{}]".format(d['model::lbit::J1_mean'], d['model::lbit::J1_wdth'],
                                         d['model::lbit::J2_mean'], d['model::lbit::J2_wdth'],
                                         d['model::lbit::J3_mean'], d['model::lbit::J3_wdth'])
    str_rL = "L" if d['model::lbit::J2_span'] == '-1' else str(d['model::lbit::J2_span'])
    str_x = str(d['model::lbit::xi_Jcls'])
    str_u = str(d['model::lbit::u_depth'])
    str_f = str(d['model::lbit::u_fmix'])
    str_circuit = f"d{str_u}_f{str_f}"

    u_tstd, u_tstds = [x.get('model::lbit::u_tstd') for x in [d, dl]]
    u_tgw8, u_tgw8s = [x.get('model::lbit::u_tgw8') for x in [d, dl]]
    u_cstd, u_cstds = [x.get('model::lbit::u_cstd') for x in [d, dl]]
    u_cgw8, u_cgw8s = [x.get('model::lbit::u_cgw8') for x in [d, dl]]
    u_bond, u_bonds = [x.get('flbit::cls::mpo_circuit_svd_bondlim') for x in [d, dl]]

    str_circuit += f"_tw{u_tstd}" if u_tstds is not None and (len(u_tgw8s) > 1 or len(u_tstds) > 1) else ''
    str_circuit += f"{u_tgw8[:2]}" if u_tgw8s is not None and (len(u_tgw8s) > 1 or len(u_tstds) > 1) else ''
    str_circuit += f"_cw{u_cstd}" if u_cstds is not None and (len(u_cgw8s) > 1 or len(u_cstds) > 1) else ''
    str_circuit += f"{u_cgw8[:2]}" if u_cgw8s is not None and (len(u_cgw8s) > 1 or len(u_cstds) > 1) else ''
    str_circuit += f"_bond{u_bond}" if u_bonds is not None and len(u_bonds) > 1 else ''
    return f"L{str_L}{delim}J{str_J}{delim}x{str_x}{delim}r{str_rL}{delim}u[{str_circuit}]"

def get_config_filename(d: dict, dl: dict, p: dict):
    unique_path = get_unique_config_string(d,dl, '_')
    return f"{p['config_dir']}/{p['output_stem']}_{unique_path}.cfg"


def get_output_filepath(d: dict, dl: dict, p: dict):
    unique_path = get_unique_config_string(d,dl, '/')
    # config_filename += f"_u[{str_circuit}].cfg"
    return f"{p['output_dir']}/{unique_path}/{p['output_stem']}.h5"

@cache
def get_config_linenumber(config,key):
    with open(config, "r") as file:
        for idx, line in enumerate(file):
            if line.startswith(key):
                return idx+1
def get_config_value(config, key):
    lineno = get_config_linenumber(config,key)
    if lineno is None:
        raise LookupError(f"Could not find [{key=}] in  [{config=}]")
    return linecache.getline(config, lineno=lineno).split()[2]


def get_max_time(d: dict, dl: dict, p: dict):
    L  = float(d['model::model_size'])
    w1 = float(d['model::lbit::J1_wdth'])  # The width of distribution for on-site field.
    w2 = float(d['model::lbit::J2_wdth'])  # The width of distribution for pairwise interactions. The distribution is either U(J2_mean-w,J2_mean+w) or N(J2_mean,w)
    w3 = float(d['model::lbit::J3_wdth'])  # The width of distribution for three-body interactions.
    x  = float(d['model::lbit::xi_Jcls'])
    r = int(d['model::lbit::J2_span'])

    if r == -1:
        r = L

    tmax1 = 1.0 / w1 if w1 != 0 else 0.0

    r2max = np.min([r, L])  # Number of sites from the center site to the edge site, max(|i-j|)/2
    Jmin2 = np.exp(-r2max / x) * (w2 if w2 != 0 else 1.0) * np.sqrt(2 / np.pi)  # Order of magnitude of the smallest 2-body terms (furthest neighbor, up to L/2)
    tmax2 = 1.0 / Jmin2  # (0.5 to improve fits) Time that it takes for the most remote site to interact with the middle
    tmax3 = 1.0 / w3 if w3 != 0 else 0.0
    tmax = np.max([tmax1, tmax2, tmax3])
    tmax = 10 ** np.ceil(np.log10(tmax))
    print(f"{L=} {tmax=}")
    if L == 12 and tmax != 1e6:
        raise AssertionError(f"{L=} is supposed to have tmax=1e6. Got {tmax=:.1e}")
    if L == 16 and tmax != 1e8:
        raise AssertionError(f"{L=} is supposed to have tmax=1e8. Got {tmax=:.1e}")
    if L == 20 and tmax != 1e9:
        raise AssertionError(f"{L=} is supposed to have tmax=1e9. Got {tmax=:.1e}")
    if L == 24 and tmax != 1e11:
        raise AssertionError(f"{L=} is supposed to have tmax=1e11. Got {tmax=:.1e}")
    if L == 28 and tmax != 1e13:
        raise AssertionError(f"{L=} is supposed to have tmax=1e13. Got {tmax=:.1e}")
    if L == 32 and tmax != 1e14:
        raise AssertionError(f"{L=} is supposed to have tmax=1e14. Got {tmax=:.1e}")
    return '{:.1e}'.format(tmax)

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
def get_entry_indices(dset, expected_num, kind = 256):
    return np.where(dset['event'] == kind)[0][:expected_num]


    # return np.where(dset['event'] == kind)[0]
def get_num_entries(dsets, expected_num, kind = 256): # 256 is the enum value for ITER_STATE
    indices = {}
    skip = ['/common/finished_all', '/fLBIT/state_real/initial_pattern','/fLBIT/state_real/number_probabilities']
    for dset in dsets:
        # Check if there is an 'iter' field
        if dset is not None:
            if dset.name in skip:
                indices[dset.name] = None
            elif dset.dtype.fields is not None and 'event' in dset.dtype.fields.keys():
                # idx = get_entry_indices(dset)
                # It can happen that some entries are repeated (due to the old way of resuming simulations)
                indices[dset.name] = len(get_entry_indices(dset,expected_num, kind))
            else:
                indices[dset.name] = np.shape(dset)[-1]
    return indices
def get_h5_status(filename, batch):
    if os.path.isfile(filename):
        config = f'{batch["config_file"]}'
        expected_dset_paths = [
            '/common/finished_all',
            '/fLBIT/state_real/measurements',
            '/fLBIT/state_real/status',
            '/fLBIT/state_real/mem_usage',
            '/fLBIT/state_real/number_probabilities',
        ]
        optional_dset_paths = [
            '/fLBIT/state_real/initial_pattern'
        ]
        optional_link_attrs = {
            'initial_state': '/fLBIT/state_real',
            'time_scale': '/fLBIT/state_real',
        }

        try:
            with h5py.File(filename, 'r') as h5file:
                expected_dsets = [h5file.get(path) for path in expected_dset_paths]
                optional_dsets = [h5file.get(path) for path in optional_dset_paths]
                optional_attrs = [h5file.get(path).attrs.get(attr) for attr,path in optional_link_attrs.items() if path in h5file]
                missing_dsets = [path for dset,path in zip(expected_dsets,expected_dset_paths) if dset is None]
                if len(missing_dsets) > 0:
                    return f"FAILED|missing datasets:{missing_dsets}"

                length = expected_dsets[1]['length'][0]

                if optional_dsets[0] is not None:
                    evn_neel = np.resize([0,1], int(length))
                    odd_neel = np.resize([1,0], int(length))
                    has_neel_init_pattern = np.all(optional_dsets[0][()] == evn_neel) or np.all(optional_dsets[0][()] == odd_neel)
                    should_be_neel = 'neel' in filename or 'lbit93-precision' in filename or '-lin' in filename
                    if should_be_neel and not has_neel_init_pattern:
                        return f"FAILED|initial state is not neel"
                    if not should_be_neel and has_neel_init_pattern:
                        return f"FAILED|initial state is neel:{filename}"
                if optional_attrs[0] is not None:
                    evn_neel = 'b'+''.join(np.resize(['0','1'], int(length)))
                    odd_neel = 'b'+''.join(np.resize(['1','0'], int(length)))
                    has_neel_init_pattern = np.all(optional_attrs[0][()] == evn_neel) or np.all(optional_attrs[0][()] == odd_neel)
                    should_be_neel = 'neel' in filename or 'lbit93-precision' in filename or '-lin' in filename
                    if should_be_neel and not has_neel_init_pattern:
                        return f"FAILED|initial state is not neel"
                    if not should_be_neel and has_neel_init_pattern:
                        return f"FAILED|initial state is neel:{filename}"
                should_have_linspace = '-lin' in filename or 'singlet' in filename
                if optional_attrs[1] is not None:
                    expected_timescale = "LINSPACED" if should_have_linspace else "LOGSPACED"
                    if optional_attrs[1] != expected_timescale:
                        return f"FAILED|unexpected time scale {optional_attrs[1]}. Expected {expected_timescale}"

                time_steps = batch['time_steps']
                num_entries = get_num_entries(expected_dsets, expected_num=time_steps, kind=256)
                num_are_equal = all_equal([x for x in num_entries.values() if x is not None])
                if not num_are_equal:
                    return  f"FAILED|unequal iters: {num_are_equal}: expected: {time_steps})"

                for dset, num in num_entries.items():
                    if num is not None and num < time_steps:
                        return f"TIMEOUT|found too few iters: {dset}: found {num}, expected: {time_steps})"
                    if num is not None and num > time_steps:
                        return f"FAILED|found too many iters: {dset}: found {num}, expected: {time_steps})"

                t_idx = get_entry_indices(expected_dsets[1], expected_num=time_steps)
                if not arrays_equal(t_idx, np.arange(time_steps)):
                    return f"FAILED|iter entries are not a contiguous range: {dset}: found {num}, expected: {t_idx})"

                times = expected_dsets[1]['physical_time'][t_idx].astype(float)
                expected_tmin = abs(complex(float(get_config_value(config, "flbit::time_start_real")),
                                            float(get_config_value(config, "flbit::time_start_imag"))))
                expected_tmax = abs(complex(float(get_config_value(config, "flbit::time_final_real")),
                                            float(get_config_value(config, "flbit::time_final_imag"))))

                expected_times = np.linspace(expected_tmin,expected_tmax, time_steps,
                                             endpoint=True) if should_have_linspace else np.logspace(
                                 np.log10(expected_tmin),np.log10(expected_tmax), time_steps, endpoint=True)
                if not np.allclose(times, expected_times, atol=1):
                    return f"FAILED|mismatching times"

                found_tmax = times[-1]
                has_expected_tmax  = found_tmax == expected_tmax
                if not has_expected_tmax:
                    return f"FAILED|found {found_tmax=:.1e}!={expected_tmax=:.1e}"
                has_finished_all   = expected_dsets[0][()]
                if has_finished_all:
                    return f"FINISHED"
                else:
                    return f"TIMEOUT|{time_steps=}"
        except Exception as e:
            print(traceback.format_exc())
            return f"FAILED|{e}"
    else:
        return "MISSING"


def write_batch_status_single_thread(batch):
    parser = argparse.ArgumentParser(description='SLURM batch generator')
    parser.add_argument('--update-status', action='store_true', help='Update status files', default=False)
    args = parser.parse_args()

    config_file = f'{batch["config_file"]}'
    status_file = f'{batch["status_dir"]}/{Path(config_file).stem}.status'
    seed_status = copy(batch.get('seed_status')) # Old one
    batch['seed_status'] = []
    if platform.node() != "neumann" and args.update_status:
        raise AssertionError("--update-status is only valid on neumann")

    if platform.node() == "neumann" and args.update_status:
        Path(status_file).parent.mkdir(parents=True, exist_ok=True)
        output_base = '/mnt/WDB-AN1500/mbl_transition'
        output_path = f'{output_base}/{batch["projectname"]}/{batch["output_path"]}'
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
                            h5status = seed_status[sidx] #Short circuit
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

def write_batch_status(batch_filename):
    parser = argparse.ArgumentParser(description='SLURM batch generator')
    parser.add_argument('--update-status', action='store_true', help='Update status files', default=False)
    args = parser.parse_args()
    with open(batch_filename, 'r') as fp:
        batch = json.load(fp)
        config_file = f'{batch["config_file"]}'
        status_file = f'{batch["status_dir"]}/{Path(config_file).stem}.status'
        seed_status = copy(batch.get('seed_status')) # Old one
        batch['seed_status'] = []
        if platform.node() != "neumann" and args.update_status:
            raise AssertionError("--update-status is only valid on neumann")

        if platform.node() == "neumann" and args.update_status:
            Path(status_file).parent.mkdir(parents=True, exist_ok=True)
            output_base = '/mnt/WDB-AN1500/mbl_transition'
            output_path = f'{output_base}/{batch["projectname"]}/{batch["output_path"]}'
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
                                h5status = seed_status[sidx] #Short circuit
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



def update_batch_status(config_paths):
    batch_filenames = sorted(Path(config_paths['config_dir']).rglob('*.json'))
    with Pool() as pool:
        batch_jsons = pool.map(write_batch_status, batch_filenames)
        for batch_filename, batch_json in zip(batch_filenames, batch_jsons):
            with open(batch_filename, 'w') as fp:
                json.dump(batch_json, fp, sort_keys=True, indent=4)




