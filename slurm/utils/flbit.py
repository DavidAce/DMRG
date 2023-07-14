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

def get_max_time(d: dict, dl: dict, p: dict):
    L  = float(d['model::model_size'])
    w1 = float(d['model::lbit::J1_wdth'])  # The width of distribution for on-site field.
    w2 = float(d['model::lbit::J2_wdth'])  # The width of distribution for pairwise interactions. The distribution is either U(J2_mean-w,J2_mean+w) or N(J2_mean,w)
    w3 = float(d['model::lbit::J3_wdth'])  # The width of distribution for three-body interactions.
    x  = float(d['model::lbit::xi_Jcls'])
    r = int(d['model::lbit::J2_span'])

    if r == -1:
        r = L

    tmax1 = 1.0 / w1

    r2max = np.min([r, L])  # Number of sites from the center site to the edge site, max(|i-j|)/2
    Jmin2 = np.exp(-r2max / x) * w2 * np.sqrt(2 / np.pi)  # Order of magnitude of the smallest 2-body terms (furthest neighbor, up to L/2)
    tmax2 = 1.0 / Jmin2  # (0.5 to improve fits) Time that it takes for the most remote site to interact with the middle
    tmax3 = 1.0 / w3
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
                time_steps=len(expected_dsets[1])
                has_finished_all   = expected_dsets[0][()]
                r2max=float(expected_dsets[1]['length'][0])
                Jmin2 = np.exp(-r2max / 1) * 1 * np.sqrt(2 / np.pi)  # Order of magnitude of the smallest 2-body terms (furthest neighbor, up to L/2)
                tmax2 = 1.0 / Jmin2
                tmax = 10 ** np.ceil(np.log10(tmax2))
                found_tmax = expected_dsets[1]['physical_time'][-1].astype(float)
                has_expected_tmax  = found_tmax == tmax
                has_expected_iters = time_steps >= batch['time_steps']
                has_exceeded_iters = time_steps > batch['time_steps']
                if has_expected_iters and has_finished_all:
                    if not has_expected_tmax:
                        return f"FAILED|found {found_tmax=:.1e}!={tmax=:.1e}"
                    if has_exceeded_iters:
                        return f"FINISHED|{time_steps=}"
                    else:
                        return f"FINISHED"
                if has_expected_iters and not has_finished_all:
                    return f"TIMEOUT|{time_steps=}"
                if not has_expected_iters and has_finished_all:
                    return f"FAILED|{time_steps=}"
                if not has_expected_iters and not has_finished_all:
                    return f"TIMEOUT|{time_steps=}"
                return "FAILED|unknown reason"
        except Exception as e:
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
            rclone.copy(f'neumann:{batch["output_prfx"]}/{batch["projectname"]}/{status_file}',
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
                rclone.copy(f'neumann:{batch["output_prfx"]}/{batch["projectname"]}/{status_file}',
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




