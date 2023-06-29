from itertools import product
from pathlib import Path
import json
import yaml
import os
import numpy as np
import linecache
def get_config_product(config_ranges : dict, p: dict):
    c = []
    for vals in product(*config_ranges.values()):
        d = dict(zip(config_ranges.keys(),vals))
        for key,val in d.items():
            if callable(val):
                d[key] = val(d, config_ranges, p)
        c.append(d)
    return c


def replace_value(line,pos,val):
    old_val = line.split()[pos]
    index_start = line.find(old_val)
    index_end   = index_start + len(old_val)
    len_diff = len(old_val) - len(val)
    return line[:index_start] + val + ' '*len_diff + line[index_end:]





def get_batch_status(batch_setup, configs, config_paths):
    for config in configs:
        config_filepath = Path(config['filename'])

        # Now we need to now which seed set corresponds to this config file
        batch_keys = []
        for key in batch_setup['batch'].keys():
            if all(x in config_filepath.name for x in key.split('|')):
                batch_keys.append(key)
        if len(batch_keys) == 0:
            print(f'No batch keys found for config: {config_filepath.name}')
            continue
        elif len(batch_keys) > 1:
            raise AssertionError(f'Found multiple seed keys matching config: {config_filepath.name}\n'
                                 f'\n{batch_keys=}')
        batch = batch_setup['batch'][batch_keys[0]]
        # batch['seed_status'] = []
        # Now we add the status
        status_file = f'{config_paths["status_dir"]}/{config_filepath.stem}.status'
        status_count = 0
        for offset, extent in zip(batch['seed_offset'], batch['seed_extent']):
            extent_size = len(batch['seed_extent'])
            offset_size = len(batch['seed_offset'])
            if offset_size != extent_size:
                raise ValueError(
                    f"offset:{offset_size} and extent:{extent_size} are not equal lengths")

            is_finished = True
            for idx,seed in enumerate(range(offset, offset + extent)):
                sfline = linecache.getline(status_file, idx+status_count+1).rstrip()
                print(idx,seed,sfline)
                sfseed, sfstatus = sfline.split('|',maxsplit=1)
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
        batch_setup['batch'][batch_keys[0]] = batch
    return batch_setup

def write_config_file(config, config_template, config_filename):
    Path(config_filename).parent.mkdir(parents=True, exist_ok=True)
    print(config_filename)
    print(yaml.dump(config, allow_unicode=True))

    with open(config_template, 'r') as template:
        with open(config_filename, 'w') as file:
            for line in template:
                for var,val in config.items():
                    if line.find(var) >= 0:
                        line = replace_value(line, 2, val)
                file.write(line)



def write_batch_files(batch_setup, configs, config_paths):
    for config in configs:
        config_filepath = Path(config['filename'])
        batch_filename = '{}/{}.json'.format(config_paths['config_dir'], config_filepath.stem)
        Path(batch_filename).parent.mkdir(parents=True, exist_ok=True)
        if os.path.isfile(batch_filename):
            with open(batch_filename, 'r') as fp:
                batchjson = json.load(fp)
        else:
            batchjson = {
                'config_file' : str(config_filepath),
                'output_path' : str(Path(config["storage::output_filepath"]).parent),
                'output_stem' : config_paths['output_stem'],
                'projectname': batch_setup['projectname'],
                'seed_extent': [],
                'seed_offset': [],
                # 'seed_status': [],
            }

        # Now we need to now which seed set corresponds to this config file
        seed_keys = []
        for key in batch_setup['batch'].keys():
            if all(x in batch_filename for x in key.split('|')):
                seed_keys.append(key)
        if len(seed_keys) == 0:
            print(f'No seed keys found for config: {config_filepath.name}')
            continue
        elif len(seed_keys) > 1:
            raise AssertionError(f'Found multiple seed keys matching config: {config_filepath.name}\n'
                                 f'\n{seed_keys=}')
        batch = batch_setup['batch'][seed_keys[0]]
        batchjson['time_steps'] = batch['time_steps']

        for offset, extent in zip(batch['seed_offset'], batch['seed_extent']):
            extent_size = len(batchjson['seed_extent'])
            offset_size = len(batchjson['seed_offset'])
            # status_size = len(batchjson['seed_status'])
            if offset_size != extent_size:
                raise ValueError(
                    f"offset:{offset_size} and extent:{extent_size} are not equal lengths")
            # if status_size != extent_size:
            #     raise ValueError(
            #         f"status:{status_size} and extent:{extent_size} are not equal lengths")

            # Check if this offset is already in jobdict for this config
            offset_index = batchjson['seed_offset'].index(offset) if offset in batchjson['seed_offset'] else -1
            if offset_index < 0:
                batchjson['seed_extent'].append(extent)
                batchjson['seed_offset'].append(offset)
                # batchjson['seed_status'].append(status)
            elif batchjson['seed_extent'][offset_index] != extent:
                raise ValueError(f"{batch_filename} has {offset=} at index {offset_index} with "
                                 f"extent {batchjson['seed_extent'][offset_index]}.\n"
                                 f"The new extent {extent} is incompatible")

        batchjson['seed_counts'] = int(np.sum(batchjson['seed_extent']))
        # print(json.dumps(batchjson, sort_keys=True, indent=4))

        with open(batch_filename, 'w') as fp:
            json.dump(batchjson, fp, sort_keys=True, indent=4)



