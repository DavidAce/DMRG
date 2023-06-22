from itertools import product
from pathlib import Path
import json
import yaml
import os
import numpy as np

def get_config_product(dictoflists : dict, p: dict):
    c = []
    for vals in product(*dictoflists.values()):
        d = dict(zip(dictoflists.keys(),vals))
        for key,val in d.items():
            if callable(val):
                d[key] = val(d, dictoflists, p)
        c.append(d)
    return c


def replace_value(line,pos,val):
    old_val = line.split()[pos]
    index_start = line.find(old_val)
    index_end   = index_start + len(old_val)
    len_diff = len(old_val) - len(val)
    return line[:index_start] + val + ' '*len_diff + line[index_end:]

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



def write_seed_files(seed_setup, config_paths):
    for idx, config_filepath in enumerate(sorted(Path(config_paths['config_dir']).glob('*.cfg'))):
        seed_filename = '{}/{}.json'.format(config_paths['seed_dir'], config_filepath.stem)
        Path(seed_filename).parent.mkdir(parents=True, exist_ok=True)
        if os.path.isfile(seed_filename):
            with open(seed_filename, 'r') as fp:
                seedjson = json.load(fp)
                if seedjson['config_index'] != idx:
                    raise ValueError(f"seed file {seed_filename} has {seedjson['config_index']=} "
                                     f"which does not match the current index {idx}")
        else:
            seedjson = {
                "config_index": idx,
                "config_file" : str(config_filepath),
                'projectname': seed_setup['projectname'],
                'seed_extent': [],
                'seed_offset': [],
            }

        # Now we need to now which seed set corresponds to this config file
        seed_keys = []
        for key in seed_setup['seeds'].keys():
            if all(x in seed_filename for x in key.split('|')):
                seed_keys.append(key)
        if len(seed_keys) == 0:
            print(f'No seed keys found for config: {config_filepath.name}')
            continue
        elif len(seed_keys) > 1:
            raise AssertionError(f'Found multiple seed keys matching config: {config_filepath.name}\n'
                                 f'\n{seed_keys=}')
        seeds = seed_setup['seeds'][seed_keys[0]]
        for offset, extent in zip(seeds['offset'], seeds['extent']):

            extent_size = len(seedjson['seed_extent'])
            offset_size = len(seedjson['seed_offset'])
            if offset_size != extent_size:
                raise ValueError(
                    f"offset:{offset_size} and extent:{extent_size} are not equal lengths")

            # Check if this offset is already in jobdict for this config
            offset_index = seedjson['seed_offset'].index(offset) if offset in seedjson['seed_offset'] else -1
            if offset_index < 0:
                seedjson['seed_extent'].append(extent)
                seedjson['seed_offset'].append(offset)
            elif seedjson['seed_extent'][offset_index] != extent:
                raise ValueError(f"{seed_filename} has {offset=} at index {offset_index} with "
                                 f"extent {seedjson['seed_extent'][offset_index]}.\n"
                                 f"The new extent {extent} is incompatible")

        seedjson['seed_counts'] = int(np.sum(seedjson['seed_extent']))
        print(json.dumps(seedjson, sort_keys=True, indent=4))

        with open(seed_filename, 'w') as fp:
            json.dump(seedjson, fp, sort_keys=True, indent=4)
