from src.io.h5ops import *

tgt_file = "merged/averaged_merged.h5"

src_files = {
    # 'A': "/mnt/Barracuda/Projects/mbl_transition/data114-small-interaction-larger-systems/analysis/merged/averaged.h5",
    'A': "/mnt/Barracuda/Projects/mbl_transition/data116-l0.01-larger-systems/analysis/merged/averaged.h5",
    'B': "/mnt/Barracuda/Projects/mbl_transition/data115-large-interactions/analysis/merged/averaged.h5",
    'C': "/mnt/Barracuda/Projects/mbl_transition/data112-uptoeleven-resume/analysis/merged/averaged.h5",
}

src_tgt_map = {
    'A': {
        'L': {},
        'l': {
            'l_0': 'l_0.010',
            'l_0.1': 'l_0.100',
            'l_0.01': 'l_0.010',
            'l_0.010': 'l_0.010',
        },
        'J': {'J_0': 'J_1.000'},
        'h': {'h_0': 'h_1.000'}
    },
    'B': {
        'L': {},
        'l': {
            'l_0': 'l_0.010',
            'l_0.1': 'l_0.100',
            'l_0.01': 'l_0.010',
        },
        'J': {'J_0': 'J_1.000'},
        'h': {'h_0': 'h_1.000'}
    },
    'C': {
        'L': {},
        'l': {
            'l_0': 'l_0.000',
            'l_0.1': 'l_0.100',
            'l_0.01': 'l_0.010',
        },
        'J': {'J_0': 'J_1.000'},
        'h': {'h_0': 'h_1.000'}
    }

}


def get_mapped_path(src_key, L_key, l_key, J_key, h_key):
    L_map = src_tgt_map[src_key]['L'][L_key] if L_key in src_tgt_map[src_key]['L'] else L_key
    l_map = src_tgt_map[src_key]['l'][l_key] if l_key in src_tgt_map[src_key]['l'] else l_key
    J_map = src_tgt_map[src_key]['J'][J_key] if J_key in src_tgt_map[src_key]['J'] else J_key
    h_map = src_tgt_map[src_key]['h'][h_key] if h_key in src_tgt_map[src_key]['h'] else h_key
    print("Map:", "{}/{}/{}/{}".format(L_key, l_key, J_key, h_key), " -- > ", "{}/{}/{}/{}".format(L_map, l_map, J_map, h_map))
    return "{}/{}/{}/{}".format(L_map, l_map, J_map, h_map)


tgt_h5file = h5py.File(tgt_file, 'w')

for src_key, src_file in src_files.items():
    src_h5file = h5py.File(src_file, 'r')
    for L_key, L_node in src_h5file.items():
        for l_key, l_node in L_node.items():
            for J_key, J_node in l_node.items():
                for h_key, h_node in J_node.items():
                    src_path = "{}/{}/{}/{}".format(L_key, l_key, J_key, h_key)
                    tgt_path = get_mapped_path(src_key, L_key, l_key, J_key, h_key)
                    # tgt_node = tgt_h5file.require_group(tgt_path)
                    tgt_h5file.copy(source=h_node, dest=tgt_path)
