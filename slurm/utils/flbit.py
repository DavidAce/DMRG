import numpy as np

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
    # config_filename += f"_u[{str_circuit}].cfg"
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
    Jmin2 = np.exp(-(r2max - 1) / x) * w2 * np.sqrt(2 / np.pi)  # Order of magnitude of the smallest 2-body terms (furthest neighbor, up to L/2)
    tmax2 = 1.0 / Jmin2  # (0.5 to improve fits) Time that it takes for the most remote site to interact with the middle
    tmax3 = 1.0 / w3
    tmax = np.max([tmax1, tmax2, tmax3])
    return '{:.1e}'.format(10 ** np.ceil(np.log10(tmax)))


