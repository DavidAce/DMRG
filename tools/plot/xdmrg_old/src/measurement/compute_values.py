import warnings

from ..io.h5ops import *


# from merge_settings import *

# LEGEND:
# A    single-site left-canonical MPS. Equivalent to "lambda * gamma"  in Vidal notation
# B    single-site right-canonical MPS. Equivalent to "gamma * lambda" in Vidal notation
# C    center-site bond matrix, equivalent to a center "lambda" in Vidal notation
# S    single-site bond matrix, any one on the chain (not the center one necessarily)
# H    single-site Hamiltonian MPO
# L    left environment. Each site has a corresponding left environment
# R    right environmnent. Each site has a corresponding left environment
# L2   left environment with 2 MPO-legs for variance calculations.
# R2   right environment with 2 MPO-legs for variacne calculations
#
#

def get_center_pos(hdf5_sim):
    key_sorted = sorted(hdf5_sim['chain/MPS'].keys(), key=natural_keys)
    A_keys = [s for s in key_sorted if 'A' in s]
    return len(A_keys) - 1


def get_middle_pos(hdf5_sim):
    chain_length = hdf5_sim['measurements']['2site']['length'][0]
    return int((chain_length - 1) / 2)


def mps_norm(hdf5_sim):
    # Start by getting the position of the "center" bond between A's and B's.
    center_position = get_center_pos(hdf5_sim)
    # Load objects around the center.
    A = load_component_cplx(hdf5_sim, path='chain/MPS', keypattern='A_' + str(center_position))[0]
    B = load_component_cplx(hdf5_sim, path='chain/MPS', keypattern='B_' + str(center_position + 1))[0]
    C = load_component_cplx(hdf5_sim, path='chain/MPS', keypattern='C')[0]
    # Contract them
    temp = np.tensordot(A, A.conjugate(), axes=([0, 1], [0, 1]))
    tempC = np.tensordot(temp, np.diag(C), axes=([0], [0]))
    temp = np.tensordot(tempC, np.diag(C), axes=([0], [0]))
    tempB = np.tensordot(temp, B, axes=([0], [1]))
    temp = np.tensordot(tempB, B.conjugate(), axes=([0, 1], [1, 0]))
    return np.real(np.trace(temp))


def mps_norm_full(hdf5_sim):
    A_list = load_component_cplx(hdf5_sim, path='chain/MPS', keypattern='A')
    B_list = load_component_cplx(hdf5_sim, path='chain/MPS', keypattern='B')
    C_list = load_component_cplx(hdf5_sim, path='chain/MPS', keypattern='C')
    # Compute norm of chain
    temp = np.tensordot(A_list[0], A_list[0].conjugate(), axes=([0, 1], [0, 1]))
    # np.shape(A_list)
    for A in A_list[1:]:
        tempA = np.tensordot(temp, A, axes=([0], [1]))
        temp = np.tensordot(tempA, A.conjugate(), axes=([0, 1], [1, 0]))
    tempC = np.tensordot(temp, np.diag(C_list[0]), axes=([0], [0]))
    temp = np.tensordot(tempC, np.diag(C_list[0]), axes=([0], [0]))

    for B in B_list:
        tempB = np.tensordot(temp, B, axes=([0], [1]))
        temp = np.tensordot(tempB, B.conjugate(), axes=([0, 1], [1, 0]))
    return np.real(np.trace(temp))


def mps_pauli_operator(hdf5_sim, op):
    A_list = load_component_cplx(hdf5_sim, path='chain/MPS', keypattern='A')
    B_list = load_component_cplx(hdf5_sim, path='chain/MPS', keypattern='B')
    C_list = load_component_cplx(hdf5_sim, path='chain/MPS', keypattern='C')
    # Compute norm of chain

    tempOp = np.tensordot(A_list[0], op, axes=([0], [0]))
    temp = np.tensordot(tempOp, A_list[0].conjugate(), axes=([2, 0], [0, 1]))

    for A in A_list[1:]:
        tempA = np.tensordot(temp, A, axes=([0], [1]))
        tempOp = np.tensordot(tempA, op, axes=([1], [0]))
        temp = np.tensordot(tempOp, A.conjugate(), axes=([0, 2], [1, 0]))
    tempC = np.tensordot(temp, np.diag(C_list[0]), axes=([0], [0]))
    temp = np.tensordot(tempC, np.diag(C_list[0]), axes=([0], [0]))

    for B in B_list:
        tempB = np.tensordot(temp, B, axes=([0], [1]))
        tempOp = np.tensordot(tempB, op, axes=([1], [0]))
        temp = np.tensordot(tempOp, B.conjugate(), axes=([0, 2], [1, 0]))
    # print(temp)
    return np.real(np.trace(temp))


def mps_energy(hdf5_sim):
    # Start by getting the position of the "center" bond between A's and B's.
    center_position = get_center_pos(hdf5_sim)
    # Load objects around the center.
    A = load_component_cplx(hdf5_sim, path='chain/MPS', keypattern='A_' + str(center_position))[0]
    B = load_component_cplx(hdf5_sim, path='chain/MPS', keypattern='B_' + str(center_position + 1))[0]
    C = load_component_cplx(hdf5_sim, path='chain/MPS', keypattern='C')[0]
    HA = load_component_cplx(hdf5_sim, path='chain/MPO', keypattern='H_' + str(center_position))[0]
    HB = load_component_cplx(hdf5_sim, path='chain/MPO', keypattern='H_' + str(center_position + 1))[0]
    L = load_component_cplx(hdf5_sim, path='chain/ENV', keypattern='L')[0]
    R = load_component_cplx(hdf5_sim, path='chain/ENV', keypattern='R')[0]

    # Contract them
    L = np.tensordot(L, A, axes=([0], [1]))
    L = np.tensordot(L, A.conjugate(), axes=([0], [1]))
    L = np.tensordot(L, HA, axes=([0, 1, 3], [0, 2, 3]))
    L = np.tensordot(L, np.diag(C), axes=([0], [0]))
    L = np.tensordot(L, np.diag(C), axes=([0], [0]))
    L = np.transpose(L, (2, 1, 0))
    L = np.tensordot(L, B, axes=([0], [1]))
    L = np.tensordot(L, B.conjugate(), axes=([0], [1]))
    L = np.tensordot(L, HB, axes=([0, 1, 3], [0, 2, 3]))
    L = np.tensordot(L, R, axes=([0, 1, 2], [0, 1, 2]))
    return np.real(L)


def mps_energy_full(hdf5_sim):
    # Compute energies
    chain_length = hdf5_sim.get('xDMRG')['chain_length'][-1]
    A_list = load_component_cplx(hdf5_sim, path='chain/MPS', keypattern='A')
    B_list = load_component_cplx(hdf5_sim, path='chain/MPS', keypattern='B')
    C_list = load_component_cplx(hdf5_sim, path='chain/MPS', keypattern='C')
    H_list = load_component_cplx(hdf5_sim, path='chain/MPO')

    L = "CONTRACT L0 A0 A0' "
    R = "CONTRACT L[-1] B[-1] and B[-1]'"

    # L      = load_component_cplx(hdf5_sim,path='chain/ENV',filter_name='L_0')
    # R      = load_component_cplx(hdf5_sim,path='chain/ENV',filter_name='R_' + str(chain_length-1))

    # Take the first elements of the edge lists
    # L    = L[0]
    # R    = R[0]
    for A, H in zip(A_list, H_list):
        L = np.tensordot(L, A, axes=([0], [1]))
        L = np.tensordot(L, A.conjugate(), axes=([0], [1]))
        L = np.tensordot(L, H, axes=([0, 1, 3], [0, 2, 3]))
    L = np.tensordot(L, np.diag(C_list[0]), axes=([0], [0]))
    L = np.tensordot(L, np.diag(C_list[0]), axes=([0], [0]))
    L = np.transpose(L, (2, 1, 0))
    for B, H in zip(B_list, H_list[len(A_list):]):
        L = np.tensordot(L, B, axes=([0], [1]))
        L = np.tensordot(L, B.conjugate(), axes=([0], [1]))
        L = np.tensordot(L, H, axes=([0, 1, 3], [0, 2, 3]))
    L = np.tensordot(L, R, axes=([0, 1, 2], [0, 1, 2]))
    return np.real(L)


def mps_variance(hdf5_sim):
    # Start by getting the position of the "center" bond between A's and B's.
    center_position = get_center_pos(hdf5_sim)
    # Load objects around the center.
    A = load_component_cplx(hdf5_sim, path='chain/MPS', keypattern='A_' + str(center_position))[0]
    B = load_component_cplx(hdf5_sim, path='chain/MPS', keypattern='B_' + str(center_position + 1))[0]
    C = load_component_cplx(hdf5_sim, path='chain/MPS', keypattern='C')[0]
    HA = load_component_cplx(hdf5_sim, path='chain/MPO', keypattern='H_' + str(center_position))[0]
    HB = load_component_cplx(hdf5_sim, path='chain/MPO', keypattern='H_' + str(center_position + 1))[0]
    L2 = load_component_cplx(hdf5_sim, path='chain/ENV2', keypattern='L')[0]
    R2 = load_component_cplx(hdf5_sim, path='chain/ENV2', keypattern='R')[0]

    # Compute variances
    L2 = np.tensordot(L2, A, axes=([0], [1]))
    L2 = np.tensordot(L2, A.conjugate(), axes=([0], [1]))
    L2 = np.tensordot(L2, HA, axes=([0, 2], [0, 2]))
    L2 = np.tensordot(L2, HA, axes=([0, 5, 2], [0, 2, 3]))
    L2 = np.tensordot(L2, np.diag(C), axes=([0], [0]))
    L2 = np.tensordot(L2, np.diag(C), axes=([0], [0]))
    L2 = np.transpose(L2, (2, 3, 0, 1))
    L2 = np.tensordot(L2, B, axes=([0], [1]))
    L2 = np.tensordot(L2, B.conjugate(), axes=([0], [1]))
    L2 = np.tensordot(L2, HB, axes=([0, 2], [0, 2]))
    L2 = np.tensordot(L2, HB, axes=([0, 5, 2], [0, 2, 3]))
    L2 = np.tensordot(L2, R2, axes=([0, 1, 2, 3], [0, 1, 2, 3]))
    return np.abs(L2 - mps_energy(hdf5_sim) ** 2)


def mps_variance_full(hdf5_sim):
    # Compute variances
    chain_length = hdf5_sim.get('xDMRG')['chain_length'][-1]
    A_list = load_component_cplx(hdf5_sim, path='chain/MPS', keypattern='A')
    B_list = load_component_cplx(hdf5_sim, path='chain/MPS', keypattern='B')
    C_list = load_component_cplx(hdf5_sim, path='chain/MPS', keypattern='C')
    H_list = load_component_cplx(hdf5_sim, path='chain/MPO')
    # L2     = load_component_cplx(hdf5_sim,path='chain/ENV2',filter_name='L_0')
    # R2     = load_component_cplx(hdf5_sim,path='chain/ENV2',filter_name='R_' + str(chain_length-1))
    # Take the first elements of the edge lists

    L2 = "CONTRACT L0 A0 A0' "
    R2 = "CONTRACT L[-1] L[-1] B[-1] and B[-1]'"

    # L2    = L2[0]
    # R2    = R2[0]
    for A, H in zip(A_list, H_list):
        L2 = np.tensordot(L2, A, axes=([0], [1]))
        L2 = np.tensordot(L2, A.conjugate(), axes=([0], [1]))
        L2 = np.tensordot(L2, H, axes=([0, 2], [0, 2]))
        L2 = np.tensordot(L2, H, axes=([0, 5, 2], [0, 2, 3]))

    L2 = np.tensordot(L2, np.diag(C_list[0]), axes=([0], [0]))
    L2 = np.tensordot(L2, np.diag(C_list[0]), axes=([0], [0]))
    L2 = np.transpose(L2, (2, 3, 0, 1))

    for B, H in zip(B_list, H_list[len(A_list):]):
        L2 = np.tensordot(L2, B, axes=([0], [1]))
        L2 = np.tensordot(L2, B.conjugate(), axes=([0], [1]))
        L2 = np.tensordot(L2, H, axes=([0, 2], [0, 2]))
        L2 = np.tensordot(L2, H, axes=([0, 5, 2], [0, 2, 3]))
    L2 = np.tensordot(L2, R2, axes=([0, 1, 2, 3], [0, 1, 2, 3]))
    return np.abs(L2 - mps_energy_full(hdf5_sim) ** 2)


def mps_keep_state(hdf5_sim):
    variance = hdf5_sim['measurements/full/energy_variance_per_site_mpo'][0]
    energy = hdf5_sim['measurements/full/energy_per_site_mpo'][0]
    energy_min = hdf5_sim['measurements'].get('simulation_progress')['energy_min'][-1]
    energy_max = hdf5_sim['measurements'].get('simulation_progress')['energy_max'][-1]
    energy_dens = (energy - energy_min) / (energy_max - energy_min)
    return np.abs(energy_dens - energy_window) < 0.5 and variance < variance_threshold_upper


def mps_entanglement_entropy(hdf5_sim, mode='middle', variance_threshold=0):
    # The lambda matrices for a chain of length N are named as follows:
    #
    #   L_0, L_1, ... L_10, L_11 , L_11_12_C, L_12, L_13, ..., L_N
    #
    # The L_11_12_C matrix is special; it is the one between the A and B MPS's.
    # However, it is not necessarily in the middle of the chain! The middle of the chain
    # is simply position (N-1)/2. If the middle happened to be 11 in this example, both
    # L_11 and L_11_12_C would match the string filter. Therefore, always pick the last of
    # the matched elements.
    if (mode == 'middle'):
        if mps_keep_state(hdf5_sim):
            return hdf5_sim['measurements/2site/entanglement_entropy'][0]
        else:
            return np.nan
    elif (mode == 'all'):
        if mps_keep_state(hdf5_sim):
            return np.asarray(hdf5_sim['measurements/full/entanglement_entropies'])
        else:
            return np.nan * np.asarray(hdf5_sim['measurements/full/entanglement_entropies'])
    elif (isinstance(mode, int)):
        position = mode
        try:
            if mps_keep_state(hdf5_sim):
                return hdf5_sim['measurements/full/entanglement_entropies'][position]
            else:
                return np.nan
        except:
            warnings.warn("Error when reading entanglement entropy at position: " + str(position))
            print(hdf5_sim['measurements/full/entanglement_entropies'])
            exit(1)
    else:
        print('Error: argument mode= is not \"all\", \"middle\" or <integer>.')
        exit(1)


def mps_bond_dimension(hdf5_sim, mode='middle'):
    # The lambda matrices for a chain of length N are named as follows:
    #
    #   L_0, L_1, ... L_10, L_11 , L_11_12_C, L_12, L_13, ..., L_N
    #
    # The L_11_12_C matrix is special; it is the one between the A and B MPS's.
    # However, it is not necessarily in the middle of the chain! The middle of the chain
    # is simply position (N-1)/2. If the middle happened to be 11 in this example, both
    # L_11 and L_11_12_C would match the string filter. Therefore, always pick the last of
    # the matched elements.
    if (mode == 'middle'):
        if mps_keep_state(hdf5_sim):
            return hdf5_sim['measurements']['2site']['bond_dimension'][0]
        else:
            return np.nan
    elif (mode == 'all'):
        if mps_keep_state(hdf5_sim):
            return np.asarray(hdf5_sim['measurements']['full']['bond_dimensions'])
        else:
            return np.nan * np.asarray(hdf5_sim['measurements']['full']['bond_dimensions'])
    elif (isinstance(mode, int)):
        position = mode
        if mps_keep_state(hdf5_sim):
            return hdf5_sim['measurements']['full']['bond_dimensions'][position]
        else:
            return np.nan
    else:
        print('Error: argument mode= is not \"all\", \"middle\" or <integer>.')
        exit(1)
