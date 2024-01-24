import argparse


def parse(algo, batches=None, basedir='/mnt/WDB-AN1500/mbl_transition'):
    if batches is not None:
        if not isinstance(batches, list):
            raise TypeError("batches must be a list")
    flbit_batches = [
        'lbit25',
        'lbit26',
        'lbit27',
        'lbit28',
        'lbit29', # swap-less: Tried to reduce the number of swap gates; L:24, w2:0.50,0.25, x:0.5, 1.0, f:0.15,0.25, u:2,4,6
        'lbit30', # time-less: Skip time-evo when gates are ~1. L:24, w2:0.50,0.25, x:0.5, 1.0, f:0.15,0.25, u:2,4,6
        'lbit31', # Try higher f to promote lbit movement. L:24, w2:0.25, x:0.5, f:0.30, u:5
        'lbit32', # Check new numb.entr, time-evo and xi. L:24, w2:0.25, x:0.5, f:0.25, u:5
        'lbit33', # Push to L32 with new numb.entr, time-evo and xi.  L:32, w2:0.25, x:0.5, f:0.25, u:5
        'lbit34', # Compare u4 and u4: L:24, w2:0.25, x:0.5, f:0.25, u:4,5
        'lbit35', # unwave: Check if more precision can remove waviness (NO):  L:24, w2:0.25, x:0.5, f:0.25, u:4
        'lbit36', # dearest-neighbor: Update exp decay to exp(-(r-1)/xi. L:24, w2:0.25, x:0.5, f:0.25, u:4
        'lbit37', # unwave-bump-x: Check if f ~ 0.15 can remove waviness (NO). L:24, w2:0.25, x:0.5, f:0.15, u:4
        'lbit38', # unwave-bump-x: Check if x ~ 1.00 can remove waviness (YES). L:24, w2:0.25, x:1.0, f:0.25, u:4
        'lbit39', # serious-trial-L28: Check loglog window. Predict ~5 decades of S_N growth (YES!):  L:12-28, w2:0.25, x:0.8, f:0.25, u:4
        'lbit40', # Mini test to check progress on new storage form
        'lbit41', # Mini test to check progress on new storage form
        'lbit42', # Big run, 5000 realizations L = 8..24
        'lbit43', # Small test w-dependence, L = 16
        'lbit44',  # Small test of normally distributed random circuit L = [4...16]
        'lbit45',  # Small test of squared distributed random circuit L = [4...16]
        'lbit46',  # Small test of squared distributed random circuit L = [8...16]
        'lbit47',  # Small test of big bias on the J2 and J3 parameters L = [8...16]
        'lbit48',  # Small test of choked unitaries on the J2 and J3 parameters L = [8...16]
        'lbit49',  # Small test of choked unitaries on the J2 and J3 parameters L = [8...16]
        'lbit50',  # 5000 realizations of choked unitaries, u=[4,5], f=[0.25,0.35], L = [8...16]
        'lbit51',  # FOR POSTER! exp(-2dh) constriction in unitaries, u=[5], f=[0.45], L = [8...20]
        'lbit52', # Trial increasing u and L together... not very useful
        'lbit53', # 10000 realizations squared constriction, u8, f0.25, L=[8,12,16]
        'lbit54',  # 10000 realizations squared constriction, u8, f0.50, L=[8,12,16]
        'lbit55',  # 25000 realizations constricted, u8, f=[0.1-0.5], L=[8,12,16]
        'lbit56',  # 25000 realizations constricted, random init, u8, f=[0.1-1.0], L=[8,12,16]
        'lbit57',  # 25000 realizations constricted, random init, u8, f=[2.0-5.0], L=[16]
        'lbit58',  # 10000 realizations constricted, random init, u8, f=0.5, L=[16], varying w
        'lbit59', # 1000 realizations blocked (constricted), random init, u8, f=0.5, L=12, increased c = 2
        'lbit60', # 1000 realizations u8, f=0.5, L=12, testing weights ID/EX std=[0.1,1.0,4.0]
        'lbit61', # 1000 realizations u8, f=0.5, L=14, testing weights ID/EX std=[0.25, 0.5 ...2.0]. This one motivates the choice f=1.0, cw=EX, sigma_c=1.0, and wt=EX/ID (doesn't matter), sigma_t=1.0 in the unitary gates
        'lbit62', # 6000 realizations u8, f=1.0, L=14, testing xi=[0.8, 0.9 ... 2.0]. This shows that the SN(t=inf) values are kind of random, mostly decreasing with xi. Probably we need more realizations or lower svd threshold.
        'lbit63', # 2000 realizations u8, f=1.0, L=14, testing xi=[0.8, 0.9 ... 2.0] with lower svd threshold: Didn't help...
        'lbit64',  # 2000 realizations u8, f=1.0, L=14, testing J and w to get more SN increase.
        'lbit65', # 1000 realizations u8, f=1.0, L=[8,12,16,20], checking how well SN grows with the new parameter set.
        'lbit66',   # 10000 realizations u8, f=1.0, L=12, u=8,10,12,14,16, checking how well SN grows with u depth.
        'lbit67', # 2000 realizations u[8,16...64], f=1.0, L=16, checking how well SN grows with u depth.
        'lbit68', # 300 realizations u[8,16...20], f=1.0, L=[8,...32], checking how well localization length xi.
        'lbit69', # 300 realizations u[8,16...32], f=1.0, L=[8,...32] and mpo bond dim [64...256], checking if mpo bond dim affects xi
        'lbit70', # 1500 realizations u[8,16...64], f=1.0, L=[8,...32] and mpo bond dim [128], checking what u is needed to reach xi=1
        'lbit71', # 1500 realizations u[8,16...80], f=1.0, L=[8,...40] and mpo bond dim [16] (new faster method), checking what u is needed to reach xi=1
        'lbit72', # 100 realizations u[16] f=1.0, L=[12,...32] with dynamic max-time for each L, to check if it sets the correct maximum
        'lbit73', # 3000 realizations u[16] f=1.0, L=[12,...32] with dynamic max-time for each L, serious run
        'lbit79', # 30000 realizations u[16] f=1.0, L=[12,...32] with dynamic max-time for each L, serious run
        'lbit80', # 3000 realizations u[8,16] f=0.5,1.0, L=[8,12,16,20] with dynamic max-time for each L, investigate entropy inversion
        'lbit81', # 2400 realizations u[8,16] f=0.5,1.0, L=[8,12,16,20] with dynamic max-time for each L, investigate entropy inversion with more precise svd
        'lbit82', # 6000 realizations u[4,8,16] f=0.1,0.25,0.5, L=[12,16] with dynamic max-time for each L, investigate entropy inversion with more precise svd
        'lbit83', # 1500 realizations u[4...16] f=[0.1...1.0], L=[12,16] with dynamic max-time for each L, investigate entropy inversion with more precise svd
        'lbit84', # 50000 realizations u[16] f=[0.2...0.5], L=[12,16] dynamic max-time, serious run with trnc_err 1e-6 and init shuffled neel state
        'lbit85', # 20000 realizations u[16] f=[0.2...0.5], L=[12,16] dynamic max-time, test to see if Neel initial state recovers L dependence. it did not
        'lbit86', # 10000 realizations u[16] f=[0.4], L=[12,16] check if the poster hamiltonian recovers L dependence
        'lbit90',  # 500 realizations u[16] f=[0.4], L=[12,16] fixed many bugs and increased precision with __float128: Check L dependence in early time
        'lbit91', # 20000 realizations u[16] f=[0.3], L=[12,16] fixed many bugs and increased precision with __float128: Check L dependence in early time
        'lbit92', # 20000 realizations u[16] f=[0.2-0.5], L=[12,16] fixed many bugs and increased precision with __float128: Check L dependence in early time FAILED
        'lbit93', # 200000 realizations u[16] f=[0.2-0.5], L=[12-28] fixed many bugs and increased precision with __float128:
        'lbit94', # 80000 realizations u[16] f=[0.2], L=[16] to test eps1e-8
        'lbit95',  # 80000 realizations u[16] f=[0.2], L=[16] to test eps1e-5 with the same seeds as lbit94
        'lbit96',  # filter bond_lim_2048 from lbit93
        'lbit97',  # filter bond_lim_8192 from lbit93
        'lbit98',  # L16 u8 to test weird time scales
        'lbit99',  # L16 u16 tgw8 expdecay to check entanglement entropy distribution
        'lbit100',  # L16 u16 tgw8 expdecay with random product state init to check entanglement entropy distribution SUCCESS!!!
        'lbit101',  # L16 u16 tgw8 expdecay with random product state init to check entanglement entropy distribution FINAL FULL SET OF SIMS
        'lbit103',  # L16 u16 tgw8 expdecay with neel init in the non-interacting limit
        'lbit104',  # L16 u16 tgw8 expdecay with two down one site from the middle (2d1)
        'lbit106',  # L12,L16 u16 tgw8 expdecay to check linspaced sampling, with neel initial state
        'lbit107',  # L12,L16 u16 tgw8 expdecay to check linspaced sampling, with dislocated neel initial state
        'lbit108',  # L16 u16 tgw8 r=1,2,3, expdecay, to check linspaced sampling with neel, for J2 range 1,2 and 3, with J3
        'lbit109',  # L16 u16 tgw8 r=1,2,3,4, expdecay, to check linspaced sampling with neel, for J2 ranges, without J3
        'lbit110',  # L16 u16 tgw8 r=2,3,4, expdecay, to check why there is a difference at small times
        'lbit111',  # L16 u16 tgw8 r=2,3,4, expdecay, with low svd precision to check why there is a difference at small times
        'lbit112',  # L16 u16 tgw8 , expdecay, with midchain singlet neel state to check if min/max SN evolves more symmetrically
        'lbit113',  # L12-28 to calculate lbit stats for f0.2, u16 tgw8 expdecay
        'lbit114',  # L16 to calculate the SN distribution without weights at f=0.06
        'lbit115',  # L16 xi_J=0.5,2.0 to see how it affects the saturation time of SE
        'lbit116',  # L16 anderson mode to see how it affects the saturation time of SN
        'lbit117',  # L[12-20] anderson mode to see how it affects SN(L)
        'lbit118',  # L[12-20] mbl mode to see how it affects SN(L) compared to lbit117
        'lbit119',  # FAILED (typo in the hermitian matrix in the circuit) L[8-16] vary the new lambda parameter
        'lbit120',  # L16 with MATRIX_(V1|V2|V3) varying lambda and to make sure the same f can be used in all three

    ]
    xdmrg_batches = ['data170',  #
                    ]



    parser = argparse.ArgumentParser(description='dmrg-plot')
    parser.add_argument('--clear', action='store_true', help='Remake averaged.h5')
    parser.add_argument('--basedir', type=str, help='Main directory for all mbl batch data', default=basedir)
    parser.add_argument('--algos', type=list, help='List of algorithms to plot data for', choices=['fLBIT', 'xDMRG'], default=[algo])
    parser.add_argument('--states', type=list, help='List of states to plot data for', default=['state_real'])
    parser.add_argument('--models', type=list, help='List of models to plot data for', default=['analysis', 'model'])
    parser.add_argument('--points', type=list, help='List of points to plot data for', default=['tables'])
    if algo == 'fLBIT':
        parser.add_argument('--batches', type=list, help='List of batches to plot data for',
                            default=[flbit_batches[-1]] if batches is None else batches)
    elif algo == 'xDMRG':
        parser.add_argument('--batches', type=list, help='List of batches to plot data for',
                            default=[xdmrg_batches[-1]] if batches is None else batches)

    args = parser.parse_args()
    return args
