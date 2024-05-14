def get_xdmrg_batch_setup(projectname):
    batch_setup = {
        'xdmrg1-fse': {
            'description': "First batch in 2024 for xDMRG using forward subsystem expansion and deparallelization."
                           "The seeds encode the Hamiltonian parameters Eg: 16|00|00000"
                           "The seeds start with 1 for negative delta, and 2 for positive."
                           "The next three second digit represents the g-value."
                           "The remaining digits reveal are the decimals of g.",
            'projectname': 'xdmrg1-fse',
            'batch': {
                'd-6.00|g0.00': {'seed_extent': [1000], 'seed_offset': [16000000000], },
                'd-5.00|g0.00': {'seed_extent': [1000], 'seed_offset': [15000000000], },
                'd-4.00|g0.00': {'seed_extent': [1000], 'seed_offset': [14000000000], },
                'd-3.00|g0.00': {'seed_extent': [1000], 'seed_offset': [13000000000], },
                'd-2.00|g0.00': {'seed_extent': [1000], 'seed_offset': [12000000000], },
                'd-1.00|g0.00': {'seed_extent': [1000], 'seed_offset': [11000000000], },
                'd+0.00|g0.00': {'seed_extent': [1000], 'seed_offset': [20000000000], },
                'd+1.00|g0.00': {'seed_extent': [1000], 'seed_offset': [21000000000], },
                'd+2.00|g0.00': {'seed_extent': [1000], 'seed_offset': [22000000000], },
                'd+3.00|g0.00': {'seed_extent': [1000], 'seed_offset': [23000000000], },
                'd+4.00|g0.00': {'seed_extent': [1000], 'seed_offset': [24000000000], },
                'd+5.00|g0.00': {'seed_extent': [1000], 'seed_offset': [25000000000], },
                'd+6.00|g0.00': {'seed_extent': [1000], 'seed_offset': [26000000000], },

                'd-6.00|g0.01': {'seed_extent': [1000], 'seed_offset': [16001000000], },
                'd-5.00|g0.01': {'seed_extent': [1000], 'seed_offset': [15001000000], },
                'd-4.00|g0.01': {'seed_extent': [1000], 'seed_offset': [14001000000], },
                'd-3.00|g0.01': {'seed_extent': [1000], 'seed_offset': [13001000000], },
                'd-2.00|g0.01': {'seed_extent': [1000], 'seed_offset': [12001000000], },
                'd-1.00|g0.01': {'seed_extent': [1000], 'seed_offset': [11001000000], },
                'd+0.00|g0.01': {'seed_extent': [1000], 'seed_offset': [20001000000], },
                'd+1.00|g0.01': {'seed_extent': [1000], 'seed_offset': [21001000000], },
                'd+2.00|g0.01': {'seed_extent': [1000], 'seed_offset': [22001000000], },
                'd+3.00|g0.01': {'seed_extent': [1000], 'seed_offset': [23001000000], },
                'd+4.00|g0.01': {'seed_extent': [1000], 'seed_offset': [24001000000], },
                'd+5.00|g0.01': {'seed_extent': [1000], 'seed_offset': [25001000000], },
                'd+6.00|g0.01': {'seed_extent': [1000], 'seed_offset': [26001000000], },

                'd-6.00|g0.02': {'seed_extent': [1000], 'seed_offset': [16002000000], },
                'd-5.00|g0.02': {'seed_extent': [1000], 'seed_offset': [15002000000], },
                'd-4.00|g0.02': {'seed_extent': [1000], 'seed_offset': [14002000000], },
                'd-3.00|g0.02': {'seed_extent': [1000], 'seed_offset': [13002000000], },
                'd-2.00|g0.02': {'seed_extent': [1000], 'seed_offset': [12002000000], },
                'd-1.00|g0.02': {'seed_extent': [1000], 'seed_offset': [11002000000], },
                'd+0.00|g0.02': {'seed_extent': [1000], 'seed_offset': [20002000000], },
                'd+1.00|g0.02': {'seed_extent': [1000], 'seed_offset': [21002000000], },
                'd+2.00|g0.02': {'seed_extent': [1000], 'seed_offset': [22002000000], },
                'd+3.00|g0.02': {'seed_extent': [1000], 'seed_offset': [23002000000], },
                'd+4.00|g0.02': {'seed_extent': [1000], 'seed_offset': [24002000000], },
                'd+5.00|g0.02': {'seed_extent': [1000], 'seed_offset': [25002000000], },
                'd+6.00|g0.02': {'seed_extent': [1000], 'seed_offset': [26002000000], },

                'd-6.00|g0.03': {'seed_extent': [1000], 'seed_offset': [16003000000], },
                'd-5.00|g0.03': {'seed_extent': [1000], 'seed_offset': [15003000000], },
                'd-4.00|g0.03': {'seed_extent': [1000], 'seed_offset': [14003000000], },
                'd-3.00|g0.03': {'seed_extent': [1000], 'seed_offset': [13003000000], },
                'd-2.00|g0.03': {'seed_extent': [1000], 'seed_offset': [12003000000], },
                'd-1.00|g0.03': {'seed_extent': [1000], 'seed_offset': [11003000000], },
                'd+0.00|g0.03': {'seed_extent': [1000], 'seed_offset': [20003000000], },
                'd+1.00|g0.03': {'seed_extent': [1000], 'seed_offset': [21003000000], },
                'd+2.00|g0.03': {'seed_extent': [1000], 'seed_offset': [22003000000], },
                'd+3.00|g0.03': {'seed_extent': [1000], 'seed_offset': [23003000000], },
                'd+4.00|g0.03': {'seed_extent': [1000], 'seed_offset': [24003000000], },
                'd+5.00|g0.03': {'seed_extent': [1000], 'seed_offset': [25003000000], },
                'd+6.00|g0.03': {'seed_extent': [1000], 'seed_offset': [26003000000], },

            }
        },
        'xdmrg2-tuned': {
            'description': "Second batch in 2024 for xDMRG using forward subsystem expansion and deparallelization."
                           "The first batch had several types of convergence issues on different rare events that "
                           "should now be addressed."
                           "The seeds encode the Hamiltonian parameters Eg: 16|00|00000"
                           "The seeds start with 1 for negative delta, and 2 for positive."
                           "The next three second digit represents the g-value."
                           "The remaining digits reveal are the decimals of g.",
            'projectname': 'xdmrg2-tuned',
            'batch': {
                'd-6.00|g0.00': {'seed_extent': [1000], 'seed_offset': [16000000000], },
                'd-5.00|g0.00': {'seed_extent': [1000], 'seed_offset': [15000000000], },
                'd-4.00|g0.00': {'seed_extent': [1000], 'seed_offset': [14000000000], },
                'd-3.00|g0.00': {'seed_extent': [1000], 'seed_offset': [13000000000], },
                'd-2.00|g0.00': {'seed_extent': [1000], 'seed_offset': [12000000000], },
                'd-1.00|g0.00': {'seed_extent': [1000], 'seed_offset': [11000000000], },
                'd+0.00|g0.00': {'seed_extent': [1000], 'seed_offset': [20000000000], },
                'd+1.00|g0.00': {'seed_extent': [1000], 'seed_offset': [21000000000], },
                'd+2.00|g0.00': {'seed_extent': [1000], 'seed_offset': [22000000000], },
                'd+3.00|g0.00': {'seed_extent': [1000], 'seed_offset': [23000000000], },
                'd+4.00|g0.00': {'seed_extent': [1000], 'seed_offset': [24000000000], },
                'd+5.00|g0.00': {'seed_extent': [1000], 'seed_offset': [25000000000], },
                'd+6.00|g0.00': {'seed_extent': [1000], 'seed_offset': [26000000000], },

                'd-6.00|g0.01': {'seed_extent': [1000], 'seed_offset': [16001000000], },
                'd-5.00|g0.01': {'seed_extent': [1000], 'seed_offset': [15001000000], },
                'd-4.00|g0.01': {'seed_extent': [1000], 'seed_offset': [14001000000], },
                'd-3.00|g0.01': {'seed_extent': [1000], 'seed_offset': [13001000000], },
                'd-2.00|g0.01': {'seed_extent': [1000], 'seed_offset': [12001000000], },
                'd-1.00|g0.01': {'seed_extent': [1000], 'seed_offset': [11001000000], },
                'd+0.00|g0.01': {'seed_extent': [1000], 'seed_offset': [20001000000], },
                'd+1.00|g0.01': {'seed_extent': [1000], 'seed_offset': [21001000000], },
                'd+2.00|g0.01': {'seed_extent': [1000], 'seed_offset': [22001000000], },
                'd+3.00|g0.01': {'seed_extent': [1000], 'seed_offset': [23001000000], },
                'd+4.00|g0.01': {'seed_extent': [1000], 'seed_offset': [24001000000], },
                'd+5.00|g0.01': {'seed_extent': [1000], 'seed_offset': [25001000000], },
                'd+6.00|g0.01': {'seed_extent': [1000], 'seed_offset': [26001000000], },

                'd-6.00|g0.02': {'seed_extent': [1000], 'seed_offset': [16002000000], },
                'd-5.00|g0.02': {'seed_extent': [1000], 'seed_offset': [15002000000], },
                'd-4.00|g0.02': {'seed_extent': [1000], 'seed_offset': [14002000000], },
                'd-3.00|g0.02': {'seed_extent': [1000], 'seed_offset': [13002000000], },
                'd-2.00|g0.02': {'seed_extent': [1000], 'seed_offset': [12002000000], },
                'd-1.00|g0.02': {'seed_extent': [1000], 'seed_offset': [11002000000], },
                'd+0.00|g0.02': {'seed_extent': [1000], 'seed_offset': [20002000000], },
                'd+1.00|g0.02': {'seed_extent': [1000], 'seed_offset': [21002000000], },
                'd+2.00|g0.02': {'seed_extent': [1000], 'seed_offset': [22002000000], },
                'd+3.00|g0.02': {'seed_extent': [1000], 'seed_offset': [23002000000], },
                'd+4.00|g0.02': {'seed_extent': [1000], 'seed_offset': [24002000000], },
                'd+5.00|g0.02': {'seed_extent': [1000], 'seed_offset': [25002000000], },
                'd+6.00|g0.02': {'seed_extent': [1000], 'seed_offset': [26002000000], },

                'd-6.00|g0.03': {'seed_extent': [1000], 'seed_offset': [16003000000], },
                'd-5.00|g0.03': {'seed_extent': [1000], 'seed_offset': [15003000000], },
                'd-4.00|g0.03': {'seed_extent': [1000], 'seed_offset': [14003000000], },
                'd-3.00|g0.03': {'seed_extent': [1000], 'seed_offset': [13003000000], },
                'd-2.00|g0.03': {'seed_extent': [1000], 'seed_offset': [12003000000], },
                'd-1.00|g0.03': {'seed_extent': [1000], 'seed_offset': [11003000000], },
                'd+0.00|g0.03': {'seed_extent': [1000], 'seed_offset': [20003000000], },
                'd+1.00|g0.03': {'seed_extent': [1000], 'seed_offset': [21003000000], },
                'd+2.00|g0.03': {'seed_extent': [1000], 'seed_offset': [22003000000], },
                'd+3.00|g0.03': {'seed_extent': [1000], 'seed_offset': [23003000000], },
                'd+4.00|g0.03': {'seed_extent': [1000], 'seed_offset': [24003000000], },
                'd+5.00|g0.03': {'seed_extent': [1000], 'seed_offset': [25003000000], },
                'd+6.00|g0.03': {'seed_extent': [1000], 'seed_offset': [26003000000], },

            }
        }
    }

    return batch_setup[projectname]
