def get_fdmrg_batch_setup(projectname):
    batch_setup = {
        'fdmrg-isingmajorana1-init': {
            # Calculate subsystem entanglement entropies for  ground states
            'projectname': 'fdmrg-isingmajorana1-init',
            'batch': {
                'L16|g0.5': {
                    'seed_extent': [100],
                    'seed_offset': [160500000],
                },
            }
        },
        'fdmrg-see-test1': {
            # Calculate subsystem entanglement entropies for  ground states
            'projectname': 'fdmrg-see-test1',
            'batch': {
                'L16|g0.5|d-8.00': {
                    'seed_extent': [1000],
                    'seed_offset': [160500000],
                },
                'L16|g0.5|d0.00': {
                    'seed_extent': [1000],
                    'seed_offset': [160510000],
                },
                'L16|g0.5|d8.00': {
                    'seed_extent': [1000],
                    'seed_offset': [160520000],
                },
            }
        },
        'fdmrg-see-test2': {
            # Calculate subsystem entanglement entropies for  ground states
            'projectname': 'fdmrg-see-test2',
            'batch': {
                'L16|g0.5|d-8.00': {
                    'seed_extent': [1000],
                    'seed_offset': [160500000],
                },
                'L16|g0.5|d0.00': {
                    'seed_extent': [1000],
                    'seed_offset': [160510000],
                },
                'L16|g0.5|d8.00': {
                    'seed_extent': [1000],
                    'seed_offset': [160520000],
                },
                'L17|g0.5|d-8.00': {
                    'seed_extent': [1000],
                    'seed_offset': [170500000],
                },
                'L17|g0.5|d0.00': {
                    'seed_extent': [1000],
                    'seed_offset': [170510000],
                },
                'L17|g0.5|d8.00': {
                    'seed_extent': [1000],
                    'seed_offset': [170520000],
                },
            }
        },
        'fdmrg-see-test3': {
            # Calculate subsystem entanglement entropies for  ground states
            'projectname': 'fdmrg-see-test3',
            'batch': {
                'L25|g0.5|d-8.00': {
                    'seed_extent': [1000],
                    'seed_offset': [251000000],
                },
                'L25|g0.5|d0.00': {
                    'seed_extent': [1000],
                    'seed_offset': [252000000],
                },
                'L25|g0.5|d8.00': {
                    'seed_extent': [1000],
                    'seed_offset': [253000000],
                },
            }
        },
        'fdmrg-see-test4': {
            # Calculate subsystem entanglement entropies for  ground states
            'projectname': 'fdmrg-see-test4',
            'batch': {
                'L25|g0.5|d-8.00': {
                    'seed_extent': [1000],
                    'seed_offset': [251000000],
                },
                'L25|g0.5|d0.00': {
                    'seed_extent': [1000],
                    'seed_offset': [252000000],
                },
                'L25|g0.5|d8.00': {
                    'seed_extent': [1000],
                    'seed_offset': [253000000],
                },
            }
        },
        'fdmrg-see-test5': {
            # Calculate subsystem entanglement entropies for  ground states
            'projectname': 'fdmrg-see-test5',
            'batch': {
                'L31|g0.5|d-8.00': {
                    'seed_extent': [1000],
                    'seed_offset': [311000000],
                },
                'L31|g0.5|d0.00': {
                    'seed_extent': [1000],
                    'seed_offset': [312000000],
                },
                'L31|g0.5|d8.00': {
                    'seed_extent': [1000],
                    'seed_offset': [313000000],
                },
            }
        },
        'fdmrg6-see': {
            'description': "Production run with L=35 to calculate subsystem entanglements with truncation error 1e-6."
                           "The seeds encode the Hamiltonian parameters"
                           " Eg:  1|600|0000|00000"
                           "     +-| d | g  | extra ",
            'projectname': 'fdmrg6-see',
            'batch': {
                'd+3.00|g0.500': {'seed_extent': [5000], 'seed_offset': [2_300_0500_000000], },
                'd+6.00|g0.500': {'seed_extent': [5000], 'seed_offset': [2_600_0500_000000], },
                'd+9.00|g0.500': {'seed_extent': [5000], 'seed_offset': [2_900_0500_000000], },
            }
        },
    }

    return batch_setup[projectname]
