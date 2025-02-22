def get_batch_setup(projectname):
    batch_setup = {
        'lbit93-precision' : {
            'projectname': 'lbit93-precision',
            'batch': {  # Number of seeds to run
                'L8|f0.2': {
                    'seed_extent': [100],
                    'seed_offset': [20000],
                    'time_steps': 100,
                },
                'L8|f0.3': {
                    'seed_extent': [100],
                    'seed_offset': [30000],
                    'time_steps': 100,
                },
                'L8|f0.4': {
                    'seed_extent': [100],
                    'seed_offset': [40000],
                    'time_steps': 100,
                },
                'L12|f0.2': {
                    'seed_extent': [10000, 16000, 26000, 8000, 10000, 4000, 6000, 120000],
                    'seed_offset': [10030000, 15045000, 16050000, 20022500, 25030000, 30007500, 35015000, 12200000],
                    'time_steps': 200,
                },
                'L12|f0.3': {
                    'seed_extent': [10000, 16000, 26000, 8000, 10000, 4000, 6000],
                    'seed_offset': [10010000, 15015000, 16025000, 20007500, 25010000, 30002500, 35005000],
                    'time_steps': 200,
                },
                'L12|f0.4': {
                    'seed_extent': [10000, 16000, 26000, 8000, 10000, 4000, 6000, 120000],
                    'seed_offset': [10070000, 15105000, 16125000, 20052500, 25070000, 30017500, 35035000,12400000],
                    'time_steps': 200,
                },
                'L16|f0.2': {
                    'seed_extent': [10000, 16000, 26000, 8000, 10000, 4000, 6000, 120000],
                    'seed_offset': [10040000, 15060000, 16075000, 20030000, 25040000, 30010000, 35020000,16200000],
                    'time_steps': 200,
                },
                'L16|f0.3': {
                    'seed_extent': [10000, 16000, 26000, 8000, 10000, 4000, 6000],
                    'seed_offset': [10060000, 15090000, 16100000, 20045000, 25060000, 30015000, 35030000],
                    'time_steps': 200,
                },
                'L16|f0.4': {
                    'seed_extent': [10000, 16000, 26000, 8000, 10000, 4000, 6000, 120000],
                    'seed_offset': [10000000, 15000000, 16000000, 20000000, 25000000, 30000000, 35000000,16400000],
                    'time_steps': 200,
                },
                'L20|f0.2': {
                    'seed_extent': [15000, 15000, 10000, 10000, 30000,120000],
                    'seed_offset': [16030000, 16115000, 26020000, 26110000, 60600000,20200000],
                    'time_steps': 200,
                },
                'L20|f0.3': {
                    'seed_extent': [15000, 15000, 10000, 10000, 30000],
                    'seed_offset': [16045000, 16130000, 26030000, 26120000, 60700000],
                    'time_steps': 200,
                },
                'L20|f0.4': {
                    'seed_extent': [15000, 15000, 10000, 10000, 30000, 120000],
                    'seed_offset': [16000000, 16100000, 26000000, 26100000, 60800000,20400000],
                    'time_steps': 200,
                },
                'L24|f0.2': {
                    'seed_extent': [5000, 5000, 5000, 10000, 1000, 1000, 1000, 1000, 1000, 100, 9900, 10000, 50000],
                    'seed_offset': [16310000, 26310000, 46000000, 46100000, 46300000, 46400000, 46500000, 46600000, 46700000,
                               57000000, 57000300, 61000000, 24020000],
                    'time_steps': 100,
                },
                'L24|f0.3': {
                    'seed_extent': [5000, 5000, 5000, 10000, 1000, 1000, 1000, 1000, 1000, 100, 9900, 10000],
                    'seed_offset': [16305000, 26305000, 46015000, 46130000, 46302000, 46402000, 46502000, 46602000, 46702000,
                               57000200, 57020100, 61100000],
                    'time_steps': 100,
                },
                'L24|f0.4': {
                    'seed_extent': [5000, 5000, 5000, 10000, 1000, 1000, 1000, 1000, 1000, 100, 9900, 10000],
                    'seed_offset': [16300000, 26300000, 46005000, 46110000, 46301000, 46401000, 46501000, 46601000, 46701000,
                               57000100, 57010200, 61200000],
                    'time_steps': 100,
                },
                'L28|f0.2': {
                    'seed_extent': [1000, 4000, 5000, 5000, 5000, 5000, 5000, 70000],
                    'seed_offset': [47000000, 47001000, 47005000, 47010000, 47015000, 47020000, 47025000, 28200000],
                    'time_steps': 100,
                },
                'L32|f0.2': {
                    'seed_extent': [25000],
                    'seed_offset': [32200000],
                    'time_steps': 100,
                },
                # 'L28|f0.2': {
                #     'seed_extent': [1000, 4000, 45000],
                #     'seed_offset': [47000000, 47100000, 62000000],
                #     'time_steps': 100,
                # },
                # 'L28|f0.3': {
                #     'seed_extent': [1000, 4000, 45000],
                #     'seed_offset': [47002000, 47108000, 62100000],
                #     'time_steps': 100,
                # },
                # 'L28|f0.4': {
                #     'seed_extent': [1000, 4000, 45000],
                #     'seed_offset': [47001000, 47104000, 62200000],
                #     'time_steps': 100,
                # },
                # 'L32|f0.2': {
                #     'seed_extent': [1000, 24000],
                #     'seed_offset': [48000000, 63000000],
                #     'time_steps': 100,
                # },
                # 'L32|f0.3': {
                #     'seed_extent': [1000, 24000],
                #     'seed_offset': [48001000, 63100000],
                #     'time_steps': 100,
                # },
                # 'L32|f0.4': {
                #     'seed_extent': [1000, 24000],
                #     'seed_offset': [48002000, 63200000],
                #     'time_steps': 100,
                # },
            },
        },
        'lbit94-eps1e-8': {
            'projectname': 'lbit94-eps1e-8',
            'batch': {  # Number of seeds to run
                'L16|f0.2': {
                    'seed_extent': [80000],
                    'seed_offset': [10000000],
                    'time_steps': 100,
                },
            },
        },
        'lbit95-eps1e-5': {
            'projectname': 'lbit95-eps1e-5',
            'batch': {  # Number of seeds to run
                'L16|f0.2': {
                    'seed_extent': [80000],
                    'seed_offset': [10000000],
                    'time_steps': 100,
                },
            },
        },
        'lbit98-L16-u8': {
            'projectname': 'lbit98-L16-u8',
            'batch': {  # Number of seeds to run
                'L16|f0.2': {
                    'seed_extent': [10000],
                    'seed_offset': [10000000],
                    'time_steps': 100,
                },
            },
        },
        'lbit99-tgw8-ex': {
            'projectname': 'lbit99-tgw8-ex',
            'batch': {  # Number of seeds to run
                'L16|f0.2': {
                    'seed_extent': [10000],
                    'seed_offset': [10000000],
                    'time_steps': 100,
                },
            },
        },
        'lbit100-rps': {
            'projectname': 'lbit100-rps',
            'batch': {  # Number of seeds to run
                'L16|f0.2': {
                    'seed_extent': [10000],
                    'seed_offset': [10000000],
                    'time_steps': 100,
                },
            },
        },
        'lbit101-rps': {
            'projectname': 'lbit101-rps',
            'batch': {  # Number of seeds to run
                'L12|f0.2': {
                    'seed_extent': [100000, 100000],
                    'seed_offset': [12200000,12300000],
                    'time_steps': 200,
                },
                'L12|f0.4': {
                    'seed_extent': [100000,100000],
                    'seed_offset': [12400000,12500000],
                    'time_steps': 200,
                },
                'L14|f0.2': {
                    'seed_extent': [80000],
                    'seed_offset': [14200000],
                    'time_steps': 200,
                },
                'L14|f0.4': {
                    'seed_extent': [80000],
                    'seed_offset': [14400000],
                    'time_steps': 200,
                },
                'L16|f0.2': {
                    'seed_extent': [100000, 100000],
                    'seed_offset': [16200000,16300000],
                    'time_steps': 200,
                },
                'L16|f0.4': {
                    'seed_extent': [100000, 100000],
                    'seed_offset': [16400000,16500000],
                    'time_steps': 200,
                },
                'L18|f0.2': {
                    'seed_extent': [80000],
                    'seed_offset': [18200000],
                    'time_steps': 200,
                },
                'L18|f0.4': {
                    'seed_extent': [80000],
                    'seed_offset': [18400000],
                    'time_steps': 200,
                },
                'L20|f0.2': {
                    'seed_extent': [100000, 100000],
                    'seed_offset': [20200000, 20300000],
                    'time_steps': 200,
                },
                'L20|f0.4': {
                    'seed_extent': [100000, 100000],
                    'seed_offset': [20400000, 20500000],
                    'time_steps': 200,
                },
                'L24|f0.2': {
                    'seed_extent': [50000],
                    'seed_offset': [24200000],
                    'time_steps': 100,
                },
                'L24|f0.4': {
                    'seed_extent': [50000],
                    'seed_offset': [24400000],
                    'time_steps': 100,
                },
                'L28|f0.2': {
                    'seed_extent': [30000, 20000],
                    'seed_offset': [28200000,28230000],
                    'time_steps': 100,
                },
                'L28|f0.4': {
                    'seed_extent': [30000, 20000],
                    'seed_offset': [28400000, 28430000],
                    'time_steps': 100,
                },
                'L32|f0.2': {
                    'seed_extent': [10000, 5000],
                    'seed_offset': [32200000, 32210000],
                    'time_steps': 100,
                },

            },
        },
        'lbit103-nil': {
            'projectname': 'lbit103-nil',
            'batch': {  # Number of seeds to run
                'L12|f0.2': {
                    'seed_extent': [200000],
                    'seed_offset': [12200000],
                    'time_steps': 200,
                },
                'L12|f0.4': {
                    'seed_extent': [200000],
                    'seed_offset': [12400000],
                    'time_steps': 200,
                },
                'L16|f0.2': {
                    'seed_extent': [200000],
                    'seed_offset': [16200000],
                    'time_steps': 200,
                },
                'L16|f0.4': {
                    'seed_extent': [200000],
                    'seed_offset': [16400000],
                    'time_steps': 200,
                },
                'L20|f0.2': {
                    'seed_extent': [200000],
                    'seed_offset': [20200000],
                    'time_steps': 200,
                },
                'L20|f0.4': {
                    'seed_extent': [200000],
                    'seed_offset': [20400000],
                    'time_steps': 200,
                },
                'L24|f0.2': {
                    'seed_extent': [50000],
                    'seed_offset': [24200000],
                    'time_steps': 100,
                },
                'L24|f0.4': {
                    'seed_extent': [50000],
                    'seed_offset': [24400000],
                    'time_steps': 100,
                },
                'L28|f0.2': {
                    'seed_extent': [30000],
                    'seed_offset': [28200000],
                    'time_steps': 100,
                },
                'L28|f0.4': {
                    'seed_extent': [30000],
                    'seed_offset': [28400000],
                    'time_steps': 100,
                },
                'L32|f0.2': {
                    'seed_extent': [10000],
                    'seed_offset': [32200000],
                    'time_steps': 100,
                },

            },
        },
        'lbit104-2d1': {
            'projectname': 'lbit104-2d1',
            'batch': {  # Number of seeds to run
                'L16|f0.2': {
                    'seed_extent': [50000],
                    'seed_offset': [16200000],
                    'time_steps': 200,
                },
                'L16|f0.4': {
                    'seed_extent': [50000],
                    'seed_offset': [16400000],
                    'time_steps': 200,
                },
            },
        },
        'lbit105-2d2': {
            'projectname': 'lbit105-2d2',
            'batch': {  # Number of seeds to run
                'L16|f0.2': {
                    'seed_extent': [50000],
                    'seed_offset': [16200000],
                    'time_steps': 200,
                },
                'L16|f0.4': {
                    'seed_extent': [50000],
                    'seed_offset': [16400000],
                    'time_steps': 200,
                },
            },
        },
        'lbit106-lin':{
            'projectname' : 'lbit106-lin',
            'batch': {
                'L12|f0.2': {
                    'seed_extent': [10000],
                    'seed_offset': [12200000],
                    'time_steps': 4001,
                },
                'L16|f0.2': {
                    'seed_extent': [10000],
                    'seed_offset': [16200000],
                    'time_steps': 4001,
                },
            }
        },
        'lbit107-lind': {
            'projectname': 'lbit107-lind',
            'batch': {
                'L12|f0.2': {
                    'seed_extent': [10000],
                    'seed_offset': [12200000],
                    'time_steps': 4001,
                },
                # 'L16|f0.2': {
                #     'seed_extent': [10000],
                #     'seed_offset': [16200000],
                #     'time_steps': 4001,
                # },
            }
        },
        'lbit108-linr': {
            'projectname': 'lbit108-linr',
            'batch': {
                'L16|f0.2': {
                    'seed_extent': [10000],
                    'seed_offset': [16200000],
                    'time_steps': 4001,
                },
            }
        },
        'lbit109-linr': {
            # Disabled J3 for this one.. otherwise it's the same as lbit108
            'projectname': 'lbit109-linr',
            'batch': {
                'L16|f0.2': {
                    'seed_extent': [10000],
                    'seed_offset': [16200000],
                    'time_steps': 4001,
                },
            }
        },
        'lbit110-tsmall': {
            # Investigate at what timescale different J2_span diverge (answer, at around 1 second in disorder average)
            'projectname': 'lbit110-tsmall',
            'batch': {
                'L16|f0.2': {
                    'seed_extent': [10000],
                    'seed_offset': [16200000],
                    'time_steps': 100,
                },
            }
        },
        'lbit111-lowprec': {
            # Investigate at what timescale different J2_span diverge when using svd_threshold=1e-5
            'projectname': 'lbit111-lowprec',
            'batch': {
                'L16|f0.2': {
                    'seed_extent': [10000],
                    'seed_offset': [16200000],
                    'time_steps': 100,
                },
            }
        },
        'lbit112-singlet': {
            # Investigate at what timescale different J2_span diverge when using svd_threshold=1e-5
            'projectname': 'lbit112-singlet',
            'batch': {
                'L16|f0.2': {
                    'seed_extent': [10000],
                    'seed_offset': [16200000],
                    'time_steps': 4001,
                },
            }
        },
        'lbit113-lbit': {
            # Calculate the lbit supports for many system sizes
            'projectname': 'lbit113-lbit',
            'batch': {
                'L12': {
                    'seed_extent': [10000],
                    'seed_offset': [12200000],
                    'time_steps': 0,
                },
                'L16': {
                    'seed_extent': [10000],
                    'seed_offset': [16200000],
                    'time_steps': 0,
                },
                'L20': {
                    'seed_extent': [10000],
                    'seed_offset': [20200000],
                    'time_steps': 0,
                },
                'L24': {
                    'seed_extent': [10000],
                    'seed_offset': [24200000],
                    'time_steps': 0,
                },
                'L28': {
                    'seed_extent': [10000],
                    'seed_offset': [28200000],
                    'time_steps': 0,
                },
                'L32': {
                    'seed_extent': [10000],
                    'seed_offset': [32200000],
                    'time_steps': 0,
                },
            }
        },
        'lbit114-now8': {
            # Calculate the lbit supports for many system sizes
            'projectname': 'lbit114-now8',
            'batch': {
                'L16': {
                    'seed_extent': [10000],
                    'seed_offset': [16200000],
                    'time_steps': 0,
                },
            }
        },
        'lbit115-xi': {
            # Calculate lbit93 at xi=0.5 and 2.0 to see how the slope of SE is affected
            'projectname': 'lbit115-xi',
            'batch': {
                'L16': {
                    'seed_extent': [1000],
                    'seed_offset': [16200000],
                    'time_steps': 200,
                },
            }
        },
        'lbit116-anderson': {
            # Calculate lbit93 at with the anderson unitary circuit to see how SN is affected
            'projectname': 'lbit116-anderson',
            'batch': {
                'L16': {
                    'seed_extent': [50000],
                    'seed_offset': [16200000],
                    'time_steps': 200,
                },
            }
        },
        'lbit117-anderson': {
            # Calculate lbit93 with the anderson circuit to check SN(L)
            'projectname': 'lbit117-anderson',
            'batch': {
                'L12': {
                    'seed_extent': [50000],
                    'seed_offset': [12000000],
                    'time_steps': 200,
                },
                'L16': {
                    'seed_extent': [50000],
                    'seed_offset': [16000000],
                    'time_steps': 200,
                },
                'L20': {
                    'seed_extent': [50000],
                    'seed_offset': [20000000],
                    'time_steps': 200,
                },
            }
        },
        'lbit118-mbl': {
            # Calculate lbit93 with the mbl circuit to check SN(L)
            'projectname': 'lbit118-mbl',
            'batch': {
                'L12': {
                    'seed_extent': [50000],
                    'seed_offset': [12000000],
                    'time_steps': 200,
                },
                'L16': {
                    'seed_extent': [50000],
                    'seed_offset': [16000000],
                    'time_steps': 200,
                },
                'L20': {
                    'seed_extent': [50000],
                    'seed_offset': [20000000],
                    'time_steps': 200,
                },
            }
        },
        'lbit119-m3': {
            # Calculate SN(L) varying the new lambda parameter
            'projectname': 'lbit119-m3',
            'batch': {
                'L8': {
                    'seed_extent': [10000],
                    'seed_offset': [8000000],
                    'time_steps': 200,
                },
                'L12': {
                    'seed_extent': [10000],
                    'seed_offset': [12000000],
                    'time_steps': 200,
                },
                'L16': {
                    'seed_extent': [10000],
                    'seed_offset': [16000000],
                    'time_steps': 200,
                },
                'L20': {
                    'seed_extent': [10000],
                    'seed_offset': [20000000],
                    'time_steps': 200,
                },
            }
        },
        'lbit120-m123': {
            # Calculate SN(L) varying the new lambda parameter
            'projectname': 'lbit120-m123',
            'batch': {
                'L16': {
                    'seed_extent': [10000],
                    'seed_offset': [16000000],
                    'time_steps': 200,
                },
            }
        },
        'lbit121-opdm': {
            # Calculate SN(L) varying the new lambda parameter
            'projectname': 'lbit121-opdm',
            'batch': {
                'L16': {
                    'seed_extent': [10000],
                    'seed_offset': [16000000],
                    'time_steps': 200,
                },
            }
        },
        'lbit122-fleff': {
            # Calculate SN(L) varying the new lambda parameter
            'projectname': 'lbit122-fleff',
            'batch': {
                'L8': {
                    'seed_extent': [1000],
                    'seed_offset': [8000000],
                    'time_steps': 200,
                },
                'L10': {
                    'seed_extent': [1000],
                    'seed_offset': [10000000],
                    'time_steps': 200,
                },
                'L12': {
                    'seed_extent': [1000],
                    'seed_offset': [12000000],
                    'time_steps': 200,
                },
                'L14': {
                    'seed_extent': [1000],
                    'seed_offset': [14000000],
                    'time_steps': 200,
                },
            }
        },
        'lbit123-fleff': {
            # Calculate SN(L) varying the new lambda parameter
            'projectname': 'lbit123-fleff',
            'batch': {
                'L12|l0.1': {
                    'seed_extent': [10000],
                    'seed_offset': [12000000],
                    'time_steps': 200,
                },
                'L12|l1.0': {
                    'seed_extent': [10000],
                    'seed_offset': [12100000],
                    'time_steps': 200,
                },
                'L12|l5.0': {
                    'seed_extent': [10000],
                    'seed_offset': [12500000],
                    'time_steps': 200,
                },
                'L16|l0.1': {
                    'seed_extent': [10000],
                    'seed_offset': [16000000],
                    'time_steps': 200,
                },
                'L16|l1.0': {
                    'seed_extent': [10000],
                    'seed_offset': [16100000],
                    'time_steps': 200,
                },
                'L16|l5.0': {
                    'seed_extent': [10000],
                    'seed_offset': [16500000],
                    'time_steps': 200,
                },
                'L20|l0.1': {
                    'seed_extent': [10000],
                    'seed_offset': [20000000],
                    'time_steps': 200,
                },
                'L20|l1.0': {
                    'seed_extent': [10000],
                    'seed_offset': [20100000],
                    'time_steps': 200,
                },
                'L20|l5.0': {
                    'seed_extent': [10000],
                    'seed_offset': [20500000],
                    'time_steps': 200,
                },
            }
        },
        'lbit124-linspaced': {
            # Time evolve to linearly spaced time points after saturation to calculate min/max SN(L)
            'projectname': 'lbit124-linspaced',
            'batch': {
                'L12': {
                    'seed_extent': [10000],
                    'seed_offset': [12000000],
                    'time_steps': 41,
                },
                'L16': {
                    'seed_extent': [10000],
                    'seed_offset': [16000000],
                    'time_steps': 41,
                },
                'L20': {
                    'seed_extent': [10000],
                    'seed_offset': [20000000],
                    'time_steps': 41,
                },
                'L24': {
                    'seed_extent': [10000],
                    'seed_offset': [24000000],
                    'time_steps': 41,
                },
                'L28': {
                    'seed_extent': [10000],
                    'seed_offset': [28000000],
                    'time_steps': 41,
                },
                'L32': {
                    'seed_extent': [1000],
                    'seed_offset': [32000000],
                    'time_steps': 41,
                },
            }
        },
        'lbit125-linspaced': {
            # Time evolve to linearly spaced time points after saturation to calculate min/max SN(L)
            'projectname': 'lbit125-linspaced',
            'batch': {
                'L12': {
                    'seed_extent': [10000],
                    'seed_offset': [12000000],
                    'time_steps': 41,
                },
                'L16': {
                    'seed_extent': [10000],
                    'seed_offset': [16000000],
                    'time_steps': 41,
                },
                'L20': {
                    'seed_extent': [10000],
                    'seed_offset': [20000000],
                    'time_steps': 41,
                },
                'L24': {
                    'seed_extent': [10000],
                    'seed_offset': [24000000],
                    'time_steps': 41,
                },
                'L28': {
                    'seed_extent': [10000],
                    'seed_offset': [28000000],
                    'time_steps': 41,
                },
                'L32': {
                    'seed_extent': [1000],
                    'seed_offset': [32000000],
                    'time_steps': 41,
                },
            }
        },
        'lbit126-fleff': {
            # Calculate SN(L) varying the new lambda parameter
            'projectname': 'lbit126-fleff',
            'batch': {
                'L12|l0.1': {
                    'seed_extent': [10000],
                    'seed_offset': [12000000],
                    'time_steps': 200,
                },
                'L12|l1.0': {
                    'seed_extent': [10000],
                    'seed_offset': [12100000],
                    'time_steps': 200,
                },
                'L12|l5.0': {
                    'seed_extent': [10000],
                    'seed_offset': [12500000],
                    'time_steps': 200,
                },
                'L16|l0.1': {
                    'seed_extent': [10000],
                    'seed_offset': [16000000],
                    'time_steps': 200,
                },
                'L16|l1.0': {
                    'seed_extent': [10000],
                    'seed_offset': [16100000],
                    'time_steps': 200,
                },
                'L16|l5.0': {
                    'seed_extent': [10000],
                    'seed_offset': [16500000],
                    'time_steps': 200,
                },
                'L20|l0.1': {
                    'seed_extent': [10000],
                    'seed_offset': [20000000],
                    'time_steps': 200,
                },
                'L20|l1.0': {
                    'seed_extent': [10000],
                    'seed_offset': [20100000],
                    'time_steps': 200,
                },
                'L20|l5.0': {
                    'seed_extent': [10000],
                    'seed_offset': [20500000],
                    'time_steps': 200,
                },
            }
        },
        'lbit127-linr': {
            # Same as lbit124 but taking time step 0.25 between times 0.25 and 1600
            'projectname': 'lbit127-linr',
            'batch': {
                'L16|f0.2': {
                    'seed_extent': [10000],
                    'seed_offset': [16200000],
                    'time_steps': 6400,
                },
            }
        },
        'lbit128-linspaced': {
            # Time evolve to linearly spaced time points after saturation to calculate min/max SN(L)
            'projectname': 'lbit128-linspaced',
            'batch': {
                'L12': {
                    'seed_extent': [10000],
                    'seed_offset': [12000000],
                    'time_steps': 201,
                },
                'L16': {
                    'seed_extent': [10000],
                    'seed_offset': [16000000],
                    'time_steps': 201,
                },
                'L20': {
                    'seed_extent': [10000],
                    'seed_offset': [20000000],
                    'time_steps': 201,
                },
                'L24': {
                    'seed_extent': [10000],
                    'seed_offset': [24000000],
                    'time_steps': 201,
                },
                'L28': {
                    'seed_extent': [10000],
                    'seed_offset': [28000000],
                    'time_steps': 201,
                },
                'L32': {
                    'seed_extent': [1000],
                    'seed_offset': [32000000],
                    'time_steps': 201,
                },
            }
        },
        'lbit129-linr': {
            # Same as lbit124 but taking time step 0.25 between times 0.25 and 1600
            'projectname': 'lbit129-linr',
            'batch': {
                'L16|f0.2': {
                    'seed_extent': [10000],
                    'seed_offset': [16200000],
                    'time_steps': 6400,
                },
            }
        },
    }

    return batch_setup[projectname]