# variance_threshold_upper = 5e-11
# variance_threshold_lower = 5e-13
from .plotting.style import *
import numpy as np

# variance_window_limits = [[1e-20, 10   ],
#                           [1e-18, 1e-6]]
# variance_window_names  = [['$\sigma_E^2 \in $['+str(variance_window_limits[0][0])  + ',' + str(variance_window_limits[0][1]) + ']', None],
#                           ['$\sigma_E^2 \in $['+str(variance_window_limits[1][0])  + ',' + str(variance_window_limits[1][1]) + ']', None]]
variance_window_limits = [[None, 1e+4]]
variance_window_names = [['$\sigma_E^2 \in $[' + str(variance_window_limits[0][0]) + ',' + str(variance_window_limits[0][1]) + ']', None]]

energy_window_limits = [[0.0, 1.0]]
energy_window_names = [['$E \in $[' + str(energy_window_limits[0][0]) + ',' + str(energy_window_limits[0][1]) + ']', None]]

# current_palette = itertools.cycle(sns.color_palette())
# variance_colors = [current_palette[0], current_palette[1]]
# variance_colors = itertools.islice(sns.color_palette(),0, None)
entanglement_threshold_upper = 100
entanglement_threshold_lower = 0
energy_window = 1.0
energy_target = 0.5
max_bond_dimension = 512
max_realizations = 1000000

eps_array = np.linspace(-16, 1, 25)
eps_array = 10 ** eps_array

zero2one = np.linspace(0, 1, 25)
winsize = 0.1
ewin_array = [[x - winsize / 2.0, x + winsize / 2.0] for x in zero2one]
