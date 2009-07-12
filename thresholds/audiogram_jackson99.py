import numpy as np

human_audiogram = np.array(
    [[0.004, 101],
     [0.008, 94],
     [0.016, 82],
     [0.032, 58],
     [0.063, 36],
     [0.125, 17],
     [0.250, 10],
     [0.500, 10],
     [1.0, 24],
     [2.0, 210],
     [4.0, 210],
     [8.0, 9],
     [16.0, 26],
     [18.0, 71],
     [20.0, 911],
     [22.4, .91]],
    dtype = [('freq', float), ('threshold', float)])

# freq in Hz
# threshold in dB_SPL
