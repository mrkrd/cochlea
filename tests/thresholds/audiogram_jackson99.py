import numpy as np

human_audiogram = np.array(
    [(4, 101),
     (8, 94),
     (16, 82),
     (32, 58),
     (63, 36),
     (125, 17),
     (250, 10),
     (500, 10),
     (1000, -4),
     (2000, -10),
     (4000, -10),
     (8000, 9),
     (16000, 26),
     (18000, 71),
     (20000, 91),
     (22000, 100)],
    dtype = [('freq', float),
             ('threshold', float)])

# freq in Hz
# threshold in dB_SPL
