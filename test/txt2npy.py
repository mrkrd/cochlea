#!/usr/bin/env python

import sys
import os.path
import numpy as np


def main():
    for file_name in sys.argv[1:]:
        print file_name

        base_name, ext_name = os.path.splitext(file_name)

        if ext_name != '.txt':
            assert False


        tab = np.loadtxt(file_name)

        np.save(base_name, tab)




if __name__ == "__main__":
    main()
