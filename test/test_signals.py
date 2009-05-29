import numpy as np
import matplotlib.pyplot as plt

def main():
    time_function = np.load("time_function.npy")
    forward = np.load("forward.npy")
    ihcrp = np.load("ihcrp.npy")

    # plt.plot(time_function)
    plt.plot(ihcrp[:,70])
    plt.show()



if __name__ == "__main__":
    main()
