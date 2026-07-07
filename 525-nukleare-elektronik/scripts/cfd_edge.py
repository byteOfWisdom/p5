import std
import propeller as p
import numpy as np

def load_data(side):
    fname = f"../data/fits/{side}/na_cfd.txt"
    data = std.load_csv(fname, skiprows=1)

    filter = np.array(energy[side][element]) == np.array(energy[side][element])
    return (data[1])[filter], np.array(energy[side][element])[filter]

    


if __name__ == "__main__":
    
