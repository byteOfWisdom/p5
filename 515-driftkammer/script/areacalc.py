import propeller as p
import numpy as np

def main():
    a = p.ev(28.5, 0.5)
    b = p.ev(5.0, 0.5)
    a = a*1e-3
    b = b*1e-3
    area = a * b
    rate = 1/(1e6*1) #1 per m^2 per min
    print(rate)
    count =  area / rate # counts m^2 per m^2 per min = counts/min
    count_s = count / 60 # counts/s
    print(count_s.format())


if __name__ == "__main__":
    main()
