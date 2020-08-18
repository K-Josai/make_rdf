import numpy as np
import pandas as pd


def pbc_dist(LCEL, xlength):
    L2 = LCEL/2.0
    if xlength >= L2:
        result = xlength-LCEL
    elif xlength < -L2:
        result = xlength+LCEL
    else:
        result = xlength
    return result


def make_g_of_r():
    df = pd.read_csv("final_config.txt",
                     delim_whitespace=True, header=None, skiprows=2)
    df2 = pd.read_csv("final_config.txt", delim_whitespace=True,
                      header=None, skiprows=lambda x: x not in [1])
    NATOM, NRHO, NCYCLE, NPRINT = df2.values[0]
    NATOM = int(NATOM)
    NCYCLE = int(NCYCLE)
    NPRINT = int(NPRINT)
    LCEL = np.sqrt(float(NATOM/NRHO))  # Box Size
    x = np.array(df.values[0])
    y = np.array(df.values[1])

    r_list = []  # list of distances
    for i in range(NATOM):
        for j in range(i+1, NATOM):
            xij = pbc_dist(LCEL, x[i]-x[j])
            yij = pbc_dist(LCEL, y[i]-y[j])
            rij = np.sqrt(xij**2 + yij**2)
            r_list.append(rij)

    k = (1+np.log2(len(r_list)))*400.0  # Sturges' formula
    width = (max(r_list)-min(r_list))/k
    r_floored_list = np.floor(np.array(r_list)/width)
    freq_r = np.array((pd.Series(r_floored_list).value_counts()).sort_index())
    yokojiku = (np.array(range(len(freq_r)))+0.50)*width
    g_of_r = np.zeros(len(freq_r))
    for i in range(len(freq_r)):
        g_of_r[i] = freq_r[i]/(np.pi*yokojiku[i]*width*NRHO*NATOM)

    with open("./g_of_r2.txt", "w") as f4:
        for i in range(len(freq_r)):
            print("{} {}".format(yokojiku[i], g_of_r[i]), file=f4)

    return 0


if __name__ == "__main__":
    make_g_of_r()
