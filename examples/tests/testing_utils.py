import numpy as np

def output(verbose, error_threshold, solution, result, name):
    if verbose:
        print("{:<40} {} with error {}".format(name+" solution", solution, result-solution))
        #print(name, 'solution\t\t', solution, 'with error', result-solution)
    if max(abs((result - solution)/solution)) > error_threshold:
        print("error using", name, "larger than expected", abs((result-solution)/solution))

def plot(reference_fname, fname_labels, title=None, show=True):
    from matplotlib import pyplot as plt
    plt.figure()
    data = np.loadtxt(reference_fname,delimiter=',')
    plt.plot(data[:,0], data[:,1:], label = 'Exact', linestyle='solid')

    linestyles=['dashed', 'dotted', 'dashdot']
    linestyle_idx = 0

    for fname in fname_labels:
        try:
            data = np.loadtxt(fname,delimiter=',')
        except:
            data = np.loadtxt(fname,delimiter=',', dtype=np.complex128)
        plt.plot(data[:,0], data[:,1:], label=fname_labels[fname], linestyle=linestyles[linestyle_idx])
        linestyle_idx += 1
        linestyle_idx = linestyle_idx % len(linestyles)
    plt.legend()
    if title is not None:
        plt.title(title)
    if show:
        plt.show()
