import numpy as np

# compute an error measure that considers the full time scale.
# if there is an early failure, the error measure will still be returned
# if no progress has been made (or the solution hasn't reached the first
#    reference time step) the returned error will be 100
def mrms_error_V(fname, reference_fname):
    test = np.loadtxt(fname,delimiter=',')
    ref = np.loadtxt(reference_fname,delimiter=',')
    if len(test.shape) == 1:
        return 100
    ref_idx = 0
    err = 0
    err_n = 0
    for line in test:
        # advance until the reference time step is at least as big as the test time step
        while line[0]-1e-6 > ref[ref_idx,0]:
            ref_idx += 1
        # if the time steps are the same, include in the error measurement
        if abs(line[0] - ref[ref_idx][0]) < 1e-6:
            line2 = line[1:len(ref[ref_idx][1:])+1]
            err += sum(((line2 - ref[ref_idx][1:]) / (1 + abs(ref[ref_idx][1:])))**2)
            err_n += len(line2)
    if err_n == 0:
        return 100
    err = (err / err_n)**0.5
    return err

# Computes the error at a time step.
# If the desired time is not in the reference file, returns 100
def error_V(result, reference_fname, t):
    ref = np.loadtxt(reference_fname,delimiter=',')
    ref_idx = 0
    while t - 1e-6 > ref[ref_idx,0]:
        ref_idx += 1
    if abs(t - ref[ref_idx][0]) < 1e-6:
        line = result[1:len(ref[ref_idx][1:])+1]
        err = sum(((line - ref[ref_idx][1:]) / (1 + abs(ref[ref_idx][1:])))**2)
        err = (err / len(line)) ** 0.5
        return err
    return 100

def plot_V_slices(result, Nx, Ny, Nz=7):
    result = result.reshape((-1, Nx, Ny, Nz))[0,:,:,:]
    import matplotlib.pyplot as plt
    from matplotlib import colors

    norm = colors.Normalize(vmin = np.min(result), vmax=np.max(result))
    
    fix, axs = plt.subplots(2,4)
    for i, ax in enumerate(axs.flatten()):
        im=ax.contourf(result[:,:,i],norm=norm,levels=100)
    fig.colorbar(im, ax=axes,orientation='horizontal',fraction=.1)
    
 
