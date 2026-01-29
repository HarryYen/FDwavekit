def output(arr, specify_data):
    f = open(f"3d_model_{specify_data}.dat", "w+")
    for ii in range(len(arr)):
        f.write(f"{arr[ii]:.3f}\n")
    f.close()
