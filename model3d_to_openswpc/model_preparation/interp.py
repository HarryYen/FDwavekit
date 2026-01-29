def parallel_create_final_arr(interp_func, x_arr, y_arr, z_arr, params):
    ii, jj, kk = params[0], params[1], params[2]
    value = interp_func([x_arr[ii], y_arr[jj], z_arr[kk]])[0]

    return value
