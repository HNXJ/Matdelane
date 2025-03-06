function flipObj = jVFLIP(x)

    data = permute(x, [2, 3, 1]);
    flipObj = vFLIP2(data(:, :, :), 'DataType', 'raw_cut', 'fsample', 1000, 'intdist', 0.04, 'plot_result', true);

end