function mi = MutualInformation(X, Y, num_bins)
% CALCULATE_MUTUAL_INFORMATION Calculates the mutual information between two time series.
%
%   MI = CALCULATE_MUTUAL_INFORMATION(X, Y, NUM_BINS)
%
%   Inputs:
%       X: A time series vector (column or row vector).
%       Y: A second time series vector (must be the same size as X).
%       NUM_BINS: The number of bins to use for discretization (scalar integer).
%
%   Output:
%       MI: The mutual information in bits.
%
%   Note: This implementation uses a histogram-based approach, which is
%   sensitive to the choice of NUM_BINS. For more advanced or robust
%   methods, consider k-nearest neighbor (k-NN) estimators.

    % 1. Input Validation and Reshaping
    if length(X) ~= length(Y)
        error('Time series X and Y must have the same length.');
    end
    X = X(:); % Ensure X is a column vector
    Y = Y(:); % Ensure Y is a column vector
    
    % 2. Calculate Joint Histogram (2D)
    % histcounts2 returns N, the count matrix. We also get the bin edges.
    [N_joint, ~, ~] = histcounts2(X, Y, num_bins);
    
    % 3. Calculate Joint Probability Distribution P(x,y)
    P_joint = N_joint / numel(X);
    
    % 4. Calculate Marginal Probability Distributions P(x) and P(y)
    P_x = sum(P_joint, 2); % Sum across columns to get marginal for X
    P_y = sum(P_joint, 1); % Sum across rows to get marginal for Y
    
    % 5. Calculate Mutual Information
    
    % Find indices of non-zero joint probabilities to avoid log(0)
    non_zero_indices = P_joint > 0;
    
    % Initialize MI
    mi = 0;
    
    % Calculate MI contribution for each bin
    for i = 1:num_bins
        for j = 1:num_bins
            p_xy = P_joint(i, j);
            if p_xy > 0
                p_x = P_x(i);
                p_y = P_y(j);
                
                % Contribution to MI: p(x,y) * log2( p(x,y) / (p(x)*p(y)) )
                mi = mi + p_xy * log2( p_xy / (p_x * p_y) );
            end
        end
    end
    
    % Optional: Using array operations for a more compact (but maybe less readable) MI calculation:
    % P_x_grid = repmat(P_x, 1, num_bins);
    % P_y_grid = repmat(P_y, num_bins, 1);
    % P_prod = P_x_grid .* P_y_grid;
    % mi_array = P_joint .* log2(P_joint ./ P_prod);
    % mi = sum(mi_array(non_zero_indices));
end