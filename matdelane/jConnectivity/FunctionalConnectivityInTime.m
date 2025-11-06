function [functional_connectivity_matrix] = FunctionalConnectivityInTime(X_signal, Method, WindowSize, Overlap)
% FUNCTIONALCONNECTIVITYINTIME Calculates time-resolved functional connectivity.
%
% This function applies a sliding window approach to calculate pairwise
% functional connectivity matrices across multiple time bins.
%
% Inputs:
%   X_signal : Time series array/matrix. Expected size M x N x T, where
%              M is the number of observations/trials, N is the number of
%              neurons/channels/ROIs, and T is the total number of time points.
%              (Dimensions: $M \times N \times T$)
%   Method   : Connectivity method. Currently supported: 'Correlation'.
%              Placeholders included for 'Granger' and 'MutualInfo'.
%   WindowSize: Moving window's length (in time points, T).
%   Overlap  : Overlap between the moving windows (in time points).
%
% Outputs:
%   functional_connectivity_matrix : Connectivity matrices in time bins.
%              Size is N x N x TimeBins.
%
% Note on 'Granger' and 'MutualInfo': Robust implementation of these methods
% typically requires external MATLAB toolboxes (e.g., MVGC toolbox for
% Granger Causality). The case statements below serve as structural
% placeholders.

% --- 1. Input Validation and Setup ---
[M, N, T] = size(X_signal);

if N < 2
    error('X_signal must contain at least 2 neurons (N > 1) for connectivity analysis.');
end
if WindowSize > T
    error('WindowSize cannot be larger than the total number of time points T.');
end

% Ensure the method input is case-insensitive
Method = lower(Method);

% --- 2. Window Index Calculation ---
StepSize = WindowSize - Overlap;
if StepSize <= 0
    error('Overlap must be less than WindowSize to ensure progression.');
end

% Calculate starting points for each window
StartPoints = 1:StepSize:(T - WindowSize + 1);

% Calculate the total number of time bins
TimeBins = length(StartPoints);

% Pre-allocate the output matrix
functional_connectivity_matrix = zeros(N, N, TimeBins);

fprintf('Calculating functional connectivity for %d time bins using %s method.\n', TimeBins, Method);
fprintf('Window Size: %d, Overlap: %d, Step Size: %d\n', WindowSize, Overlap, StepSize);

% --- 3. Main Sliding Window Loop ---
for k = 1:TimeBins
    % Determine the start and end time indices for the current window
    start_idx = StartPoints(k);
    end_idx = start_idx + WindowSize - 1;

    % Extract the data for the current window
    % Size is M (observations) x N (neurons) x WindowSize (time points)
    X_window = X_signal(:, :, start_idx:end_idx);

    % Flatten the windowed data: Combine M and WindowSize into one effective
    % time dimension. This creates a matrix of size (M * WindowSize) x N.
    % This is crucial for calculating a single correlation coefficient
    % across all available data points (M observations * W timepoints) for each pair of neurons.
    X_flat = reshape(X_window, M * WindowSize, N);

    % Calculate the connectivity matrix based on the selected method
    switch Method
        case 'correlation'
            % The corrcoef function computes the Pearson correlation matrix.
            % R is an N x N matrix, where R(i, j) is the correlation
            % between neuron i and neuron j across all M*WindowSize data points.
            R = corrcoef(X_flat);
            
            % Store the result
            functional_connectivity_matrix(:, :, k) = R;

        case 'granger'
            % --- Placeholder for Granger Causality ---
            % Requires specialized code, typically VAR model fitting and
            % application of the MVGC toolbox (or similar).
            % Example steps would involve:
            % 1. Fit a Vector Autoregressive (VAR) model to X_flat.
            % 2. Calculate the Pairwise Conditional Granger Causality.
            
            % Since external toolboxes are required, we skip the calculation
            % and store a placeholder matrix (e.g., NaNs) to prevent errors.
            % If you integrate a toolbox, replace this block with the
            % actual Granger calculation.
            warning('Granger Causality selected. Requires external toolboxes (e.g., MVGC) for robust calculation. Storing NaN placeholder.');
            functional_connectivity_matrix(:, :, k) = NaN(N, N);

        case 'mutualinfo'
            % --- Placeholder for Mutual Information ---
            % Requires estimation of joint and marginal probability density
            % functions, typically using k-nearest neighbor methods (e.g., Kraskov estimator).
            % The code requires an Information Theory toolbox.
            warning('Mutual Information selected. Requires external toolboxes (e.g., Kraskov estimator) for robust calculation. Storing NaN placeholder.');
            functional_connectivity_matrix(:, :, k) = NaN(N, N);

        otherwise
            error('Unsupported connectivity method: %s. Supported methods are ''Correlation'', ''Granger'', ''MutualInfo''.', Method);
    end
end

% Optional: Display a success message
disp('Functional connectivity calculation complete.');

end
