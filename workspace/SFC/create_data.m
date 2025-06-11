function [lfp, dataspike] = create_data()

% Define parameters
f = 40; % Oscillation frequency in Hz
eigenvalue = 0.94; % Eigenvalue for damping - stable until eigenvalue=1
Fs = 1000; % Sampling rate in Hz
T  = 100; % Number of seconds
N  = T*Fs; 

% Convert the frequency to radians per sample
omega = 2 * pi * f / Fs;

% Calculate the poles using the damping factor and frequency
r = eigenvalue; % The magnitude of the poles
theta = omega; % The angle of the poles

% Calculate AR coefficients using the pole representation
a2 = -r^2; % modulus 
a1 = 4*a2*cos(theta) / (a2-1); 

% Generate a white noise sequence
w = randn(N, 1); % Standard normal white noise

% Pre-allocate the signal array
x = zeros(N, 1);

b = 1;  % 
a = [1, -a1, -a2];  % MATLAB expects the feedback coefficients with opposite signs

% innovation
x = randn(1,N); % Random initial condition for x[1]
x = filter(b,a,x); 
x = x-min(x); 
x = x./max(x);  

% restructure the data 
nTrials = 100; 
xdata = reshape(x,[Fs 100])'; 
for iTrial = 1:nTrials
    lfp.trial{iTrial} = xdata(iTrial,:); 
    lfp.time{iTrial} = [0:Fs-1]/Fs;
    lfp.label{1} = 'chan1'; 
    lfp.fsample = Fs; 
end
lfp = ft_checkdata(lfp,'datatype', 'raw', 'feedback', 'yes');

% simulate 100 neurons
nNeurons = 100; 
firing_rates = exprnd(5,[1 nNeurons]);

% generate spikes 
for iNeuron = 1:nNeurons
    Fr = firing_rates(iNeuron); 
    p = (1./Fs)*Fr*x;
    s = 0; 
    while s<3 % make sure the neuron has at least 3 spikes
        r = rand(1,N); 
        S = double(p>r); 
        s = sum(S(:));
    end

    % create a similar spike structure
    nTrials = 100; 
    S = reshape(S,[Fs 100])'; 
    for iTrial = 1:nTrials
        dataspike.trial{iTrial}(iNeuron,:) = S(iTrial,:); 
        dataspike.time{iTrial} = [0:Fs-1]/Fs;
        dataspike.label{iNeuron} = ['neuron' num2str(iNeuron)]; 
        dataspike.fsample = Fs; 
    end
end