%% Clear workspace, close figures, and clear command window
clear; close all; clc;

%% Step 1: Generate Synthetic Data for Electrical Loads
numDataPoints = 101;  % Total days (data points)

% Electrical load parameters
elec_baseline = 1.49;        % Initial electrical load
elec_slope = 0.0012;         % Daily increase for electrical load
elec_noise_amp = 0.0005;     % Noise amplitude for electrical load

% Time vector (days 0 to 100)
time = (0:numDataPoints-1)';

% Generate synthetic data with a linear trend and noise
elec_data = elec_baseline + elec_slope * time + elec_noise_amp * randn(numDataPoints,1);

%% Step 2: Prepare Data for Training (1D input: only electrical load)
% Use first 100 days of electrical load as input
inputSequence = elec_data(1:100)';   % Size: [1 x 100]

% The target is the electrical load on day 101
targetValue = elec_data(101);

% MATLAB's trainNetwork expects a cell array of sequences
XTrain = {inputSequence};  % Each cell contains a [1 x 100] sequence
YTrain = targetValue;

%% Step 3: Define the LSTM Network Architecture (1 input feature)
layers = [ ...
    sequenceInputLayer(1)               % Now only 1 feature per time step
    lstmLayer(50, 'OutputMode', 'last') % LSTM with 50 hidden units
    dropoutLayer(0.2)                   % Dropout for regularization
    fullyConnectedLayer(1)             % Output a single value
    regressionLayer                    % For continuous output
];

%% Step 4: Specify Training Options
options = trainingOptions('adam', ...
    'MaxEpochs', 250, ...
    'GradientThreshold', 1, ...
    'InitialLearnRate', 0.005, ...
    'Verbose', 0, ...
    'Plots', 'training-progress');

%% Step 5: Train the LSTM Network
net = trainNetwork(XTrain, YTrain, layers, options);

%% Step 6: Predict the 101st Electrical Load Using the Trained Model
YPred = predict(net, XTrain);

%% Display the Results
fprintf('Predicted 101st Electrical Load: %f\n', YPred);
fprintf('True 101st Electrical Load: %f\n', YTrain);
