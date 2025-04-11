%% Clear workspace, close figures, and clear command window
clear; close all; clc;

%% Step 1: Generate Synthetic Data
% We want 101 data points over 100 days where:
%   - The first (day 1) value is 5.38.
%   - The values increase roughly linearly over time.
%   - Some random noise is added to simulate error.
numDataPoints = 101;
baseline = 5.38;       % Initial value at day 1
slope = 0.001;           % Linear increase per day (adjust as needed)
noise_amp = 0.0005;      % Noise amplitude

% Generate the data: 
% For day index 0 to 100, the ideal value is baseline + slope*day, plus noise.
data = baseline + slope * (0:numDataPoints-1)' + noise_amp * randn(numDataPoints, 1);

%% Step 2: Prepare Data for Training
% Use the first 100 values as the input sequence (one feature per time step)
% and the 101st value as the target.
inputSequence = data(1:100)';    % Transpose to get a [1 x 100] row vector
targetValue = data(101);         % The target is the 101st valuet

% MATLAB’s trainNetwork requires the inputs as a cell array of sequences.
XTrain = {inputSequence};        % Cell array with one sequence
YTrain = targetValue;            % A scalar target

%% Step 3: Define the LSTM Network Architecture
% The network will process a sequence (100 time steps, 1 feature each) and output
% the prediction for the next time step.
layers = [ ...
    sequenceInputLayer(1)               % One feature per time step
    lstmLayer(50, 'OutputMode', 'last')   % LSTM layer with 50 hidden units; output only last time step
    dropoutLayer(0.2)                   % Dropout layer with 20% dropout for regularization
    fullyConnectedLayer(1)              % Fully connected layer to map LSTM output to one prediction
    regressionLayer                     % Regression layer for continuous output
];

%% Step 4: Specify Training Options
options = trainingOptions('adam', ...
    'MaxEpochs', 250, ...                % Maximum number of training epochs
    'GradientThreshold', 1, ...          % Gradient threshold to prevent exploding gradients
    'InitialLearnRate', 0.005, ...       % Initial learning rate
    'Verbose', 0, ...                    % Suppress verbose output
    'Plots', 'training-progress');       % Display training progress plot

%% Step 5: Train the LSTM Network
% Note: Training with a single sample is only for demonstration.
net = trainNetwork(XTrain, YTrain, layers, options);

%% Step 6: Predict the 101st Value Using the Trained Model
YPred = predict(net, XTrain);

%% Display the Results
fprintf('Predicted 101st value: %f\n', YPred);
fprintf('True 101st value: %f\n', YTrain);
