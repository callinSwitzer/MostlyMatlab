%% Sine wave output for shaker
amplitude = 0.3; % Amplitude for shaking in Volts (1 is Max)
frequency = 280; % in HZ
seconds = .5; % length of signal in seconds



% Output frequency signal to shaker

s = daq.createSession('ni');
s.addAnalogOutputChannel('Dev2','ao0','Voltage');
s.Rate = 100000;

% Generate output signal
outputSignal =  amplitude * sin(linspace(0, pi*2* frequency, s.Rate)');
outputSignal = outputSignal(1:seconds*length(outputSignal));

% plot(outputSignal); xlabel('Time'); ylabel('Voltage');

%
s.queueOutputData([outputSignal]);
s.startForeground();

%% amplitude sweep at 280 Hz
samps = randsample(3, 20, true);

amplitude = 0.406; 

for ii = 1:20
    
    
    amplitude = amplitude ; % Amplitude for shaking in Volts (1 is Max)
    ([amplitude, ii])
    frequency = 280; % in HZ
    seconds = 1; % length of signal in seconds
    
    % Output frequency signal to shaker
    s = daq.createSession('ni');
    s.addAnalogOutputChannel('Dev2','ao0','Voltage');
    s.Rate = 100000;

    % Generate output signal
    outputSignal =  amplitude * sin(linspace(0, pi*2* frequency, s.Rate)');
    outputSignal = outputSignal(1:seconds*length(outputSignal));
    s.queueOutputData([outputSignal]);
    s.startForeground();
end
