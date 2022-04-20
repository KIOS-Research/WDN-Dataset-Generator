function pat = demandGeneratorModel(patType)
%% Demand Generator
% Creates a unique demand pattern based on real data

%% Load fourier coeff based on pattern type:
if strcmp(patType,'Residential')
    load('extractedCoefficients/fourier122.mat')
elseif strcmp(patType,'Commercial')
    load('extractedCoefficients/fourier129.mat')
else
    error('Unknown demand pattern type')
end

%% Define times:
lengthDays=7;
tdaySteps = (24*60/tstep);

%% Create yearly component
%%% Create one year pattern with "tstep" (minute) time step:
periodDays=365;
T=periodDays*tdaySteps;
k=1:lengthDays*tdaySteps;

w=2*pi/T; % angular frequency
n=ny; % number of fourier coefficients
Hy=ones(length(k),1);
for i =1:n
    Hy=[Hy sin(i*w.*k)' cos(i*w.*k)'];
end
AyR = Ay; % randomize fourier coefficients
yearOffset = Hy*AyR;

%% Create weekly component
%%% Create one week pattern with "tstep" (minute) time step:
periodDays=7;
T=periodDays*tdaySteps;
k=1:lengthDays*tdaySteps;

w=2*pi/T;
n=nw; % number of fourier coefficients
Hw=ones(length(k),1);
for i =1:n
    Hw=[Hw sin(i*w.*k)' cos(i*w.*k)'];
end
AwR = Aw; % randomize fourier coefficients
weekYearPat = Hw*AwR;

%% Create random component
random = 0; %normally distributed random numbers

%% Create demand pattern
yearOffset = 0;
pat = (yearOffset+1) .* (weekYearPat+1) .* (random+1);

%% Plot demand
% plot(yearOffset)
% hold all
% plot(weekYearPat)
% plot(random)
% 
% figure
% hold all
% plot(pat)
end