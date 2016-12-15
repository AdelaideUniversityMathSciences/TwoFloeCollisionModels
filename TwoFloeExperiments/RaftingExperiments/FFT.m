% function [a1] = FFT(ydata,sampfreq,f_actual)
%
% LJ YIEW
% Created on  Feb 2014
% Last edited Oct 2016
%
% Calculates the corresponding wave properties of a wave field for a given
% input parameter, according to the dispersion relation.
%
% INPUTS:
%  ydata    = Y data
%  sampfreq = sampling frequency (for Y data)
%  f_actual = estimated actual frequency
%
% OUTPUTS:
%  a1 = [amp1 freq1] % amplitudes and frequencies from FFT

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

function [a1] = FFT(ydata,sampfreq,f_actual)

X1 = ydata;


Fs = sampfreq;
T = 1/Fs;         % Sampling period
L = length(X1);   % Length of signal
t = (1:L)*T;      % Time vector

Y1 = fft(X1);
P = abs(Y1/L);
P1 = P(1:L/2+1);
P1(2:end-1) = 2*P1(2:end-1);

f = Fs*(0:(L/2))/L;

h = figure;
hold on
plot(f,P1,'b')
title('Single-Sided Amplitude Spectrum of X(t)')
xlabel('f (Hz)')
ylabel('|P1(f)|')
xlim([0 10])
  
  
%%
% EXTRACT AMPLITUDE AND FREQUENCY
% clear all 
f = gcf;
axesObjs = get(f, 'Children');
dataObjs = get(axesObjs, 'Children'); 
objTypes = get(dataObjs, 'Type');

frequency = get(dataObjs, 'XData');
amplitude = get(dataObjs, 'YData');

adjustedamp1 = amplitude./abs(frequency - f_actual).^2;

[~,j] = max(adjustedamp1(2:end));

j = j + 1;

amp1  = amplitude;
freq1 = frequency;

a1 = [amp1(j) freq1(j)];

close(h)



