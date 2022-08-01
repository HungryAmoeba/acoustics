%%% BOXCAR SIGNAL
% boxcar in the time domain -> frequency domain -> time domain

fs = 1000;
% time betweeen -10 and 10 seconds
t = -10:(1/fs):(10-1/fs);

% let endpoints be [-.5, .5]
left_endpoint = -.5;

%make signal a boxcar of length 1 and given left endpoint
S = heaviside(t - left_endpoint) - heaviside(t-left_endpoint-1);
figure(1); clf;
plot(t, S)
title('Original boxcar signal')
xlabel('Time')
ylim([0 1.2])
n = length(S);
X = fft(S);
f = (0:n-1)*(fs/n);
%figure(2); clf;
%plot(f, real(X));
%figure(2); clf;
%plot(t, sinc(t));
Y = fftshift(X);
yshift = (-n/2:n/2 - 1)*(fs/n);
figure(3); clf;
plot(yshift, real(Y)/fs);
hold on;
plot(yshift, sinc(yshift));
hold on;
plot(yshift, -sinc(yshift));
xlim([-10 10])
title('FFT of boxcar signal')
xlabel("Frequency")
legend('FFT of signal', 'Upper sinc envelope','lower sinc envelope')

figure(4); clf;
S_reverse = ifft(Y);
t_transformed = yshift  * (20/fs); % scaling factor of 20 
% not sure why it's 20, or if it's even correct 
plot(t_transformed, abs(S_reverse))
xlabel("Time")
title('IFFT of transformed signal')

%%

% Linear chirp

t = 0:1/1e3:2;
y = chirp(t,0,1,250);

figure(1); clf;
pspectrum(y,1e3,'spectrogram','TimeResolution',0.1, ...
    'OverlapPercent',99,'Leakage',0.85);

figure(2); clf;
plot(t, y)

figure(3); clf;
plot(t, abs(fft(y)))

f = zeros(500,1);
for ind = 100:200 
    f(ind) = 1;
end

figure(4); clf;
plot(f)

figure(5); clf;
plot(abs(fftshift(fft(f))));

figure(6); clf;
plot(angle(ifft(f)));


%% 

% signal with some sines

Fs = 1000;
T = 1/Fs;
L = 1500;
t = (0:L-1) * T;

S = 0.7*sin(2*pi*50*t) + sin(2*pi*120*t);
X = S + 2*randn(size(t));

figure(1); clf;
plot(1000*t(1:50), X(1:50))

Y = fft(X);

P2 = abs(Y/L);

fs = 100; % sampling freq
t = 0:(1/fs):(10-1/fs);
S = cos(2*pi*15*t);
n = length(S);
X = fft(S);
f = (0:n-1)*(fs/n);
%figure(2); clf;
%plot(t, S);
figure(3); clf;
plot(f, abs(X));

Y = fftshift(X);
fshift = (-n/2:n/2 - 1)*(fs/n);
figure(4); clf;
plot(fshift, abs(Y))



