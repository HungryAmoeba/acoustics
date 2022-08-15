%%% BOXCAR SIGNAL

% do the check with the point source in homogeneous media

% zero order hankel function of the first kind 

close all; clear all;
path(path, '/Users/charlesxu/Documents/WHOI_2022/reflection_coeff/subroutines')

% boxcar in the time domain -> frequency domain -> time domain
fs = 1000;
n = 20000;
T = n/fs;
df = 1/(n/fs);
t = [0:n-1]/fs; 
f = [0:ceil(n/2)-1 ceil(-n/2):-1]*df;

% source function as ricker wavelet 

% convolve with data
% sometimes impulsive sources, if speakers then any arbitrary func



%%  frequency to time domain

mode = 3; % 1 for 5 layered, 2 for tunneling, 3 for 3 layered

f0 = 100;
f1 = 300;
IDF = find(f>=f0 & f<=f1);

spect = zeros(size(f)); 
window_1_endpoints = find(f==f0 | f == f1);
window_2_endpoints = find(f == -f0 | f == -f1);
window = ones(size(f));
window(window_1_endpoints(1):window_1_endpoints(2)) = tukeywin(diff(window_1_endpoints) + 1, 1);
window(window_2_endpoints(1):window_2_endpoints(2)) = tukeywin(diff(window_2_endpoints) + 1, 1);

% define the medium properties

% Attenuation by introducing complex sound speed

% lambda = c/(2 * pi * f)
% alpha = 20 * log10(exp(.3))./ (c/ (2 * pi * freq) )


if mode == 1
    rho = [1;1.5;1.7;1.9;2];
    %c = [1500 - 2i; 1500 - 1i; 1575- 3i; 1650- 4i; 1800- 1i];
    c = [1500; 1500; 1575; 1650; 1800];
    h = [200;300;80];
end

if mode == 2
    %tunneling
    c_0 = 1000; c_1 = 2000;
    rho_0 = 1; rho_1 = 1.5;

    c = [c_0, c_1, c_0]; rho = [rho_0, rho_1, rho_0];

    h=10;
end

if mode == 3
    c = [1500, 1700,  1800];
    rho = [1000, 1500, 2000];
    h_arr = 5:160:310;
    h = 200;
end


%just pick some random angle for now
angle_arr = 3 * pi/10;%0:pi/10:pi/2;

%pick some random x_0
x0_arr = 100; %0:100:200;
counter = 1;

%Refl_arr = zeros(length(angle_arr), length(IDF));
for ind = 1%:length(h_arr)
    for angle = angle_arr
        for x0 = x0_arr
            for idf = IDF
                %h = h_arr(ind);
                freq = f(idf);    
                
                % if attenuation is desired
                lambda = c/(freq);
                beta = .1; % typical is .1 - .3
                alpha = beta./lambda * log(10)/20; % imaginary part of k
                k = 2*pi*freq./c + 1i .* alpha;
                
                % Ref Coef Computation at freq
                Refl = rayleigh_strat(rho, c, h, angle, freq);
                Refl = Refl(1);
                k = 2*pi*freq/c(1);

                % multiply by exp(ikx_0) to get p_r
                spect(idf) = Refl*exp(1i*k*x0);

                idf_neg = find(f == -freq);

                spect(idf_neg) = conj(spect(idf)); 
            end
            windowed_spect = spect .* window;
            Ref_time = fft(windowed_spect);
            figure(2 * counter); clf;
            plot(f, real(windowed_spect))
            title(x0)
            figure(2* ind + 2 * counter + 1); clf;
            %figure(counter); 
            plot(t-max(t)/2, abs(fftshift(Ref_time)))
            title("x0 = " + x0 + ", angle = " + angle * 180/pi + ", h = " + h)
            xlim([0 1.4])
            counter = counter + 1;
            %hold on;
        end
    end
end

%% test a chirp

y = chirp(t,10,15,20);
fft_y = fft(y);
figure(82); clf
ax(1) = subplot(211);
plot(f,abs(fft_y)); xlabel('freq (hz)'); ylabel('magnitude')
ax(2) = subplot(212);

plot(f,angle(fft_y)/pi*180); xlabel('freq (hz)'); ylabel('phase (deg)')

linkaxes(ax,'x')
%% 

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



