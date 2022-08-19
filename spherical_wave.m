addpath subroutines

r_arr = 1:1000;
z = 20;
z_0 = 8;
p_arr  =zeros(size(r_arr));

for ind = 1:length(r_arr)
    p_arr(ind) = p_hankel(r_arr(ind), z, z_0);
end

figure(1); clf;
plot(r_arr, real(p_arr));
grid on;
xlabel("r")
ylabel("real p(r)")

figure(2); clf;
subplot(2,1,1);
plot(r_arr, abs(p_arr));
grid on;
xlabel("r")
ylabel("abs(p(r))")

subplot(2,1,2);
plot(r_arr, angle(p_arr))
xlabel("r")
ylabel("angle(p(r))")

%%

% do a similar experiment as the fft practice file

z = 20; z_0 = 8; %fix these for now

fs = 1000;
n = 20000;
T = n/fs;
df = 1/(n/fs);
t = [0:n-1]/fs; 
f = [0:ceil(n/2)-1 ceil(-n/2):-1]*df;

f0 = 100;
f1 = 300;
IDF = find(f>=f0 & f<=f1);

spect = zeros(size(f)); 
window_1_endpoints = find(f==f0 | f == f1);
window_2_endpoints = find(f == -f0 | f == -f1);
window = ones(size(f));
window(window_1_endpoints(1):window_1_endpoints(2)) = tukeywin(diff(window_1_endpoints) + 1, 1);
window(window_2_endpoints(1):window_2_endpoints(2)) = tukeywin(diff(window_2_endpoints) + 1, 1);
beta = [.0001,.1, .2];
%just pick some random angle for now
angle_arr = 3 * pi/10;%0:pi/10:pi/2;

%pick some random x_0
x0_arr = 100; %0:100:200;

for ind = 1%:length(h_arr)
    for theta = angle_arr
        for x0 = x0_arr
            for idf = IDF
                %h = h_arr(ind);
                freq = f(idf);    
                
%                 % if attenuation is desired
%                 lambda = c/(freq);
%                 beta = [.0001,.1, .2]; % typical is .1 - .5
%                 alpha = beta./lambda * log(10)/20; % imaginary part of k
%                 k = 2*pi*freq./c + 1i .* alpha;
                
%                 % Ref Coef Computation at freq
%                 Refl = rayleigh_strat(rho, c, h, angle, freq,k); %rayleigh_strat(rho,c,h,angle,freq) if no attenuation
%                 Refl = Refl(1);
%                 k = 2*pi*freq/c(1);

                % multiply by exp(ikx_0) to get p_r
                h = 200;
                spect(idf) = p_hankel(x0, z, z_0, freq);

                idf_neg = find(f == -freq);

                spect(idf_neg) = conj(spect(idf)); 
            end
            windowed_spect = spect .* window;
            Ref_time = fft(windowed_spect);
            figure(1); clf;
            plot(f, real(windowed_spect))
            title(x0)
            figure(2); clf;
            %figure(counter); 
            plot(t-max(t)/2, abs(fftshift(Ref_time)))
            title("x0 = " + x0 + ", angle = " + theta * 180/pi + ", h = " + h + ", beta = [" + beta(1) +"," + beta(2) + "," + beta(3) + "]")
            %fix the title later lol
            xlim([0 1.4])
            %hold on;
        end
    end
end


%%

%UAV simulation
z = 5; z_0 = 5;
r = sqrt(5^2 + 10^2);
h_arr = [100; 150; 100; 80; 120];
c_arr = [1500, 1500, 1575;
    1500, 1550, 1650;
    1500, 1600, 1700;
    1480, 1590, 1550;
    1520, 1580, 1690];
rho_arr = [1000, 1500, 2000; 
    1000, 1600, 2000;
    1000, 1500, 2000;
    1000, 1400, 1900;
    1000, 1500, 2000];

num_positions = length(h_arr);

% set up signal
for d=1
    fs = 1000;
    n = 20000;
    T = n/fs;
    df = 1/(n/fs);
    t = [0:n-1]/fs; 
    f = [0:ceil(n/2)-1 ceil(-n/2):-1]*df;

    f0 = 100;
    f1 = 300;
    IDF = find(f>=f0 & f<=f1);

    spect = zeros(num_positions, length(f)); 
    window_1_endpoints = find(f==f0 | f == f1);
    window_2_endpoints = find(f == -f0 | f == -f1);
    window = ones(size(f));
    window(window_1_endpoints(1):window_1_endpoints(2)) = tukeywin(diff(window_1_endpoints) + 1, 1);
    window(window_2_endpoints(1):window_2_endpoints(2)) = tukeywin(diff(window_2_endpoints) + 1, 1);
    beta = [.0001,.1, .2];
end

fig_counter = 1;
% simulate returns at each position
for pos = 1:num_positions
    h = h_arr(pos);
    rho = rho_arr(pos, :);
    c = c_arr(pos, :);
    for idf = IDF
        %h = h_arr(ind);
        freq = f(idf); 
        spect(pos, idf) = p_hankel(r, z, z_0, freq, c, rho, h, beta);
        idf_neg = find(f == -freq);
        spect(pos, idf_neg) = conj(spect(idf)); 
    end
    windowed_spect = spect(pos, :) .* window;
    Ref_time = fft(windowed_spect);
    figure(fig_counter); clf;
    fig_counter = fig_counter + 1;
    plot(f, real(windowed_spect))
    title("UAV position " + pos)
    figure(fig_counter); clf;
    fig_counter = fig_counter + 1;
    plot(t-max(t)/2, abs(fftshift(Ref_time)))
    title("UAV position " + pos)
    %fix the title later lol
    xlim([-.05 .4])
end


%%
% compare this with the plane wave solution.
theta = atan(2);

% set up signal
for d=1
    fs = 1000;
    n = 20000;
    T = n/fs;
    df = 1/(n/fs);
    t = [0:n-1]/fs; 
    f = [0:ceil(n/2)-1 ceil(-n/2):-1]*df;

    f0 = 100;
    f1 = 400;
    IDF = find(f>=f0 & f<=f1);

    spect_plane = zeros(size(f)); 
    spect_spherical = zeros(size(f));
    window_1_endpoints = find(f==f0 | f == f1);
    window_2_endpoints = find(f == -f0 | f == -f1);
    window = ones(size(f));
    window(window_1_endpoints(1):window_1_endpoints(2)) = tukeywin(diff(window_1_endpoints) + 1, 1);
    window(window_2_endpoints(1):window_2_endpoints(2)) = tukeywin(diff(window_2_endpoints) + 1, 1);
end

%define medium properties
mode = 3;
if mode == 1
    rho = [1,1.5,1.7,1.9,2];
    %c = [1500 - 2i; 1500 - 1i; 1575- 3i; 1650- 4i; 1800- 1i];
    c = [1500, 1500, 1575, 1650, 1800];
    h = [200,300,80];
    beta = [.0001, .1, .15, .2, .1];
end

if mode == 2
    %tunneling
    c_0 = 1000; c_1 = 2000;
    rho_0 = 1; rho_1 = 1.5;

    c = [c_0, c_1, c_0]; rho = [rho_0, rho_1, rho_0];
    beta = [.0001, .1, .0001];

    h=10;
end

if mode == 3
    c = [1500, 1700,  1800];
    rho = [1000, 1500, 2000];
    h = 200;
    beta = [.0001, .1, .15];
end

% c = [1500, 1700,  1800];
% rho = [1000, 1500, 2000];
% h = 500;
z = 50; z_0 = 50;
r = sqrt(50^2 + 100^2);

fig_counter = 1;


for idf = IDF
    %h = h_arr(ind);
    freq = f(idf); 
    spect_spherical(idf) = p_hankel(r, z, z_0, freq, c, rho, h, beta);
    lambda = c./(freq);
    % make sure they are going in the same direction!!!
    alpha = beta./lambda * log(10)/20; % imaginary part of k
    k = 2*pi*freq./c + 1i .* alpha;
    
    %Refl = rayleigh_strat(rho, c, h, theta, freq,k); %rayleigh_strat(rho,c,h,angle,freq) if no attenuation
    kx = k .* cos(theta);
    Refl = reflection_coeff(kx, rho, k, h);
    Refl = Refl(1);

    spect_plane(idf) = Refl*exp(1i*k(1)*r);
    idf_neg = find(f == -freq);
    spect_spherical(idf_neg) = conj(spect_spherical(idf)); 
    spect_plane(idf_neg) = conj(spect_plane(idf));
end

%%
fig_counter = 1;

% Generate plots for spherical wave
windowed_spect = spect_spherical .* window;
Ref_time = fft(windowed_spect);
figure(fig_counter); clf;
fig_counter = fig_counter + 1;
plot(f, real(windowed_spect))
title("Real part of windowed spectrum spherical, r = " + r)
xlabel("freq")
ylabel("p(freq)")
figure(fig_counter); clf;
fig_counter = fig_counter + 1;
%figure(counter); 
plot(t-max(t)/2, abs(fftshift(Ref_time)))
title("Time domain solution spherical, r = " + r)
xlabel("t")
ylabel("p(t)")
xlim([0 1.4])

% Generate plots for plane wave
windowed_spect = spect_plane .* window;
Ref_time = fft(windowed_spect);
figure(fig_counter); clf;
fig_counter = fig_counter + 1;
plot(f, real(windowed_spect))
title("Real part of windowed spectrum plane, r = " + r)
xlabel("freq")
ylabel("p(freq)")
figure(fig_counter); clf;
fig_counter = fig_counter + 1;
plot(t-max(t)/2, abs(fftshift(Ref_time)))
title("Time domain solution plane, r = " + r)
xlabel("t")
ylabel("p(t)")
xlim([0 1.4])
fig_counter = fig_counter + 1;
