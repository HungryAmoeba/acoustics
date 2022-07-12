mud_speed = 1480;
mud_density = 1.5; % relative to water
sand_speed = 1750;
sand_density = 1.7; % relative to water
water_speed = 1500;
water_density = 1; % reference point

%% 

% mud over sand base

close all;

c = [water_speed, mud_speed, sand_speed];
rho = [water_density, mud_density, sand_density];
h_arr = 2:20:62; % vary this parameter

angle_arr = 0:.01:pi/2;
freq_arr = 5:10:500;
R_arr = zeros(length(h_arr), length(angle_arr), length(freq_arr)); 


for freq_index = 1:length(freq_arr)
    for angle_index = 1:length(angle_arr)
        for h_index = 1:length(h_arr)
            R_array = rayleigh_strat(rho,c,h_arr(h_index),angle_arr(angle_index),freq_arr(freq_index));
            R_arr(h_index, angle_index, freq_index) = R_array(1);
        end
    end
end

[X,Y] = meshgrid(freq_arr, angle_arr);
% figure(18); clf;
for ind = 1:length(h_arr)
    figure(ind)
    %surf(X, Y * 180/pi, abs(squeeze(R_arr(ind,:,:))));
    imagesc([0 90], [5 500], abs(squeeze(R_arr(ind,:,:)))');
    title("|R| for mud h = " + h_arr(ind) + " over sand")
    ylabel('Frequency')
    xlabel('Angle')
    % saveas(gcf, "mud_sand_" + h_arr(ind) + ".png")
end
% 
% surf(X, Y * 180/pi, abs(squeeze(R_arr(2,:,:))))
% 
% title('|R| for mud-sand')
% xlabel('Frequency')
% ylabel('Angle')

% now fix some angle and 

%% 

% sand-mud-sand base

close all;

n_1 = water_speed/sand_speed;
m_1 = sand_density/water_density;
brewster_angle = asin(sqrt((m_1^2 - n_1^2)/(m_1^2 - 1)));


c = [water_speed, sand_speed, mud_speed, sand_speed];
rho = [water_density, sand_density, mud_density, sand_density];
h = [600,3]; % vary this parameter
mud_thickness = 5; % also vary this parameter
h_arr = ones(5, 2) * mud_thickness;
max_sand_thickness = 500;
for ind = 1:length(h_arr)
    h_arr(ind,1) = 1 + (ind-1) * max_sand_thickness/length(h_arr);
end

angle_arr = 0:.02:pi/2;
freq_arr = 50:10:500;
R_arr = zeros(length(angle_arr), length(freq_arr),length(h_arr)); 

for freq_index = 1:length(freq_arr)
    for angle_index = 1:length(angle_arr)
        for h_index = 1:length(h_arr)
            R_array = rayleigh_strat(rho,c,h_arr(h_index, :),angle_arr(angle_index),freq_arr(freq_index));
            R_arr(angle_index, freq_index, h_index) = R_array(1);   
        end
    end
end


for h_index = 1:length(h_arr)
    %[X,Y] = meshgrid(freq_arr, angle_arr);
    figure(h_index); clf;
    %surf(X, Y * 180/pi, abs(R_arr(:,:,h_index)));
    imagesc([0 90], [5 500], abs(R_arr(:,:,h_index))');
    title("|R| for sand-mud-sand, top sand h = " + h_arr(h_index) +", mud h = " + mud_thickness)
    ylabel("Frequency")
    xlabel("Angle")
    colorbar
    saveas(gcf, "sand_mud_" + h_arr(h_index,1) + "sand.png")
end
% 
% [X,Y] = meshgrid(freq_arr, angle_arr);
% figure(19); clf;
% 
% surf(X, Y * 180/pi, abs(R_arr))
% 
% title('|R| for sand-mud-sand')
% xlabel('Frequency')
% ylabel('Angle')
% 
% figure(20); clf;
% 
% surf(X, Y * 180/pi, angle(R_arr))
% 
% title("Phase for sand-mud-sand")
% xlabel("frequency")
% ylabel("Angle")

%% 

%water-sand-mud

mud_speed = 1480;
mud_density = 1.5; % relative to water
sand_speed = 1750;
sand_density = 1.7; % relative to water
water_speed = 1500;
water_density = 1; % reference point

n_1 = water_speed/sand_speed;
m_1 = sand_density/water_density;
brewster_angle = asin(sqrt((m_1^2 - n_1^2)/(m_1^2 - 1)));

% mud over sand base
c = [water_speed, sand_speed, mud_speed];
rho = [water_density, sand_density, mud_density];
h = 20; % vary this parameter

angle_arr = 0:.01:pi/2;
freq_arr = 50:10:500;
R_arr = zeros(length(angle_arr), length(freq_arr)); 

for freq_index = 1:length(freq_arr)
    for angle_index = 1:length(angle_arr)
        R_array = rayleigh_strat(rho,c,h,angle_arr(angle_index),freq_arr(freq_index));
        R_arr(angle_index, freq_index) = R_array(1);      
    end
end

[X,Y] = meshgrid(freq_arr, angle_arr);
figure(21); clf;

surf(X, Y * 180/pi, abs(R_arr))

title('|R| for sand-mud')
xlabel('Frequency')
ylabel('Angle')

figure(22); clf;

surf(X, Y * 180/pi, angle(R_arr))

title("Phase for sand-mud")
xlabel("frequency")
ylabel("Angle")

%% 

mud_speed = 1480;
mud_density = 1.5; % relative to water
sand_speed = 1750;
sand_density = 1.7; % relative to water
water_speed = 1500;
water_density = 1; % reference point

% mud over sand base
c = [water_speed, sand_speed, water_speed];
rho = [water_density, sand_density, water_density];
h_arr = 5:5:20; % vary this parameter

angle_arr = 0:.01:pi/2;
freq_arr = 50:10:500;
R_arr = zeros(length(h_arr), length(angle_arr), length(freq_arr)); 


for freq_index = 1:length(freq_arr)
    for angle_index = 1:length(angle_arr)
        for h_index = 1:length(h_arr)
            R_array = rayleigh_strat(rho,c,h_arr(h_index),angle_arr(angle_index),freq_arr(freq_index));
            R_arr(h_index, angle_index, freq_index) = R_array(1);
        end
    end
end

[X,Y] = meshgrid(freq_arr, angle_arr);
% figure(18); clf;
for ind = 1:length(h_arr)
    figure(ind)
    surf(X, Y * 180/pi, abs(squeeze(R_arr(ind,:,:))));
    title(h_arr(ind))
    xlabel('Frequency')
    ylabel('Angle')

end
%surf(X, Y * 180/pi, abs(squeeze(R_arr(2,:,:))))
