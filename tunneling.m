% experiments with wave tunneling
close all; clear all;
path(path, '/Users/charlesxu/Documents/WHOI_2022/reflection_coeff/subroutines')

% 3 layered medium
%top layer and bottom half space are same material (rho, c)
%middle layer is rho_1 > rho, c_1 > c

%choose some arbitrary c, c_1, rho, rho_1 satisfying rules
c = 1000; c_1 = 2000;
rho = 1; rho_1 = 1.5;

theta_c = asin(c/c_1) * 180/pi;

c_arr = [c, c_1, c]; rho_arr = [rho, rho_1, rho];
h_arr = [2:2:20];%1:5:11;
freq_arr = 500:500:2000; % vary this paramter
angle_arr = 0:pi/100:pi/2; % change depending on run time

R_arr = zeros(length(angle_arr), length(h_arr), length(freq_arr));
%R_arr_2 = zeros(length(angle_arr), length(h_arr), length(freq_arr));
T_arr = zeros(length(angle_arr), length(h_arr), length(freq_arr));

for angle_index = 1:length(angle_arr)
    for h_index = 1:length(h_arr)
        for freq_index = 1:length(freq_arr)
            R_array = rayleigh_strat(rho_arr,c_arr,h_arr(h_index),angle_arr(angle_index),freq_arr(freq_index));
            T_total = multi_layer_transmission(rho_arr, c_arr, h_arr(h_index), angle_arr(angle_index), freq_arr(freq_index));
            R_arr(angle_index, h_index, freq_index) = R_array(1);
            T_arr(angle_index, h_index, freq_index) = T_total;
            %R_arr_2(angle_index, h_index, freq_index) = R_array(2);
        end
    end
end

%% 
%section for plotting

close all;

for freq_index = 1:length(freq_arr)
    figure(freq_index); clf;
    for h_index = 1:length(h_arr)
        plot(angle_arr * 180/pi, abs(R_arr(:,h_index, freq_index)));
        hold on;
    end
    xline(theta_c,'--r',{'Critical','angle'},'LineWidth',2,'HandleVisibility','off');
    legendCell = cellstr(num2str(h_arr', 'h=%-d'));
    %legendCell(length(legendCell) + 1) = {"Critical angle"};
    legend(legendCell);
    hold off;
    title("|R| for " +freq_arr(freq_index) +  " Hz");
    xlabel("angle");
    saveas(gcf, "R_Tunneling_"+freq_arr(freq_index) + "_hz.png");

end

for freq_index = 1:length(freq_arr)
    figure(length(freq_arr) + freq_index); clf; 
    for h_index = 1:length(h_arr)
        
        plot(angle_arr * 180/pi, abs(T_arr(:,h_index, freq_index)) );
        hold on;
    end
    xline(theta_c,'--r',{'Critical','angle'},'LineWidth',2,'HandleVisibility','off');
    legendCell = cellstr(num2str(h_arr', 'h=%-d'));
    %legendCell(length(legendCell) + 1) = {"Critical angle"};
    %legend('1','5','10');
    legend(legendCell);
    hold off;
    title("|T| for " +freq_arr(freq_index) +  " Hz");
    xlabel("angle");
    saveas(gcf, "T_Tunneling_"+freq_arr(freq_index) + "_hz.png");

end

%plot(angle_arr * 180/pi, abs_R);

% saveas(gcf, 'Figure 11.png')