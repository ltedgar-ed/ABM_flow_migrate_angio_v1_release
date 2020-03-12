function [] = generate_ideal_cap_bed2(num_hc_long, num_hc_high, Ncell_cap)

close all			  

num_hc_long = double(num_hc_long);
num_hc_high = double(num_hc_high);

a = 25;

b = a/sqrt(2);
Lhc = (1 + sqrt(2))*a;



hc_v1 = [sind(45) cosd(45)];
hc_v2 = [sind(45) -cosd(45)];
hc_v3 = [1 0];
hc_v4 = [1 0];
hc_v5 = [cosd(45) -sind(45)];
hc_v6 = [cosd(45) sind(45)];

Lseg = 10;

Nseg_a_hc = a/Lseg;


%figure(1), hold on

vess_seg = [];

for j = 1:num_hc_long

    origin = [];
    
    if (mod(j,2) ~= 0)
        for i = 1:2:2*num_hc_high
            if (j == 1)
                next_origin = [a/2 i*b];
            else
                next_origin = [floor(j/2)*(Lhc+a)+(a/2) i*b];
            end
            
            origin = [origin; next_origin];
        end
        
        [len wid] = size(origin);
        
        for i = 1:len
            if (j == 1)
                vess_seg = [vess_seg; 0 a/2 origin(i,2) origin(i,2) Ncell_cap];
                %plot([0 a/2], [origin(i,2) origin(i,2)], 'r.-', 'MarkerSize', 8, 'LineWidth', d_cap/2)
            end
            
            vess_seg = [vess_seg; origin(i,1) origin(i,1)+a*hc_v1(1) origin(i,2) origin(i,2)+a*hc_v1(2) Ncell_cap];
            %plot([origin(i,1) origin(i,1)+a*hc_v1(1)], [origin(i,2) origin(i,2)+a*hc_v1(2)], 'r.-', 'MarkerSize', 8, 'LineWidth', d_cap/2) 
            
            vess_seg = [vess_seg; origin(i,1) origin(i,1)+a*hc_v2(1) origin(i,2) origin(i,2)+a*hc_v2(2) Ncell_cap];
            %plot([origin(i,1) origin(i,1)+a*hc_v2(1)], [origin(i,2) origin(i,2)+a*hc_v2(2)], 'r.-', 'MarkerSize', 8, 'LineWidth', d_cap/2)
            
            vess_seg = [vess_seg; origin(i,1)+a*hc_v1(1) origin(i,1)+a*hc_v1(1)+a*hc_v3(1) origin(i,2)+a*hc_v1(2) origin(i,2)+a*hc_v1(2)+a*hc_v3(2) Ncell_cap];
            %plot([origin(i,1)+a*hc_v1(1) origin(i,1)+a*hc_v1(1)+a*hc_v3(1)], [origin(i,2)+a*hc_v1(2) origin(i,2)+a*hc_v1(2)+a*hc_v3(2)], 'r.-', 'MarkerSize', 8, 'LineWidth', d_cap/2)
            
            vess_seg = [vess_seg; origin(i,1)+a*hc_v2(1) origin(i,1)+a*hc_v2(1)+a*hc_v4(1) origin(i,2)+a*hc_v2(2) origin(i,2)+a*hc_v2(2)+a*hc_v4(2) Ncell_cap];
            %plot([origin(i,1)+a*hc_v2(1) origin(i,1)+a*hc_v2(1)+a*hc_v4(1)], [origin(i,2)+a*hc_v2(2) origin(i,2)+a*hc_v2(2)+a*hc_v4(2)], 'r.-', 'MarkerSize', 8, 'LineWidth', d_cap/2)
            
            vess_seg = [vess_seg; origin(i,1)+a*hc_v1(1)+a*hc_v3(1) origin(i,1)+a*hc_v1(1)+a*hc_v3(1)+a*hc_v5(1) origin(i,2)+a*hc_v1(2)+a*hc_v3(2) origin(i,2)+a*hc_v1(2)+a*hc_v3(2)+a*hc_v5(2) Ncell_cap];
            %plot([origin(i,1)+a*hc_v1(1)+a*hc_v3(1) origin(i,1)+a*hc_v1(1)+a*hc_v3(1)+a*hc_v5(1)], [origin(i,2)+a*hc_v1(2)+a*hc_v3(2) origin(i,2)+a*hc_v1(2)+a*hc_v3(2)+a*hc_v5(2)], 'r.-', 'MarkerSize', 8, 'LineWidth', d_cap/2)
            
            vess_seg = [vess_seg; origin(i,1)+a*hc_v2(1)+a*hc_v4(1) origin(i,1)+a*hc_v2(1)+a*hc_v4(1)+a*hc_v6(1) origin(i,2)+a*hc_v2(2)+a*hc_v4(2) origin(i,2)+a*hc_v2(2)+a*hc_v4(2)+a*hc_v6(2) Ncell_cap];
            %plot([origin(i,1)+a*hc_v2(1)+a*hc_v4(1) origin(i,1)+a*hc_v2(1)+a*hc_v4(1)+a*hc_v6(1)], [origin(i,2)+a*hc_v2(2)+a*hc_v4(2) origin(i,2)+a*hc_v2(2)+a*hc_v4(2)+a*hc_v6(2)], 'r.-', 'MarkerSize', 8, 'LineWidth', d_cap/2)
        
            if (j == num_hc_long)
                vess_seg = [vess_seg; origin(i,1)+a*hc_v2(1)+a*hc_v4(1)+a*hc_v6(1) origin(i,1)+a*hc_v2(1)+a*hc_v4(1)+a*hc_v6(1)+(a/2) origin(i,2)+a*hc_v2(2)+a*hc_v4(2)+a*hc_v6(2) origin(i,2)+a*hc_v2(2)+a*hc_v4(2)+a*hc_v6(2) Ncell_cap]; 
                %plot([origin(i,1)+a*hc_v2(1)+a*hc_v4(1)+a*hc_v6(1) origin(i,1)+a*hc_v2(1)+a*hc_v4(1)+a*hc_v6(1)+(a/2)], [origin(i,2)+a*hc_v2(2)+a*hc_v4(2)+a*hc_v6(2) origin(i,2)+a*hc_v2(2)+a*hc_v4(2)+a*hc_v6(2)], 'r.-', 'MarkerSize', 8, 'LineWidth', d_cap/2)
            end
        end
    else
        for i = 1:num_hc_high-1
            next_origin = [(Lhc-b)+((j/2)-1)*(a+Lhc)+(a/2) 2*i*b];
            origin = [origin; next_origin];
        end
        
        [len wid] = size(origin);
        
        for i = 1:len
            vess_seg = [vess_seg; origin(i,1) origin(i,1)+a*hc_v1(1) origin(i,2) origin(i,2)+a*hc_v1(2) Ncell_cap];
            %plot([origin(i,1) origin(i,1)+a*hc_v1(1)], [origin(i,2) origin(i,2)+a*hc_v1(2)], 'r.-', 'MarkerSize', 8, 'LineWidth', d_cap/2) 
            
            vess_seg = [vess_seg; origin(i,1) origin(i,1)+a*hc_v2(1) origin(i,2) origin(i,2)+a*hc_v2(2) Ncell_cap];
            %plot([origin(i,1) origin(i,1)+a*hc_v2(1)], [origin(i,2) origin(i,2)+a*hc_v2(2)], 'r.-', 'MarkerSize', 8, 'LineWidth', d_cap/2)
            
            vess_seg = [vess_seg; origin(i,1)+a*hc_v1(1) origin(i,1)+a*hc_v1(1)+a*hc_v3(1) origin(i,2)+a*hc_v1(2) origin(i,2)+a*hc_v1(2)+a*hc_v3(2) Ncell_cap];
            %plot([origin(i,1)+a*hc_v1(1) origin(i,1)+a*hc_v1(1)+a*hc_v3(1)], [origin(i,2)+a*hc_v1(2) origin(i,2)+a*hc_v1(2)+a*hc_v3(2)], 'r.-', 'MarkerSize', 8, 'LineWidth', d_cap/2)
            
            vess_seg = [vess_seg; origin(i,1)+a*hc_v2(1) origin(i,1)+a*hc_v2(1)+a*hc_v4(1) origin(i,2)+a*hc_v2(2) origin(i,2)+a*hc_v2(2)+a*hc_v4(2) Ncell_cap];
            %plot([origin(i,1)+a*hc_v2(1) origin(i,1)+a*hc_v2(1)+a*hc_v4(1)], [origin(i,2)+a*hc_v2(2) origin(i,2)+a*hc_v2(2)+a*hc_v4(2)], 'r.-', 'MarkerSize', 8, 'LineWidth', d_cap/2)
            
            vess_seg = [vess_seg; origin(i,1)+a*hc_v1(1)+a*hc_v3(1) origin(i,1)+a*hc_v1(1)+a*hc_v3(1)+a*hc_v5(1) origin(i,2)+a*hc_v1(2)+a*hc_v3(2) origin(i,2)+a*hc_v1(2)+a*hc_v3(2)+a*hc_v5(2) Ncell_cap];
            %plot([origin(i,1)+a*hc_v1(1)+a*hc_v3(1) origin(i,1)+a*hc_v1(1)+a*hc_v3(1)+a*hc_v5(1)], [origin(i,2)+a*hc_v1(2)+a*hc_v3(2) origin(i,2)+a*hc_v1(2)+a*hc_v3(2)+a*hc_v5(2)], 'r.-', 'MarkerSize', 8, 'LineWidth', d_cap/2)
            
            vess_seg = [vess_seg; origin(i,1)+a*hc_v2(1)+a*hc_v4(1) origin(i,1)+a*hc_v2(1)+a*hc_v4(1)+a*hc_v6(1) origin(i,2)+a*hc_v2(2)+a*hc_v4(2) origin(i,2)+a*hc_v2(2)+a*hc_v4(2)+a*hc_v6(2) Ncell_cap];
            %plot([origin(i,1)+a*hc_v2(1)+a*hc_v4(1) origin(i,1)+a*hc_v2(1)+a*hc_v4(1)+a*hc_v6(1)], [origin(i,2)+a*hc_v2(2)+a*hc_v4(2) origin(i,2)+a*hc_v2(2)+a*hc_v4(2)+a*hc_v6(2)], 'r.-', 'MarkerSize', 8, 'LineWidth', d_cap/2)
            
            if (j == num_hc_long)
                vess_seg = [vess_seg; origin(i,1)+a*hc_v2(1)+a*hc_v4(1)+a*hc_v6(1) origin(i,1)+a*hc_v2(1)+a*hc_v4(1)+a*hc_v6(1)+(a/2) origin(i,2)+a*hc_v2(2)+a*hc_v4(2)+a*hc_v6(2) origin(i,2)+a*hc_v2(2)+a*hc_v4(2)+a*hc_v6(2) Ncell_cap];
                %plot([origin(i,1)+a*hc_v2(1)+a*hc_v4(1)+a*hc_v6(1) origin(i,1)+a*hc_v2(1)+a*hc_v4(1)+a*hc_v6(1)+(a/2)], [origin(i,2)+a*hc_v2(2)+a*hc_v4(2)+a*hc_v6(2) origin(i,2)+a*hc_v2(2)+a*hc_v4(2)+a*hc_v6(2)], 'r.-', 'MarkerSize', 8, 'LineWidth', d_cap/2)
            end
        end
    end
    
%     if (j == 1)
%         y_max = 2*num_hc_high*b;
%         
%         origin(:,1) = zeros(len,1);
%         
%         artery = [0 0; origin; 0 y_max];
%         
%         artery_seg = [];
%         
%         for i = 1:length(artery)-1
%             artery_seg = [artery_seg; artery(i,1) artery(i+1,1) artery(i,2) artery(i+1,2) Ncell_art];
%             
%             %plot([artery(i,1) artery(i+1,1)], [artery(i,2) artery(i+1,2)], 'r.-', 'MarkerSize', 8, 'LineWidth', d_art/2)
%         end
%     end
%     
%     vess_seg = [artery_seg; vess_seg];
%     
%     if (j == num_hc_long)
%         x_max = floor(num_hc_long/2)*(Lhc+a)+(a/2);
%         
%         if (mod(num_hc_long,2) ~= 0)
%             x_max = x_max + Lhc + (a/2);
%         else
%             x_max = x_max + b + (a/2);
%         end
%         
%         origin(:,1) = x_max*ones(len,1);
%         
%         %vein = [x_max y_max; flipud(origin); x_max 0];
%         vein = [x_max 0; origin; x_max y_max];
%         
%         for i = 1:length(vein)-1
%             vess_seg = [vess_seg; vein(i,1) vein(i+1,1) vein(i,2) vein(i+1,2) Ncell_vein];
%             %plot([vein(i,1) vein(i+1,1)], [vein(i,2) vein(i+1,2)], 'r.-', 'MarkerSize', 8, 'LineWidth', d_vein/2)
%         end
%     end
end

% vess_diam = zeros(length(vess_seg), 1);
% 
% for i = 1:length(vess_seg)
%     vess_diam(i,1) = vess_seg(i,5)*cell_size/pi;
% end
% 
% for i = 1:length(vess_seg)
%     plot([vess_seg(i,1) vess_seg(i,2)], [vess_seg(i,3) vess_seg(i,4)], 'r.-', 'LineWidth', vess_diam(i,1)/2)
% end

% axis([-0.1*x_max 1.1*x_max -0.1*y_max 1.1*y_max])
% grid on

vess_seg = [vess_seg(:,1) vess_seg(:,3) vess_seg(:,2) vess_seg(:,4) vess_seg(:,5)];



save('ideal_cap_bed2.mat', 'vess_seg')
