function [] = plot_ABM_output_ideal_cap2(file_name)

close all
clc

save_directory = [pwd '\Simulation Data'];
cd(save_directory)
load([file_name '.mat'])

% file_name = uigetfile;
% load([file_name])

vidObj = VideoWriter(file_name,'MPEG-4');
vidObj.FrameRate = 4;
vidObj.Quality = 100;
open(vidObj);

num_nodes = length(nodes);
[num_vess num_timesteps] = size(vess_diameter);

nodes = nodes*1e6;
vess_conn = vess_conn + ones(num_vess, 2);

%nodal_pressures = nodal_pressures/98;

vess_diameter = vess_diameter*1e6;
vess_flow = vess_flow*3.6e12;

% Provide information on each honeycomb
a = 40;
b = a/sqrt(2);
Lhc = (1 + sqrt(2))*a;

num_hc_long = round(max(nodes(:,1))/((1 + sqrt(2))*a));
num_hc_high = round(max(nodes(:,2))/(2*b));

Lseg = 5;
cellspeed = 3;
t_time = Lseg/cellspeed;

time = linspace(1,num_timesteps,num_timesteps) - 1;
time = time*t_time/24;

% Set font of labels
pf = 10;
qf = 8;
lf = 10;

for i = 1:length(PBCs)
    PBCs(i,1) = PBCs(i,1) + 1;
end

figure(1)

% vess_type = strings(num_vess,1);
% 
% for v = 1:num_vess
%     if (vess_num_cells(v,1) == 10)
%         vess_type(v,1) = "A";
%     end
%     
%     if (vess_num_cells(v,1) == 5)
%         vess_type(v,1) = "C";
%     end
%     
%     if (vess_num_cells(v,1) == 20)
%         vess_type(v,1) = "V";
%     end
% end

vess_DP_flow = zeros(num_vess, 2);

for v = 1:num_vess
    node1 = vess_conn(v,1);
    node2 = vess_conn(v,2);
    
    vess_DP_flow(v,1) = nodal_pressures(node2,1) - nodal_pressures(node1,1);
    vess_DP_flow(v,2) = vess_flow(v,1);
end

vess_unit_vect = zeros(num_vess, 2);

for v = 1:num_vess
    node1 = vess_conn(v,1);
    node2 = vess_conn(v,2);
    
    vess_vect = nodes(node2,:) - nodes(node1,:);
    vess_vect = vess_vect/norm(vess_vect);
    
    vess_unit_vect(v,:) = vess_vect;
end

flow_max = max(max(vess_flow));
flow_min = min(min(vess_flow));
flow_zero = 1e-4;

red_max = [1 0.5 0.5];
red_min = [0.9 0 0];
blue_max = [0.5 0.5 1];
blue_min = [0 0 0.9];
zero_gray = [0.8 0.8 0.8];

x_range = max(nodes(:,1)) - min(nodes(:,1));
y_range = max(nodes(:,2)) - min(nodes(:,2));

x_min = min(nodes(:,1)) -(x_range/4);
x_max = max(nodes(:,1)) +(x_range/4);
y_min = min(nodes(:,2)) -(y_range/4);
y_max = max(nodes(:,2)) +(y_range/4);

for t = 1:num_timesteps
close (1)
    figure(1), hold on, grid on, axis image, axis([x_min x_max y_min y_max])
    
    if (x_range > y_range)
        set(gcf, 'Position', [50 50 (x_range/y_range)*600 600])
    else
        set(gcf, 'Position', [50 50 600 (y_range/x_range)*600])
    end
    set(gcf, 'Color', 'w')
    %axis off
    
    vess_colors = zeros(num_vess,3);
    
    for v = 1:num_vess
        x1 = nodes(vess_conn(v, 1), 1);
        y1 = nodes(vess_conn(v, 1), 2);
        x2 = nodes(vess_conn(v, 2), 1);
        y2 =  nodes(vess_conn(v, 2), 2);
        
        if (vess_flow(v,t) > flow_zero)
            xi = vess_flow(v,t)/flow_max;
            
            vess_colors(v,:) = (1 - xi)*red_min + xi*red_max;
        end
        
        if (vess_flow(v,t) < -flow_zero)
                xi = vess_flow(v,t)/flow_min;
                
                vess_colors(v,:) = (1 - xi)*blue_min + xi*blue_max;
        end
        
        if (abs(vess_flow(v,t)) < flow_zero)
            vess_colors(v,:) = zero_gray;
        end
        
                
        if (vess_diameter(v,t) ~= 0)
            plot([x1 x2], [y1 y2], 'Color', vess_colors(v,:), 'MarkerSize',14,'LineWidth', vess_diameter(v, t)/3)
        end
        
        arrow3([x1 y1], [x2 y2], 'd', 0.6, 0.6)
    end
    
    flow_range = [linspace(flow_min, -flow_zero, 5) linspace(flow_zero, flow_max, 5)]';
    
    leg_plot_length = (max(nodes(:,2)) - min(nodes(:,1)))/length(flow_range);
    leg_colormap = [];
    
    for i = 1:length(flow_range)
        if (flow_range(i) < -flow_zero)
            xi = flow_range(i)/flow_min;
            
            leg_plot_color = (1 - xi)*blue_min + xi*blue_max;
            
            %plot([x_max-x_range/10 x_max-x_range/10], [min(nodes(:,2))+(i-1)*leg_plot_length min(nodes(:,2))+i*leg_plot_length], 'Color', leg_plot_color, 'LineWidth', 10)
        end
        
        if (abs(flow_range(i)) <= flow_zero)
            leg_plot_color = zero_gray;
            
            %plot([x_max-x_range/10 x_max-x_range/10], [min(nodes(:,2))+(i-1)*leg_plot_length min(nodes(:,2))+i*leg_plot_length], 'Color', leg_plot_color, 'LineWidth', 10)
        end
        
        if (flow_range(i) > flow_zero)
            xi = flow_range(i)/flow_max;
            
            leg_plot_color = (1 - xi)*red_min + xi*red_max;
            
            %plot([x_max-x_range/10 x_max-x_range/10], [min(nodes(:,2))+(i-1)*leg_plot_length min(nodes(:,2))+i*leg_plot_length], 'Color', leg_plot_color, 'LineWidth', 10)
        end
        
        leg_colormap = [leg_colormap; leg_plot_color];
    end
    
    leg_tick_labels = cell(length(flow_range)+1,1);
    
    for i = 1:length(flow_range)/2
        leg_tick_labels{i} = num2str(flow_range(i));
    end
    
    leg_tick_labels{length(flow_range)/2 + 1} = num2str(0);
    
    for i = length(flow_range)/2 + 2:length(flow_range)+1
        leg_tick_labels{i} = num2str(flow_range(i-1));
    end
    
    colormap(leg_colormap);
    c = colorbar;  
    set(c, 'Ticks', linspace(0,1,length(flow_range)+1));
    set(c, 'TickLabels', leg_tick_labels);
    c.Label.String = ' Flow (\muL/hr)';
    c.Label.FontSize = 12;
    
    title([' day ' num2str(time(t))])
    
%     text(min(nodes(:,1))-30, min(nodes(:,2))-8, ['  P_i_n_1 = ' num2str(round(nodal_pressures(PBCs(1,1), t)*0.0075)) ' mmHg'],'Color','k','FontSize', pf)
%     text(min(nodes(:,1))-30, max(nodes(:,2))+8, ['P_o_u_t_1 = ' num2str(round(nodal_pressures(PBCs(2,1), t)*0.0075)) ' mmHg'],'Color','k','FontSize', pf)
%     text(max(nodes(:,1))-30, min(nodes(:,2))-8, ['  P_i_n_2 = ' num2str(round(nodal_pressures(PBCs(3,1), t)*0.0075)) ' mmHg'],'Color','k','FontSize', pf)
%     text(max(nodes(:,1))-30, max(nodes(:,2))+8, ['P_o_u_t_2 = ' num2str(round(nodal_pressures(PBCs(4,1), t)*0.0075)) ' mmHg'],'Color','k','FontSize', pf)
%     
%     for i = 1:length(bifur_nodes)
%         bif_node = bifur_nodes(i)+1;
%         
%         if (bifur(i,t) == 1)
%             plot(nodes(bif_node,1), nodes(bif_node,2), 'ro', 'MarkerSize', 20)
%         else if (bifur(i,t) == -1)
%                 plot(nodes(bif_node,1), nodes(bif_node,2), 'kd', 'MarkerSize', 20)
%             end
%         end   
%     end
    
    axis off
    
    % Write to the video file
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame)
    
    %pause(0.01)
end

close(vidObj);
