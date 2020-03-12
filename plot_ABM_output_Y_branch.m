function [] = plot_ABM_output_Y_branch(file_name)

close all
clc

% file_name = uigetfile;
% load(file_name)

save_directory = [pwd '\Simulation Data'];
cd(save_directory)
load([file_name '.mat'])

vidObj = VideoWriter(file_name,'MPEG-4');
vidObj.FrameRate = 4;
vidObj.Quality = 100;
open(vidObj);

num_nodes = length(nodes);
[num_vess num_timesteps] = size(vess_diameter);

nodes = nodes*1e6;
vess_conn = vess_conn + ones(num_vess, 2);

nodal_pressures = nodal_pressures;

vess_diameter = vess_diameter*1e6;
vess_flow = vess_flow*3.6e12;

Lseg = 10;
cellspeed = 3;
t_time = Lseg/cellspeed;

time = linspace(1,num_timesteps,num_timesteps) - 1;
time = time*t_time/24;
figure(1)

for t = 1:num_timesteps
    close (1)
    figure(1), hold on, grid on, axis square, axis([-105 105 -10 200])
    
    for v = 1:num_vess
        x1 = nodes(vess_conn(v, 1), 1);
        y1 = nodes(vess_conn(v, 1), 2);
        x2 = nodes(vess_conn(v, 2), 1);
        y2 =  nodes(vess_conn(v, 2), 2);
        
        if (vess_diameter(v,t) ~= 0)
            plot([x1 x2], [y1 y2], 'r-','MarkerSize',14,'LineWidth', vess_diameter(v, t)/3)
        end
    end
    
    for v = 1:num_vess
        x1 = nodes(vess_conn(v, 1), 1);
        y1 = nodes(vess_conn(v, 1), 2);
        x2 = nodes(vess_conn(v, 2), 1);
        y2 =  nodes(vess_conn(v, 2), 2);
                
        plot([x1 x2], [y1 y2], 'k.:', 'MarkerSize', 5.0)
    end
    
    title([' day ' num2str(time(t))])
    
    text(-10,-5,['Pout = ' num2str(round(nodal_pressures(PBCs(1,1)+1, t),3))],'Color','k','FontSize',10)
    text(75,180,['Pright = ' num2str(round(nodal_pressures(PBCs(2,1)+1, t),3))],'Color','k','FontSize',10)
    text(-95,180,['Pleft = ' num2str(round(nodal_pressures(PBCs(3,1)+1, t),3))],'Color','k','FontSize',10)

    text(5,100,[num2str(round(nodal_pressures(11, t),3))],'Color','k','FontSize',10)
    
    text(5,50,[num2str(round(mean(vess_flow(1:10, t)),3))],'Color','b','FontSize',10)
    text(35,125,[num2str(round(mean(vess_flow(11:20, t)),3))],'Color','b','FontSize',10)
    text(-55,125,[num2str(round(mean(vess_flow(21:30, t)),3))],'Color','b','FontSize',10)
    
    text(-25,50,[num2str(round(mean(vess_WSS(1:10, t)),3))],'Color','r','FontSize',10)
    text(15,145,[num2str(round(mean(vess_WSS(11:20, t)),3))],'Color','r','FontSize',10)
    text(-35,145,[num2str(round(mean(vess_WSS(21:30, t)),3))],'Color','r','FontSize',10)
    
    text(-48,-25,['Pressure, P (Pa)'],'Color','k','FontSize',10)
    text(30,-25,['Flow, Q (\muL/hr)'],'Color','b','FontSize',10)
    text(110,-25,['WSS (Pa)'],'Color','r','FontSize',10)
    
    box on
    set(gcf, 'Color', 'w')
    axis off
    
    % Write to the video file
    currFrame = getframe(gcf);
    writeVideo(vidObj,currFrame)
    
    pause(0.01)
end

close(vidObj);
