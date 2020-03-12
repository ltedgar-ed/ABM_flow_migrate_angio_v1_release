close all
clc

file_name = uigetfile
load([file_name])

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

% Set font of labels
pf = 10;
qf = 8;
lf = 10;

for i = 1:length(PBCs)
    PBCs(i,1) = PBCs(i,1) + 1;
end

figure(1)

vess_type = strings(num_vess,1);

for v = 1:num_vess
    if (vess_num_cells(v,1) == 10)
        vess_type(v,1) = "A";
    end
    
    if (vess_num_cells(v,1) == 5)
        vess_type(v,1) = "C";
    end
    
    if (vess_num_cells(v,1) == 20)
        vess_type(v,1) = "V";
    end
end

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

tot_vascular_length = zeros(1,num_timesteps);
tot_vascular_volume = zeros(1,num_timesteps);

for t = 1:num_timesteps
    
    for v = 1:num_vess
        x1 = nodes(vess_conn(v, 1), 1);
        y1 = nodes(vess_conn(v, 1), 2);
        x2 = nodes(vess_conn(v, 2), 1);
        y2 =  nodes(vess_conn(v, 2), 2);
        
        vess_length = sqrt((x2 - x1)^2 + (y2 - y1)^2);
        
        if (vess_num_cells(v,t) ~= 0)
            tot_vascular_length(t) = tot_vascular_length(t) + vess_length;
            tot_vascular_volume(t) = tot_vascular_volume(t) + pi*((vess_diameter(v,t)/2)^2)*vess_length;
        end
    end
end

figure(1), plot(tot_vascular_length, 'Color', [0.64 0.08 0.18], 'LineWidth', 3), grid on, xlabel(' Time Steps '), ylabel(' Total Vascular Length (\mum) '), axis([1 num_timesteps min(tot_vascular_length) max(tot_vascular_length)])
figure(2), plot(tot_vascular_volume, 'LineWidth', 3), grid on, xlabel(' Time Steps '), ylabel(' Total Vascular Volume (\mum^3) '), axis([1 num_timesteps min(tot_vascular_volume) max(tot_vascular_volume)])

start_vascular_length = tot_vascular_length(1)
end_vascular_length = tot_vascular_length(num_timesteps)
start_vascular_volume = tot_vascular_volume(1)
end_vascular_volume = tot_vascular_volume(num_timesteps)
start_end_ration = tot_vascular_volume(1)/tot_vascular_volume(num_timesteps)
end_start_ratio = tot_vascular_volume(num_timesteps)/tot_vascular_volume(1)

node_degree = zeros(num_nodes,1);

for v = 1:num_vess
    node1 = vess_conn(v,1);
    node2 = vess_conn(v,2);
    
    node_degree(node1) = node_degree(node1) + 1;
    node_degree(node2) = node_degree(node2) + 1;
end

num_deg3_nodes = 0;

for i = 1:num_nodes
    if (node_degree(i) == 3)
        num_deg3_nodes = num_deg3_nodes + 1;
    end
end

bifurcations = [find(node_degree ==3) zeros(num_deg3_nodes,3)];

for i = 1:num_deg3_nodes
    node = bifurcations(i,1);
    
    for v = 1:num_vess
        if (vess_conn(v,1) == node) || (vess_conn(v,2) == node)
            if (bifurcations(i,2) == 0)
                bifurcations(i,2) = v;
            else if (bifurcations(i,3) == 0)
                    bifurcations(i,3) = v;
                else if (bifurcations(i,4) == 0)
                        bifurcations(i,4) = v;
                    end
                end
            end
        end
    end
end

num_cap_bifurcations_start = 0;

num_cap_bifurcations_remaining = 0;
figure(2), hold on

for i = 1:num_deg3_nodes
    node = bifurcations(i,1);
    v1 = bifurcations(i,2);
    v2 = bifurcations(i,3);
    v3 = bifurcations(i,4);
    
    if (vess_type(v1) == "C") && (vess_type(v2) == "C") && (vess_type(v3) == "C")
        if (vess_num_cells(v1,1) ~= 0) && (vess_num_cells(v2,1) ~= 0) && (vess_num_cells(v3,1) ~= 0)
            %figure(2), plot(nodes(node,1), nodes(node,2), 'bo', 'MarkerSize',10)
            num_cap_bifurcations_start = num_cap_bifurcations_start + 1;
        end
    end
    
    if (vess_type(v1) == "C") && (vess_type(v2) == "C") && (vess_type(v3) == "C")
        if (vess_num_cells(v1,t) ~= 0) && (vess_num_cells(v2,t) ~= 0) && (vess_num_cells(v3,t) ~= 0)
            %figure(2), plot(nodes(node,1), nodes(node,2), 'r*', 'MarkerSize',10)
            num_cap_bifurcations_remaining = num_cap_bifurcations_remaining + 1;
        end
    end
end

num_cap_bifurcations_start
num_cap_bifurcations_remaining

