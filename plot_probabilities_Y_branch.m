clc
clear all
close all

load('ABM_output_Y_branch_Nt_36_Ncell_8_bifrule_7_alpha_1.000_run10.mat')

[num_vess, Nt] = size(vess_WSS);

v1 = 11;
v2 = 21;

for t = 1:Nt
    Pt1(t) = vess_WSS(v1,t)/(vess_WSS(v1,t) + vess_WSS(v2,t));
    Pn1(t) = vess_num_cells(v1,t)/(vess_num_cells(v1,t) + vess_num_cells(v2,t));
    
    Pt2(t) = vess_WSS(v2,t)/(vess_WSS(v1,t) + vess_WSS(v2,t));
    Pn2(t) = vess_num_cells(v2,t)/(vess_num_cells(v1,t) + vess_num_cells(v2,t));

    n1 = vess_conn(v1,1)+1;
    n2 = vess_conn(v1,2)+1;
    dp1(t) = nodal_pressures(n2,t) - nodal_pressures(n1,t);
    
    n1 = vess_conn(v2,1)+1;
    n2 = vess_conn(v2,2)+1;
    dp2(t) = nodal_pressures(n2,t) - nodal_pressures(n1,t);
end

figure(1), hold on
plot(Pt1,'b-', 'LineWidth', 3)
plot(Pn1,'r-', 'LineWidth', 3)
%plot((Pt1 + Pn1)/2, 'g-', 'LineWidth', 3)

plot(Pt2,'b-.', 'LineWidth', 3)
plot(Pn2,'r-.', 'LineWidth', 3)
%plot((Pt2 + Pn2)/2, 'g-.', 'LineWidth', 3)

legend('Pt1', 'Pn1', 'P1', 'Pt2', 'Pn2', 'P2', 'location', 'northeastoutside')

axis([1 Nt-1 0 1])
xlabel(' time steps ')
ylabel(' probability ')
box on
set(gca, 'FontName', 'Cambria Math')
set(gca, 'FontSize', 26)
set(gca, 'LineWidth', 2)
set(gcf, 'Color', 'w')

% figure(2), hold on
% set(gca, 'FontName', 'Cambria Math')
% set(gca, 'FontSize', 26)
% set(gca, 'LineWidth', 2)
% set(gcf, 'Color', 'w')
% xlabel(' time steps ')
% box on
% yyaxis left
% plot(-dp1, 'k', 'LineWidth', 3)
% ylabel(' |dp| (Pa) ')
% yyaxis right
% plot(vess_num_cells(v1,:), 'r', 'LineWidth', 3)
% ylabel(' cell number ')
% %title(' branch 1 (high flow) ')
% ax = gca;
% ax.XAxis.Limits = [1 Nt-1];
% ax.YAxis(1).Limits = [0 100];
% ax.YAxis(2).Limits = [0 15];
% ax.YAxis(1).Color = 'k';
% ax.YAxis(2).Color = 'r';
% 
% 
% figure(3), hold on
% set(gca, 'FontName', 'Cambria Math')
% set(gca, 'FontSize', 26)
% set(gca, 'LineWidth', 2)
% set(gcf, 'Color', 'w')
% xlabel(' time steps ')
% box on
% yyaxis left
% plot(-dp2, 'k', 'LineWidth', 3)
% ylabel(' |dp| (Pa) ')
% yyaxis right
% plot(vess_num_cells(v2,:), 'r', 'LineWidth', 3)
% ylabel(' cell number ')
% %title(' branch 2 (low flow) ')
% ax = gca;
% ax.XAxis.Limits = [1 Nt-1];
% ax.YAxis(1).Limits = [0 100];
% ax.YAxis(2).Limits = [0 15];
% ax.YAxis(1).Color = 'k';
% ax.YAxis(2).Color = 'r';
