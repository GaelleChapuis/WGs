% Depends on IBLLIB
cd /home/olivier/MatlabProjects/WGs/Physiology
close all


%[V, h, cs] = brainatlas.get_allen(100);
%  [V, h, cs] = brainatlas.get_allen(50);
[V, h, cs, labels] = brainatlas.get_whs12;

% X: ML (pitch), 2nd_dim (right to left)
% Y: AP (roll), 3d_dim 
% Z: DV (yaw), 1st_dim

P01_make_all_electrodes;
P02_prune_electrodes;
% P03_make_sections;

% 
% 
% set(h_.fig_volume, 'Position', [647         123        1080         773])
% view(h_.ax, -116.9542, 41.6299);
% axis(h_.ax, [-0.006 0.006 -0.012 0.008 -0.006 0.004])

