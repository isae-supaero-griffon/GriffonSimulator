% Design output code post-processing script. The following script pretends
% to nucleate the output files in order to compare rocket performance under
% different configurations.
% Author: Jose Felix Zapata Usandivaras
% Date: 6/2/2019
% ISAE-SUPAERO Space Section / Griffon Project

clc
clear
close all
addpath(genpath(pwd))

%% Get the file names and import their data

directory = [pwd '\ThreePortGeometry\'];                        % Set the directory
csv_list = dir(directory);                                              % Get files in directory
csv_list([1 2]) = [];                                                   % Remove initial non-representative values
is_csv = cellfun(@(str) strcmp(str(end-3:end),'.csv'), ...              % Get only those files who have a csv extension
                {csv_list.name}', 'UniformOutput', false);
csv_list = csv_list(cell2mat(is_csv));                                  % Filter those files out
IX = regexp({csv_list.name}', '(\d+)', 'tokens');
IX = arrayfun(@(c) str2double(c{1}{1}{1}), IX);
[~, IX] = sort(IX);                                                     % Sort the structure array

csv_list = csv_list(IX);

% Parse the csv_files
results = struct('title', {}, ...
                 'geometry', {}, ...
                 'rocket_mass', {}, ...
                 'average_values', {}, ...
                 'drag', {}, ...
                 'trajectory', {}, ...
                 'flight_table', {});
             
r = cellfun(@(file_name) parse_csv_file([directory file_name]), {csv_list.name}', 'UniformOutput', false);
for i = 1:1:numel(r), results(end+1) = r{i} ; end %#ok

%% Load results already stored if they exist

% load('run_results_compilation_three_port_geometry.mat');

%% Make plots 

close all
% Plot3s: variables to be plotted: [thrust, burn_time, max_altitude, max_velocity] 

% ------------------------ Start with Figure 1 --------------------------
figure('name', 'Global Magnitudes', 'units', 'normalized', ...
       'outerposition', [0 0 1 1])
   
xdata = arrayfun(@(str) str.geometry.r_ext, results);
ydata = arrayfun(@(str) str.average_values.ox_flow, results);

%%%%%%%% Thrust Plot
subplot(2, 2, 1);

zdata = arrayfun(@(str) str.average_values.thrust, results);

plot3(xdata, ydata, zdata, 'color', 'b', 'DisplayName', 'Mean Thrust [N]', ...
     'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', ...
     'b', 'MarkerEdgeColor', 'b');
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on');
lg = legend(gca);
set(lg, 'box', 'off', 'location', 'southoutside', 'orientation', 'horizontal');
xlabel(gca, 'Chamber Ext Radius [m]', 'fontsize', 11);
ylabel(gca, 'Ox Flow [kg/sec]', 'fontsize', 11);
zlabel(gca, 'Mean Thrust [N]', 'fontsize', 11);
title(gca, 'Mean Thrust', 'fontsize', 12);

%%%%%%%%% ISP

subplot(2, 2, 2);

zdata = arrayfun(@(str) str.average_values.isp, results);
plot3(xdata, ydata, zdata, 'color', 'r', 'DisplayName', 'Mean Thrust [N]', ...
     'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', ...
     'r', 'MarkerEdgeColor', 'r')
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on');
lg = legend(gca);
set(lg, 'box', 'off', 'location', 'southoutside', 'orientation', 'horizontal');
xlabel(gca, 'Chamber Ext Radius [m]', 'fontsize', 11);
ylabel(gca, 'Ox Flow [kg/sec]', 'fontsize', 11);
zlabel(gca, 'Isp [secs]', 'fontsize', 11);
title(gca, 'Isp', 'fontsize', 12);

%%%%%%%%%%% Max Altitude

subplot(2, 2, 3);

zdata = arrayfun(@(str) str.average_values.burn_time, results);
plot3(xdata, ydata, zdata, 'color', 'g', 'DisplayName', 'Mean Thrust [N]', ...
     'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', ...
     'g', 'MarkerEdgeColor', 'g')
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on');
lg = legend(gca);
set(lg, 'box', 'off', 'location', 'southoutside', 'orientation', 'horizontal');
xlabel(gca, 'Chamber Ext Radius [m]', 'fontsize', 11);
ylabel(gca, 'Ox Flow [kg/sec]', 'fontsize', 11);
zlabel(gca, 'Burn Time [secs]', 'fontsize', 11);
title(gca, 'Burn Time', 'fontsize', 12);

%%%%%%%%%% Max Velocity

subplot(2, 2, 4);

zdata = arrayfun(@(str) str.average_values.of, results);

plot3(xdata, ydata, zdata, 'color', 'k', 'DisplayName', 'Mean Thrust [N]', ...
     'LineStyle', 'none', 'Marker', 'o', 'MarkerSize', 6, 'MarkerFaceColor', ...
     'k', 'MarkerEdgeColor', 'k')
set(gca, 'XGrid', 'on', 'YGrid', 'on', 'ZGrid', 'on');
lg = legend(gca);
set(lg, 'box', 'off', 'location', 'southoutside', 'orientation', 'horizontal');
xlabel(gca, 'Chamber Ext Radius [m]', 'fontsize', 11);
ylabel(gca, 'Ox Flow [kg/sec]', 'fontsize', 11);
zlabel(gca, 'OF', 'fontsize', 11);
title(gca, 'OF', 'fontsize', 12);


% --------------- Plot the flight envelopes

%%%%%%%% ALTITUDE PLOTS

figure('name', 'Flight Envelopes', 'units', 'normalized', 'outerposition', [0 0 1 1]);

for str = results
    plot_name = sprintf('L : %5.2f, r_ext : %5.2f', str.geometry.length, str.geometry.r_ext);
    plot(str.flight_table.time, str.flight_table.altitude, 'LineStyle', '--', ...
         'color', 'b', 'DisplayName', plot_name);
    hold(gca, 'on')
end

set(gca, 'XGrid', 'on', 'YGrid', 'on')
lg = legend(gca);
% set(lg, 'box', 'off', 'location', 'southoutside', 'orientation', 'horizontal')
xlabel(gca, 'Time [secs]', 'fontsize', 12);
ylabel(gca, 'Altitude [m]', 'fontsize', 12);
title(gca, 'Flight Path Rocket', 'fontsize', 12);

%%%%%%%% SPEED PLOTS

figure('name', 'Speed Envelopes', 'units', 'normalized', 'outerposition', [0 0 1 1]);

for str = results
    plot_name = sprintf('L : %5.2f, r_ext : %5.2f', str.geometry.length, str.geometry.r_ext);
    plot(str.flight_table.time, str.flight_table.velocity, 'LineStyle', '--', ...
         'color', 'r', 'DisplayName', plot_name);
    hold(gca, 'on')
end

set(gca, 'XGrid', 'on', 'YGrid', 'on')
lg = legend(gca);
% set(lg, 'box', 'off', 'location', 'southoutside', 'orientation', 'horizontal')
xlabel(gca, 'Time [secs]', 'fontsize', 12);
ylabel(gca, 'Speed [m/seg]', 'fontsize', 12);
title(gca, 'Speed Rocket', 'fontsize', 12);

%%%%%%%% ALTITUDE vs. SPEED PLOTS

figure('name', 'Flight Curve', 'units', 'normalized', 'outerposition', [0 0 1 1]);

for str = results
    plot_name = sprintf('L : %5.2f, r_ext : %5.2f', str.geometry.length, str.geometry.r_ext);
    plot(str.flight_table.altitude, str.flight_table.velocity, 'LineStyle', '--', ...
         'color', 'k', 'DisplayName', plot_name);
    hold(gca, 'on')
end

set(gca, 'XGrid', 'on', 'YGrid', 'on', 'XLim', [0, 20], 'YLim', [-50, 50]);
lg = legend(gca);
% set(lg, 'box', 'off', 'location', 'westoutside', 'orientation', 'vertical')
xlabel(gca, 'Altitude [m]', 'fontsize', 12);
ylabel(gca, 'Speed [m/seg]', 'fontsize', 12);
title(gca, 'Speed vs. Altitude', 'fontsize', 12);

%%%%%%
