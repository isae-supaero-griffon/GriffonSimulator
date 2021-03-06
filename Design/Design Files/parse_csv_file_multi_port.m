function results = parse_csv_file_multi_port(file_name)
% parse_csv_file parses the structure of the file so as to treat the
% results

% Load the csv file into memory
fID = fopen(file_name);
data = textscan(fID, '%s', 'delimiter', {',', '\n'});
fclose(fID);

% Get the simulation title
title = data{1}{1};

% Get the geometry parameters of interest (length and radius
geometry.type = data{1}{157};
geometry.length = str2double(data{1}{159});
geometry.ports_number = str2double(data{1}{161});
geometry.r_ext = str2double(data{1}{163});
geometry.r_int = str2double(data{1}{173});
geometry.minimum_thickness = str2double(data{1}{56});
geometry.max_regression_rate = str2double(data{1}{147});

% Get the rockets total wet mass
mass = data{1}{130};
mass = strrep(mass, ' kgs', '');
rocket_mass = str2double(mass);

% Get the average values results
average_values.thrust = str2double(data{1}{203});
average_values.isp = str2double(data{1}{205});
average_values.pressure = str2double(data{1}{207});
average_values.temperature = str2double(data{1}{209});
average_values.of = str2double(data{1}{213});
average_values.impulse = str2double(data{1}{215});
average_values.burn_time = str2double(data{1}{217});
average_values.ox_flow = str2double(data{1}{63});

% Print the values
fprintf('L: %5.2f, Ox Flow: %5.2f \n', geometry.length, average_values.ox_flow);

% Get the drag data
drag.cd = str2double(data{1}{226});
drag.S = str2double(data{1}{228});

% Get the trajectory results
trajectory.max_altitude = str2double(data{1}{233});
trajectory.max_velocity = str2double(data{1}{235});
trajectory.max_acceleration = str2double(data{1}{237});

% Flight path
headers = data{1}(239:248);
ncols = length(headers);
values = str2double(data{1}(249:end));
nrows = length(values)/ncols;
matrix = reshape(str2double(data{1}(249:end)), ncols, nrows)';
flight_table = array2table(matrix, 'VariableNames', headers');

% Allocate everything into the results structure
results = struct('title', title, ...
                 'geometry', geometry, ...
                 'rocket_mass', rocket_mass, ...
                 'average_values', average_values, ...
                 'drag', drag, ...
                 'trajectory', trajectory, ...
                 'flight_table', flight_table);

end % parse_csv_file