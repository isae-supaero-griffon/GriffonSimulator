function results = parse_csv_file(file_name)
% parse_csv_file parses the structure of the file so as to treat the
% results

% Load the csv file into memory
fID = fopen(file_name);
data = textscan(fID, '%s', 'delimiter', {',', '\n'});
fclose(fID);

% Get the simulation title
title = data{1}{1};

% Get the geometry parameters of interest (length and radius
geometry.type = data{1}{145};
geometry.length = str2double(data{1}{147});
geometry.ports_number = str2double(data{1}{149});
geometry.r_ext = str2double(data{1}{151});
geometry.r_int = str2double(data{1}{161});
geometry.minimum_thickness = str2double(data{1}{44});
geometry.max_regression_rate = str2double(data{1}{135});

% Get the rockets total wet mass
mass = data{1}{118};
mass = strrep(mass, ' kgs', '');
rocket_mass = str2double(mass);

% Get the average values results
average_values.thrust = str2double(data{1}{179});
average_values.isp = str2double(data{1}{181});
average_values.pressure = str2double(data{1}{183});
average_values.temperature = str2double(data{1}{185});
average_values.of = str2double(data{1}{189});
average_values.impulse = str2double(data{1}{191});
average_values.burn_time = str2double(data{1}{193});
average_values.ox_flow = str2double(data{1}{42});

% Get the drag data
drag.cd = str2double(data{1}{202});
drag.S = str2double(data{1}{204});

% Get the trajectory results
trajectory.max_altitude = str2double(data{1}{209});
trajectory.max_velocity = str2double(data{1}{211});
trajectory.max_acceleration = str2double(data{1}{213});

% Flight path
headers = data{1}(215:224);
ncols = length(headers);
values = str2double(data{1}(225:end));
nrows = length(values)/ncols;
matrix = reshape(str2double(data{1}(225:end)), ncols, nrows)';
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