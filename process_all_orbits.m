processed_data = [];

output_excel = 'C:\Users\DOSIMETER-PC\Downloads\out_full.xlsx';
output_statistics_excel = 'C:\Users\DOSIMETER-PC\Downloads\out_statistcs_full_1.xlsx';
jsonFile = '\\10.131.30.162\MicroSat-2C\On Orbit Data\Data Processing GUI\v6\Aditya\processed_files.json';

fit_params  = [3400,    175,    1100,   1800,   75,     5400,   200];
 
SUNLIT_DURATION = fit_params(1);
SUNRISE_LEFT    = fit_params(2);
SUNRISE_RIGHT   = fit_params(3);
LINEAR_DURATION = fit_params(4);
SUNSET_RIGHT    = fit_params(5);
ORBIT_DURATION  = fit_params(6);
TRANSITION      = fit_params(7);
  
all_orbit_wise_statistics_UVA = [];
all_orbit_wise_statistics_UVB = [];
all_orbit_wise_statistics_UVC = [];
all_orbit_wise_statistics_temp = [];
all_orbit_timestamp = [];

% Step 1: Choose files using uigetfile with multi-select

[fileNames, filePath] = uigetfile('data/*.*', 'Select One or More Files', 'MultiSelect', 'on');
    
if isequal(fileNames, 0)
    disp('No file selected');
    return;
end
    
% Full file paths (combine with the directory)
fullFilePaths = fullfile(filePath, fileNames);
    
% Step 2: Load the existing JSON file if it exists

if isfile(jsonFile)
    % If the JSON file exists, load its content
    fileStatus = jsondecode(fileread(jsonFile));
else
    % If no JSON file, initialize an empty struct
    fileStatus = struct();
end
    
% Step 3: Process the selected files
    for i_files = 1:length(fullFilePaths)
        fileName = fileNames{i_files};
        
        % Sanitize the file name to make it a valid field name
        validFieldName = matlab.lang.makeValidName(fileName);
        
        % Check if this file has been processed
        if isfield(fileStatus, validFieldName) && fileStatus.(validFieldName)
            disp(['Skipping processed file: ', fileName]);
            continue;  % Skip to the next file
        end
        
        % Process the file (replace with your actual processing code)
        disp(['Processing file: ', fileName]);
        
        % Example of file processing (reading contents, analyzing, etc.)
        % Here, you would insert your actual processing logic

        fileID = fopen(fullFilePaths{i_files}, 'r', 'n', 'UTF-8');
        raw_dat = fread(fileID, '*ubit8');
        fclose(fileID);
        len_space_packet = 130+252;
        num_space_packets = floor(length(raw_dat)/len_space_packet);
        num_data_packets = num_space_packets * 12;
        raw_dat = raw_dat(1:num_space_packets*len_space_packet);
        raw_dat = reshape(raw_dat, [len_space_packet, num_space_packets])';
        raw_timestamp = (single(raw_dat(:,76:103)));
        raw_dat = reshape(raw_dat(:,131:end)', [21 num_data_packets])';
        raw_dat = (single(raw_dat));

        % Processing DN channel wise
        uva = (256*raw_dat(:,3)) + raw_dat(:,4);
        uva(uva==0) = NaN;
        uva = fillmissing(uva, 'nearest');

        uvb = (256*raw_dat(:,7)) + raw_dat(:,8);
        uvb(uvb==0) = NaN;
        uvb = fillmissing(uvb, 'nearest');

        uvc = [ 256*raw_dat(:,1) + raw_dat(:,2), 256*raw_dat(:,5) + raw_dat(:,6), 256*raw_dat(:,9) + raw_dat(:,10), ...
        256*raw_dat(:,11) + raw_dat(:,12), 256*raw_dat(:,13) + raw_dat(:,14), 256*raw_dat(:,17) + raw_dat(:,18)];   
        uvc = reshape(uvc', [6*num_data_packets, 1]);
        uvc(uvc==0) = NaN;
        uvc = fillmissing(uvc, 'nearest');

        uvc_1hz = 256*raw_dat(:,1) + raw_dat(:,2);
        uvc_1hz(uvc_1hz==0) = NaN;
        uvc_1hz = fillmissing(uvc_1hz, 'nearest');

        therm = (256*raw_dat(:,15)) + raw_dat(:,16);
        therm(therm==0) = NaN;
        therm = fillmissing(therm, 'nearest');

        therm_c = [therm therm therm therm therm therm];
        therm_c = reshape(therm_c', [6*num_data_packets, 1]);

        num_space_packets = floor(size(raw_timestamp, 1));
        num_data_packets = num_space_packets * 12;

        year = (raw_timestamp(:,1) + raw_timestamp(:,2)*256);
        month = (raw_timestamp(:,5));
        day = (raw_timestamp(:,9));
        hour = (raw_timestamp(:,13));
        minute = (raw_timestamp(:,17));
        second = (raw_timestamp(:,21));
        millisecond = (17 * raw_timestamp(:,28));

        year = [year year year year year year year year year year year year];
        month = [month month month month month month month month month month month month];
        day = [day day day day day day day day day day day day];
        hour = [hour hour hour hour hour hour hour hour hour hour hour hour];
        minute = [minute minute minute minute minute minute minute minute minute minute minute minute];
        second = [second second+1 second+2 second+3 second+4 second+5 second+6 second+7 second+8 second+9 second+10 second+11];
        millisecond = [millisecond millisecond millisecond millisecond millisecond millisecond millisecond millisecond millisecond millisecond millisecond millisecond];

        year = reshape(year', [num_data_packets,1]);
        month = reshape(month', [num_data_packets,1]);
        day = reshape(day', [num_data_packets,1]);
        hour = reshape(hour', [num_data_packets,1]);
        minute = reshape(minute', [num_data_packets,1]);
        second = reshape(second', [num_data_packets,1]);
        millisecond = reshape(millisecond', [num_data_packets,1]);

        year_c = [year year year year year year];
        month_c = [month month month month month month];
        day_c = [day day day day day day];
        hour_c = [hour hour hour hour hour hour];
        minute_c = [minute minute minute minute minute minute];
        second_c = [second second second second second second];
        millisecond_c = ([millisecond millisecond+200 millisecond+400 millisecond+500 millisecond+600 millisecond+800]);

        year_c = reshape(year_c', [6*num_data_packets,1]);
        month_c = reshape(month_c', [6*num_data_packets,1]);
        day_c = reshape(day_c', [6*num_data_packets,1]);
        hour_c = reshape(hour_c', [6*num_data_packets,1]);
        minute_c = reshape(minute_c', [6*num_data_packets,1]);
        second_c = reshape(second_c', [6*num_data_packets,1]);
        millisecond_c = reshape(millisecond_c', [6*num_data_packets, 1]);

        timestamp = (datetime(year, month, day, hour, minute, second));
        timestamp_c = (datetime(year_c, month_c, day_c, hour_c, minute_c, second_c, millisecond_c, 'Format', 'dd MMMM yyyy, HH:mm:ss.SSS'));
        processsed_data = [];
        % Calculating irradiances
        % Final calibraiton values in mW/cm2
        uva_res = 240*1E-6 * 1.21;
        irrad_uva = uva_res*(uva-6);

        uvb_res = 144*1E-6 * 1.40;
        irrad_uvb = uvb_res*(uvb-3);

        uvc_res = 66.6*1E-6 * 0.90;
        irrad_uvc_1hz = uvc_res*(uvc_1hz-5);
        irrad_uvc = uvc_res*(uvc-5);

        countsToTemperature = @(T_counts) 1./(1.48E-3 + 2.36E-4 * log(T_counts - 29.1)) - 273.15;
        temperature = real(countsToTemperature(therm));
        temperature_c = real(countsToTemperature(therm_c));
        % end
        % Creating tables
        processed_data = timetable(timestamp, uva, uvb, uvc_1hz, therm, irrad_uva, irrad_uvb, irrad_uvc_1hz, temperature, ...
                'VariableNames', {'UVA', 'UVB', 'UVC', 'Thermistor', 'Irradiance_UVA', 'Irradiance_UVB', 'Irradiance_UVC', 'Temperature'});
        processed_data = sortrows(processed_data);
        unique_times = unique(processed_data.timestamp);
        processed_data = retime(processed_data, unique_times);


        processed_data_c = timetable(timestamp_c, uvc, therm_c, irrad_uvc, temperature_c, ...
                'VariableNames', {'UVC', 'Thermistor', 'Irradiance_UVC', 'Temperature'});
        processed_data_c = sortrows(processed_data_c);
        unique_times_c = unique(processed_data_c.timestamp_c);
        processed_data_c = retime(processed_data_c, unique_times_c);

sunrise_id      = [1; find((diff(processed_data.UVA))>10); length(processed_data.UVA)];
sunrise_id      = [1; sunrise_id(diff(sunrise_id) > TRANSITION)];
sunset_id       = [1; find((diff(processed_data.UVA))<-10); length(processed_data.UVA)];
sunset_id       = [sunset_id(diff(sunset_id) > TRANSITION); length(processed_data.UVA)];

% getting valid pairs of sunrise and sunset indicies
valid_sunrise_id    = [];
valid_sunset_id     = [];

for rise_id = 1:length(sunrise_id)
    for set_id = 1:length(sunset_id)
        if ( (sunset_id(set_id) - sunrise_id(rise_id)) < ORBIT_DURATION )
            if ( (sunset_id(set_id) - sunrise_id(rise_id)) > SUNLIT_DURATION )
                valid_sunrise_id    = [valid_sunrise_id sunrise_id(rise_id)];
                valid_sunset_id     = [valid_sunset_id sunset_id(set_id)];
            end
        end
    end
end

orbit_flag  = zeros(size(processed_data.timestamp)); % used to flag irregular orbits

status_flag = zeros(size(processed_data.timestamp)); % Initialize with zeros

for i = 1:length(valid_sunrise_id)
    sunrise_start = max(1, valid_sunrise_id(i) - SUNRISE_LEFT);
    sunrise_stop = sunrise_start;
    linear_stop = sunrise_stop + LINEAR_DURATION;
    sunset_start = linear_stop;
    sunset_stop = valid_sunset_id(i) + SUNSET_RIGHT;
    
    linear_stop = min(linear_stop, length(status_flag));
    sunset_stop = min(sunset_stop, length(status_flag));

    status_flag(sunrise_start:sunrise_stop) = 1;
    status_flag(sunrise_stop:linear_stop) = 2;
    status_flag(sunset_start:sunset_stop) = -1;
end

processed_data.status_flag = status_flag;

% method to check for irregular orbits
for i = 1:length(valid_sunrise_id)
    sunrise_time = valid_sunrise_id(i);
    end_time = sunrise_time + SUNLIT_DURATION;

    if any(sunrise_id > sunrise_time & sunrise_id < end_time)
        orbit_flag(i) = 1;
    end
end

    % Constants for extra values before sunrise and after sunset
    PRE_SUNRISE_WINDOW = 2000;  % Number of values before the sunrise index
    POST_SUNSET_WINDOW = 2000;  % Number of values after the sunset index

 % **NEW: Initialize the orbit_numbers field to NaN (no orbit assigned)**
    orbit_numbers = NaN(size(processed_data.timestamp)); % Initially NaN (no orbit assigned)
    current_orbit_number = 1; % Start orbit numbering from 1

    % **NEW: Assign orbit numbers to each data point based on extended windows**
    for i = 1:length(valid_sunrise_id)
        sunrise_time = valid_sunrise_id(i);
        sunset_time  = valid_sunset_id(i);

        % **Extend the start and end indices to include extra points before and after**
        orbit_start = max(1, sunrise_time - PRE_SUNRISE_WINDOW);   % Ensure start is within bounds
        orbit_end   = min(length(processed_data.timestamp), sunset_time + POST_SUNSET_WINDOW);  % Ensure end is within bounds

        % All indices between the extended range will get the current orbit number
        orbit_numbers(orbit_start:orbit_end) = current_orbit_number;

        % Increment orbit number for the next orbit
        current_orbit_number = current_orbit_number + 1;
    end

    % **NEW: Add the orbit_numbers field to the processed_data struct**
    processed_data.orbit_numbers = orbit_numbers;

    assignin('base', "processed_data", processed_data);

counts_uva = processed_data.UVA;
counts_uvb = processed_data.UVB;
counts_uvc = processed_data.UVC;
t_stamps = processed_data.timestamp;

orbit_wise_statistics_UVA = []; % stores some statistcs of each valid orbit
orbit_wise_statistics_UVB = [];
orbit_wise_statistics_UVC = [];
orbit_timestamp           = []; % stores the start time of each valid orbit

    % Calculate orbit statistics for each file
    for j = 1:length(valid_sunrise_id)
        if orbit_flag(j) == 0 % orbits with no flags (irregularity) is plotted
            valid_sunrise_id_start = valid_sunrise_id(j) + SUNRISE_RIGHT;
            valid_sunrise_end = valid_sunrise_id(j) + LINEAR_DURATION;
            
            duty_perOrbit_UVA = counts_uva(valid_sunrise_id_start:valid_sunrise_end);
            duty_perOrbit_UVB = counts_uvb(valid_sunrise_id_start:valid_sunrise_end);
            duty_perOrbit_UVC = counts_uvc(valid_sunrise_id_start:valid_sunrise_end);

            temp_perOrbit = temperature(valid_sunrise_id_start:valid_sunrise_end);
            
            % Compute the statistics for UVA, UVB, UVC, and Temperature
            x = (1:length(duty_perOrbit_UVA))';
            p = polyfit(x, duty_perOrbit_UVA, 1);
            y_int = polyval(p, x);
            shift_value = p(1) * x + p(2);
            y_shifted = round(duty_perOrbit_UVA - shift_value + y_int(1));
            mean_UVA = mean(y_shifted);
            min_UVA = min(y_shifted);
            max_UVA = max(y_shifted);
            std_UVA = std(y_shifted);
            
            % Repeat for UVB, UVC, and Temperature
            x = (1:length(duty_perOrbit_UVB))';
            p = polyfit(x, duty_perOrbit_UVB, 1);
            y_int = polyval(p, x);
            shift_value = p(1) * x + p(2);
            y_shifted = round(duty_perOrbit_UVB - shift_value + y_int(1));
            mean_UVB = mean(y_shifted);
            min_UVB = min(y_shifted);
            max_UVB = max(y_shifted);
            std_UVB = std(y_shifted);
            
            x = (1:length(duty_perOrbit_UVC))';
            p = polyfit(x, duty_perOrbit_UVC, 1);
            y_int = polyval(p, x);
            shift_value = p(1) * x + p(2);
            y_shifted = round(duty_perOrbit_UVC - shift_value + y_int(1));
            mean_UVC = mean(y_shifted);
            min_UVC = min(y_shifted);
            max_UVC = max(y_shifted);
            std_UVC = std(y_shifted);
            
            mean_temp = mean(temp_perOrbit);
            min_temp = min(temp_perOrbit);
            max_temp = max(temp_perOrbit);
            std_temp = std(temp_perOrbit);
            
            % Store the statistics for this orbit
            all_orbit_wise_statistics_UVA = [all_orbit_wise_statistics_UVA; mean_UVA, min_UVA, max_UVA, std_UVA];
            all_orbit_wise_statistics_UVB = [all_orbit_wise_statistics_UVB; mean_UVB, min_UVB, max_UVB, std_UVB];
            all_orbit_wise_statistics_UVC = [all_orbit_wise_statistics_UVC; mean_UVC, min_UVC, max_UVC, std_UVC];
            all_orbit_wise_statistics_temp = [all_orbit_wise_statistics_temp; mean_temp, min_temp, max_temp, std_temp];
            all_orbit_timestamp = [all_orbit_timestamp; processed_data.timestamp(valid_sunrise_id_start)];
            
        end
    end

orbit_wise_statistics = table(all_orbit_timestamp, all_orbit_wise_statistics_UVA, ...
                              all_orbit_wise_statistics_UVB, all_orbit_wise_statistics_UVC, ...
                              all_orbit_wise_statistics_temp, ...
                              'VariableNames', {'orbit_timestamp', 'UVA_Statistics', 'UVB_Statistics', 'UVC_Statistics', 'Temperature_Statistics'});
                              'VariableNames', {'orbit_timestamp', 'UVA_Statistics', 'UVB_Statistics', 'UVC_Statistics', 'Temperature_Statistics'});


writetable(orbit_wise_statistics, 'C:\Users\DOSIMETER-PC\Downloads\out_statistcs_full_1.xlsx');

        
        % Mark this file as processed
        fileStatus.(validFieldName) = true;
    end
    
    % Step 4: Save the updated status to the JSON file
    % Convert fileStatus struct to JSON and save it (without overwriting, just update)
    fid = fopen(jsonFile, 'w');  % Open in write mode ('w') to overwrite the file
    if fid == -1
        error('Cannot open JSON file for writing.');
    end
    fwrite(fid, jsonencode(fileStatus), 'char');
    fclose(fid);
    
    disp('Processing completed.');



