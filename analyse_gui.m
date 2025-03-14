 function varargout = analyse_gui(varargin)
% ANALYSE_GUI MATLAB code for analyse_gui.fig
%      ANALYSE_GUI, by itself, creates a new ANALYSE_GUI or raises the existing
%      singleton*.
%
%      H = ANALYSE_GUI returns the handle to a new ANALYSE_GUI or the handle to
%      the existing singleton*.
%
%      ANALYSE_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in ANALYSE_GUI.M with the given input arguments.
%
%      ANALYSE_GUI('Property','Value',...) creates a new ANALYSE_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before analyse_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to analyse_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help analyse_gui

% Last Modified by GUIDE v2.5 31-Dec-2024 15:11:30

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @analyse_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @analyse_gui_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end
% End initialization code - DO NOT EDIT


% --- Executes just before analyse_gui is made visible.
function analyse_gui_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to analyse_gui (see VARARGIN)

% Choose default command line output for analyse_gui
handles.output = hObject;

% Program variables
handles.filepath = [];
handles.filepath_sunvec = [];

handles.processed_data = [];
handles.processed_data_c = [];

warning('off','all')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes analyse_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);

% --- Outputs from this function are returned to the command line.
function varargout = analyse_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in popup_good_orbits.
function popup_good_orbits_Callback(hObject, eventdata, handles)
% hObject    handle to popup_good_orbits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_good_orbits contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_good_orbits
contents = cellstr(get(hObject,'String'));
orbit_to_display = contents{get(hObject,'Value')};
handles.good_orbit_display_number = str2double(orbit_to_display);

guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function popup_good_orbits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_good_orbits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in popup_irregular_orbits.
function popup_irregular_orbits_Callback(hObject, eventdata, handles)
% hObject    handle to popup_irregular_orbits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popup_irregular_orbits contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popup_irregular_orbits

contents = cellstr(get(hObject,'String'));
orbit_to_display = contents{get(hObject,'Value')};
handles.irregular_orbit_display_number = str2double(orbit_to_display);
guidata(hObject, handles);
% --- Executes during object creation, after setting all properties.
function popup_irregular_orbits_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popup_irregular_orbits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in re_button.
function analyze_button_Callback(hObject, eventdata, handles)
% hObject    handle to analyze_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

if handles.button_regular_orbits.Value == 1
    
    orbit_flag_condition = 0;
    orbit_number_chosen = handles.popup_good_orbits.Value; % this should be equal to regular_orbit_flag in table
    
    % Given condition should look like this
    
%     SELECT timestamp, UVA, UVB, UVC
%     FROM processed_data
%     WHERE processed_data.orbit_flag == orbit_flag_condition AND processed_data.regular_orbit_flag == orbit_number_chosen

% preprocessed_processed_data(preprocessed_processed_data.orbit_flag == 0 & preprocessed_processed_data.regular_orbit_flag == 2, :)
        
    axes(handles.axes1);
    inst = handles.processed_data(handles.processed_data.orbit_numbers == orbit_number_chosen, :);
    furtinst = handles.processed_data(handles.processed_data.status_flag == 2 & handles.processed_data.orbit_numbers == orbit_number_chosen, :);
    assignin('base', "furtinst", furtinst);
    if handles.radio_UVA.Value == 1
        plot(inst.timestamp, inst.UVA, 'LineWidth',2)
        xlabel("Time", 'FontSize', 14, 'Interpreter', 'latex');
        ylabel("Irradiance (counts)", 'FontSize', 14, 'Interpreter', 'latex')

%         plot(furtinst.timestamp, furtinst.UVA, 'LineWidth',2);
        
        handles.label_uv_mean.String = 'Mean: ' + string(handles.stats_UVA(orbit_number_chosen, 1));
        handles.label_uv_std.String = 'Std: ' + string(handles.stats_UVA(orbit_number_chosen, 4));
        
        handles.label_temp_mean.String = 'Mean: ' + string(handles.stats_temp(orbit_number_chosen, 1));
        handles.label_temp_std.String = 'Std: ' + string(handles.stats_temp(orbit_number_chosen, 4));
        
    elseif handles.radio_UVB.Value == 1
        
        plot(inst.timestamp, inst.UVB, 'LineWidth',2)
        xlabel("Time", 'FontSize', 14, 'Interpreter', 'latex');
        ylabel("Irradiance (counts)", 'FontSize', 14, 'Interpreter', 'latex');

%         plot(furtinst.timestamp, furtinst.UVB, 'LineWidth',2);
        
        handles.label_uv_mean.String = 'Mean: ' + string(handles.stats_UVB(orbit_number_chosen, 1));
        handles.label_uv_std.String = 'Std: ' + string(handles.stats_UVB(orbit_number_chosen, 4));
        
        handles.label_temp_mean.String = 'Mean: ' + string(handles.stats_temp(orbit_number_chosen, 1));
        handles.label_temp_std.String = 'Std: ' + string(handles.stats_temp(orbit_number_chosen, 4));
        
    elseif handles.radio_UVC.Value == 1

        plot(inst.timestamp, inst.UVC, 'LineWidth',2)
        xlabel("Time", 'FontSize', 14, 'Interpreter', 'latex');
        ylabel("Irradiance (counts)", 'FontSize', 14, 'Interpreter', 'latex')
        
%         plot(furtinst.timestamp, furtinst.UVC, 'LineWidth',2);
        
        handles.label_uv_mean.String = 'Mean: ' + string(handles.stats_UVC(orbit_number_chosen, 1));
        handles.label_uv_std.String = 'Std: ' + string(handles.stats_UVC(orbit_number_chosen, 4));
        
        handles.label_temp_mean.String = 'Mean: ' + string(handles.stats_temp(orbit_number_chosen, 1));
        handles.label_temp_std.String = 'Std: ' + string(handles.stats_temp(orbit_number_chosen, 4));
    end
        
elseif handles.radio_irregular_orbits.Value == 1
    orbit_flag_condition = 1;
    orbit_number_chosen = handles.popup_irregular_orbits.Value;
    
    inst = handles.processed_data(handles.processed_data.irregular_orbit_flag == orbit_number_chosen, :);
    if handles.radio_UVA.Value == 1
        plot(inst.timestamp, inst.UVA, 'LineWidth',2);
    elseif handles.radio_UVB.Value == 1
        plot(inst.timestamp, inst.UVB, 'LineWidth',2);
    elseif handles.radio_UVC.Value == 1
        plot(inst.timestamp, inst.radio_UVC, 'LineWidth',2);
    end
end

% --- Executes on button press in load_button.
function load_button_Callback(hObject, eventdata, handles)
% hObject    handle to load_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file, location] = uigetfile('data/*.*');
if file == 0
    disp('File selection cancelled');
else
    handles.filepath= fullfile(location, file);
    handles.label_fileName.String  = file;
end

[handles.processed_data, handles.processed_data_c] = read_dosimeter_data(handles);

guidata(hObject, handles);



% --- Executes on button press in clear_button.
function clear_button_Callback(hObject, eventdata, handles)
% hObject    handle to clear_button (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

handles.filepath = [];
handles.filepath_sunvec = [];

handles.processed_data = [];
handles.processed_data_c = [];

handles.label_fileName.String = '<Data FileName>';

handles.label_uv_mean.String = 'Mean: ';
handles.label_uv_std.String = 'Standard Deviation: ';

handles.label_temp_mean.String = 'Mean: ';
handles.label_temp_std.String = 'Standard Deviation: ';

guidata(hObject, handles);

% --- Executes on button press in button_preprocess.
function button_preprocess_Callback(hObject, eventdata, handles)
% hObject    handle to button_preprocess (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
fit_params  = [3400,    175,    1100,   1800,   75,     5400,   200];
handles = analyse_orbitwise(handles.processed_data, fit_params, handles);
guidata(hObject, handles);

function [processed_data, processed_data_c] = read_dosimeter_data(handles)

%loading the selected file
fileID = fopen(handles.filepath, 'r', 'n', 'UTF-8');
raw_dat = fread(fileID, '*ubit8');

% Reshaping the data according to the space packet structure provided by BDH
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

[timestamp, timestamp_c] = read_timestamp(raw_timestamp);

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


% Processing timestamp
function [timestamp, timestamp_c] = read_timestamp(raw_timestamp)

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

function handles = analyse_orbitwise(~, fit_params, handles)

fit_params  = [3400,    175,    1100,   1800,   75,     5400,   200];

SUNLIT_DURATION = fit_params(1);
SUNRISE_LEFT    = fit_params(2);
SUNRISE_RIGHT   = fit_params(3);
LINEAR_DURATION = fit_params(4);
SUNSET_RIGHT    = fit_params(5);
ORBIT_DURATION  = fit_params(6);
TRANSITION      = fit_params(7);


% initializing indicies for UVA
sunrise_id      = [1; find((diff(handles.processed_data.UVA))>10); length(handles.processed_data.UVA)];
sunrise_id      = [1; sunrise_id(diff(sunrise_id) > TRANSITION)]; %1 at start was there
sunset_id       = [1; find((diff(handles.processed_data.UVA))<-10); length(handles.processed_data.UVA)];
sunset_id       = [sunset_id(diff(sunset_id) > TRANSITION); length(handles.processed_data.UVA)];

sunrise_id = unique(sunrise_id);
sunset_id = unique(sunset_id);

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


status_flag = zeros(size(handles.processed_data.timestamp)); % Initialize with zeros

for i = 1:length(valid_sunrise_id)
    sunrise_start = max(1, valid_sunrise_id(i) - SUNRISE_LEFT);
    sunrise_stop = valid_sunrise_id(i) + SUNRISE_RIGHT; %sunrise_start;
    linear_stop = sunrise_stop + LINEAR_DURATION;
    sunset_start = linear_stop;
    sunset_stop = valid_sunset_id(i) + SUNSET_RIGHT;
    
    linear_stop = min(linear_stop, length(status_flag));
    sunset_stop = min(sunset_stop, length(status_flag));

    status_flag(sunrise_start:sunrise_stop) = 1;
    status_flag(sunrise_stop:linear_stop) = 2;
    status_flag(sunset_start:sunset_stop) = -1;
end

handles.processed_data.status_flag = status_flag;

% method to check for irregular orbits
orbit_flag  = zeros(size(handles.processed_data.timestamp)); % used to flag irregular orbits
orbit_flag_stats = 0 * valid_sunrise_id;

regular_orbit_flag = NaN(size(handles.processed_data.timestamp));
irregular_orbit_flag = NaN(size(handles.processed_data.timestamp));

regular_orbit_number = 0;
irregular_orbit_number = 0;
for i = 1:length(valid_sunrise_id)
    sunrise_time = valid_sunrise_id(i);
    end_time = sunrise_time + SUNLIT_DURATION;

    if any(sunrise_id > sunrise_time & sunrise_id < end_time)
        orbit_flag(sunrise_time:end_time) = 1;
        orbit_flag_stats(i) = 1;
  
        irregular_orbit_number = irregular_orbit_number + 1;
        irregular_orbit_flag(sunrise_time:end_time) = irregular_orbit_number;
    else
        regular_orbit_number = regular_orbit_number + 1;
        regular_orbit_flag(sunrise_time:end_time) = regular_orbit_number;
    end
end

handles.processed_data.orbit_flag = orbit_flag;
handles.processed_data.regular_orbit_flag = regular_orbit_flag;
handles.processed_data.irregular_orbit_flag = irregular_orbit_flag;

    % Constants for extra values before sunrise and after sunset
    PRE_SUNRISE_WINDOW = SUNRISE_LEFT;  % Number of values before the sunrise index
    POST_SUNSET_WINDOW = SUNSET_RIGHT;  % Number of values after the sunset index

 % **NEW: Initialize the orbit_numbers field to NaN (no orbit assigned)**
    orbit_numbers = NaN(size(handles.processed_data.timestamp)); % Initially NaN (no orbit assigned)
    current_orbit_number = 1; % Start orbit numbering from 1
    
    % **NEW: Assign orbit numbers to each data point based on extended windows**
    for i = 1:length(valid_sunrise_id)
        sunrise_time = valid_sunrise_id(i);
        sunset_time  = valid_sunset_id(i);

        % **Extend the start and end indices to include extra points before and after**
        orbit_start = max(1, sunrise_time - PRE_SUNRISE_WINDOW);   % Ensure start is within bounds
        orbit_end   = min(length(handles.processed_data.timestamp), sunset_time + POST_SUNSET_WINDOW);  % Ensure end is within bounds
        
%         if any(sunrise_time > orbit_start & sunset_time < orbit_end)
%             continue;
%         end
        
        % All indices between the extended range will get the current orbit number
        orbit_numbers(orbit_start:orbit_end) = current_orbit_number;

        % Increment orbit number for the next orbit
        current_orbit_number = current_orbit_number + 1;
    end

    % **NEW: Add the orbit_numbers field to the processed_data struct**
    handles.processed_data.orbit_numbers = orbit_numbers;

    %adding statistics

all_orbit_wise_statistics_UVA = [];
all_orbit_wise_statistics_UVB = [];
all_orbit_wise_statistics_UVC = [];
all_orbit_wise_statistics_temp = [];
all_orbit_timestamp = [];

counts_uva = handles.processed_data.UVA;
counts_uvb = handles.processed_data.UVB;
counts_uvc = handles.processed_data.UVC;
temperature   = handles.processed_data.Temperature;

orbit_wise_statistics_UVA = []; % stores some statistcs of each valid orbit
orbit_wise_statistics_UVB = [];
orbit_wise_statistics_UVC = [];
orbit_timestamp           = []; % stores the start time of each valid orbit

    % Calculate orbit statistics for each file
    for j = 1:length(valid_sunrise_id)
        if orbit_flag_stats(j) == 0 % orbits with no flags (irregularity) is plotted
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
        end
    end
    
    handles.stats_UVA = all_orbit_wise_statistics_UVA;
    handles.stats_UVB = all_orbit_wise_statistics_UVB;
    handles.stats_UVC = all_orbit_wise_statistics_UVC;
    handles.stats_temp = all_orbit_wise_statistics_temp;

    %thing to do with gui 
    reg_dropdown_array = [1:regular_orbit_number];
    irreg_dropdown_array = [1:irregular_orbit_number];
    
    reg_options = arrayfun(@num2str, reg_dropdown_array, 'UniformOutput', false);
    irreg_options = arrayfun(@num2str, irreg_dropdown_array, 'UniformOutput', false);
    
    set(handles.popup_good_orbits, 'String', reg_options);
    set(handles.popup_irregular_orbits, 'String', irreg_options);
    
    set(handles.popup_good_orbits, 'Value', 1);
    set(handles.popup_irregular_orbits, 'Value', 1);
    
    assignin('base', 'processed_data', handles.processed_data)
    assignin('base', 'handles', handles); 
    

function label_fileName_Callback(hObject, eventdata, handles)
% hObject    handle to label_fileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of label_fileName as text
%        str2double(get(hObject,'String')) returns contents of label_fileName as a double


% --- Executes during object creation, after setting all properties.
function label_fileName_CreateFcn(hObject, eventdata, handles)
% hObject    handle to label_fileName (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in radio_UVA.
function radio_UVA_Callback(hObject, eventdata, handles)
% hObject    handle to radio_UVA (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_UVA

guidata(hObject, handles);


% --- Executes on button press in radio_UVB.
function radio_UVB_Callback(hObject, eventdata, handles)
% hObject    handle to radio_UVB (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_UVB
guidata(hObject, handles);

% --- Executes on button press in radio_UVC.
function radio_UVC_Callback(hObject, eventdata, handles)
% hObject    handle to radio_UVC (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of radio_UVC
guidata(hObject, handles);


% --- Executes on button press in button_regular_orbits.
function button_regular_orbits_Callback(hObject, eventdata, handles)
% hObject    handle to button_regular_orbits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
isChecked = get(hObject, 'Value');
if isChecked
    set(handles.popup_good_orbits, 'Enable', 'on');
    set(handles.popup_irregular_orbits, 'Enable', 'off');
end
% Hint: get(hObject,'Value') returns toggle state of button_regular_orbits
guidata(hObject, handles);


% --- Executes on button press in radio_irregular_orbits.
function radio_irregular_orbits_Callback(hObject, eventdata, handles)
% hObject    handle to radio_irregular_orbits (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
isChecked = get(hObject, 'Value');
if isChecked
    set(handles.popup_good_orbits, 'Enable', 'off');
    set(handles.popup_irregular_orbits, 'Enable', 'on');
end
% Hint: get(hObject,'Value') returns toggle state of radio_irregular_orbits
guidata(hObject, handles);



function label_temp_mean_Callback(hObject, eventdata, handles)
% hObject    handle to label_temp_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of label_temp_mean as text
%        str2double(get(hObject,'String')) returns contents of label_temp_mean as a double


% --- Executes during object creation, after setting all properties.
function label_temp_mean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to label_temp_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function label_temp_std_Callback(hObject, eventdata, handles)
% hObject    handle to label_temp_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of label_temp_std as text
%        str2double(get(hObject,'String')) returns contents of label_temp_std as a double


% --- Executes during object creation, after setting all properties.
function label_temp_std_CreateFcn(hObject, eventdata, handles)
% hObject    handle to label_temp_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function label_uv_mean_Callback(hObject, eventdata, handles)
% hObject    handle to label_uv_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of label_uv_mean as text
%        str2double(get(hObject,'String')) returns contents of label_uv_mean as a double


% --- Executes during object creation, after setting all properties.
function label_uv_mean_CreateFcn(hObject, eventdata, handles)
% hObject    handle to label_uv_mean (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function label_uv_std_Callback(hObject, eventdata, handles)
% hObject    handle to label_uv_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of label_uv_std as text
%        str2double(get(hObject,'String')) returns contents of label_uv_std as a double


% --- Executes during object creation, after setting all properties.
function label_uv_std_CreateFcn(hObject, eventdata, handles)
% hObject    handle to label_uv_std (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
