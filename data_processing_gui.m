
% v0 : 01-08-2024 : Baseline program
% v1 : 12-08-2024 : Timestamping implemented
% v2 : 14-06-2024 : Sun Vector computation
% v3 : 17-08-2024 : Changes made at URSC for data representation
% v4 : 23-08-2024 : Timestamping scheme updated. Analyse functionality
% added. Temperature Compensation added.
% v5 : 10-10-2024 : Added functionality to handle sunvec computation.
% Processing time is improved.
% v6 : 04-11-2024 : Orbitwise statistics and further computation has been
% added
% v7 : 10-11-2024 : Mandated the need for both DOSM and SUNVEC files
% Exporting of data is improved.

% Pending activities - Temperature compensation, Elevation angle
% compensation, Offset removal

% Save the data file as is
% Save the sunvec file with .dat extension

function varargout = data_processing_gui(varargin)
% DATA_PROCESSING_GUI MATLAB code for data_processing_gui.fig
%      DATA_PROCESSING_GUI, by itself, creates a new DATA_PROCESSING_GUI or raises the existing
%      singleton*.
%
%      H = DATA_PROCESSING_GUI returns the handle to a new DATA_PROCESSING_GUI or the handle to
%      the existing singleton*.
%
%      DATA_PROCESSING_GUI('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DATA_PROCESSING_GUI.M with the given input arguments.
%
%      DATA_PROCESSING_GUI('Property','Value',...) creates a new DATA_PROCESSING_GUI or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before data_processing_gui_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to data_processing_gui_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help data_processing_gui

% Last Modified by GUIDE v2.5 23-Aug-2024 15:24:15

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @data_processing_gui_OpeningFcn, ...
                   'gui_OutputFcn',  @data_processing_gui_OutputFcn, ...
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



% --- Executes just before data_processing_gui is made visible.
function data_processing_gui_OpeningFcn(hObject, eventdata, handles, varargin) %#ok<*INUSL>
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to data_processing_gui (see VARARGIN)

% Choose default command line output for data_processing_gui
handles.output = hObject;

% Program variables
handles.filepath = [];
handles.filepath_sunvec = [];

handles.processed_data = [];
handles.processed_data_c = [];

warning('off','all')

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes data_processing_gui wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = data_processing_gui_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in button_open_data.
function button_open_data_Callback(hObject, eventdata, handles) %#ok<*DEFNU>
% hObject    handle to button_open_data (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file, location] = uigetfile('data/*.*');
if file == 0
    disp('File selection cancelled');
else
    handles.filepath= fullfile(location, file);
    handles.label_filename.String  = file;
    guidata(hObject, handles);
end


% --- Executes on button press in button_open_sunvec.
function button_open_sunvec_Callback(hObject, eventdata, handles)
% hObject    handle to button_open_sunvec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[file, location] = uigetfile('*.*');
if file == 0
    disp('File selection cancelled');
else
    handles.filepath_sunvec = fullfile(location, file);
    handles.label_filename_sunvec.String  = file;
    guidata(hObject, handles);
end


% --- Executes on button press in button_load.
function button_load_Callback(hObject, eventdata, handles)
% hObject    handle to button_load (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% if sunvec file is not selected during initializaiton.
if isempty(handles.filepath_sunvec)
    disp('Sunvec file not loaded');
    return;

else
    
    % takes the regular expression of the strings of filenames to check if
    % it is the same orbit number.

    % note that dosm file has six integer values and sunvec has five.

    if isequal(str2double(regexp(handles.label_filename.String, '\d{6}', 'match')), (str2double(regexp(handles.label_filename_sunvec.String, '\d{5}', 'match'))))
        
        [handles.processed_data, handles.processed_data_c] = read_dosimeter_data(handles);
    else
        disp('Sunvec file is not of same orbit')
        return;
    end
end

handles = angle_compensation(handles);

% if handles.toggle_temperature_compensation.Value == 1 
%     
%     uva_res_tc = 240E-6 * 1.21;
%     irrad_uva_tc =  uva_res_tc * round((handles.processed_data.UVA - 7) ./ (1 + ((handles.processed_data.Temperature - 22)*0.00282)));
%     
%     uvb_res_tc = 144E-6 * 1.40;
%     irrad_uvb_tc = uvb_res_tc * round((handles.processed_data.UVB - 3) ./ (1 + ((handles.processed_data.Temperature - 22)*0.00055)));
% 
%     uvc_res_tc = 66.6E-6 * 0.90;
%     irrad_uvc_tc =  uvc_res_tc * round((handles.processed_data_c.UVC - 5) ./ (1 - (handles.processed_data_c.Temperature - 19)*0.0008));
%     irrad_uvc_tc_1hz =  uvc_res_tc * round((handles.processed_data.UVC - 5) ./ (1 - (handles.processed_data.Temperature - 19)*0.0008));
%     
%     handles.processed_data.UVA_T_comp = irrad_uva_tc;
%     handles.processed_data.UVB_T_comp = irrad_uvb_tc;
%     handles.processed_data.UVC_T_comp = irrad_uvc_tc_1hz;
%     
%     handles.processed_data_c.UVC_T_comp = irrad_uvc_tc;
% end


% Plotting the data
plot_data(hObject, handles);

assignin('base', 'processed_data', handles.processed_data);
assignin('base', 'processed_data_c', handles.processed_data_c);

guidata(hObject, handles);



% --- Executes on button press in button_clear.
function button_clear_Callback(hObject, eventdata, handles)
% hObject    handle to button_clear (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

axes(handles.axes_uva);
cla reset
axes(handles.axes_uvb);
cla reset
axes(handles.axes_uvc);
cla reset
axes(handles.axes_therm);
cla reset

handles.filepath = [];
handles.filepath_sunvec = [];
handles.label_filename.String  = "<Filename>";
handles.label_filename_sunvec.String  = "<SunVec Filename>";
handles.text_stage.String  = [];

handles.counts_uva = [];
handles.counts_uvb = [];
handles.counts_uvc = [];
handles.counts_therm = [];
handles.tm_data = [];

handles.irrad_uva = [];
handles.irrad_uvb = [];
handles.irrad_uvc = [];
handles.temperature_data = [];


handles.sunvec_a = [];
handles.sunvec_b = [];
handles.sunvec_c = [];

handles.elevation_a = [];
handles.elevation_b = [];
handles.elevation_c = [];

handles.timestamp_a = [];
handles.timestamp_b = [];
handles.timestamp_c = [];
handles.timestamp_d = [];


guidata(hObject, handles);



function text_stage_Callback(hObject, eventdata, handles)

% --- Executes during object creation, after setting all properties.
function text_stage_CreateFcn(hObject, eventdata, handles)
% hObject    handle to text_stage (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in button_export.
function button_export_Callback(hObject, eventdata, handles)
% hObject    handle to button_export (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


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
sunrise_id      = [1; sunrise_id(diff(sunrise_id) > TRANSITION)];
sunset_id       = [1; find((diff(handles.processed_data.UVA))<-10); length(handles.processed_data.UVA)];
sunset_id       = [sunset_id(diff(sunset_id) > TRANSITION); length(handles.processed_data.UVA)];

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

orbit_flag  = zeros(size(handles.processed_data.timestamp)); % used to flag irregular orbits

status_flag = zeros(size(handles.processed_data.timestamp)); % Initialize with zeros

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

handles.processed_data.status_flag = status_flag;

% method to check for irregular orbits
for i = 1:length(valid_sunrise_id)
    sunrise_time = valid_sunrise_id(i);
    end_time = sunrise_time + SUNLIT_DURATION;

    if any(sunrise_id > sunrise_time & sunrise_id < end_time)
        orbit_flag(sunrise_time:end_time) = 1;
    end
end
handles.processed_data.orbit_flag = orbit_flag;
    % Constants for extra values before sunrise and after sunset
    PRE_SUNRISE_WINDOW = 2000;  % Number of values before the sunrise index
    POST_SUNSET_WINDOW = 2000;  % Number of values after the sunset index

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

        % All indices between the extended range will get the current orbit number
        orbit_numbers(orbit_start:orbit_end) = current_orbit_number;

        % Increment orbit number for the next orbit
        current_orbit_number = current_orbit_number + 1;
    end

    % **NEW: Add the orbit_numbers field to the processed_data struct**
    handles.processed_data.orbit_numbers = orbit_numbers;

    assignin('base', 'processed_data_latest', handles.processed_data);
    

% originally meant for the counts, using irradiance values instead.    
counts_uva = handles.processed_data.Irradiance_UVA;
counts_uvb = handles.processed_data.Irradiance_UVB;
counts_uvc = handles.processed_data.Irradiance_UVC;
t_stamps = handles.processed_data.timestamp;
temperature = handles.processed_data.Temperature;

orbit_wise_statistics_UVA = []; % stores some statistcs of each valid orbit
orbit_wise_statistics_UVB = [];
orbit_wise_statistics_UVC = [];
orbit_wise_statistics_temp = [];
orbit_timestamp           = []; % stores the start time of each valid orbit

all_orbit_wise_statistics_UVA = [];
all_orbit_wise_statistics_UVB = [];
all_orbit_wise_statistics_UVC = [];
all_orbit_wise_statistics_temp = [];
all_orbit_timestamp = [];


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
            y_shifted = duty_perOrbit_UVA - shift_value + y_int(1);
            mean_UVA = mean(y_shifted);
            min_UVA = min(y_shifted);
            max_UVA = max(y_shifted);
            std_UVA = std(y_shifted);
            
            % Repeat for UVB, UVC, and Temperature
            x = (1:length(duty_perOrbit_UVB))';
            p = polyfit(x, duty_perOrbit_UVB, 1);
            y_int = polyval(p, x);
            shift_value = p(1) * x + p(2);
            y_shifted = duty_perOrbit_UVB - shift_value + y_int(1);
            mean_UVB = mean(y_shifted);
            min_UVB = min(y_shifted);
            max_UVB = max(y_shifted);
            std_UVB = std(y_shifted);
            
            x = (1:length(duty_perOrbit_UVC))';
            p = polyfit(x, duty_perOrbit_UVC, 1);
            y_int = polyval(p, x);
            shift_value = p(1) * x + p(2);
            y_shifted = duty_perOrbit_UVC - shift_value + y_int(1);
            mean_UVC = mean(y_shifted);
            min_UVC = min(y_shifted);
            max_UVC = max(y_shifted);
            std_UVC = std(y_shifted);
            
            mean_temp = mean(temp_perOrbit);
            min_temp = min(temp_perOrbit);
            max_temp = max(temp_perOrbit);
            std_temp = std(temp_perOrbit);
            
            % Store the statistics for this orbit
            orbit_wise_statistics_UVA = [orbit_wise_statistics_UVA; mean_UVA, min_UVA, max_UVA, std_UVA];
            orbit_wise_statistics_UVB = [orbit_wise_statistics_UVB; mean_UVB, min_UVB, max_UVB, std_UVB];
            orbit_wise_statistics_UVC = [orbit_wise_statistics_UVC; mean_UVC, min_UVC, max_UVC, std_UVC];
            orbit_wise_statistics_temp = [orbit_wise_statistics_temp; mean_temp, min_temp, max_temp, std_temp];
            orbit_timestamp = [orbit_timestamp; handles.processed_data.timestamp(valid_sunrise_id_start)];
            
        end
    end

orbit_wise_statistics = table(orbit_timestamp, orbit_wise_statistics_UVA, orbit_wise_statistics_UVB,...
    orbit_wise_statistics_UVC, orbit_wise_statistics_temp,...
    'VariableNames', {'orbit_timestamp', 'UVA_Statistics', 'UVB_Statistics', 'UVC_Statistics', 'Temperature_Statistics'});


% writetable(orbit_wise_statistics, 'C:\Users\DOSIMETER-PC\Downloads\out_statistcs_full_1.xlsx');

    
    guidata(hObject, handles);
   
% this file stores the files that have been already exported. 
jsonFile = pwd + "\exported_files.json";

if isfile(jsonFile)
    fileStatus = jsondecode(fileread(jsonFile));
else
    fileStatus = struct();
end

% filename = regexp(handles.label_filename.String, '\d{6}', 'match');
filename = handles.label_filename.String;

validFieldName = matlab.lang.makeValidName(filename);

if isfield(fileStatus, validFieldName) && fileStatus.(validFieldName)
    disp('File has already been exported.');
    return;
else

compressed_filename = regexp(handles.label_filename.String, '\d{6}', 'match');

table_1hz = timetable2table(handles.processed_data);
table_6hz = timetable2table(handles.processed_data_c);
table_stats = (orbit_wise_statistics);

writefile = pwd + "\processed_data\" + string(compressed_filename);
writestats = pwd + "\statistics\" + string(compressed_filename);

writetable(table_1hz, writefile + "_1hz.csv");
writetable(table_6hz, writefile + "_6hz.csv");
writetable(table_stats, writestats + "_stat.csv");
fileStatus.(validFieldName) = true;

end

fid = fopen(jsonFile, 'w');
if fid == -1
    error('Cannot open JSON file for writing.');
end
fwrite(fid, jsonencode(fileStatus), 'char');
fclose(fid);

guidata(hObject, handles);


function label_filename_Callback(hObject, eventdata, handles)
% hObject    handle to label_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of label_filename as text
%        str2double(get(hObject,'String')) returns contents of label_filename as a double


% --- Executes during object creation, after setting all properties.
function label_filename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to label_filename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --------------------------------------------------------------------

% --- Executes during object creation, after setting all properties.
function label_filename_sunvec_CreateFcn(hObject, eventdata, handles)
% hObject    handle to label_filename_sunvec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function label_filename_sunvec_Callback(hObject, eventdata, handles)
% hObject    handle to label_filename_sunvec (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of label_filename_sunvec as text
%        str2double(get(hObject,'String')) returns contents of label_filename_sunvec as a double
    
   
% --- Executes on button press in toggle_angle_compensation.
function toggle_angle_compensation_Callback(hObject, eventdata, handles)
% hObject    handle to toggle_angle_compensation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of toggle_angle_compensation

function separatePlot(handles)

% Plotting the data

uva = figure;
color_uva = [0.6350 0.0780 0.1840];
uva.Position(3:4) = [800 400];
plot(handles.timestamp_a, handles.irrad_uva, 'Linewidth', 2, 'Color', color_uva );
title("UV A", 'FontSize', 14, 'Interpreter', 'latex')
xlabel("Time", 'FontSize', 14, 'Interpreter', 'latex')
% ylim((uva_lim-7)*uva_res)
% yticks((uva_ticks-7)*uva_res)
ylim([0, 5])
% yticks()
grid on
text(0.9, 0.7, sprintf("Mean = %0.3f $mW/cm^2$", mean(handles.irrad_uva)), 'Units', 'normalized', 'FontSize', 14, 'HorizontalAlignment', 'right', 'Interpreter', 'latex')
% text(0.9, 0.6, sprintf("3$\sigma$ = %0.3f $mW/cm^2$", 3*std(handles.irrad_uva)), 'Units', 'normalized', 'FontSize', 14, 'HorizontalAlignment', 'right', 'Interpreter', 'latex')
ylabel("Irradiance $(mW/cm^{2})$", 'FontSize', 14, 'Interpreter', 'latex')

uvb = figure;
color_uvb = [0.4660 0.6740 0.1880];
uvb.Position(3:4) = [800 400];
plot(handles.timestamp_b, handles.irrad_uvb, 'Linewidth', 2, 'Color', color_uvb);
title("UV B", 'FontSize', 14, 'Interpreter', 'latex')
xlabel("Time", 'FontSize', 14, 'Interpreter', 'latex')
ylim([0 2.5])
grid on
text(0.9, 0.6, sprintf("Mean = %0.3f $mW/cm^2$", mean(handles.irrad_uvb)), 'Units', 'normalized', 'FontSize', 14, 'HorizontalAlignment', 'right', 'Interpreter', 'latex')
ylabel("Irradiance $(mW/cm^{2})$", 'FontSize', 14, 'Interpreter', 'latex')

uvc = figure;
color_uvc = [0.4940 0.1840 0.5560];
uvc.Position(3:4) = [800 400];
plot(handles.timestamp_c, handles.irrad_uvc, 'Linewidth', 2, 'Color', color_uvc );
title("UV C", 'FontSize', 14, 'Interpreter', 'latex', 'FontSize', 14)
xlabel("Time", 'FontSize', 14, 'Interpreter', 'latex')
ylim([0 0.15])
grid on
text(0.9, 0.7, sprintf("Mean = %0.3f $mW/cm^2$", mean(handles.irrad_uvc)), 'Units', 'normalized', 'FontSize', 14, 'HorizontalAlignment', 'right', 'Interpreter', 'latex')
ylabel("Irradiance $(mW/cm^{2})$", 'FontSize', 14, 'Interpreter', 'latex')

uvd = figure;
color_uvd = [0.9290 0.6940 0.1250];
uvd.Position(3:4) = [800 400];
plot(handles.timestamp_d, handles.temperature_data, 'Linewidth', 2, 'Color', color_uvd );
title("Thermistor", 'FontSize', 14, 'Interpreter', 'latex')
xlabel("Time", 'FontSize', 14, 'Interpreter', 'latex')
ylim([-10 50])
yticks(-10:10:50)
grid on
text(0.9, 0.7, sprintf("Mean = %0.3f °C", mean(handles.temperature_data)), 'Units', 'normalized', 'FontSize', 14, 'HorizontalAlignment', 'right', 'Interpreter', 'latex')
ylabel("Temperature $(C)$", 'FontSize', 14, 'Interpreter', 'latex')


% --- Executes on button press in pushbutton8.
function pushbutton8_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton8 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

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

%loading sunvec directly within the same load as data


% Processing timestamp
function [timestamp, timestamp_c] = read_timestamp(raw_timestamp, handles)

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


function plot_data(hObject, handles)

axes(handles.axes_uva);
cla reset
color_uva = [0.6350 0.0780 0.1840];
plot(handles.processed_data.timestamp, handles.processed_data.Irradiance_UVA, 'Linewidth', 2, 'Color', color_uva );
title("UV A", 'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'tex')
xlabel("Time", 'FontWeight', 'bold', 'FontSize', 16, 'Interpreter', 'tex')
ylim([0, 6])
grid on
ylabel("{Irradiance (mW/cm^{2})}", 'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'tex')
set(gca, 'FontWeight', 'bold')
set(gca, 'FontSize', 13)

axes(handles.axes_uvb);
cla reset
color_uvb = [0.4660 0.6740 0.1880];
plot(handles.processed_data.timestamp, handles.processed_data.Irradiance_UVB, 'Linewidth', 2, 'Color', color_uvb);
title("UV B", 'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'tex')
xlabel("Time", 'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'tex')
ylim([0 3])
grid on
ylabel("{Irradiance (mW/cm^{2})}", 'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'tex')
set(gca, 'FontWeight', 'bold')
set(gca, 'FontSize', 13)

axes(handles.axes_uvc);
cla reset
color_uvc = [0.4940 0.1840 0.5560];
plot(handles.processed_data_c.timestamp_c, handles.processed_data_c.Irradiance_UVC, 'Linewidth', 2, 'Color', color_uvc );
title("UV C", 'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'tex', 'FontSize', 14)
xlabel("Time", 'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'tex')
ylim([0 0.1])
grid on
ylabel("{Irradiance (mW/cm^{2})}", 'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'tex')
set(gca, 'FontWeight', 'bold')
set(gca, 'FontSize', 13)

axes(handles.axes_therm);
cla reset
color_uvd = [0.9290 0.6940 0.1250];
plot(handles.processed_data.timestamp, handles.processed_data.Temperature, 'Linewidth', 2, 'Color', color_uvd );
title("Thermistor", 'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'tex')
xlabel("Time", 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'tex')
ylim([-10 50])
yticks(-10:10:50)
grid on
text(0.9, 0.65, sprintf("Mean = %0.3f °C", mean(handles.processed_data.Temperature)), 'Units', 'normalized', 'FontSize', 14, 'HorizontalAlignment', 'right', 'Interpreter', 'tex')
ylabel("{Temperature (°C)}", 'FontWeight', 'bold', 'FontSize', 14, 'Interpreter', 'tex')
set(gca, 'FontWeight', 'bold')
set(gca, 'FontSize', 13)

% if handles.toggle_temperature_compensation.Value == 1
%     axes(handles.axes_uva);
%     hold on
%     plot(handles.processed_data.timestamp, handles.processed_data.UVA_T_comp, '--', 'Linewidth', 2, 'Color', color_uva);
%     
%     axes(handles.axes_uvb);
%     hold on
%     plot(handles.processed_data.timestamp, handles.processed_data.UVB_T_comp, '--', 'Linewidth', 2, 'Color', color_uvb);
%     
%     axes(handles.axes_uvc);
%     hold on
%     plot(handles.processed_data_c.timestamp_c, handles.processed_data_c.UVC_T_comp, '--', 'Linewidth', 2, 'Color', color_uvc);
% end

    axes(handles.axes_therm);
    hold on
    yyaxis('right')
    ylim([-5 180]);
    yticks(0:30:180);
    title("Thermistor / Sun Angle", 'FontSize', 14, 'FontWeight', 'bold', 'Interpreter', 'tex')
    ylabel("Angle (°)")
    
    plot(handles.processed_data.timestamp, 180/pi*acos(handles.processed_data.Pitch), 'Linewidth', 2, 'Color', 'cyan');
    legend(["Temperature", "Sun Zenith Angle"])

guidata(hObject, handles);


function handles = angle_compensation(handles)
% raw_sunvec = datastore(handles.filepath_sunvec, 'Type', 'tabulartext', ...
%     'VariableNames', {'Year', 'Month', 'Day', 'Hour', 'Minute', 'Second', 'Nanosecond', 'Yaw', 'Roll', 'Pitch'});
% raw_sunvec = tall(raw_sunvec);

% raw_sunvec = readtable(handles.filepath_sunvec, 'FileType', 'delimitedtext');

raw_sunvec = readtable(handles.filepath_sunvec);
raw_sunvec.Properties.VariableNames = {'Year', 'Month', 'Day', 'Hour', 'Minute', 'Second', 'Nanosecond', 'Yaw', 'Roll', 'Pitch'};
raw_sunvec.Time = datetime(raw_sunvec.Year, raw_sunvec.Month, raw_sunvec.Day, raw_sunvec.Hour, raw_sunvec.Minute, raw_sunvec.Second);  

raw_sunvec(:,1:9) = [];
raw_sunvec.Pitch = -raw_sunvec.Pitch;
id = (raw_sunvec.Pitch > 1) | (raw_sunvec.Pitch < -1) ;
raw_sunvec.Pitch(id) =  NaN;

raw_sunvec = table2timetable(raw_sunvec);
raw_sunvec = sortrows(raw_sunvec);

tt_sunvec = retime(raw_sunvec, handles.processed_data.timestamp, 'nearest');
tt_sunvec_c = retime(raw_sunvec, handles.processed_data_c.timestamp_c, 'nearest');

assignin('base', "tt_sunvec", tt_sunvec);

handles.processed_data = [handles.processed_data, tt_sunvec];
handles.processed_data_c = [handles.processed_data_c, tt_sunvec_c];

assignin('base', "processed_data", handles.processed_data);
assignin('base', 'handles', handles);

% --- Executes on button press in toggle_temperature_compensation.
% function toggle_temperature_compensation_Callback(hObject, eventdata, handles)
