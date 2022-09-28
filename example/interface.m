function varargout = interface(varargin)
% INTERFACE MATLAB code for interface.axis
%      INTERFACE, by itself, creates a new INTERFACE or raises the existing
%      singleton*.
%
%      H = INTERFACE returns the handle to a new INTERFACE or the handle to
%      the existing singleton*.
%
%      INTERFACE('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in INTERFACE.M with the given input arguments.
%
%      INTERFACE('Property','Value',...) creates a new INTERFACE or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before interface_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to interface_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help interface

% Last Modified by GUIDE v2.5 14-Sep-2022 22:02:55

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
    'gui_Singleton',  gui_Singleton, ...
    'gui_OpeningFcn', @interface_OpeningFcn, ...
    'gui_OutputFcn',  @interface_OutputFcn, ...
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
end

% --- Executes just before interface is made visible.
function interface_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to interface (see VARARGIN)

% Choose default command line output for interface
handles.output = hObject;

% Update handles structure
guidata(hObject, handles);

% UIWAIT makes interface wait for user response (see UIRESUME)
% uiwait(handles.figure1);
end

% --- Outputs from this function are returned to the command line.
function varargout = interface_OutputFcn(hObject, eventdata, handles)
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;
end

function resp = runCommand(txt, hObject, eventdata, handles)
write(handles.userdata.cli, uint8(sprintf('%s\n', txt)));
resp = [];
while true
    t = read(handles.userdata.cli, 1);
    if t == 0
        break;
    end
    resp = [resp, t];
end
resp = char(resp);
end

function gettmpfile(hObject, eventdata, handles)
% disp('send command');
write(handles.userdata.cli, uint8(sprintf('%s\n', 'gettmp')));
pause(1);
fsize = [];
while true
    t = read(handles.userdata.cli, 1);
%     disp(t);
    if t == 0
        break;
    end
    fsize = [fsize, t];
end
% disp(fsize);
fsize = char(fsize);
% disp(fsize);
fn = str2double(fsize);
% fprintf('filesize: %d\n', fn);
% buffer = zeros(fn, 1, 'uint8');
% for i = 1:fn
%     buffer(i) = read(handles.userdata.cli, 1);
% end
buffer = read(handles.userdata.cli, fn, 'uint8');
% disp('tmp file downloaded');
handles.userdata.figpath = '.tmpplot.mat';
guidata(hObject, handles);
fid = fopen(handles.userdata.figpath, 'w');
fwrite(fid, buffer, 'uint8');
fclose(fid);
read(handles.userdata.cli, 1);
end

function l = autoscale(x)
c = (max(x)+min(x))/2;
r = (max(x)-min(x))/2;
l = [-1, 1]*1.2*r + c;
end

function plotstationscatter(hObject, eventdata, handles)
runCommand('plot station', hObject, eventdata, handles);
if ~strcmp(handles.ServerIP.String, '127.0.0.1') && ~strcmp(handles.ServerIP.String, 'localhost')
    gettmpfile(hObject, eventdata, handles);
end
load(handles.userdata.figpath, 'lon1', 'lat1', 'tag1', 'lon2', 'lat2', 'tag2', ...
    'elon', 'elat', 'type');
flag = true;
for t = {'lon1', 'lat1', 'tag1', 'lon2', 'lat2', 'tag2', 'elon', 'elat', 'type'}
    flag = flag && exist(t{1}, 'var');
end
if ~flag
    disp('variable not enough');
    return;
end
if strcmp(type, 'station')
    cla(handles.axis);
    scatter(lon1, lat1, 'k^', 'markerfacecolor', 'k', 'markeredgecolor', 'k', 'sizedata', 40); hold on;
    scatter(lon2, lat2, 'k^', 'markerfacecolor', 'w', 'markeredgecolor', 'k', 'sizedata', 40);
    scatter(elon, elat, 'rp', 'markerfacecolor', 'r', 'markeredgecolor', 'r', 'sizedata', 60);
    text([lon1; lon2], [lat1; lat2], [tag1; tag2], ...
        'HorizontalAlignment', 'center',...
        'VerticalAlignment', 'top');
    set(handles.axis, 'XLim', autoscale([lon1;lon2]), 'YLim', autoscale([lat1;lat2]), ...
        'XScale', 'linear', 'Box', 'on', 'FontSize', 14, 'FontName', 'Georgia');
    handles.userdata.curfigtype = type;
    guidata(hObject, handles);
end
end

function [x, y] = rect(x1, x2, y1, y2)
x = [x1, x1, x2, x2, x1];
y = [y1, y2, y2, y1, y1];
end

function plotwave(hObject, eventdata, handles)
csta = handles.stationlistbox.String{handles.stationlistbox.Value};
runCommand(sprintf('csta %s', csta), hObject, eventdata, handles);
runCommand('plot wave', hObject, eventdata, handles);
if ~strcmp(handles.ServerIP.String, '127.0.0.1') && ~strcmp(handles.ServerIP.String, 'localhost')
    gettmpfile(hObject, eventdata, handles);
end
load(handles.userdata.figpath, 'tl', 'ww', 'tg', 'wg', 'sw', 'sg', 'type', 'bt', 'win');
flag = true;
for t = {'tl', 'ww', 'tg', 'wg', 'sw', 'sg', 'type', 'bt', 'win'}
    flag = flag && exist(t{1}, 'var');
end
if ~flag
    disp('variable not enough');
    return;
end
phasetype = {};
phaseat = [];
phasett = [];
trims = {};
st = 0.02;
if strcmp(type, 'wave')
    cla(handles.axis);
    btg = tg{1}(1);
    for i = 1:length(win)
        t = win{i};
        phasetype = [phasetype, t.type];
        phaseat = [phaseat, t.at];
        phasett = [phasett, t.tt + btg];
        for j = 1:3
            fn = fieldnames(t);
            for k = 1:length(fn)
                if contains(fn{k}, 'trim')
                    trims = [trims, fn{k}];
                    ttag = strrep(fn{k}, '_trim', '');
                    w = t.(fn{k});
                    [wx, wy] = rect(t.at+w(1), t.at+w(2), 0.5, 3.5+st*k);
                    text(t.at+w(1), 3.5+st*k, ['<-', ttag], 'FontName', 'Georgia', 'FontSize', 10, 'Rotation', 90, ...
                        'Interpreter', 'none');
                    text(t.at+w(2), 3.5+st*k, ['<-', ttag], 'FontName', 'Georgia', 'FontSize', 10, 'Rotation', 90, ...
                        'Interpreter', 'none');
                    fill(wx, wy, 'k', 'EdgeColor', 'w', 'FaceAlpha', 0.02, 'EdgeAlpha', 0.02);
                    [wx, wy] = rect(t.tt+w(1)+btg, t.tt+w(2)+btg, 0.5-st*k, 3.5);
                    fill(wx, wy, 'k', 'EdgeColor', 'w', 'FaceAlpha', 0.02, 'EdgeAlpha', 0.02);
                    text(t.tt+w(1)+btg, 0.5-st*k, [ttag, '->'], 'FontName', 'Georgia', 'FontSize', 10, 'Rotation', 90, ...
                        'HorizontalAlignment', 'right', 'Interpreter', 'none');
                    text(t.tt+w(2)+btg, 0.5-st*k, [ttag, '->'], 'FontName', 'Georgia', 'FontSize', 10, 'Rotation', 90, ...
                        'HorizontalAlignment', 'right', 'Interpreter', 'none');
                end
            end
        end
    end
    for i = 1:length(win)
        t = win{i};
        for j = 1:3
            line([1, 1]*t.at, [-1, 1]*0.2+sw(j), 'Color', 'r');
            line([1, 1]*t.tt+btg, [-1, 1]*0.2+sg(j), 'Color', 'r');
        end
    end
    mint = 0;
    maxt = 0;
    for i = 1:length(tl)
        line(tl{i}, ww{i}, 'Color', [1, 1, 1]*0.5);
        line(tg{i}, wg{i}, 'Color', [0.4, 0.4, 0.7]);
        if mint > min(tl{i})
            mint = min(tl{i});
        end
        if maxt < max(tl{i})
            maxt = max(tl{i});
        end
        if mint > min(tg{i})
            mint = min(tg{i});
        end
        if maxt < max(tg{i})
            maxt = max(tg{i});
        end
    end
    for i = 1:length(win)
        t = win{i};
        for j = 1:3
            line([1, 1]*t.at, [-1, 1]*0.3+sw(j), 'Color', 'r');
            line([1, 1]*t.tt+btg, [-1, 1]*0.3+sg(j), 'Color', 'r');
        end
    end
    set(handles.axis, 'XLim', [mint, maxt], 'YLim', [0, 4], 'XScale', 'linear', 'Box', 'on', 'FontSize', 14, ...
        'FontName', 'Georgia');
    v1 = str2double(handles.xlim1.String);
    v2 = str2double(handles.xlim2.String);
    set(handles.axis, 'XLim', [v1, v2]);
    handles.userdata.curfigtype = type;
    handles.userdata.gbtime = btg;
    handles.userdata.btime = datetime(bt, 'InputFormat', 'yyyy-MM-dd''T''HH:mm:ss.SSS');
    handles.userdata.attime = phaseat;
    handles.userdata.tttime = phasett;
    handles.pickphasetype.String = phasetype;
    handles.listbox_trimtype.String = unique(trims);
    guidata(hObject, handles);
end
end

function plotspec(hObject, eventdata, handles)
runCommand('plot spec', hObject, eventdata, handles);
if ~strcmp(handles.ServerIP.String, '127.0.0.1') && ~strcmp(handles.ServerIP.String, 'localhost')
    gettmpfile(hObject, eventdata, handles);
end
load(handles.userdata.figpath, 'band', 'fs', 'spec', 'type');
flag = true;
for t = {'band', 'fs', 'spec', 'type'}
    flag = flag && exist(t{1}, 'var');
end
if ~flag
    disp('variable not enough');
    return;
end
if strcmp(type, 'spec')
    cla(handles.axis);
    maxfs = 0.0;
    minfs = Inf;
    bandtype = {};
    for i = 1:length(band)
        t = band{i};
        [wx, wy] = rect(t.band(1), t.band(2), 1-0.02*i, 4);
        % line(wx, wy, 'Color', 'r');
        fill(wx, wy, 'r', 'EdgeColor', 'w', 'FaceAlpha', 0.02);
        text(t.band(1), 1-0.02*i, sprintf('%s/%s ->', t.phase, t.key), 'HorizontalAlignment', 'right', ...
            'FontName', 'Georgia', 'FontSize', 10, 'Rotation', 30, 'Interpreter', 'none');
        text(t.band(2), 1-0.02*i, sprintf('%s/%s ->', t.phase, t.key), 'HorizontalAlignment', 'right', ...
            'FontName', 'Georgia', 'FontSize', 10, 'Rotation', 30, 'Interpreter', 'none');
        bandtype = [bandtype, sprintf('%s %s', t.phase, t.key)];
    end
    for i = 1:length(fs)
        line(fs{i}, spec{i}, 'Color', 'k');
        if maxfs < max(fs{i})
            maxfs = max(fs{i});
        end
        if minfs > min(fs{i})
            minfs = min(fs{i});
        end
    end
    set(handles.axis, 'XLim', [minfs, maxfs], 'YLim', [0.5, 4], ...
        'XScale', 'log', 'Box', 'on', 'FontSize', 14, 'FontName', 'Georgia');
    handles.userdata.curfigtype = type;
    handles.bandtype.String = bandtype;
    guidata(hObject, handles);
end
end

function Command_Callback(hObject, eventdata, handles)
% hObject    handle to Command (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of Command as text
%        str2double(get(hObject,'String')) returns contents of Command as a double
end

% --- Executes during object creation, after setting all properties.
function Command_CreateFcn(hObject, eventdata, handles)
% hObject    handle to Command (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in pushbuttonplotfigure.
function pushbuttonplotfigure_Callback(hObject, eventdata, handles)
% hObject    handle to pushbuttonplotfigure (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
t = handles.userdata.willplot;
if strcmp(t, 'spec')
    plotspec(hObject, eventdata, handles);
elseif strcmp(t, 'station')
    plotstationscatter(hObject, eventdata, handles);
elseif strcmp(t, 'wave')
    plotwave(hObject, eventdata, handles);
end
end

% --- Executes on button press in buttonplotstation.
function buttonplotstation_Callback(hObject, eventdata, handles)
% hObject    handle to buttonplotstation (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of buttonplotstation
handles.userdata.willplot = 'station';
guidata(hObject, handles);
end

% --- Executes on button press in buttonplotwave.
function buttonplotwave_Callback(hObject, eventdata, handles)
% hObject    handle to buttonplotwave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of buttonplotwave
handles.userdata.willplot = 'wave';
guidata(hObject, handles);
end

% --- Executes on button press in buttonplotspectrum.
function buttonplotspectrum_Callback(hObject, eventdata, handles)
% hObject    handle to buttonplotspectrum (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of buttonplotspectrum
handles.userdata.willplot = 'spec';
guidata(hObject, handles);
end

% --- Executes on selection change in stationlistbox.
function stationlistbox_Callback(hObject, eventdata, handles)
% hObject    handle to stationlistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns stationlistbox contents as cell array
%        contents{get(hObject,'Value')} returns selected item from stationlistbox
station = handles.stationlistbox.String{handles.stationlistbox.Value};
runCommand(['csta ', station], hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function stationlistbox_CreateFcn(hObject, eventdata, handles)
% hObject    handle to stationlistbox (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --------------------------------------------------------------------
function initbutton_ClickedCallback(hObject, eventdata, handles)
% hObject    handle to initbutton (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if isfield(handles, 'userdata')
    if isfield(handles.userdata, 'cli')
        msgbox('Connection has already been established.');
        return;
    end
end
for i = 1:10
    if isfield(handles, 'userdata')
        if isfield(handles.userdata, 'cli')
            break;
        end
    end
    try
        handles.userdata.cli = tcpclient(handles.ServerIP.String, str2double(handles.serverPort.String), 'Timeout', 60);
    catch
        pause(1);
    end
end
if ~isfield(handles, 'userdata') || ~isfield(handles.userdata, 'cli')
    errordlg('Cannot connect to server');
    return
end
handles.userdata.willplot = 'station';
runCommand('plot station', hObject, eventdata, handles);
if ~strcmp(handles.ServerIP.String, '127.0.0.1') && ~strcmp(handles.ServerIP.String, 'localhost')
    gettmpfile(hObject, eventdata, handles);
    handles.userdata.figpath = '.tmpplot.mat';
else
    handles.userdata.figpath = runCommand('tmpfile', hObject, eventdata, handles);
end
load(handles.userdata.figpath, 'tag1', 'tag2');
if exist('tag1', 'var') && exist('tag2', 'var')
    handles.stationlistbox.String = sort([tag1; tag2]);
end
guidata(hObject, handles);
end

% --- Executes on selection change in pickphasetype.
function pickphasetype_Callback(hObject, eventdata, handles)
% hObject    handle to pickphasetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pickphasetype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pickphasetype
handles.userdata.currentphasetype = handles.pickphasetype.String{handles.pickphasetype.Value};
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function pickphasetype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pickphasetype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in buttonPickphase.
function buttonPickphase_Callback(hObject, eventdata, handles)
% hObject    handle to buttonPickphase (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axis);
[x, ~] = ginput(1);
% fprintf('pick arrival time at: %fs\n', x);
t = handles.userdata.btime + seconds(x);
currentStation = handles.stationlistbox.String{handles.stationlistbox.Value};
currentPhase = handles.pickphasetype.String{handles.pickphasetype.Value};
runCommand(['csta ', currentStation], hObject, eventdata, handles);
runCommand(sprintf('set phase %s at %s', currentPhase, datestr(t, 'yyyy-mm-ddTHH:MM:SS.FFF')), ...
    hObject, eventdata, handles);
end

% --- Executes on button press in togglebuttonfilter.
function togglebuttonfilter_Callback(hObject, eventdata, handles)
% hObject    handle to togglebuttonfilter (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of togglebuttonfilter
od = str2double(handles.edit_filterOrder.String);
b1 = str2double(handles.edit_filterband1.String);
b2 = str2double(handles.edit_filterband2.String);
if handles.togglebuttonfilter.Value
    runCommand('set status filtwave true', hObject, eventdata, handles);
    runCommand(sprintf('set status filterorder %d', od), hObject, eventdata, handles);
    runCommand(sprintf('set status filterband %f %f', b1, b2), hObject, eventdata, handles);
    handles.togglebuttonfilter.String = 'ON';
    handles.togglebuttonfilter.ForegroundColor = [0, 0, 0];
else
    runCommand('set status filtwave false', hObject, eventdata, handles);
    runCommand(sprintf('set status filterorder %d', od), hObject, eventdata, handles);
    runCommand(sprintf('set status filterband %f %f', b1, b2), hObject, eventdata, handles);
    handles.togglebuttonfilter.String = 'OFF';
    handles.togglebuttonfilter.ForegroundColor = [0.8, 0.8, 0.8];
end
guidata(hObject, handles);
end

function edit_filterOrder_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filterOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_filterOrder as text
%        str2double(get(hObject,'String')) returns contents of edit_filterOrder as a double
t = str2double(get(hObject, 'String'));
t = round(t);
if t > 0
    runCommand(sprintf('set status filterorder %d', t), hObject, eventdata, handles);
end
end

% --- Executes during object creation, after setting all properties.
function edit_filterOrder_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filterOrder (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit_filterband1_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filterband1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_filterband1 as text
%        str2double(get(hObject,'String')) returns contents of edit_filterband1 as a double
b1 = str2double(handles.edit_filterband1.String);
b2 = str2double(handles.edit_filterband2.String);
if b1 > 0 && b2 > 0
    runCommand(sprintf('set status filterband %f %f', b1, b2), hObject, eventdata, handles);
end
end

% --- Executes during object creation, after setting all properties.
function edit_filterband1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filterband1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

function edit_filterband2_Callback(hObject, eventdata, handles)
% hObject    handle to edit_filterband2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_filterband2 as text
%        str2double(get(hObject,'String')) returns contents of edit_filterband2 as a double
b1 = str2double(handles.edit_filterband1.String);
b2 = str2double(handles.edit_filterband2.String);
if b1 > 0 && b2 > 0
    runCommand(sprintf('set status filterband %f %f', b1, b2), hObject, eventdata, handles);
end
end

% --- Executes during object creation, after setting all properties.
function edit_filterband2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_filterband2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end

% --- Executes on button press in pushbutton_pickwindow.
function pushbutton_pickwindow_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_pickwindow (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
axes(handles.axis);
[x, ~] = ginput(2);
phaseAT = handles.userdata.attime(handles.pickphasetype.Value);
phaseTT = handles.userdata.tttime(handles.pickphasetype.Value);
currentStation = handles.stationlistbox.String{handles.stationlistbox.Value};
currentPhase = handles.pickphasetype.String{handles.pickphasetype.Value};
currentTrim = handles.listbox_trimtype.String{handles.listbox_trimtype.Value};
if strcmp(currentTrim, 'polarity_trim')
    t = x - phaseTT;
else
    t = x - phaseAT;
end
runCommand(['csta ', currentStation], hObject, eventdata, handles);
runCommand(sprintf('set phase %s %s %f %f', currentPhase, currentTrim, t(1), t(2)), ...
    hObject, eventdata, handles);
end

% --- Executes on selection change in listbox_trimtype.
function listbox_trimtype_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_trimtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_trimtype contents as cell array

%        contents{get(hObject,'Value')} returns selected item from listbox_trimtype
handles.userdata.currentTrimType = handles.listbox_trimtype.String{handles.listbox_trimtype.Value};
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function listbox_trimtype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_trimtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in pushbutton_picksyntt.
function pushbutton_picksyntt_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_picksyntt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)axes(handles.axis);
axes(handles.axis);
[x, ~] = ginput(1);
t = x - handles.userdata.gbtime;
currentStation = handles.stationlistbox.String{handles.stationlistbox.Value};
currentPhase = handles.pickphasetype.String{handles.pickphasetype.Value};
runCommand(['csta ', currentStation], hObject, eventdata, handles);
runCommand(sprintf('set phase %s tt %f', currentPhase, t), hObject, eventdata, handles);
end



function xlim1_Callback(hObject, eventdata, handles)
% hObject    handle to xlim1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xlim1 as text
%        str2double(get(hObject,'String')) returns contents of xlim1 as a double
v1 = str2double(handles.xlim1.String);
v2 = str2double(handles.xlim2.String);
set(handles.axis, 'XLim', [v1, v2]);
end

% --- Executes during object creation, after setting all properties.
function xlim1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xlim1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function xlim2_Callback(hObject, eventdata, handles)
% hObject    handle to xlim2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of xlim2 as text
%        str2double(get(hObject,'String')) returns contents of xlim2 as a double
v1 = str2double(handles.xlim1.String);
v2 = str2double(handles.xlim2.String);
set(handles.axis, 'XLim', [v1, v2]);
end

% --- Executes during object creation, after setting all properties.
function xlim2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to xlim2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


function waveamp_Callback(hObject, eventdata, handles)
% hObject    handle to waveamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of waveamp as text
%        str2double(get(hObject,'String')) returns contents of waveamp as a double
v = str2double(get(hObject,'String'));
runCommand(sprintf('set status waveamp %f', v), hObject, eventdata, handles);
end

% --- Executes during object creation, after setting all properties.
function waveamp_CreateFcn(hObject, eventdata, handles)
% hObject    handle to waveamp (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in pushbutton6.
function pushbutton6_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton6 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
currentStation = handles.stationlistbox.String{handles.stationlistbox.Value};
runCommand(['set status select ' currentStation], hObject, eventdata, handles);
end


% --- Executes on selection change in bandtype.
function bandtype_Callback(hObject, eventdata, handles)
% hObject    handle to bandtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns bandtype contents as cell array
%        contents{get(hObject,'Value')} returns selected item from bandtype
end

% --- Executes during object creation, after setting all properties.
function bandtype_CreateFcn(hObject, eventdata, handles)
% hObject    handle to bandtype (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in pushbutton_setband.
function pushbutton_setband_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_setband (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if handles.bandtype.Value > 0
    bd = handles.bandtype.String{handles.bandtype.Value};
    b1 = str2double(handles.edit_filterband1.String);
    b2 = str2double(handles.edit_filterband2.String);
    runCommand(sprintf('set phase %s %f %f', bd, b1, b2), hObject, eventdata, handles);
end
end


% --- Executes when user attempts to close figure1.
function figure1_CloseRequestFcn(hObject, eventdata, handles)
% hObject    handle to figure1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: delete(hObject) closes the figure
if isfield(handles, 'userdata')
    try
        runCommand('save .data.bak.jld2', hObject, eventdata, handles);
        runCommand('exit', hObject, eventdata, handles);
    catch
        warning('Cannot connect to server. Quit directly');
    end
end
delete(hObject);
end

% --- Executes on selection change in listbox_maxlag_type.
function listbox_maxlag_type_Callback(hObject, eventdata, handles)
% hObject    handle to listbox_maxlag_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox_maxlag_type contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox_maxlag_type
end

% --- Executes during object creation, after setting all properties.
function listbox_maxlag_type_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox_maxlag_type (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in pushbutton_pickmaxlag.
function pushbutton_pickmaxlag_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_pickmaxlag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
lagtype = handles.listbox_maxlag_type.String{handles.listbox_maxlag_type.Value};
currentStation = handles.stationlistbox.String{handles.stationlistbox.Value};
currentPhase = handles.pickphasetype.String{handles.pickphasetype.Value};
phaseAT = handles.userdata.attime(handles.pickphasetype.Value);
axes(handles.axis);
[t, ~] = ginput(1);
lag = abs(t-phaseAT);
runCommand(['csta ' currentStation], hObject, eventdata, handles);
runCommand(sprintf('set phase %s %s_maxlag %f', currentPhase, lagtype, lag), ...
    hObject, eventdata, handles);
end


% --- Executes on selection change in popupmenu_commands.
function popupmenu_commands_Callback(hObject, eventdata, handles)
% hObject    handle to popupmenu_commands (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns popupmenu_commands contents as cell array
%        contents{get(hObject,'Value')} returns selected item from popupmenu_commands
handles.Command.String = handles.popupmenu_commands.String{handles.popupmenu_commands.Value};
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function popupmenu_commands_CreateFcn(hObject, eventdata, handles)
% hObject    handle to popupmenu_commands (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in pushbutton_clearcmd.
function pushbutton_clearcmd_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_clearcmd (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if strcmp(handles.Command.String, 'Input Command')
    return;
end
if isempty(handles.historylist.String)
    handles.historylist.String = {handles.Command.String};
else
    idxs = [];
    for i = 1:length(handles.historylist.String)
        if ~strcmp(handles.Command.String, handles.historylist.String{i})
            idxs = [idxs, i];
        end
    end
    if length(handles.historylist.String) > 10
        idxs = idxs(1:end-1);
    end
    handles.historylist.String = [handles.Command.String; handles.historylist.String(idxs)];
end
handles.historylist.Value = 1;
runCommand(handles.Command.String, hObject, eventdata, handles);
handles.Command.String = 'Input Command';
guidata(hObject, handles);
end


% --- Executes on button press in pushbutton_cmpE.
function pushbutton_cmpE_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cmpE (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
strs = [handles.Command.String ' -c E'];
handles.Command.String = strs;
guidata(hObject, handles);
end


% --- Executes on button press in pushbutton_cmpN.
function pushbutton_cmpN_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cmpN (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
strs = [handles.Command.String ' -c N'];
handles.Command.String = strs;
guidata(hObject, handles);
end


% --- Executes on button press in pushbutto_cmpZ.
function pushbutto_cmpZ_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutto_cmpZ (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
strs = [handles.Command.String ' -c Z'];
handles.Command.String = strs;
guidata(hObject, handles);
end


% --- Executes on button press in pushbutton_save.
function pushbutton_save_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_save (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
cmd = ['save ' handles.edit_savefilename.String '.jld2'];
runCommand(cmd, hObject, eventdata, handles);
end


function edit_savefilename_Callback(hObject, eventdata, handles)
% hObject    handle to edit_savefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of edit_savefilename as text
%        str2double(get(hObject,'String')) returns contents of edit_savefilename as a double
end

% --- Executes during object creation, after setting all properties.
function edit_savefilename_CreateFcn(hObject, eventdata, handles)
% hObject    handle to edit_savefilename (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end


% --- Executes on button press in pushbutton_cutwave.
function pushbutton_cutwave_Callback(hObject, eventdata, handles)
% hObject    handle to pushbutton_cutwave (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
if ~isfield(handles, 'userdata') || ~isfield(handles.userdata, 'btime')
    return;
end
axes(handles.axis);
[x, ~] = ginput(2);
x = sort(x);
t = handles.userdata.btime + seconds(x);
currentStation = handles.stationlistbox.String{handles.stationlistbox.Value};
runCommand(['csta ', currentStation], hObject, eventdata, handles);
runCommand(sprintf('set station base_trim %s %s', ...
    datestr(t(1), 'yyyy-mm-ddTHH:MM:SS.FFF'), ...
    datestr(t(2), 'yyyy-mm-ddTHH:MM:SS.FFF')), hObject, eventdata, handles);
end


% --- Executes on selection change in historylist.
function historylist_Callback(hObject, eventdata, handles)
% hObject    handle to historylist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns historylist contents as cell array
%        contents{get(hObject,'Value')} returns selected item from historylist
handles.Command.String = handles.historylist.String{handles.historylist.Value};
guidata(hObject, handles);
end

% --- Executes during object creation, after setting all properties.
function historylist_CreateFcn(hObject, eventdata, handles)
% hObject    handle to historylist (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function ServerIP_Callback(hObject, eventdata, handles)
% hObject    handle to ServerIP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ServerIP as text
%        str2double(get(hObject,'String')) returns contents of ServerIP as a double
end

% --- Executes during object creation, after setting all properties.
function ServerIP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ServerIP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end



function serverPort_Callback(hObject, eventdata, handles)
% hObject    handle to serverPort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of serverPort as text
%        str2double(get(hObject,'String')) returns contents of serverPort as a double
end

% --- Executes during object creation, after setting all properties.
function serverPort_CreateFcn(hObject, eventdata, handles)
% hObject    handle to serverPort (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
end
