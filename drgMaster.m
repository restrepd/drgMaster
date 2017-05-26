function varargout = drgMaster(varargin)
% DRGMASTER MATLAB code for drgMaster.fig
%      DRGMASTER, by itself, creates a new DRGMASTER or raises the existing
%      singleton*.
%
%      H = DRGMASTER returns the handle to a new DRGMASTER or the handle to
%      the existing singleton*.
%
%      DRGMASTER('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DRGMASTER.M with the given input arguments.
%
%      DRGMASTER('Property','Value',...) creates a new DRGMASTER or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before drgMaster_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to drgMaster_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help drgMaster

% Last Modified by GUIDE v2.5 30-Apr-2017 10:59:21

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @drgMaster_OpeningFcn, ...
                   'gui_OutputFcn',  @drgMaster_OutputFcn, ...
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


% --- Executes just before drgMaster is made visible.
function drgMaster_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to drgMaster (see VARARGIN)

% Choose default command line output for drgMaster
handles.output = hObject;
handles.analysisNoOsc=1;
handles.analysisNoSpikes=1;
handles.peakLowF=6;
handles.peakHighF=10;
handles.burstLowF=65;
handles.burstHighF=95;
handles.peakLFPNo=1; %3 for Anan, 2 for Justin
handles.burstLFPNo=1;
handles.evTypeNo=1;
handles.evTypeNo2=1;
handles.unitNo=1;
handles.unitNo2=1;
handles.time_start=-0.2;
handles.time_end=2.2; %2.2 for Shane
handles.time_pad=0.2;
handles.startRef=-2.2;
handles.endRef=0.2;
handles.n_phase_bins=50;
handles.sessionNo=1;
handles.trialNo=1;
handles.lastTrialNo=40;
handles.firstEvNo=1;
handles.lastEvNo=40;
handles.data_vs_simulate=0;
handles.subtractRef=0;
handles.window=1; %This is the FFT window in sec  %Tort is 1 sec, old DR 0.37
handles.noverlap=handles.window*0.9; 
handles.deltaLowF_PAC=1;  %2
handles.deltaHighF_PAC=5; %30 
handles.bandwidth_lowF=3;
handles.bandwidth_highF=25;
handles.which_method=1;
handles.amplitudeHz=80;
handles.phaseHz=10;
handles.unitNo=1;
handles.analysisNoBeh=1;
handles.notch60=0;
handles.displayData=1;
handles.save_drgb=0;
handles.save_events=0;
handles.autoscale=1;
handles.maxLogP=0.2;
handles.minLogP=-0.3;
handles.perTetrode=0;
handles.max_dt_between_events=3.5; %Maximum time between events (odorOn to reinforcement)
handles.read_entire_file=0;
handles.corr_window=0.0004;
handles.analysisSync=1;


% Update handles structure
guidata(hObject, handles);

%Populate string for drgMaster.m modification date
if strcmp('diegorereposmbp.ucdenver.pvt',getComputerName)
    p = mfilename('fullpath');
    fileDirectory = fileparts(p);
    moddatedir = dir(fileDirectory);
    moddatetag = moddatedir.date;
    set(handles.moddate,'String',['ver ',moddatetag]);
end

% UIWAIT makes drgMaster wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = drgMaster_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on selection change in whichOscillatoryAnalysis.
function whichOscillatoryAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to whichOscillatoryAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns whichOscillatoryAnalysis contents as cell array
%        contents{get(hObject,'Value')} returns selected item from whichOscillatoryAnalysis
handles.analysisNoOsc=get(handles.whichOscillatoryAnalysis,'Value');
drgDoOscillatoryAnal(handles)

% --- Executes during object creation, after setting all properties.
function whichOscillatoryAnalysis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichOscillatoryAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in whichEvent.
function whichEvent_Callback(hObject, eventdata, handles)
% hObject    handle to whichEvent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns whichEvent contents as cell array
%        contents{get(hObject,'Value')} returns selected item from whichEvent
handles.evTypeNo=get(handles.whichEvent,'Value');
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichEvent_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichEvent (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in openDrg.
function openDrg_Callback(hObject, eventdata, handles)
% hObject    handle to openDrg (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

[FileName,PathName] = uigetfile('jt_times*.*','Select . drg file to open');
handles.fullName=[PathName,FileName];
handles.FileName=FileName;
handles.PathName=PathName;
set(handles.whichFile,'String',FileName);
load(handles.fullName);
handles.drg=drg;
set(handles.whichEvent,'String',drg.draq_d.eventlabels)
set(handles.whichEvent2,'String',drg.draq_d.eventlabels)

if handles.read_entire_file==1
    handles=drgReadAllDraOrDg(handles);
end


%Set the last trial to the last trial in the session
handles.lastTrialNo=handles.drg.session(handles.sessionNo).noTrials; 
set(handles.whichLastTrial,'String',num2str(handles.lastTrialNo))

for nUn=1:handles.drg.session(handles.sessionNo).noUnits
   String{nUn}=handles.drg.unit(nUn).ch_un; 
end
set(handles.whichUnit,'String',String);
set(handles.whichUnitNo2,'String',String);
% Update the handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function whichFile_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called




function whichTrial_Callback(hObject, eventdata, handles)
% hObject    handle to whichTrial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of whichTrial as text
%        str2double(get(hObject,'String')) returns contents of whichTrial as a double
trialNo=str2double(get(hObject,'String'));
if trialNo>handles.drg.session(handles.sessionNo).noTrials
   trialNo=handles.drg.session(handles.sessionNo).noTrials; 
end
if trialNo<1
   trialNo=1; 
end

%Find the closest event No
sessionNo=handles.sessionNo;
[mint,eventNo]=min(abs(handles.drg.session(sessionNo).events(handles.evTypeNo).times-handles.drg.session(sessionNo).trial_start(trialNo)));
handles.firstEvNo=eventNo(1);
set(handles.firstEventNo,'String',num2str(handles.firstEvNo));

%Now set the trial number
trialNo=find((handles.drg.session(sessionNo).events(handles.evTypeNo).times(handles.firstEvNo)>=handles.drg.session(sessionNo).trial_start)&...
    (handles.drg.session(sessionNo).events(handles.evTypeNo).times(handles.firstEvNo)<=handles.drg.session(sessionNo).trial_start+handles.drg.session(sessionNo).sec_per_trial));

handles.trialNo=trialNo;
set(handles.whichTrial,'String',num2str(handles.trialNo));

% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichTrial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichTrial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


function whichLowF1_Callback(hObject, eventdata, handles)
% hObject    handle to whichLowF1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of whichLowF1 as text
%        str2double(get(hObject,'String')) returns contents of whichLowF1 as a double
handles.peakLowF=str2double(get(hObject,'String'));
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichLowF1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichLowF1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function whichLowF2_Callback(hObject, eventdata, handles)
% hObject    handle to whichLowF2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of whichLowF2 as text
%        str2double(get(hObject,'String')) returns contents of whichLowF2 as a double
handles.peakHighF=str2double(get(hObject,'String'));
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichLowF2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichLowF2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function whichHighF1_Callback(hObject, eventdata, handles)
% hObject    handle to whichHighF1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of whichHighF1 as text
%        str2double(get(hObject,'String')) returns contents of whichHighF1 as a double
handles.burstLowF=str2double(get(hObject,'String'));
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichHighF1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichHighF1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function whichHighF2_Callback(hObject, eventdata, handles)
% hObject    handle to whichHighF2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of whichHighF2 as text
%        str2double(get(hObject,'String')) returns contents of whichHighF2 as a double
handles.burstHighF=str2double(get(hObject,'String'));
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichHighF2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichHighF2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function whichStartTime_Callback(hObject, eventdata, handles)
% hObject    handle to whichStartTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of whichStartTime as text
%        str2double(get(hObject,'String')) returns contents of whichStartTime as a double
handles.time_start=str2double(get(hObject,'String'))-handles.time_pad;
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichStartTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichStartTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function whichEndTime_Callback(hObject, eventdata, handles)
% hObject    handle to whichEndTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of whichEndTime as text
%        str2double(get(hObject,'String')) returns contents of whichEndTime as a double
handles.time_end=str2double(get(hObject,'String'))+handles.time_pad;
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichEndTime_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichEndTime (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function whichPeakLFPNo_Callback(hObject, eventdata, handles)
% hObject    handle to whichPeakLFPNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of whichPeakLFPNo as text
%        str2double(get(hObject,'String')) returns contents of whichPeakLFPNo as a double
handles.peakLFPNo=str2double(get(hObject,'String'));
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichPeakLFPNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichPeakLFPNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function whichBurstLFPNo_Callback(hObject, eventdata, handles)
% hObject    handle to whichBurstLFPNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of whichBurstLFPNo as text
%        str2double(get(hObject,'String')) returns contents of whichBurstLFPNo as a double
handles.burstLFPNo=str2double(get(hObject,'String'));
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichBurstLFPNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichBurstLFPNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function whichStartRef_Callback(hObject, eventdata, handles)
% hObject    handle to whichStartRef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of whichStartRef as text
%        str2double(get(hObject,'String')) returns contents of whichStartRef as a double
handles.startRef=str2double(get(hObject,'String'))-handles.time_pad;
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichStartRef_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichStartRef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function whichEndRef_Callback(hObject, eventdata, handles)
% hObject    handle to whichEndRef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of whichEndRef as text
%        str2double(get(hObject,'String')) returns contents of whichEndRef as a double
handles.endRef=str2double(get(hObject,'String'))+handles.time_pad;
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichEndRef_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichEndRef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in subtractRef.
function subtractRef_Callback(hObject, eventdata, handles)
% hObject    handle to subtractRef (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of subtractRef

handles.subtractRef=get(hObject,'Value');
% Update the handles structure
guidata(hObject, handles);


function whichLastTrial_Callback(hObject, eventdata, handles)
% hObject    handle to whichLastTrial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of whichLastTrial as text
%        str2double(get(hObject,'String')) returns contents of whichLastTrial as a double
handles.lastTrialNo=str2double(get(hObject,'String'));
if handles.lastTrialNo>handles.drg.session(handles.sessionNo).noTrials
   handles.lastTrialNo=handles.drg.session(handles.sessionNo).noTrials; 
end
if handles.lastTrialNo<handles.trialNo
   handles.lastTrialNo=handles.trialNo; 
end


%Now set the event number
sessionNo=handles.sessionNo;

trNo=handles.lastTrialNo+1;
evNo=-1
while (evNo==-1)&(trNo>handles.trialNo)
    trNo=trNo-1;
    evNo = drgFindEvNo(handles,trNo,sessionNo);
end
handles.lastTrialNo=trNo;
set(hObject,'String',num2str(handles.lastTrialNo))
handles.lastEvNo=evNo;
set(handles.lastEventNo,'String',num2str(evNo));

% Update the handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function whichLastTrial_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichLastTrial (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in originalOrSimulate.
function originalOrSimulate_Callback(hObject, eventdata, handles)
% hObject    handle to originalOrSimulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns originalOrSimulate contents as cell array
%        contents{get(hObject,'Value')} returns selected item from originalOrSimulate
handles.data_vs_simulate=get(hObject,'Value')-1;
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function originalOrSimulate_CreateFcn(hObject, eventdata, handles)
% hObject    handle to originalOrSimulate (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in whichMethod.
function whichMethod_Callback(hObject, eventdata, handles)
% hObject    handle to whichMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns whichMethod contents as cell array
%        contents{get(hObject,'Value')} returns selected item from whichMethod
handles.which_method=get(hObject,'Value');
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichMethod_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichMethod (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function whichLowBand_Callback(hObject, eventdata, handles)
% hObject    handle to whichLowBand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of whichLowBand as text
%        str2double(get(hObject,'String')) returns contents of whichLowBand as a double
handles.bandwidth_lowF=str2double(get(hObject,'String'));
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichLowBand_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichLowBand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function whichHighBand_Callback(hObject, eventdata, handles)
% hObject    handle to whichHighBand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of whichHighBand as text
%        str2double(get(hObject,'String')) returns contents of whichHighBand as a double
handles.bandwidth_highF=str2double(get(hObject,'String'));
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichHighBand_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichHighBand (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function whichPhaseHz_Callback(hObject, eventdata, handles)
% hObject    handle to whichPhaseHz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of whichPhaseHz as text
%        str2double(get(hObject,'String')) returns contents of whichPhaseHz as a double
handles.phaseHz=str2double(get(hObject,'String'));
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichPhaseHz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichPhaseHz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function whichAmplitudeHz_Callback(hObject, eventdata, handles)
% hObject    handle to whichAmplitudeHz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of whichAmplitudeHz as text
%        str2double(get(hObject,'String')) returns contents of whichAmplitudeHz as a double
handles.amplitudeHz=str2double(get(hObject,'String'));
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichAmplitudeHz_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichAmplitudeHz (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function whichUnit_Callback(hObject, eventdata, handles)
% hObject    handle to whichUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns originalOrSimulate contents as cell array
%        contents{get(hObject,'Value')} returns selected item from originalOrSimulate
handles.unitNo=get(hObject,'Value');

if handles.unitNo<1
    handles.unitNo=1;
end
if handles.unitNo>handles.drg.session(handles.sessionNo).noUnits
    handles.unitNo=handles.drg.session(handles.sessionNo).noUnits;
end
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichUnit_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichUnit (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in whichSpikeAnalysis.
function whichSpikeAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to whichSpikeAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns whichSpikeAnalysis contents as cell array
%        contents{get(hObject,'Value')} returns selected item from whichSpikeAnalysis
handles.analysisNoSpikes=get(handles.whichSpikeAnalysis,'Value');
drgDoSpikePlotsMaster(handles)

% --- Executes during object creation, after setting all properties.
function whichSpikeAnalysis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichSpikeAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in whichBehaviorAnalysis.
function whichBehaviorAnalysis_Callback(hObject, eventdata, handles)
% hObject    handle to whichBehaviorAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns whichBehaviorAnalysis contents as cell array
%        contents{get(hObject,'Value')} returns selected item from whichBehaviorAnalysis
handles.analysisNoBeh=get(handles.whichBehaviorAnalysis,'Value');
drgDoBehavioralAnal(handles)

% --- Executes during object creation, after setting all properties.
function whichBehaviorAnalysis_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichBehaviorAnalysis (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in whichSpikeLFP.
function whichSpikeLFP_Callback(hObject, eventdata, handles)
% hObject    handle to whichSpikeLFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns whichSpikeLFP contents as cell array
%        contents{get(hObject,'Value')} returns selected item from whichSpikeLFP
handles.analysisNoOsc=get(handles.whichSpikeLFP,'Value');
drgDoSpikeLFPAnal(handles)

% --- Executes during object creation, after setting all properties.
function whichSpikeLFP_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichSpikeLFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in openJt.
function openJt_Callback(hObject, eventdata, handles)
% hObject    handle to openJt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


[jtFileName,jtPathName] = uigetfile('*.*','Select jt_times file to open');
handles.jtfullName=[jtPathName,jtFileName];
handles.jtFileName=jtFileName;
handles.jtPathName=jtPathName;
drgRead_jt_times(jtPathName,jtFileName)
FileName=[jtFileName(10:end-4) '_drg.mat'];
set(handles.whichFile,'String',FileName);
handles.fullName=[jtPathName,FileName];
load(handles.fullName);
handles.drg=drg;
handles.drg.drta_p.fullName=[jtPathName,handles.drg.drta_p.FileName];
set(handles.whichEvent,'String',drg.draq_d.eventlabels)
set(handles.whichEvent2,'String',drg.draq_d.eventlabels)

if handles.read_entire_file==1
    handles=drgReadAllDraOrDg(handles); 
end


%Set the last trial to the last trial in the session
handles.lastTrialNo=handles.drg.session(handles.sessionNo).noTrials; 
set(handles.whichLastTrial,'String',num2str(handles.lastTrialNo))

String{1}='None';
for nUn=1:handles.drg.session(handles.sessionNo).noUnits
    if handles.drg.unit(nUn).SingleUnit==1
        String{nUn}=[handles.drg.unit(nUn).ch_un ' SU']; 
    else
        String{nUn}=[handles.drg.unit(nUn).ch_un ' MU']; 
    end
end
set(handles.whichUnit,'String',String);
set(handles.whichUnitNo2,'String',String);
% Update the handles structure
guidata(hObject, handles);


% --- Executes on button press in is60HzNotchOn.
function is60HzNotchOn_Callback(hObject, eventdata, handles)
% hObject    handle to is60HzNotchOn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of is60HzNotchOn
handles.notch60=get(hObject,'Value');
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function is60HzNotchOn_CreateFcn(hObject, eventdata, handles)
% hObject    handle to is60HzNotchOn (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called


% --- Executes on button press in drgAutoscale.
function drgAutoscale_Callback(hObject, eventdata, handles)
% hObject    handle to drgAutoscale (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of drgAutoscale

handles.autoscale=get(hObject,'Value');
% Update the handles structure
guidata(hObject, handles);

function drgLogPowerMin_Callback(hObject, eventdata, handles)
% hObject    handle to drgLogPowerMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of drgLogPowerMin as text
%        str2double(get(hObject,'String')) returns contents of drgLogPowerMin as a double

handles.minLogP=str2double(get(hObject,'String'));
% Update the handles structure
guidata(hObject, handles);


% --- Executes during object creation, after setting all properties.
function drgLogPowerMin_CreateFcn(hObject, eventdata, handles)
% hObject    handle to drgLogPowerMin (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function drgLogPowerMax_Callback(hObject, eventdata, handles)
% hObject    handle to drgLogPowerMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of drgLogPowerMax as text
%        str2double(get(hObject,'String')) returns contents of drgLogPowerMax as a double

handles.maxLogP=str2double(get(hObject,'String'));
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function drgLogPowerMax_CreateFcn(hObject, eventdata, handles)
% hObject    handle to drgLogPowerMax (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on selection change in whichEvent2.
function whichEvent2_Callback(hObject, eventdata, handles)
% hObject    handle to whichEvent2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns whichEvent2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from whichEvent2
handles.evTypeNo2=get(handles.whichEvent2,'Value');
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichEvent2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichEvent2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function firstEventNo_Callback(hObject, eventdata, handles)
% hObject    handle to firstEventNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of firstEventNo as text
%        str2double(get(hObject,'String')) returns contents of firstEventNo as a double
handles.firstEvNo=str2double(get(hObject,'String'));
sessionNo=handles.sessionNo;
if handles.firstEvNo>handles.drg.session(sessionNo).events(handles.evTypeNo).noTimes
   handles.firstEvNo=handles.drg.session(sessionNo).events(handles.evTypeNo).noTimes; 
end
if handles.firstEvNo<1
   handles.firstEvNo=1; 
end
set(hObject,'String',num2str(handles.firstEvNo))

%Now set the trial number
trialNo=find((handles.drg.session(sessionNo).events(handles.evTypeNo).times(handles.firstEvNo)>=handles.drg.session(sessionNo).trial_start)&...
    (handles.drg.session(sessionNo).events(handles.evTypeNo).times(handles.firstEvNo)<=handles.drg.session(sessionNo).trial_start+handles.drg.session(sessionNo).sec_per_trial));
handles.trialNo=trialNo;
set(handles.whichTrial,'String',num2str(handles.trialNo));

% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function firstEventNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to firstEventNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lastEventNo_Callback(hObject, eventdata, handles)
% hObject    handle to lastEventNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lastEventNo as text
%        str2double(get(hObject,'String')) returns contents of lastEventNo as a double
handles.lastEvNo=str2double(get(hObject,'String'));
sessionNo=handles.sessionNo;
if handles.lastEvNo>handles.drg.session(sessionNo).events(handles.evTypeNo).noTimes
   handles.lastEvNo=handles.drg.session(sessionNo).events(handles.evTypeNo).noTimes; 
end
if handles.lastEvNo<=handles.firstEvNo
   handles.lastEvNo=handles.firstEvNo; 
end
set(hObject,'String',num2str(handles.lastEvNo))

%Now set the event number
trialNo=find((handles.drg.session(sessionNo).events(handles.evTypeNo).times(handles.lastEvNo)>=handles.drg.session(sessionNo).trial_start)&...
    (handles.drg.session(sessionNo).events(handles.evTypeNo).times(handles.lastEvNo)<=handles.drg.session(sessionNo).trial_start+handles.drg.session(sessionNo).sec_per_trial));
handles.lastTrialNo=trialNo;
set(handles.whichLastTrial,'String',num2str(handles.lastTrialNo));
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function lastEventNo_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lastEventNo (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in whichUnitNo2.
function whichUnitNo2_Callback(hObject, eventdata, handles)
% hObject    handle to whichUnitNo2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns whichUnitNo2 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from whichUnitNo2
handles.unitNo2=get(hObject,'Value');

if handles.unitNo2<1
    handles.unitNo2=1;
end
if handles.unitNo2>handles.drg.session(handles.sessionNo).noUnits
    handles.unitNo2=handles.drg.session(handles.sessionNo).noUnits;
end
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichUnitNo2_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichUnitNo2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in drgDoSynch.
function drgDoSynch_Callback(hObject, eventdata, handles)
% hObject    handle to drgDoSynch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns drgDoSynch contents as cell array
%        contents{get(hObject,'Value')} returns selected item from drgDoSynch
handles.analysisSync=get(handles.drgDoSynch,'Value');
drgDoSynch(handles)

% --- Executes during object creation, after setting all properties.
function drgDoSynch_CreateFcn(hObject, eventdata, handles)
% hObject    handle to drgDoSynch (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function corr_window_Callback(hObject, eventdata, handles)
% hObject    handle to corr_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of corr_window as text
%        str2double(get(hObject,'String')) returns contents of corr_window as a double
handles.corr_window=str2double(get(hObject,'String'));
% Update the handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function corr_window_CreateFcn(hObject, eventdata, handles)
% hObject    handle to corr_window (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
