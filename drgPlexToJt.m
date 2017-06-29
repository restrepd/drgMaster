function varargout = drgPlexToJt(varargin)
% DRGPLEXTOJT MATLAB code for drgPlexToJt.fig
%      DRGPLEXTOJT, by itself, creates a new DRGPLEXTOJT or raises the existing
%      singleton*.
%
%      H = DRGPLEXTOJT returns the handle to a new DRGPLEXTOJT or the handle to
%      the existing singleton*.
%
%      DRGPLEXTOJT('CALLBACK',hObject,eventData,handles,...) calls the local
%      function named CALLBACK in DRGPLEXTOJT.M with the given input arguments.
%
%      DRGPLEXTOJT('Property','Value',...) creates a new DRGPLEXTOJT or raises the
%      existing singleton*.  Starting from the left, property value pairs are
%      applied to the GUI before drgPlexToJt_OpeningFcn gets called.  An
%      unrecognized property name or invalid value makes property application
%      stop.  All inputs are passed to drgPlexToJt_OpeningFcn via varargin.
%
%      *See GUI Options on GUIDE's Tools menu.  Choose "GUI allows only one
%      instance to run (singleton)".
%
% See also: GUIDE, GUIDATA, GUIHANDLES

% Edit the above text to modify the response to help drgPlexToJt

% Last Modified by GUIDE v2.5 29-Jun-2017 06:29:26

% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @drgPlexToJt_OpeningFcn, ...
                   'gui_OutputFcn',  @drgPlexToJt_OutputFcn, ...
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


% --- Executes just before drgPlexToJt is made visible.
function drgPlexToJt_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to drgPlexToJt (see VARARGIN)

% Choose default command line output for drgPlexToJt
handles.output = hObject;
handles.which_experiment=1;
handles.draq_p.plx.FileNameLFP='no_file.mat';
handles.threshold=0.5;
% Update handles structure
guidata(hObject, handles);

% UIWAIT makes drgPlexToJt wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = drgPlexToJt_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;


% --- Executes on button press in openSpikes1.
function openSpikes1_Callback(hObject, eventdata, handles)
% hObject    handle to openSpikes1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.draq_p.plx.FileNameSpikes1,handles.draq_p.plx.PathNameSpikes1] = uigetfile({'*.mat'},'Select spikes file to open');
set(handles.spikeFile1Label,'String',handles.draq_p.plx.FileNameSpikes1);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in openEventsFile.
function openEventsFile_Callback(hObject, eventdata, handles)
% hObject    handle to openEventsFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.draq_p.plx.FileNameEvents,handles.draq_p.plx.PathNameEvents] = uigetfile({'*.txt'},'Select events file to open');
set(handles.eventsFIleLabel,'String',handles.draq_p.plx.FileNameEvents);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in odorOnFile.
function odorOnFile_Callback(hObject, eventdata, handles)
% hObject    handle to odorOnFile (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.draq_p.plx.FileNameOdorOn,handles.draq_p.plx.PathNameOdorOn] = uigetfile({'*.mat'},'Select odor on file to open');
set(handles.odorOnLabel,'String',handles.draq_p.plx.FileNameOdorOn);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in generateJt.
function generateJt_Callback(hObject, eventdata, handles)
% hObject    handle to generateJt (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

fprintf(1, 'Reading plexon files...\n');

%Initialize draq_p
draq_p.sec_before_trigger=4;
draq_p.sec_per_trigger=10;
draq_p.pre_gain=0;
draq_p.scaling=1;
draq_p.offset=0;
draq_p.plx=handles.draq_p.plx;
draq_p.no_spike_ch=4;

drta_p.FileName=handles.draq_p.plx.FileNameLFP;
%Initialize draq_d
draq_d.noTrials=0;

%Load the odor on and setup tstart
load([handles.draq_p.plx.PathNameOdorOn handles.draq_p.plx.FileNameOdorOn])

try
    close 1
catch
end

hFig1 = figure(1);
set(hFig1, 'units','normalized','position',[.05 .75 .35 .15])

eval(['draq_p.plx.srate_odorOn=1/' handles.draq_p.plx.FileNameOdorOn(1:end-4) '_ts_step;'])
eval(['odorOn=' handles.draq_p.plx.FileNameOdorOn(1:end-4) ';'])
plot(odorOn)
threshold=handles.threshold*(max(odorOn)-min(odorOn))+min(odorOn);
hold on
plot([1 length(odorOn)],[threshold threshold],'-r')

%Find odor on trials
at_end=0;
ii=1;
while at_end==0
    delta_ii_first=find(odorOn(ii:end)>=threshold,1,'first');
    
    %Found odorOn
    if ~isempty(delta_ii_first)
        
        %Find the interval
        ii=ii+delta_ii_first;
        delta_ii_last=find(odorOn(ii:end)<threshold,1,'first');
        if ~isempty(delta_ii_last)
            ii=ii+delta_ii_last;
            
                %Found a trial
                %Full trial
                draq_d.noTrials=draq_d.noTrials+1;
                %Trial start time goes in column 1
                trials_to_sort(draq_d.noTrials,1)=(ii/draq_p.plx.srate_odorOn)-(draq_p.sec_before_trigger+1);
                trials_to_sort(draq_d.noTrials,2)=1;
            
        else
            at_end=1;
        end
    else
        at_end=1;
    end
end

full_trials=draq_d.noTrials;

fprintf(1, 'Found %d full trials...\n',full_trials);



%Find inter trials
last_trial=draq_d.noTrials;



fprintf(1, 'Finding inter trials...\n');

%Empty trial before the first trial
if trials_to_sort(1,1)>draq_p.sec_per_trigger+draq_p.sec_before_trigger
    draq_d.noTrials=draq_d.noTrials+1;
    trials_to_sort(draq_d.noTrials,1)= trials_to_sort(1,1)-draq_p.sec_per_trigger;
    trials_to_sort(draq_d.noTrials,2)=2;
end

%Empty trials between tirals
for ii=2:last_trial
    if (trials_to_sort(ii,1)-trials_to_sort(ii-1,1)-draq_p.sec_per_trigger)>draq_p.sec_per_trigger
        
        %Trial after the last one
        draq_d.noTrials=draq_d.noTrials+1;
        trials_to_sort(draq_d.noTrials,1)= trials_to_sort(ii-1,1)+draq_p.sec_per_trigger;
        trials_to_sort(draq_d.noTrials,2)=2;
        
        %Trial before the next one
        draq_d.noTrials=draq_d.noTrials+1;
        trials_to_sort(draq_d.noTrials,1)= trials_to_sort(ii,1)-draq_p.sec_per_trigger;
        trials_to_sort(draq_d.noTrials,2)=2;
    end
end

fprintf(1, 'Found %d inter trials ...\n',draq_d.noTrials-full_trials);

fprintf(1, 'Sorting trials...\n');
sorted_trials=sortrows(trials_to_sort);
draq_d.t_trial=sorted_trials(:,1)';
draq_d.odor_inter=sorted_trials(:,2);


draq_d.t_end=draq_d.t_trial(end)+draq_p.sec_per_trigger;
fprintf(1, 'Done reading Plexon files...\n');

%Generate events

%Read the text file
fileID=fopen( [handles.draq_p.plx.PathNameEvents handles.draq_p.plx.FileNameEvents],'r');
for trialNo=1:full_trials
    events{trialNo}=fscanf(fileID,'%s',1);
end
fclose(fileID);

%Which experiment?
draq_d.noEvents=0;
switch handles.which_experiment
    case 1
        %go-no go
        draq_d.nEvPerType=zeros(1,17);
        draq_d.nEventTypes=17;
        draq_d.eventlabels=cell(1,17);
        draq_d.eventlabels{1}='TStart';
        draq_d.eventlabels{2}='OdorOn';
        draq_d.eventlabels{3}='Hit';
        draq_d.eventlabels{4}='HitE';
        draq_d.eventlabels{5}='S+';
        draq_d.eventlabels{6}='S+E';
        draq_d.eventlabels{7}='Miss';
        draq_d.eventlabels{8}='MissE';
        draq_d.eventlabels{9}='CR';
        draq_d.eventlabels{10}='CRE';
        draq_d.eventlabels{11}='S-';
        draq_d.eventlabels{12}='S-E';
        draq_d.eventlabels{13}='FA';
        draq_d.eventlabels{14}='FAE';
        draq_d.eventlabels{15}='Reinf';
        draq_d.eventlabels{16}='Short';
        draq_d.eventlabels{17}='Inter';
        
        drta_p.which_c_program=2;
 
        
        odorNo=0;
        for trialNo=1:draq_d.noTrials
            drta_p.trial_ch_processed(1:16,trialNo)=zeros(16,1);
            drta_p.trial_allch_processed(trialNo)=0;
            
            
            
            
            if draq_d.odor_inter(trialNo)==1
                %This is an odor trial
                %tstart
                draq_d.noEvents=draq_d.noEvents+1;
                draq_d.events(draq_d.noEvents)=draq_d.t_trial(trialNo)+draq_p.sec_before_trigger-1.5;
                draq_d.eventType(draq_d.noEvents)=1;
                draq_d.nEvPerType(1)=draq_d.nEvPerType(1)+1;
                
                %odorOn
                draq_d.noEvents=draq_d.noEvents+1;
                draq_d.events(draq_d.noEvents)=draq_d.t_trial(trialNo)+draq_p.sec_before_trigger;
                draq_d.eventType(draq_d.noEvents)=2;
                draq_d.nEvPerType(2)=draq_d.nEvPerType(2)+1;
                
                odorNo=odorNo+1;
                
                %Hit
                if strcmp(events{odorNo},'Hit')==1
                    %Hit
                    draq_d.noEvents=draq_d.noEvents+1;
                    draq_d.events(draq_d.noEvents)=draq_d.t_trial(trialNo)+draq_p.sec_before_trigger;
                    draq_d.eventType(draq_d.noEvents)=3;
                    draq_d.nEvPerType(3)=draq_d.nEvPerType(3)+1;
                    
                    %S+
                    draq_d.noEvents=draq_d.noEvents+1;
                    draq_d.events(draq_d.noEvents)=draq_d.t_trial(trialNo)+draq_p.sec_before_trigger;
                    draq_d.eventType(draq_d.noEvents)=5;
                    draq_d.nEvPerType(5)=draq_d.nEvPerType(5)+1;
                    
                    %Reinf
                    draq_d.noEvents=draq_d.noEvents+1;
                    draq_d.events(draq_d.noEvents)=draq_d.t_trial(trialNo)+draq_p.sec_before_trigger+2;
                    draq_d.eventType(draq_d.noEvents)=15;
                    draq_d.nEvPerType(15)=draq_d.nEvPerType(15)+1;
                end
                
                %Miss
                if strcmp(events{odorNo},'Miss')==1
                    %Miss
                    draq_d.noEvents=draq_d.noEvents+1;
                    draq_d.events(draq_d.noEvents)=draq_d.t_trial(trialNo)+draq_p.sec_before_trigger;
                    draq_d.eventType(draq_d.noEvents)=7;
                    draq_d.nEvPerType(7)=draq_d.nEvPerType(7)+1;
                    
                    %S+
                    draq_d.noEvents=draq_d.noEvents+1;
                    draq_d.events(draq_d.noEvents)=draq_d.t_trial(trialNo)+draq_p.sec_before_trigger;
                    draq_d.eventType(draq_d.noEvents)=5;
                    draq_d.nEvPerType(5)=draq_d.nEvPerType(5)+1;
                end
                
                %CR
                if strcmp(events{odorNo},'CR')==1
                    %CR
                    draq_d.noEvents=draq_d.noEvents+1;
                    draq_d.events(draq_d.noEvents)=draq_d.t_trial(trialNo)+draq_p.sec_before_trigger;
                    draq_d.eventType(draq_d.noEvents)=9;
                    draq_d.nEvPerType(9)=draq_d.nEvPerType(9)+1;
                    
                    %S-
                    draq_d.noEvents=draq_d.noEvents+1;
                    draq_d.events(draq_d.noEvents)=draq_d.t_trial(trialNo)+draq_p.sec_before_trigger;
                    draq_d.eventType(draq_d.noEvents)=11;
                    draq_d.nEvPerType(11)=draq_d.nEvPerType(11)+1;
                end
                
                %FA
                if strcmp(events{odorNo},'FA')==1
                    %FA
                    draq_d.noEvents=draq_d.noEvents+1;
                    draq_d.events(draq_d.noEvents)=draq_d.t_trial(trialNo)+draq_p.sec_before_trigger;
                    draq_d.eventType(draq_d.noEvents)=13;
                    draq_d.nEvPerType(13)=draq_d.nEvPerType(13)+1;
                    
                    %S-
                    draq_d.noEvents=draq_d.noEvents+1;
                    draq_d.events(draq_d.noEvents)=draq_d.t_trial(trialNo)+draq_p.sec_before_trigger;
                    draq_d.eventType(draq_d.noEvents)=11;
                    draq_d.nEvPerType(11)=draq_d.nEvPerType(11)+1;
                end
                
                
            else
                %This is an inter trial
                draq_d.noEvents=draq_d.noEvents+1;
                draq_d.events(draq_d.noEvents)=draq_d.t_trial(trialNo)+draq_p.sec_before_trigger;
                draq_d.eventType(draq_d.noEvents)=17;
                draq_d.nEvPerType(17)=draq_d.nEvPerType(17)+1;
            end
            
        end
        
        %dropcspm
        indxOdorOn=find(strcmp('OdorOn',draq_d.eventlabels));
        evenTypeIndxOdorOn=find(draq_d.eventType==indxOdorOn);
        
        szev=size(evenTypeIndxOdorOn);
        numBlocks=ceil(szev(2)/20);
        draq_d.blocks=zeros(numBlocks,2);
        draq_d.blocks(1,1)=min(draq_d.events)-0.001;
        draq_d.blocks(numBlocks,2)=max(draq_d.events)+0.001;
        
        for blockNo=2:numBlocks
            draq_d.blocks(blockNo,1)=((draq_d.events(evenTypeIndxOdorOn((blockNo-1)*20))...
                +draq_d.events(evenTypeIndxOdorOn((blockNo-1)*20+1)))/2)-0.001;
        end
        for blockNo=1:numBlocks-1
            draq_d.blocks(blockNo,2)=((draq_d.events(evenTypeIndxOdorOn(blockNo*20))...
                +draq_d.events(evenTypeIndxOdorOn(blockNo*20+1)))/2)+0.001;
        end
    case 2
        %go-go
        
end

for trialNo=1:draq_d.noTrials
    drta_p.trial_allch_processed(trialNo)=1;
    drta_p.trial_ch_processed(1:16,trialNo)=ones(16,1);
end

%Now do the spikes
current_offset=0;
cluster_class_per_file=[];
noSpikes=zeros(1,4);
all_timestamp_per_file=[];
drta_p.ch_processed=zeros(1,4);
for tetNo=1:4
    %Is there a file with data for this tetrode?
    eval(['does_file_exist= isfield(handles.draq_p.plx,''FileNameSpikes' num2str(tetNo) ''');' ])
    offset_for_chan(tetNo)=current_offset;
    if does_file_exist==1
        %Load spikes for this tetrode
        eval([ 'load([handles.draq_p.plx.PathNameSpikes' num2str(tetNo) ' handles.draq_p.plx.FileNameSpikes' num2str(tetNo) '])'])
        switch tetNo
            case 1
                cluster_class_per_file=[cluster_class_per_file TETSPKC01(:,2)'];
                all_timestamp_per_file=[all_timestamp_per_file TETSPKC01(:,3)'];
                noSpikes(tetNo)=length(TETSPKC01(:,3));
                current_offset=current_offset+noSpikes(tetNo);
                units_per_tet(tetNo)=max(TETSPKC01(:,2));
                drta_p.ch_processed(tetNo)=1;
            case 2
                cluster_class_per_file=[cluster_class_per_file TETSPKC05(:,2)'];
                all_timestamp_per_file=[all_timestamp_per_file TETSPKC05(:,3)'];
                noSpikes(tetNo)=length(TETSPKC05(:,3));
                current_offset=current_offset+noSpikes(tetNo);
                units_per_tet(tetNo)=max(TETSPKC05(:,2));
                drta_p.ch_processed(tetNo)=1;
            case 3
                cluster_class_per_file=[cluster_class_per_file TETSPKC09(:,2)'];
                all_timestamp_per_file=[all_timestamp_per_file TETSPKC09(:,3)'];
                noSpikes(tetNo)=length(TETSPKC09(:,3));
                current_offset=current_offset+noSpikes(tetNo);
                units_per_tet(tetNo)=max(TETSPKC09(:,2));
                drta_p.ch_processed(tetNo)=1;
            case 4
                cluster_class_per_file=[cluster_class_per_file TETSPKC13(:,2)'];
                all_timestamp_per_file=[all_timestamp_per_file TETSPKC13(:,3)'];
                noSpikes(tetNo)=length(TETSPKC13(:,3));
                current_offset=current_offset+noSpikes(tetNo);
                units_per_tet(tetNo)=max(TETSPKC13(:,2));
                drta_p.ch_processed(tetNo)=1;
        end
    else
        noSpikes(tetNo)=0;
    end
end

par.doBehavior=0;

jt_times_file=[handles.draq_p.plx.PathNameEvents 'jt_times_' handles.draq_p.plx.FileNameEvents(1:end-3) 'mat'];

save(jt_times_file, 'cluster_class_per_file', 'par', 'offset_for_chan','noSpikes', 'all_timestamp_per_file','drta_p', 'draq_p', 'draq_d','units_per_tet');

msgbox('Saved jt_times file');
pffft=1


% --- Executes on selection change in whichExperiment.
function whichExperiment_Callback(hObject, eventdata, handles)
% hObject    handle to whichExperiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns whichExperiment contents as cell array
%        contents{get(hObject,'Value')} returns selected item from whichExperiment

handles.which_experiment=get(hObject,'Value');
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function whichExperiment_CreateFcn(hObject, eventdata, handles)
% hObject    handle to whichExperiment (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in openSpikes2.
function openSpikes2_Callback(hObject, eventdata, handles)
% hObject    handle to openSpikes2 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.draq_p.plx.FileNameSpikes2,handles.draq_p.plx.PathNameSpikes2] = uigetfile({'*.mat'},'Select spikes file to open');
set(handles.spikeFile2Label,'String',handles.draq_p.plx.FileNameSpikes2);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in openSpikes3.
function openSpikes3_Callback(hObject, eventdata, handles)
% hObject    handle to openSpikes3 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.draq_p.plx.FileNameSpikes3,handles.draq_p.plx.PathNameSpikes3] = uigetfile({'*.mat'},'Select spikes file to open');
set(handles.spikeFile3Label,'String',handles.draq_p.plx.FileNameSpikes3);
% Update handles structure
guidata(hObject, handles);

% --- Executes on button press in openSpikes4.
function openSpikes4_Callback(hObject, eventdata, handles)
% hObject    handle to openSpikes4 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.draq_p.plx.FileNameSpikes4,handles.draq_p.plx.PathNameSpikes4] = uigetfile({'*.mat'},'Select spikes file to open');
set(handles.spikeFile4Label,'String',handles.draq_p.plx.FileNameSpikes4);
% Update handles structure
guidata(hObject, handles);


% --- Executes on button press in openLFP.
function openLFP_Callback(hObject, eventdata, handles)
% hObject    handle to openLFP (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
[handles.draq_p.plx.FileNameLFP,handles.draq_p.plx.PathNameLFP] = uigetfile({'*.mat'},'Select spikes file to open');
set(handles.LFPfile,'String',handles.draq_p.plx.FileNameLFP);
% Update handles structure
guidata(hObject, handles);



function threshold_Callback(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of threshold as text
%        str2double(get(hObject,'String')) returns contents of threshold as a double
handles.threshold=str2double(get(hObject,'String'));
% Update handles structure
guidata(hObject, handles);

% --- Executes during object creation, after setting all properties.
function threshold_CreateFcn(hObject, eventdata, handles)
% hObject    handle to threshold (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end
