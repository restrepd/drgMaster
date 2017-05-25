function eventID=drgGetEventID(handles,trialNo)

Hit=3;
Miss=7;
CR=9;
FA=13;
odorOn=2;

eventID='N/A';

if sum(handles.drg.session(handles.sessionNo).events(Hit).times==handles.drg.session(handles.sessionNo).events(odorOn).times(trialNo))>0
    eventID='Hit';
end

if sum(handles.drg.session(handles.sessionNo).events(Miss).times==handles.drg.session(handles.sessionNo).events(odorOn).times(trialNo))>0
    eventID='Miss';
end

if sum(handles.drg.session(handles.sessionNo).events(CR).times==handles.drg.session(handles.sessionNo).events(odorOn).times(trialNo))>0
    eventID='CR';
end

if sum(handles.drg.session(handles.sessionNo).events(FA).times==handles.drg.session(handles.sessionNo).events(odorOn).times(trialNo))>0
    eventID='FA';
end

