function drgDoBehavioralAnal(handles)

switch handles.analysisNoBeh
    case 1
        drgBehaviorbyBlock(handles)
    case 2
        drgBehaviorbyTrial(handles)
    case 3
        drgProbabilisticTimecourse(handles)
    case 4 
        drgPercentLick(handles)
    case 5
        drgLickPerTrialTimecourse(handles);
    case 6
        drgITI(handles)
    case 7
        drgLickPerConc(handles)
    case 8
        drgLickTimecourse(handles);
    case 9
        drgLickTraces(handles)
end