function drgDoOscillatoryAnal(handles)

switch handles.analysisNoOsc
    case 1
        drgThetaAmpPhaseTrialRange(handles)
    case 2
        drgMIvsHiLoTrialRange(handles)
    case 3
        drgLFPspect(handles)
    case 4
        drgLFPPerTrial(handles)
    case 5 
        drgLFPandLicksPerTrial(handles)
    case 6
        drgLFPspectTimecourse(handles)
    case 7
        drgPhaseTimecourse(handles)
    case 8
        drgThetaAmpPhaseLick(handles)
    case 9
        drgLFPwaveTimecourse(handles)
    case 10 
        drgComparePhases(handles)
    case 11
        drgEventRelatedAnalysis(handles)
    case 12
        drgLFProc(handles)
end