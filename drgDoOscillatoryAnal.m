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
        drgLFPspectTimecourse(handles)
    case 6
        drgPhaseTimecourse(handles)
    case 7
        drgThetaAmpPhaseLick(handles)
    case 8
        drgLFPwaveTimecourse(handles)
end