function drgDoOscillatoryAnal(handles)

switch handles.analysisNoOsc
    case 1
        drgThetaAmpPhaseTrialRange(handles);
    case 2
        drgMIvsHiLoTrialRange(handles);
    case 3
        drgLFPspect(handles);
    case 4
        drgLFPPerTrial(handles);
    case 5 
        drgLFPandLicksPerTrial(handles);
    case 6
        drgLFPspectTimecourse(handles);
    case 7
        drgPhaseTimecourse(handles);
    case 8
        drgThetaAmpPhaseLick(handles);
    case 9
        drgLFPwaveTimecourse(handles);
    case 10 
        drgComparePhases(handles);
    case 11
        drgEventRelatedAnalysis(handles);
    case 12
        drgLFProc(handles);
    case 13
        drgEventRelatedROCAnalysis(handles);
    case 14
%         drgEventRelatedWaveletAnalysis(handles);
        drgLFP_ERWA(handles);
    case 15
         drgERWAtimecourse(handles);
    case 16
        drgLFPcohspectTimecourse(handles);
    case 17
         drgLFPxspectTimecourse(handles);
    case 18 
        drgLFPwavePerTrialPower(handles);
    case 19
        drgLFPspectrogramPerTrialPower(handles);
    case 20
        drgLFPspectPerceptron(handles);
    case 21
       drgThetaAmpPhaseTrialRangeConc(handles); 
end