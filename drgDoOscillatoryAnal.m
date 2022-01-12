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
        drgLFPimgcohspectTimecourse(handles);
    case 18
         drgLFPxspectTimecourse(handles);
    case 19 
        drgLFPwavePerTrialPower(handles);
    case 20
        drgLFPspectrogramPerTrialPower(handles);
    case 21
        drgLFPspectPerceptron(handles);
    case 22
       drgThetaAmpPhaseTrialRangeConc(handles); 
    case 23
        drgLFPwaveTimecourseCont(handles);
    case 24
        drgLFPCorrTrialRange(handles);
    case 25
%         drgLFPCorrAllPairs(handles);
        drgLFPCorrTimecourseAllTets(handles);
    case 26
        drgLFPCorrTimecourse(handles);
    case 27
        drgPLVTimecourse(handles);
    case 28
        drgPLVTimecourseAllElecs(handles);
    case 29
        drgClosedLoopPhase(handles);
    case 30
        drgLFPimgcohspectShTimecourse(handles);
end