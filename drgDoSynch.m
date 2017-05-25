function drgDoSynch(handles)

switch handles.analysisSync
    case 1
        drgPlotCrossCorr(handles)
    case 2
        drgPlotCorrPSTHEv1vs2Fct(handles)
end