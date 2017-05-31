function drgDoSpikePlotsMaster(handles)

switch handles.analysisNoSpikes
    case 1
        drgPlotPSTHBlockFct(handles.drg, handles.unitNo, handles.evTypeNo)
    case 2
        drgPlotPSTHEncRetrFct(handles)
    case 3
        drgPlotPSTHFct(handles)
    case 4
        drgScatterPlotFct(handles)
    case 5
        drgScatterPlotFctSorted(handles)
    case 6
        drgPlotAutoCorr(handles)
    case 7
        drgdFvsBFREncRetr(handles)
    case 8
        drgLick_vs_PreFR(handles);
    case 9
        drgPlotPSTHEv1vs2Fct(handles)
    case 10
        drgPlotPSTH_in_phase_Ev1vs2Fct(handles)
end