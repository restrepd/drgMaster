function handles=drgDoSpikeLFPAnal(handles)

switch handles.analysisNoOsc
    case 1
        drgLFPspikeTrigPerTr(handles)
    case 2
        handles=drgSpikePhase(handles);
end