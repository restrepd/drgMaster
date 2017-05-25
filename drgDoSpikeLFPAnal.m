function drgDoSpikeLFPAnal(handles)

switch handles.analysisNoOsc
    case 1
        drgLFPspikeTrigPerTr(handles)
    case 2
        drgSpikePhase(handles)
end