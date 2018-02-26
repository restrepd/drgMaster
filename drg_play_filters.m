            
handles.draq_p.ActualRate=24000;

Fstop1 = 800;
Fpass1 = 1000;
Fpass2 = 2800;
Fstop2 = 3000;
Astop1 = 100;
Apass  = 0.5;
Astop2 = 100;
Fs = handles.draq_p.ActualRate;

d = designfilt('bandpassfir', ...
    'StopbandFrequency1',Fstop1,'PassbandFrequency1', Fpass1, ...
    'PassbandFrequency2',Fpass2,'StopbandFrequency2', Fstop2, ...
    'StopbandAttenuation1',Astop1,'PassbandRipple', Apass, ...
    'StopbandAttenuation2',Astop2, ...
    'DesignMethod','equiripple','SampleRate',Fs);

fvtool(d)

d = designfilt('bandpassiir','FilterOrder',20, ...
    'HalfPowerFrequency1',1000,'HalfPowerFrequency2',3000, ...
    'SampleRate',handles.draq_p.ActualRate);
fvtool(d)

