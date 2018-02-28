<<<<<<< HEAD
function Probability =drgProbLick(spm,animal_licked)
% This function was written by Dr. Nicolas Busquet
%For S+ B=1, for S- B=0
%If the mouse licked in all 0.5 sec segments spm=1, if mouse did not lick in all segments spm=0

%These are test data
%Uncomment if your are testing the function
%Animal licks all the time
% spm= [ 1 1 1 1 0 1 0 1 1 0 1 0 1 1 1 1 1 0 0 1 1 1 0 0 1 0 1 0 1 0 1 1 0 1 0 0 1 ];
% animal_licked= [ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ];
% 
% %Animal licks only for S+
% spm=           [ 1 1 1 1 0 1 0 1 1 0 1 0 1 1 1 1 1 0 0 1 1 1 0 0 1 0 1 0 1 0 1 1 0 1 0 0 1 ];
% animal_licked= [ 1 1 1 1 0 1 0 1 1 0 1 0 1 1 1 1 1 0 0 1 1 1 0 0 1 0 1 0 1 0 1 1 0 1 0 0 1 ];


% For a moving window of length w
w= 19; % Outputs w+1 elements
lena=length(spm);
Probability = ones(1,w);

for i=1:lena-w
    
    %counts the number of available S+ and S- in the window
    Splus = animal_licked(i:i+w);
    AvailSplus = sum(Splus);
    Smin = ~animal_licked(i:i+w);
    AvailSminus = sum(Smin);
    Lks = spm(i:i+w);
    Nlicks = sum(Lks);

    % calculates the basic probability of licking based on the total number
    % of licks in the window
    
    % Output the mask of when spm and animal_licked are together
    HIT = spm(i:i+w)&animal_licked(i:i+w);
    N_Hits = sum(HIT);
    
    %     % Output the mask when spm and animal_licked are both zero
    %     CORRECT_REJECTION = ~spm(i:i+w)&~animal_licked(i:i+w);
    %     N_CorrRej = sum(CORRECT_REJECTION);
    
    % Output the mask when spm is 1, but animal_licked is zero
    FALSE_ALARM = spm(i:i+w)&~animal_licked(i:i+w);
    N_FalseAlarm = sum(FALSE_ALARM);
    
    
    % calculates probability (new method)
    Totalproto =[];
    protobehav=[];
    % calculates the number of ways to allocate Nlicks in 20 trials
    Totalproto(end+1)= nchoosek(20,Nlicks);
    
    Nhhh=N_Hits;
    Nfff=N_FalseAlarm;
    
    for r=1:20
        if Nfff<0 | Nhhh==AvailSplus
            break
        end
        
        protobehav(end+1) =nchoosek(AvailSplus,Nhhh)*nchoosek(AvailSminus,Nfff);
        Nhhh=Nhhh+1;
        Nfff= Nfff - 1;
        
    end
    
    Probability(i+9)=sum(protobehav)/Totalproto;
    if i==1
        Probability(1:9)=Probability(10);
    end
end

Probability(end+1:lena)=Probability(end);

=======
function Probability =drgProbLick(spm,animal_licked)
% This function was written by Dr. Nicolas Busquet
%For S+ B=1, for S- B=0
%If the mouse licked in all 0.5 sec segments spm=1, if mouse did not lick in all segments spm=0

%These are test data
%Uncomment if your are testing the function
%Animal licks all the time
% spm= [ 1 1 1 1 0 1 0 1 1 0 1 0 1 1 1 1 1 0 0 1 1 1 0 0 1 0 1 0 1 0 1 1 0 1 0 0 1 ];
% animal_licked= [ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ];
% 
% %Animal licks only for S+
% spm=           [ 1 1 1 1 0 1 0 1 1 0 1 0 1 1 1 1 1 0 0 1 1 1 0 0 1 0 1 0 1 0 1 1 0 1 0 0 1 ];
% animal_licked= [ 1 1 1 1 0 1 0 1 1 0 1 0 1 1 1 1 1 0 0 1 1 1 0 0 1 0 1 0 1 0 1 1 0 1 0 0 1 ];


% For a moving window of length w
w= 19; % Outputs w+1 elements
lena=length(spm);
Probability = ones(1,w);

for i=1:lena-w
    
    %counts the number of available S+ and S- in the window
    Splus = animal_licked(i:i+w);
    AvailSplus = sum(Splus);
    Smin = ~animal_licked(i:i+w);
    AvailSminus = sum(Smin);
    Lks = spm(i:i+w);
    Nlicks = sum(Lks);

    % calculates the basic probability of licking based on the total number
    % of licks in the window
    
    % Output the mask of when spm and animal_licked are together
    HIT = spm(i:i+w)&animal_licked(i:i+w);
    N_Hits = sum(HIT);
    
    %     % Output the mask when spm and animal_licked are both zero
    %     CORRECT_REJECTION = ~spm(i:i+w)&~animal_licked(i:i+w);
    %     N_CorrRej = sum(CORRECT_REJECTION);
    
    % Output the mask when spm is 1, but animal_licked is zero
    FALSE_ALARM = spm(i:i+w)&~animal_licked(i:i+w);
    N_FalseAlarm = sum(FALSE_ALARM);
    
    
    % calculates probability (new method)
    Totalproto =[];
    protobehav=[];
    % calculates the number of ways to allocate Nlicks in 20 trials
    Totalproto(end+1)= nchoosek(20,Nlicks);
    
    Nhhh=N_Hits;
    Nfff=N_FalseAlarm;
    
    for r=1:20
        if Nfff<0 | Nhhh==AvailSplus
            break
        end
        
        protobehav(end+1) =nchoosek(AvailSplus,Nhhh)*nchoosek(AvailSminus,Nfff);
        Nhhh=Nhhh+1;
        Nfff= Nfff - 1;
        
    end
    
    Probability(i+9)=sum(protobehav)/Totalproto;
    if i==1
        Probability(1:9)=Probability(10);
    end
end

Probability(end+1:lena)=Probability(end);

>>>>>>> 362cddba4db155729b2e4c347f0e0147096a5855
