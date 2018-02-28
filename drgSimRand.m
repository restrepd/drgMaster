<<<<<<< HEAD
function [AveSim prctile5]= drgSimRand(spm,animal_licked)

%These are test data
%Uncomment if your are testing the function
%Animal licks all the time
% spm= [ 1 1 1 1 0 1 0 1 1 0 1 0 1 1 1 1 1 0 0 1 1 1 0 0 1 0 1 0 1 0 1 1 0 1 0 0 1 ];
% animal_licked= [ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ];
%
% %Animal licks only for S+
% spm=           [ 1 1 1 1 0 1 0 1 1 0 1 0 1 1 1 1 1 0 0 1 1 1 0 0 1 0 1 0 1 0 1 1 0 1 0 0 1 ];
% animal_licked= [ 1 1 1 1 0 1 0 1 1 0 1 0 1 1 1 1 1 0 0 1 1 1 0 0 1 0 1 0 1 0 1 1 0 1 0 0 1 ];

w= 19;
lena=length(spm);
SimRand = zeros(100,lena);
AveSim =[];
stdSim =[];
infstd =[];
AveP =[];
animal_licked_sim =[];


%% randomizes animal_licked and calculates Probability

for nsim = 1:100
    ProbaSim=ones(1,w);
    animal_licked_sim = permVec(animal_licked);
    
    % For a moving window of length w
    % Outputs w+1 elements
    
    
    for i=1:lena-w
        
        %counts the number of available S+ and S- in the window
        
        Splus = animal_licked_sim(i:i+w);
        AvailSplus = sum(Splus);
        Smin = ~animal_licked_sim(i:i+w);
        AvailSminus = sum(Smin);
        Lks = spm(i:i+w);
        Nlicks = sum(Lks);
        
        % Output the mask of when spm and animal_licked are together
        HIT = spm(i:i+w)&animal_licked_sim(i:i+w);
        N_Hits = sum(HIT);
        
        % Output the mask when spm is 1, but animal_licked is zero
        FALSE_ALARM = spm(i:i+w)&~animal_licked_sim(i:i+w);
        N_FalseAlarm = sum(FALSE_ALARM);
        
        
        %calculates the probability
        %% calculates probability (new method)
        totalprotsim=[];
        protoSim=[];
        % calculates the number of ways to allocate Nlicks in 20 trials
        totalprotsim = nchoosek(20,Nlicks);
        
        Nhhh=N_Hits;
        Nfff=N_FalseAlarm;
        
        
        Nhhh=N_Hits;
        Nfff=N_FalseAlarm;
        
        for r=1:20
            if Nfff<0 | Nhhh>AvailSplus
                break
            end
            
            protoSim(end+1) =nchoosek(AvailSplus,Nhhh)*nchoosek(AvailSminus,Nfff);
            Nhhh=Nhhh+1;
            Nfff= Nfff - 1;
            
        end
        bof = sum(protoSim);
        ProbaSim(i+9)= bof/totalprotsim;
        
        if i==1
            ProbaSim(1:9)=ProbaSim(10);
        end
    end
    ProbaSim(end+1:lena)=ProbaSim(end);
    SimRand(nsim,:)= ProbaSim;
end
%%

AveSim = mean(SimRand);
prctile5 = prctile(SimRand,5);

=======
function [AveSim prctile5]= drgSimRand(spm,animal_licked)

%These are test data
%Uncomment if your are testing the function
%Animal licks all the time
% spm= [ 1 1 1 1 0 1 0 1 1 0 1 0 1 1 1 1 1 0 0 1 1 1 0 0 1 0 1 0 1 0 1 1 0 1 0 0 1 ];
% animal_licked= [ 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 1 ];
%
% %Animal licks only for S+
% spm=           [ 1 1 1 1 0 1 0 1 1 0 1 0 1 1 1 1 1 0 0 1 1 1 0 0 1 0 1 0 1 0 1 1 0 1 0 0 1 ];
% animal_licked= [ 1 1 1 1 0 1 0 1 1 0 1 0 1 1 1 1 1 0 0 1 1 1 0 0 1 0 1 0 1 0 1 1 0 1 0 0 1 ];

w= 19;
lena=length(spm);
SimRand = zeros(100,lena);
AveSim =[];
stdSim =[];
infstd =[];
AveP =[];
animal_licked_sim =[];


%% randomizes animal_licked and calculates Probability

for nsim = 1:100
    ProbaSim=ones(1,w);
    animal_licked_sim = permVec(animal_licked);
    
    % For a moving window of length w
    % Outputs w+1 elements
    
    
    for i=1:lena-w
        
        %counts the number of available S+ and S- in the window
        
        Splus = animal_licked_sim(i:i+w);
        AvailSplus = sum(Splus);
        Smin = ~animal_licked_sim(i:i+w);
        AvailSminus = sum(Smin);
        Lks = spm(i:i+w);
        Nlicks = sum(Lks);
        
        % Output the mask of when spm and animal_licked are together
        HIT = spm(i:i+w)&animal_licked_sim(i:i+w);
        N_Hits = sum(HIT);
        
        % Output the mask when spm is 1, but animal_licked is zero
        FALSE_ALARM = spm(i:i+w)&~animal_licked_sim(i:i+w);
        N_FalseAlarm = sum(FALSE_ALARM);
        
        
        %calculates the probability
        %% calculates probability (new method)
        totalprotsim=[];
        protoSim=[];
        % calculates the number of ways to allocate Nlicks in 20 trials
        totalprotsim = nchoosek(20,Nlicks);
        
        Nhhh=N_Hits;
        Nfff=N_FalseAlarm;
        
        
        Nhhh=N_Hits;
        Nfff=N_FalseAlarm;
        
        for r=1:20
            if Nfff<0 | Nhhh>AvailSplus
                break
            end
            
            protoSim(end+1) =nchoosek(AvailSplus,Nhhh)*nchoosek(AvailSminus,Nfff);
            Nhhh=Nhhh+1;
            Nfff= Nfff - 1;
            
        end
        bof = sum(protoSim);
        ProbaSim(i+9)= bof/totalprotsim;
        
        if i==1
            ProbaSim(1:9)=ProbaSim(10);
        end
    end
    ProbaSim(end+1:lena)=ProbaSim(end);
    SimRand(nsim,:)= ProbaSim;
end
%%

AveSim = mean(SimRand);
prctile5 = prctile(SimRand,5);

>>>>>>> 362cddba4db155729b2e4c347f0e0147096a5855
