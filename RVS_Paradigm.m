clear all;
rand('state',sum(100*clock));
Tstart = clock;
timestamp = [num2str(Tstart(4)) 'h' num2str(Tstart(5)) 'm' num2str(ceil(Tstart(6))) 's'];

%% *******************************************************************
%                       INTRODUCTION
%*********************************************************************
% PARADIGM:
% This script implements a tactile rehearsal/refreshment task.
%
% Subjects are presented with an initial sequence of stimuli (length 3),
% which is repeated during a Stimulation phase for until 4,6,8 stimulus is
% presented. The last presentation of stimuli is accompanied with a CUE, which tells
% the subject that from then on Refreshment starts
% During the Refreshment-Period (6,8,10 stimulus) a visual Rehearsal-cue is presented
% to speed the refreshment.
% After that period, a Target-Stimulus is presented, which either matches
% the current to-be-refreshed stimulus in the sequence, or not.
% Subjects indicate via left/right button-press if it was a correct target
% or not.
%
% TECHNICAL DETAILS:

% OUTPUTS:
% Cogent-Logfiles are generated
% A MATLAB-Logfile is generated (to be used for the SPM Analysis)


%% *******************************************************************
%                      SUBJECT SETTINGS
%*********************************************************************
subject_no    = 99;
if subject_no < 10
    subject_name  = ['sj_0' num2str(subject_no)];
else
    subject_name  = ['sj_' num2str(subject_no)];
end
run           = 8;
display       = 1;

% Load Stimuli
Stimuli = RVS_Stimuli;

% RESPONSE-BUTTON RANDOMIZATION
% acccording to randomization Scheme
if mod(subject_no,2) == 0 % Even numbers
    yes = 97; %left(97)
    no =  98; %right(98)
else
    yes = 98; %right
    no =  97; %left
end

tic_delay = 700;
%% *******************************************************************
%                  COMPUTER SPECIFIC SETTINGS LOADED
%*********************************************************************
% Set parameters which are different paths, dependent on which computer the
% Script is run. (Presentation is the computer at the scanner)

cd('C:\Tang\RVS_EEG');



%% *******************************************************************
%                        INITIATE QUAEROSYS
%*********************************************************************
% Load Stimulator-Setup files
loadlibrary('C:\Jan\TacMoDo\Quaerobox\QuaeroSys\stimlib0.dll', 'C:\Jan\TacMoDo\Quaerobox\QuaeroSys\stimlibrel.h');
calllib('stimlib0', 'setProperty', 'local_buffer_size', 2000000);
calllib('stimlib0', 'initStimulator', '98J9007703J8RJHA9J27J0B053KBKR'); % License October 16


%% *******************************************************************
%                       EXPERIMENTAL SETTINGS
%*********************************************************************
%ITIs             = [2000 3000]; % in ms
Anfangspause     = 4000;
ref_rate         = 1250;
sequence_length  = 3;
stimdur          = 250;
rampdur          = 100;   % TTS: PLEASE DOUBLECHECK WITH Jans Scripts how they did it (e-mail in case) ->YH has code
dt               = 1;     % resoultion of signal in ms (QueroSys Stimulator has a maximal resolution of 0.5ms)
T_download       = 500;   % Time given for the QuaeroSys to download stimuli in ms
ParAddr=(hex2dec('378')); %Parallel Port Adress
trigdur          = 2;     % duration of the trigger

% Compute Stimuli for downloading to QuaeroSys
x= linspace(0,stimdur, (stimdur*2)+1);
ramp1 = 0:(1/rampdur/0.5):1;
ramp2 = 1:-(1/rampdur/0.5):0;
envelope = [ramp1 ones(1,(length(x)-length(ramp1)-length(ramp2))) ramp2];


%SEQUENCES
% Are determined in separate function
sequences = RVSsequences(Stimuli.Frequencies,sequence_length);

%% *******************************************************************
%                   COMPOSITION OF DESIGN-MATRIX
%*********************************************************************
% This is the composition of the design where the order of trials is
% randomized.

% Randomization of correct/incorrect Trials (Desired Responses)
correct = [ones(12,1).' zeros(12,1).'];

% Specification of Trials within one Run including the Desired Response
Trials_all =   [sequences.fre_sequences(1:24,1:3).'                 ; % Stim_seq 
                sequences.fre_sequences(1:24,4).'                   ; % Stim_dur
                sequences.fre_sequences(1:24,5:7).'                 ; % Rehearsal_seq
                sequences.fre_sequences(1:24,8).'                   ; % Delay_dur
                randsample(correct,24)                              ; % Correct/Incorrect
                1:24                                              ] ; % Sequence

% Changing the format in lines 2,3 from seconds into ms
% while we have a scale in it, duration for each frequency is 1250ms
Trials_all(4,:) = Trials_all(4,:)*ref_rate;
Trials_all(8,:) = Trials_all(8,:)*ref_rate;


% Rename Trials_all
trial_order = randperm(length(Trials_all));
for i = 1:length(trial_order)
    Trials(:,i) = Trials_all(:,trial_order(i));
end
clear i Trials_all trial_order;


% Selecting which stimulus has to be the Target-Stimulus
Targets = [];
for i = 1: size(Trials,2)
    sequence_dur = Trials(4,i) + Trials(8,i);
    num_seq = sequence_dur/ref_rate;
    cor_next = mod(num_seq,3)+1;
    Targets_real(i) = sequences.fre_sequences(Trials(10,i),cor_next);
    if Trials(9,i) == 1 % Correct Trials
        Targets(i) = sequences.fre_sequences(Trials(10,i),cor_next);
    else                % Incorrect Trials
        alt_targets  = 1:3;
        alt_targets(cor_next) = [];
        Targets(i) = sequences.fre_sequences(Trials(10,i),randsample(alt_targets,1));
    end   
end

% Adding Tragets to Trial Matrix and adding ITIs, ITI is between 2000 ~
% 4000ms
Trials = [Trials; ...
    Targets;...
    Targets_real;...
    [2000 + randi(2000,22,1)' 3999 0]];
% The last ITI is set to zero. ATTENTION: This is hardcoded, in case
% the overall number of trial changes, this does not work anymore
%randsample(repmat(ITIs,1,
%floor((size(Trials,2)/length(ITIs)))),size(Trials,2))];
clear alt_targets num_seq sequence_dur cor_next Targets;

% Composing a vector that contains the timing information
Timing = [Anfangspause];    % Start of first trial after 6s
for i = 1:(length(Trials)-1)
    %last_onset + Sequence    + Delay       + Response +  ITI
    next_onset = Timing(i)  + Trials(4,i) + Trials(8,i) + 2000     + Trials(13,i); % Adding of time
    Timing(i+1) = next_onset;
    clear next_onset;
end

% For a simplified Overview all information are rearranged and composed to
% the Design-Matrix
Design = [ Timing      ;...      % 1:Trial_onset time in ms 
    Trials(10,:) ;...            % 2:Sequence 
    Trials(1:3,:);...            % 3-5: Stim_seq 
    Trials(4,:) ;...             % 6:Sequence Duration
    Trials(5:7,:) ;...           % 7-9: Reheasal_seq
    Trials(8,:) ;...             % 10:Delay Duration
    Trials(11,:) ;...            % 11:Target as Index of Sequence
    Trials(12,:) ;...            % 12:Real Target
    Trials(9,:) ;...             % 13:Correct/Incorrect
    Trials(13,:) ];              % 14:ITI

%Clearing of unused variables
clear Timing Trials i Anfangspause ITIs correct;


%% ************************************************************************
%                             COGENT-Setup
%**************************************************************************
config_keyboard;
config_log(fullfile(pwd, '\RVS_Logs', ['RVS_log_' subject_name...
    '_' num2str(run) '_' timestamp '.txt']));
config_display(display,3,[0 0 0],[1 1 1], 'Arial', 30, 4, 0);
% mode - window mode ( 0=window, 1=full primary screen, 2= full second screen )
% resolution - screen resolution (1=640x480, 2=800x600, 3=1024x768, 4=1152x864, 5=1280x1024, 6=1600x1200)
start_cogent;
cgloadlib;
clearkeys;
lptwrite(ParAddr,0); %initial reset Parallel Port

% *******************************************
%       Definition der Sprites
%*******************************************
% 1: Fixation
% 2: Wait for Trigger
% 11: Stimulation
% 15: Rehearsal
% 20: Response ('?')
% 21: Positive Feedback
% 22: Negative Feedback

% *******************
% Fixation ITI
cgmakesprite(1,500,100,[0,0,0]);
cgsetsprite(1);
cgpencol(0.5,0.5,0.5);
cgfont('Arial',50);
cgtext('+',0,0);
cgsetsprite(0);

%*********************
%  PLAY
cgmakesprite(11,500,100,[0,0,0]);
cgsetsprite(11);
cgfont('Arial',50);
cgtext('>',0,0);
cgsetsprite(0);

%*********************
%  Rehearsal
cgmakesprite(15,500,100,[0,0,0]);
cgsetsprite(15);
cgpencol(1,1,1);
cgfont('Arial',50);
cgtext('+',0,0);
cgsetsprite(0);

%*********************
%  Response
cgmakesprite(20,500,100,[0,0,0]);
cgsetsprite(20);
cgpencol(1,1,1);
cgfont('Arial',50);
cgtext('?',0,0);
cgsetsprite(0);

%*********************
%  Positive Feedback
cgmakesprite(21,500,100,[0,0,0]);
cgsetsprite(21);
cgpencol(1,1,1);
cgfont('Arial',50);
cgtext('+ + +',0,0);
cgsetsprite(0);

%*********************
%  Negative Feedback
cgmakesprite(22,500,100,[0,0,0]);
cgsetsprite(22);
cgpencol(1,1,1);
cgfont('Arial',50);
cgtext('- + -',0,0);
cgsetsprite(0);

%*********************
%  Performance Feedback
cgmakesprite(23,500,100,[0,0,0]);
cgsetsprite(23);
cgpencol(1,1,1);
cgfont('Arial',30);
cgtext('Your Performance is',0,0);
cgsetsprite(0);

% Finished defining sprites
cgflip(0,0,0);



%% ************************************************************************
%             Initiation of Logfile
%**************************************************************************
% MATLAB-LOGFILE
log_path = fullfile(pwd, '\RVS_Logs', ['RVS_log_' subject_name '_' ...
    num2str(run) '_' timestamp '.mat']);
log_RVS.subject_no   = subject_no;
log_RVS.subject_name = subject_name;
log_RVS.run          = run;
log_RVS.date         = date;
log_RVS.time         = clock;
log_RVS.Design       = Design;
log_RVS.responses    = zeros(3,24); % ATTENTION this is hardcoded - if amount of trials change, this need adjustment
log_RVS.timing       = zeros(20,48);
log_RVS.trigger      = zeros(24,24); % 1--sequence, 2--stimulation duration, 3:10--stimulation, 11:20--rehearsal, 21--target, 22--response, 23--correct response, 24--incorrect response
save(log_path, 'log_RVS');

% COGENT-LOGFILE
logstring(strcat('Subject: ', subject_name , '  Run Number: ', num2str(run),...
    '_', timestamp));
log_clock = clock;
logstring(['Start-time: ', num2str(log_clock(4)), ':', num2str(log_clock(5)), ':', num2str(log_clock(6)), ' Date: ', date ]);

%% ************************************************************************
%                      Start of Experiment
%%*************************************************************************
Start_time = time;     % time counter
log_RVS.start_time    = Start_time;
logstring('START');
cgflip(0,0,0); cgdrawsprite(1,0,0); cgflip(0,0,0);


% It appears best practice with the QuaeroSys Stimulator to first run one
% Stimulus, befor the real experiment starts. Otherwise the first stimulus
% is often delayed. The reason for this is unknown - this is just a
% workaround, but no real fix - However, it works.
% For this stimulus, all pins are set to an elevation of 0, for 100ms
RVS_stimDownload_trig(ones(1,100), dt);
lptwrite(ParAddr,130);
wait(trigdur);
lptwrite(ParAddr,0);
calllib('stimlib0', 'stopStimulation'); % clear local buffer
calllib('stimlib0', 'resetStimulator'); % clear remote buffer

% This is a for-loop that runs column for columnd through the Design-Matrix
% and within the loop one Trial is performed

for i=1:length(Design)
    
    % To adjust the first trigger, it can not work well
    
    % SEQUENCE-LOOP
    % This loop will be performed to play the tactile stimuli together with
    % visual cues (color-determines Pattern/Frequency-Trials)
    % The timing is done in relation to the first Trigger (or start time)
    while time <= Start_time + (Design(1,i) - tic_delay)
    end
    
    % Have the sequence of the current trial available in a Matrix
    Curr_Seq = sequences.fre_sequences(Design(2,i),1:3);
    logstring(['Trial-' num2str(i) ': StimDur: ' num2str(Design(6,i)/1000) 'sec  DelayDur: ' num2str(Design(10,i)/1000) 'sec']);
    
    %ATTENTION: tic IS 700 ms before the Trial
    %All timing within the trial is done in relation to this tic
    % For this the variable 'tic_delay' is used
    tic;
    
    % Set trial triggers, sequence and stimulation duration
    trig = 50+Design(2,i);
    lptwrite(ParAddr,trig);
    wait(trigdur);
    sent_trig = lptread(ParAddr);
    log_RVS.trigger(1,i) = sent_trig;
    lptwrite(ParAddr,0);
    wait(50);
    trig = Design(6,i)/ref_rate;
    lptwrite(ParAddr,trig);
    wait(trigdur);
    sent_trig = lptread(ParAddr);
    log_RVS.trigger(2,i) = sent_trig;
    lptwrite(ParAddr,0);
    wait(50);
    trig = 80 + Design(10,i)/ref_rate;
    lptwrite(ParAddr,trig);
    wait(trigdur);
    sent_trig = lptread(ParAddr);
    log_RVS.trigger(3,i) = sent_trig;
    lptwrite(ParAddr,0);
    
    
    for s = 1:(Design(6,i)/ref_rate) % Determine that the s-th stimulus in sequence is executed
        
        % Dirty programming to determine at what index wihtin the sequence
        % we are with the current stimulus.        
        stim_indx = mod(s,sequence_length);
        if stim_indx ==0
            stim_indx = 3;
        end
        
        % Set triggers
        trig = 200 + Curr_Seq(stim_indx)*10 + s - 1; % set triggers with position information
        
        % download stimulus 500 ms prior to the onset of each trial
        while toc < (((s-1)*ref_rate)-T_download + tic_delay)/1000
        end
        
        pinhub_freq = (sin(2*pi*(Stimuli.Frequencies(Curr_Seq(stim_indx)))*x/1000)/2)+0.5;
        pinhub_freq = pinhub_freq.*envelope;
        pinhub_freq = pinhub_freq(1:(dt/0.5):length(pinhub_freq)); % downsampling
        
        RVS_stimDownload_trig(pinhub_freq, dt);
        
        %PLAY STIMULUS (onset of the trial)
        while toc < ((s-1)*ref_rate + tic_delay)/1000
        end
        
        if s == (Design(6,i)/ref_rate) % For the last stimulus in the Sequence
            cgdrawsprite(11,0,0); cgflip(0,0,0); % PLAY
        else
            cgdrawsprite(15,0,0); cgflip(0,0,0); % Rehearsal Cue
        end
        
        % Write the stimuli information to the log file
        logstring(['Stimulus : ' num2str(Curr_Seq(stim_indx))]);
      
        % Sent the triggers to the stimulator, start the QuaeroSys
        % Stimulus
        
        lptwrite(ParAddr,trig);
        sent_trig = lptread(ParAddr);
        log_RVS.trigger(s+3,i) = sent_trig;
        wait(trigdur);
        lptwrite(ParAddr,0);
        wait(stimdur-trigdur);
        calllib('stimlib0', 'stopStimulation'); % clear local buffer        
        calllib('stimlib0', 'resetStimulator'); % clear remote buffer
        
        
        % log_RVS.timing(s,i+24) = toc - log_RVS.timing(s,i);
        cgdrawsprite(1,0,0); cgflip(0,0,0);     % Flip to fixation-cross
        
    end
    
    
    % REHEARSAL-LOOP
    % During this loop the Rehearsal cues are displayed. Programming as above
    for s = 1:(Design(10,i)/ref_rate)
        % Dirty programming to determine at what index wihtin the sequence
        % we are with the current stimulus.
        stim_indx_R = mod((s+stim_indx),sequence_length);
        if stim_indx_R ==0
            stim_indx_R = 3;
        end
        
        % trigger value
        trig = Curr_Seq(stim_indx_R)*10 + s - 1;  
        
        % DISPLAY Rehearsal CUE
        while toc < (Design(6,i)+((s-1)*ref_rate) + tic_delay)/1000  % As the tic is done 1TR before the start of the Trial, one has to wait for 2s (ATTENTION: tic/toc is in s (not in ms))
        end       
        
        cgdrawsprite(15,0,0); cgflip(0,0,0); % Rehearsal Cue
        logstring(['Rehearsal: ' num2str(Curr_Seq(stim_indx_R))]);
        
        % Sent the triggers to the EEG
        log_RVS.timing(s,i)= toc;      
        lptwrite(ParAddr,trig);  
        wait(trigdur);
        sent_trig = lptread(ParAddr);
        log_RVS.trigger((s+8)+3,i) = sent_trig;
        lptwrite(ParAddr,0); 
        wait(stimdur-trigdur);
        log_RVS.timing(s,i+24) = toc - log_RVS.timing(s,i);
        
        cgdrawsprite(1,0,0); cgflip(0,0,0); % Fixation
        
    end
    
    
    % TARGET-STIMULUS
    % Stimulus Download
    while toc < (Design(6,i)+ Design(10,i) + tic_delay - T_download)/1000  % As the tic is done 1TR before the start of the Trial, one has to wait for 2s (ATTENTION: tic/toc is in s (not in ms))
    end
    
    pinhub_freq = (sin(2*pi*(Stimuli.Frequencies(Design(11,i)))*x/1000)/2)+0.5;
    pinhub_freq = pinhub_freq.*envelope;
    pinhub_freq = pinhub_freq(1:(dt/0.5):length(pinhub_freq)); % downsampling
    RVS_stimDownload_trig(pinhub_freq, dt);
    
    % trigger value
    trig = 209 + Design(11,i)*10; % set triggers with position information
    
    % Stimulus Play + Cue
    while toc < (Design(6,i)+ Design(10,i) + tic_delay )/1000
    end
    
    % t_start_stim = toc;
    cgdrawsprite(20,0,0); cgflip(0,0,0);    % RESPONSE
    logstring(['Target-Stimulus: ' num2str(Design(11,i))]);
    
    % Start the QuaeroSys Stimulus    
    lptwrite(ParAddr,trig);
    sent_trig = lptread(ParAddr);
    log_RVS.trigger(22,i) = sent_trig;
    wait(trigdur);
    lptwrite(ParAddr,0);
    wait(stimdur-trigdur);   
%     calllib('stimlib0', 'stopStimulation'); % clear local buffer
%     calllib('stimlib0', 'resetStimulator'); % clear remote buffer
    
    
    % RESPONSE
    zeit = time;
    [response_key, press] = waitkeydown(1750, [yes no]); % Reads only the yes/no buttons in for 1500ms        
    trig = 9;
    lptwrite(ParAddr,trig);
    wait(trigdur);
    sent_trig = lptread(ParAddr);
    log_RVS.trigger(23,i) = sent_trig;
    lptwrite(ParAddr,0);
    cgdrawsprite(1,0,0); cgflip(0,0,0);   
    
    if ~isempty(response_key)
        rt=press-zeit;
        if response_key == yes
            given_response = 1;
        elseif response_key == no
            given_response =0;
        end
    end
    
    % Feedback is given imediately after the response
    if isempty(response_key)
        % If there was no response, the trials are marked with 8s...
        % Why it is 8s? They visually pop out when checking the
        % results-matrix... there is no other reason...
        logstring('Response: NoResponse');
        answer = 8;
        rt     = 8888;
        response_key = 0;
        % Flickering '- + -' Display
        cgdrawsprite(22,0,0); cgflip(0,0,0); wait(100);
        cgdrawsprite(1,0,0); cgflip(0,0,0);  wait(100);
        cgdrawsprite(22,0,0); cgflip(0,0,0); wait(100);
        cgdrawsprite(1,0,0); cgflip(0,0,0);  wait(100);
        cgdrawsprite(22,0,0); cgflip(0,0,0); wait(100);
        cgdrawsprite(1,0,0); cgflip(0,0,0);
        
    elseif (Design(13,i) == 1 && response_key == yes) || (Design(13,i) == 0 && response_key == no)
        answer = 1;
        trig = 2;
        lptwrite(ParAddr,trig);
        wait(trigdur);
        sent_trig = lptread(ParAddr);
        log_RVS.trigger(24,i) = sent_trig;
        lptwrite(ParAddr,0);        
        logstring(['Response: ' num2str(response_key) ' Meaning: ' num2str(given_response) ' Answer: ' num2str(answer)]);
        cgdrawsprite(21,0,0); cgflip(0,0,0); wait(500);
        cgdrawsprite(1,0,0); cgflip(0,0,0);
        
    else
        answer = 0;
        trig = 3;
        lptwrite(ParAddr,trig);
        wait(trigdur);
        sent_trig = lptread(ParAddr);
        log_RVS.trigger(24,i) = sent_trig;
        lptwrite(ParAddr,0);
        logstring(['Response: ' num2str(response_key) ' Meaning: ' num2str(given_response) ' Answer: ' num2str(answer)]);
        cgdrawsprite(22,0,0); cgflip(0,0,0); wait(500);
        cgdrawsprite(1,0,0); cgflip(0,0,0);
    end
    
    %reset the stimulator
    calllib('stimlib0', 'stopStimulation'); % clear local buffer
    calllib('stimlib0', 'resetStimulator'); % clear remote buffer
    
    
    % lptwrite(ParAddr,(trig+answer));
    % Saving the responses
    log_RVS.responses(:,i) = [response_key;answer;rt];
    
    % Reseting the stimulator - this is best practice, as otherwise the
    % stimulator sometimes makes problems.    

    wait(100);
    RVS_stimDownload_trig(zeros(1,100), dt);
    lptwrite(ParAddr,130);
    wait(trigdur);
    lptwrite(ParAddr,0);
    calllib('stimlib0', 'stopStimulation'); % clear local buffer
    calllib('stimlib0', 'resetStimulator'); % clear remote buffer
    
    % Saving the MATLAB-Log to hard-disc
    save(log_path, 'log_RVS');
    clearkeys;   
    
    
end


%% **************************************************************
%                       QUICK STATISTIK
%****************************************************************
perf_Rehearsal = sum(log_RVS.responses(2,:)==1) /length(log_RVS.responses(2,:))
% mean_rt        = sum(log_RVS.responses(3,log_RVS.responses(3,:)~=8888))/sum(log_RVS.responses(3,:)~=8888)
% idx = sub2ind(size(Design(3:5,:)),[1:3],Design(11,:).');



%% **************************************************************
%                         THE END
%****************************************************************
logstring('The END');
clear dio;
cgshut;
stop_cogent;
calllib('stimlib0', 'closeStimulator');
clear all;