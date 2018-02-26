function data_this_trial=drgGetThisTrialRHD(handles,trialNo)

% read_Intan_RHD2000_file
%
% Version 2.0, 20 October 2016
%
%   Reads data ffor a single trial
%
% Reads Intan Technologies RHD2000 data file generated by evaluation board
% GUI or Intan Recording Controller.  Data are parsed and placed into
% variables that appear in the base MATLAB workspace.  Therefore, it is
% recommended to execute a 'clear' command before running this program to
% clear all other variables from the base workspace.
%
% Example:
% >> clear
% >> read_Intan_RHD2000_file
% >> whos
% >> amplifier_channels(1)
% >> plot(t_amplifier, amplifier_data(1,:))

sessionNo=1;

tic;
fid = fopen(handles.drg.drta_p.fullName, 'r');

s = dir(handles.drg.drta_p.fullName);
filesize = s.bytes;

% Check 'magic number' at beginning of file to make sure this is an Intan
% Technologies RHD2000 data file.
magic_number = fread(fid, 1, 'uint32');
if magic_number ~= hex2dec('c6912702')
    error('Unrecognized file type.');
end

% Read version number.
data_file_main_version_number = fread(fid, 1, 'int16');
data_file_secondary_version_number = fread(fid, 1, 'int16');


if (data_file_main_version_number == 1)
    num_samples_per_data_block = 60;
else
    num_samples_per_data_block = 128;
end
%

%
amplifier_index = 1;
board_adc_index = 1;
board_dig_in_index = 1;

num_amplifier_channels=handles.drg.session(sessionNo).draq_d.num_amplifier_channels;
num_samples_per_data_block=handles.drg.session(sessionNo).draq_d.num_samples_per_data_block;
num_board_adc_channels=handles.drg.session(sessionNo).draq_d.num_board_adc_channels;
num_board_dig_in_channels=handles.drg.session(sessionNo).draq_d.num_board_dig_in_channels;

%Now read the data
for i=handles.drg.session(sessionNo).draq_d.start_blockNo(trialNo):handles.drg.session(sessionNo).draq_d.end_blockNo(trialNo)
    % In version 1.2, we moved from saving timestamps as unsigned
    % integeters to signed integers to accomidate negative (adjusted)
    % timestamps for pretrigger data.
    
    %Read amplifier channels
    if (num_amplifier_channels > 0)
        fseek(fid,handles.drg.session(sessionNo).draq_d.offset_start_ch(i),'bof');
        amplifier_data(:, amplifier_index:(amplifier_index + num_samples_per_data_block - 1)) = fread(fid, [num_samples_per_data_block, num_amplifier_channels], 'uint16')';
    end
    
    if (num_board_adc_channels > 0)
        fseek(fid,handles.drg.session(sessionNo).draq_d.offset_start_adc(i),'bof');
        board_adc_data(:, board_adc_index:(board_adc_index + num_samples_per_data_block - 1)) = fread(fid, [num_samples_per_data_block, num_board_adc_channels], 'uint16')';
    end
    
    if (num_board_dig_in_channels > 0)
        fseek(fid,handles.drg.session(sessionNo).draq_d.offset_start_dig(i),'bof');
        board_dig_in_raw(board_dig_in_index:(board_dig_in_index + num_samples_per_data_block - 1)) = fread(fid, num_samples_per_data_block, 'uint16');
    end
    
    amplifier_index = amplifier_index + num_samples_per_data_block;
    board_adc_index = board_adc_index + num_samples_per_data_block;
    board_dig_in_index = board_dig_in_index + num_samples_per_data_block;
    
end



%end

% Close data file.
fclose(fid);



% Extract digital input channels to separate variables.
for i=1:handles.drg.session(sessionNo).draq_d.num_board_dig_in_channels
    mask = 2^(handles.drg.session(sessionNo).draq_d.board_dig_in_channels(i).native_order) * ones(size(board_dig_in_raw));
    board_dig_in_data(i, :) = (bitand(board_dig_in_raw, mask) > 0);
end

% Scale voltage levels appropriately.
amplifier_data = 0.195 * (amplifier_data - 32768); % units = microvolts

if (handles.drg.session(sessionNo).draq_d.eval_board_mode == 1)
    board_adc_data = 152.59e-6 * (board_adc_data - 32768); % units = volts
elseif (handles.drg.session(sessionNo).draq_d.eval_board_mode == 13) % Intan Recording Controller
    board_adc_data = 312.5e-6 * (board_adc_data - 32768); % units = volts
else
    board_adc_data = 50.354e-6 * board_adc_data; % units = volts
end

%Enter the digital input channel
digital_input=board_dig_in_data(1,:)+2*board_dig_in_data(2,:)+4*board_dig_in_data(3,:)...
    +8*board_dig_in_data(4,:)+16*board_dig_in_data(5,:)+32*board_dig_in_data(6,:)...
    +64*board_dig_in_data(7,:);

%Setup the output as used by drta
data_this_trial=zeros(length(digital_input),22);

%Enter the electrode recordings
data_this_trial(:,1:16)=amplifier_data';


%Enter the trigger (bit 8)
data_this_trial(:,17)=1000*board_dig_in_data(8,:);

%Enter the four votage inputes: shiff, lick, photodiode and laser trigger
data_this_trial(:,18:21)=board_adc_data(1:4,:)';

data_this_trial(:,22)=digital_input;

return




