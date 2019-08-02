function [quality_check] = cell_quality_check_speed_coding(root,cel,cut_recording)
% checks if a cell in root is stable over the whole recording length (first
% half vs. second half), and if mean firing rate is over 1Hz; cell is considered stable, if the means of the two
% halves differ less than a factor of two; cut_recording is a
% vector with 1 for all included indices and 0 for cut indices

% get spiking rate of cell in root object
spktimes = get_spktimes_of_cel(root,cel);
ts = CMBHOME.Utils.ContinuizeEpochs(root.ts);
[spkRate,~] = get_InstFR(spktimes,ts,root.fs_video,'filter_length',125,'filter_type','Gauss');

% adjust spkRate for cut_recording
if nargin == 3
    spkRate = spkRate(logical(cut_recording));
end

% check for minimum spiking rate
if mean(spkRate) <= 1
    min_spkRate_check = 0;
else
    min_spkRate_check = 1;
end

% compare spkRate of first half with spkRate of second half of the session
spkRate_first_half_of_session = spkRate(1:round(end/2));
spkRate_second_half_of_session = spkRate(round(end/2)+1:end);

test_value = mean(spkRate_first_half_of_session)/mean(spkRate_second_half_of_session);
if test_value < 0.5 || test_value > 2
    quality_check = min_spkRate_check*0;
else
    quality_check = min_spkRate_check*1;
end

