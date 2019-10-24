function [theta_grid,phaseVec,phase_time] = theta_map(phase,timeVec,sampleRate,nbins)

% since Caitlin already has the theta phase
% just take the phase we already have, don't need to find it

%%%%%% taken out by Alec and Holger
% %compute instantaneous phase
% hilb_eeg = hilbert(filt_eeg); % compute hilbert transform
% phase = atan2(imag(hilb_eeg),real(hilb_eeg)); %inverse tangent (-pi to pi)

% rescales to 0 to 2pi
ind = phase <0;
phase(ind) = phase(ind)+2*pi; % from 0 to 2*pi

% resample to 50 Hz
phase_time = resample(phase, 50, sampleRate);

theta_grid = zeros(length(timeVec),nbins);
phaseVec = 2*pi/nbins/2:2*pi/nbins:2*pi-2*pi/nbins/2;

for i = 1:numel(timeVec)

    % figure out the theta index
%     [~, idx] = min(abs(phase_time(i)-phaseVec)); % replaced by Holger
    [~, idx] = min(abs(phase_time(i)-phaseVec));
    theta_grid(i,idx) = 1;

end

return
