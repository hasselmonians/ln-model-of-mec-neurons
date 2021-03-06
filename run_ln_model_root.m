% function [Results] = run_LNmodel_HD_root(root,active_lfp,cel)
% runs the linear non-linear model (Hardcastle et al.) on CMBHOME root object loaded into the workspace

addpath(genpath('D:\Dropbox (hasselmonians)\hdannenb\Scripts\LNmodel_HD'))
%% Description of linear non-linear model published by Hardcastle et al.
% This script is segmented into several parts. First, the data (an
% example cell) is loaded. Then, 15 LN models are fit to the
% cell's spike train. Each model uses information about
% position, head direction, running speed, theta phase,
% or some combination thereof, to predict a section of the
% spike train. Model fitting and model performance is computed through
% 10-fold cross-validation, and the minimization procedure is carried out
% through fminunc. Next, a forward-search procedure is
% implemented to find the simplest 'best' model describing this spike
% train. Following this, the firing rate tuning curves are computed, and
% these - along with the model-derived response profiles and the model
% performance and results of the selection procedure are plotted.

% Code as implemented in Hardcastle, Maheswaranthan, Ganguli, Giocomo,
% Neuron 2017
% V1: Kiah Hardcastle, March 16, 2017

% rename to match notations of Hardcastle code and get
% central position
post = root.ts;
posx_c = root.sx;
posy_c = root.sy;
% let the position points start at 0
posx_c = posx_c - min(posx_c);
posy_c = posy_c - min(posy_c);

spktimes = get_spktimes_of_cel(root,cel);
% randomly delete spikes, so that the total number of
% spikes is 1800 (equivalent to 0.75Hz in a 20-min
% recording)

% number of spikes
nmbr_spikes = length(spktimes);
% nmbr_spikes = 1800;

if length(spktimes) < nmbr_spikes
    return
end
p = randperm(length(spktimes),nmbr_spikes);
spikes = spktimes(p);
spikes = sort(spikes);

% get spike train
[real_fr,~] = get_InstFR(spktimes,post,root.fs_video,'filter_length',125,'filter_type','Gauss');
[smooth_fr,spiketrain] = get_InstFR(spikes,post,root.fs_video,'filter_length',125,'filter_type','Gauss');

box_size = 100;
eeg_sample_rate = 600;%250;
theta_freq_range = [6 10];

% get theta filtered EEG
root.active_lfp = active_lfp;
eeg_4800 = root.b_lfp(root.active_lfp).signal;
% eeg_250 = resample(eeg_4800,250,4800);
eeg_600 = resample(eeg_4800,600,4800);
% filt_eeg = CMBHOME.LFP.BandpassFilter(eeg_250, eeg_sample_rate, theta_freq_range);
filt_eeg = CMBHOME.LFP.BandpassFilter(eeg_600, eeg_sample_rate, theta_freq_range);

% description of variables included:
% box_size = length (in cm) of one side of the square box
% post = vector of time (seconds) at every 20 ms time bin
% spiketrain = vector of the # of spikes in each 20 ms time bin
% posx = x-position of left LED every 20 ms
% posx2 = x-position of right LED every 20 ms
% posx_c = x-position in middle of LEDs
% posy = y-position of left LED every 20 ms
% posy2 = y-posiiton of right LED every 20 ms
% posy_c = y-position in middle of LEDs
% filt_eeg = local field potential, filtered for theta frequency (4-12 Hz)
% eeg_sample_rate = sample rate of filt_eeg (250 Hz)
% sampleRate = sampling rate of neural data and behavioral variable (50Hz)

%% fit the model
fprintf('(2/5) Fitting all linear-nonlinear (LN) models\n')
fit_all_ln_models

%% find the simplest model that best describes the spike train
fprintf('(3/5) Performing forward model selection\n')
select_best_model

%% Compute the firing-rate tuning curves
fprintf('(4/5) Computing tuning curves\n')
compute_all_tuning_curves

%% plot the results
fprintf('(5/5) Plotting performance and parameters\n')
plot_performance_and_parameters

%% read out and save the following:
tested_models = {'phst','phs','pht','pst','hst','ph','ps','pt','hs','ht','st','p','h','s','t;'};

% cell ID
% speed tuning curve
% best model, in particular if the cell's firing rate is
% modulated by speed
Results.cell_ID = strcat(root.name,'_T', num2str(cel(1)),'C', num2str(cel(2)));
Results.speed_response = speed_response; % modeled speed response tuning curve firing rate values
Results.speed_vector = speed_vector;% x-axis for speed tuning curve
if ~isnan(selected_model)
    Results.best_model = tested_models(selected_model);
else
    Results.best_model = nan;
end                % also save information about position and head direction
% tuning
Results.pos_response = pos_response; % for plotting use: imagesc(reshape(pos_response,20,20));
Results.hd_vector = hd_vector; % x-axis for head direction tuning curve
Results.hd_response = hd_response; % modeled head direction response tuning curve firing rate values
% read out the tuning curves
Results.pos_curve = pos_curve;
Results.hd_curve = hd_curve;
Results.speed_curve = speed_curve;
% read out the scale factors (alpha)
Results.scale_factor_pos = scale_factor_pos;
Results.scale_factor_hd = scale_factor_hd;
Results.scale_factor_spd = scale_factor_spd;
% save the firing rate and speed
Results.fr = smooth_fr;
Results.speed = speed;
% finally save log-likelihoods of model performance (m x n
% matrix with m = number of cross-validations, and n =
% number of models
Results.LLH_values = LLH_values;
Results.tested_models = tested_models;

rmpath('D:\Dropbox (hasselmonians)\hdannenb\Scripts\LNmodel_HD')
