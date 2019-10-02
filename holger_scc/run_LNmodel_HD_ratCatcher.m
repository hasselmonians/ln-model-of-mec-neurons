function Results = run_LNmodel_HD_ratCatcher(filename,filecode)

import CMBHOME.*

x = load(filename);

%%%%%% Apply speed filter (only for Caitlin's rat data)
x.root.b_vel = [];
x.root = x.root.AppendKalmanVel;

% set epoch to [-inf inf] for filtering out the time scales
x.root.epoch = [-inf inf];
% fix time
x.root = x.root.FixTime;
% pull out speed and time
spd = CMBHOME.Utils.ContinuizeEpochs(x.root.vel) * x.root.spatial_scale;
ts = CMBHOME.Utils.ContinuizeEpochs(x.root.ts);

% pull out sampling rate
sr = x.root.fs_video;

Results = struct;

% get spiking rate of cell in root object for whole
% recording length (will be cut later to match min_speed
% and epochs if desired)
x.root.epoch = [-inf inf];
spktimes = get_spktimes_of_cel(x.root,filecode);

% cut recording if necessary
if exist('x.cut_recording','var')
    if ~isempty(x.cut_recording{zelle})
        inds_cut = get_inds_from_epoch(x.root,x.cut_recording{zelle});
    else
        inds_cut = 1:length(spd);
    end
else
    inds_cut = 1:length(spd);
end
logical_inds_cut = zeros(size(spd));
logical_inds_cut(inds_cut) = 1;

% combine indices for min_speed, epochs, and cut_recording
logical_inds = logical_inds_cut;
% transform to logical
logical_inds = logical(logical_inds);

% run quality check on cell
quality_check = cell_quality_check_speed_coding(x.root,filecode,logical_inds_cut);
if quality_check

    % get epochs from logical_inds
    epochs_inds_for_mle = PeriodsAboveThr(double(logical_inds),1,2);
    epochs_mle = nan(length(epochs_inds_for_mle),2); % initialize matrix for epochs
    for i = 1:length(epochs_inds_for_mle)
        epochs_mle(i,1) = ts(epochs_inds_for_mle{i}(1));
        epochs_mle(i,2) = ts(epochs_inds_for_mle{i}(end));
    end
    % set epcohs
    x.root.epoch = epochs_mle;
    x.root.cel = filecode; % set cell
    spktimes = CMBHOME.Utils.ContinuizeEpochs(x.root.cel_ts); % get spike times
    post = CMBHOME.Utils.ContinuizeEpochs(x.root.ts); % get time
    speed = CMBHOME.Utils.ContinuizeEpochs(x.root.vel) * x.root.spatial_scale; % get speed

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
    post = x.root.ts;
    posx_c = x.root.sx;
    posy_c = x.root.sy;
    % let the position points start at 0
    posx_c = posx_c - min(posx_c);
    posy_c = posy_c - min(posy_c);

    %%%%%%%%%%%
    % randomly delete spikes, so that the total number of
    % spikes is 2400 (equivalent to 1Hz in a 40-min
    % recording)
    % nmbr_spikes = length(spktimes);
    nmbr_spikes = 2400;

    if length(spktimes) < nmbr_spikes
        return
    end
    p = randperm(length(spktimes),nmbr_spikes);
    spikes = spktimes(p);
    spikes = sort(spikes);
    %%%%%%%%%%%%%%

    % get spike train
    [~,real_spiketrain] = get_InstFR(spktimes,post,sr);
    [~,spiketrain] = get_InstFR(spikes,post,sr);

    box_size = 100;
    eeg_sample_rate = 600;%250;
    theta_freq_range = [6 10];

    %%% set active lfp to the other tetrode channel
    if filecode(1) == 1
        active_lfp = 2;
    elseif filecode(1) == 2
        active_lfp = 1;
    end

    % get theta filtered EEG
    x.root.active_lfp = active_lfp;
    eeg_4800 = x.root.b_lfp(x.root.active_lfp).signal;
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
    Results.cell_ID = strcat(x.root.name,'_T', num2str(filecode(1)),'C', num2str(filecode(2)));
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
    Results.theta_response = theta_response;
    % read out the tuning curves
    Results.pos_curve = pos_curve;
    Results.hd_curve = hd_curve;
    Results.speed_curve = speed_curve;
    Results.theta_curve = theta_curve;
    % save the spike train
    Results.spktrain = spiketrain;
    Results.real_spktrain = real_spiketrain;
    Results.speed = speed;
    % finally save log-likelihoods of model performance (m x n
    % matrix with m = number of cross-validations, and n =
    % number of models
    Results.LLH_values = LLH_values;
    Results.tested_models = tested_models;
    Results.nmbr_spikes = nmbr_spikes;

    %%%%%%%%%%%%%%%%%%%%%%%
    %%%%%%%%%%%%%%%%%%%%%%%

    clearvars x.root
else
    return
end
