function run_LNmodel_HD_root_SCC(pathname)

addpath(genpath('/projectnb/hasselmogrp/hdannenb_1/Scripts/LNmodel_HD_SCC'))
addpath(genpath('/projectnb/hasselmogrp/Speed Modulation/CMBHOME'))

import CMBHOME.*

% start analysis
sprintf(['Analyzing ' pathname])

% list textfiles
% Textfiles = dir(strcat(pathname,'/','spike_sessions_speedMod_TimeScales.txt')); % for baseline and inhibition conditions
Textfiles = dir(strcat(pathname,'/','spike_sessions_speedMod_TimeScales_pooled_with_LS_control.txt')); % for 'none' condition

if ~isempty(Textfiles)
    for files_i = 1:length(Textfiles)
        textfile = Textfiles(files_i).name;
        
        fileID = fopen(strcat(pathname,'/',textfile));
        
        C = textscan(fileID,'%s');
        
        for session = 1:length(C{1})
            disp(['loading session ',C{1}{session}])
            x = load(strcat(pathname,'/',C{1}{session}));
            disp('session loaded')
            disp(['analyzing ',C{1}{session}])
            name = strsplit([pathname,'/',C{1}{session}],'/'); %name{end} contains the session name
            savename = ['LNmodel_1200spikes_' name{end}];
            
            %%%%%% Apply speed filter (only for Caitlin's rat data)
            %             x.root.b_vel = [];
            %             x.root = x.root.AppendKalmanVel;
            
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
            counter = 0; % set counter for cells
            
            for zelle = 1:size(x.root.cells,1) % for loop over all cells
                cel = x.root.cells(zelle,:);
                
                % get spiking rate of cell in root object for whole
                % recording length (will be cut later to match min_speed
                % and epochs if desired)
                x.root.epoch = [-inf inf];
                spktimes = get_spktimes_of_cel(x.root,cel);
                
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
                quality_check = cell_quality_check_speed_coding(x.root,cel,logical_inds_cut);
                if quality_check
                    fprintf('Cell %d passed quality check. Start analysis.\n',zelle)
                    counter = counter + 1; % counter for all cells
                    
                    % get epochs from logical_inds
                    epochs_inds_for_mle = PeriodsAboveThr(double(logical_inds),1,2);
                    epochs_mle = nan(length(epochs_inds_for_mle),2); % initialize matrix for epochs
                    for i = 1:length(epochs_inds_for_mle)
                        epochs_mle(i,1) = ts(epochs_inds_for_mle{i}(1));
                        epochs_mle(i,2) = ts(epochs_inds_for_mle{i}(end));
                    end
                    % set epcohs
                    x.root.epoch = epochs_mle;
                    x.root.cel = cel; % set cell
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
                    % spikes is 1200 (equivalent to 0.5Hz in a 20-min
                    % recording)
                    % nmbr_spikes = length(spktimes);
                    nmbr_spikes = 1200;
                    
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
                    
                    boxSize = 100;
                    eeg_sample_rate = 600;%250;
                    theta_freq_range = [6 10];
                    
                    %%% set active lfp to the other tetrode channel
                    if cel(1) == 1
                        active_lfp = 2;
                    elseif cel(1) == 2
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
                    % boxSize = length (in cm) of one side of the square box
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
                    Results(zelle).cell_ID = strcat(x.root.name,'_T', num2str(cel(1)),'C', num2str(cel(2)));
                    Results(zelle).speed_response = speed_response; % modeled speed response tuning curve firing rate values
                    Results(zelle).speed_vector = speed_vector;% x-axis for speed tuning curve
                    if ~isnan(selected_model)
                        Results(zelle).best_model = tested_models(selected_model);
                    else
                        Results(zelle).best_model = nan;
                    end                % also save information about position and head direction
                    % tuning
                    Results(zelle).pos_response = pos_response; % for plotting use: imagesc(reshape(pos_response,20,20));
                    Results(zelle).hd_vector = hd_vector; % x-axis for head direction tuning curve
                    Results(zelle).hd_response = hd_response; % modeled head direction response tuning curve firing rate values
                    % read out the tuning curves
                    Results(zelle).pos_curve = pos_curve;
                    Results(zelle).hd_curve = hd_curve;
                    Results(zelle).speed_curve = speed_curve;
                    % read out the scale factors (alpha)
                    Results(zelle).scale_factor_pos = scale_factor_pos;
                    Results(zelle).scale_factor_hd = scale_factor_hd;
                    Results(zelle).scale_factor_spd = scale_factor_spd;
                    % save the spike train
                    Results(zelle).spktrain = spiketrain;
                    Results(zelle).real_spktrain = real_spiketrain;
                    Results(zelle).speed = speed;
                    % finally save log-likelihoods of model performance (m x n
                    % matrix with m = number of cross-validations, and n =
                    % number of models
                    Results(zelle).LLH_values = LLH_values;
                    Results(zelle).tested_models = tested_models;
                    
                    %%%%%%%%%%%%%%%%%%%%%%%
                    %%%%%%%%%%%%%%%%%%%%%%%
                    
                    % define cell_session_animal ID
                    name_export = [name{end-1},'_',name{end}(1:end-4),'_',sprintf('T%dC%d',cel(1),cel(2))];
                    
                    clearvars x.root
                    disp([savename,', Cell ', sprintf('T%dC%d',cel(1),cel(2)), ' analyzed.'])
                end
            end
            % save variables
            save(strcat(pathname,'/',savename),'Results','-v7.3')
            disp(['Results for session ', savename, ' saved.'])
        end
    end
end

