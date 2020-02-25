function [F, EMG, Mwave, Emav, Mest, tAP, APbystim] =...
    rampVol_arbStim (tincr, Uvolramp, n_frames, samples_per_frame, t_stim, Musc, smuap, plotlevel, datasave)
% modified from MNP_sim_01.m to generate ramp vol contractions, and
% modulaed stim contractions based on vol EMG

% model of force produced by motor neuron pool based on models published by ...
%  Fuglevand Winter & Patla 1993; Zhou and Rymer 2007;
% tincr is the sample interval for continuous signals, e.g., F, EMG
% Uvol contains linear ramp parameters
% n_frames is the integer number of whole frames from t=0 to t=tstop
% t_stim is the vector of possible stim times
% stimorder(1:N) is user specified order of motor unit stimulation recruitment
% P(1:N) is the peak twitch force of each motor unit
% Tc(1:N) is the contraction time of each motor unit
% RET(1:N) is the recruitment excitation threshold of each motor unit
% smuap(1:N,:) is the array of single motor unit AP waveforms 
% plotlevel (logical) controls whether computed results are plotted 
% datasave (logical) controls whether results are saved

% Ustim(1:nsamp) is the fraction of motor units recruited by the electrical stimulation
% freqstim (global variable) is the fixed stimulus rate
%
%% simulation parameters
global freqstim;
global tvolstart;
global tstop;
global nsamp;
GSF_order = 2; % order of Graham Schmidt MAV estimator
t = 0:tincr:tstop; % vector of sample times
nsamp = length(t);
[~, aplength] = size(smuap); % number of time samples per SMU AP
twitchlength = ceil(5*Musc.('Tc')(1)/tincr); % number of time samples for 5 Tc
F = zeros(1,nsamp);% initialize total force vector to zeros
EMG = zeros(1,nsamp); % initialize total EMG vector to zeros
EMG_noise = zeros(1, nsamp);
Mwave = zeros(1,nsamp); % initialize total Mwave vector to zeros
Mest = zeros(1,nsamp); % initialize estimated Mwave vector to zeros
EMG_est = zeros(1,nsamp);
Emav = zeros(1,n_frames); % initialize mean absolute EMG vector to zeros
Emav_est = zeros(1,n_frames);
tAP = cell(Musc.('N'), 1); % this will store the firing times of smus

% force and emg summation arrays and scaling relationships
smuforce = zeros(Musc.('N'), nsamp); % initialize single motor unit force summing arrays to zeros
smuemg = zeros(Musc.('N'), nsamp); %initialize single motor unit EMG summing array to zeros
smumwave = zeros(Musc.('N'), nsamp); % initialize single motor unit Mwave summing array to zeros
EMG_frame = zeros(n_frames,samples_per_frame);
%% calculate R(k,t) for neurally activated units
% calculate MN rates for excitation ramp
% note: this only calculates the values of the voluntary rates
volmin = Uvolramp(1);
excit_slope = Uvolramp(2);
tvolstart = Uvolramp(3);
tvolstop = Uvolramp(4);
graph = false;
ivolstart = round(tvolstart/tincr,0);
ivolstop = round(tvolstop/tincr,0);
U_vol(1:ivolstart) = volmin;
for i = ivolstart+1:ivolstop
    U_vol(i) = volmin + excit_slope*(i-ivolstart)/(ivolstop-ivolstart)*(tvolstop-tvolstart);
end
U_vol(ivolstop:nsamp) = U_vol(ivolstop);

[~, R] = Contessa_ramp(volmin, excit_slope, tvolstart, tvolstop, tincr, tstop, graph);

%% calculate force vs. time and EMG vs. time vectors for each active motor unit
%   Time frames are calculated sequentially
% start by creating vectors of times for motor unit action potentials
% then calculate force responses to excitations
gnorm = (1-exp(-2*(0.4)^3))/0.4; % for MU force gain normalization

% save individual smu data in cell arrays
tAP = cell(1,Musc.('N'));
APbystim = cell(1,Musc.('N'));
integral = zeros(1,Musc.('N')); % for integrate and fire model of each neuron
threshold = zeros(1,Musc.('N'));
threshold(1:Musc.('N')) = 1; % thresholds for integrate and fire model
tstartframe = zeros(1,n_frames);
istartframe = zeros(1,n_frames);
tendframe = zeros(1,n_frames);
iendframe = zeros(1,n_frames);
U_stim = zeros(1, n_frames);
tn_open = cell(1,Musc.('N'));
ts_open = cell(1,Musc.('N'));
tAP_prior = zeros(1,Musc.('N'));
EMG_filter_frames = zeros(GSF_order+1, samples_per_frame);
EMG_frame_est = zeros(n_frames, samples_per_frame);
m_frame_est = zeros(n_frames, samples_per_frame);
for j = 1:n_frames
    istartframe(j) = (j-1)*samples_per_frame + 1;
    iendframe(j) = istartframe(j) + samples_per_frame - 1;
    tstartframe(j) = t_stim(j); % t at start of frame
    tendframe(j) = t_stim(j+1) - tincr;
end
% Create a noise background as a fraction of MAV MVC EMG
MAV_MVC_EMG = 1.217e-4; %1.217017696096749e-04
noise_to_signal_ratio = .1;
low_cutoff = 15; % Hz
high_cutoff = 400; % Hz
noise = signalnoise(nsamp, 1/tincr, low_cutoff, high_cutoff);
MAV_noise = mean(abs(noise));
EMG_noise = noise_to_signal_ratio*MAV_MVC_EMG*noise'/MAV_noise;
for j = 1:n_frames-1
% ******* DEBUGGING
% U_stim(j) = 0; % to turn off stimulation, switch with next if ... end
    if j > GSF_order+1 && ~isnan(Emav_est(j-1))
            U_stim(j) = stimlevel(Emav_est(j-1));
    end
    % fill logical vector of which axons are stimulated
    t_frame = tstartframe(j):tincr:tendframe(j); 
    stimactive = zeros(1,Musc.('N'));
    N_stim = floor(max(U_stim(j))*Musc.('N'));
    if N_stim >0
        stimactive(Musc.('stimorder')(1:N_stim)) = true;
    end
    for k=1:Musc.('N')  % calculate the times and sources of APs in all frames
        R_frame = R(k, istartframe(j):iendframe(j));
        [ tAP_frame, APbystim_frame, integral_next_APn, threshold_next_APn, tn_open_frame, ts_open_frame ]...
            = APbyFrame(tn_open{k}, R_frame, integral(k), threshold(k), ts_open{k}, stimactive(k), tincr, t_frame, ...
            Musc.('tic')(k), Musc.('tep')(k), Musc.('tr'), Musc.('ipivmin'), Musc.('ipivcv') );
        % append new APs in this frame
        tAP{k} = [tAP{k} tAP_frame]; 
        APbystim{k} = [APbystim{k} APbystim_frame]; 
        integral(k) = integral_next_APn;
        threshold(k) = threshold_next_APn;
        tn_open{k} = tn_open_frame;
        ts_open{k} = ts_open_frame;
        % calculate individual motor unit EMGs and add to EMG and Mwave totals
        if ~isempty(tAP_frame)
            N_prior_frames = 0;
            if (tAP_frame(1) < tstartframe(j))
                N_prior_frames = ceil((tstartframe(j)-tAP_frame(1))/samples_per_frame);
            end
            N_SMUAP_frames = 1 + ceil(aplength/samples_per_frame); % +1 since tAP could occur at end of frame
            t_SMUAP_frames = tstartframe(j-N_prior_frames):tincr:tendframe(min(j+1+N_SMUAP_frames, n_frames));
            [ smuemg_frames, smumwave_frames ] = EMGbyFrame( tAP_frame, APbystim_frame,...
                smuap(k,:), tincr, t_SMUAP_frames );
            EMG_indexrange = istartframe(j-N_prior_frames):(iendframe(j+N_SMUAP_frames-1));
            smuemg(k,EMG_indexrange) = smuemg(k,EMG_indexrange) + smuemg_frames(1:length(EMG_indexrange));
            smumwave(k,EMG_indexrange) = smumwave(k,EMG_indexrange) + smumwave_frames(1:length(EMG_indexrange));
            % sum voluntary EMG level in each frame
            EMG_frame(j,:) = EMG_frame(j,:) + smuemg(k,istartframe(j):iendframe(j));
        end
        % calculate individual motor unit forces and add to force totals
        if ~isempty(tAP_frame)
            N_prior_frames = 0;
            if (tAP_frame(1) < tstartframe(j))
                N_prior_frames = ceil((tstartframe(j)-tAP_frame(1))/samples_per_frame);
            end
            N_SMUforce_frames = 1 + ceil(twitchlength/samples_per_frame); % +1 since tAP could occur at end of frame
            t_SMUforce_frames = tstartframe(j-N_prior_frames):tincr:tendframe(min(j+1+N_SMUforce_frames, n_frames));
            [ SMUforce_frames, tAP_last ] = ForcebyFrame( tAP_frame, tAP_prior(k), Musc.('P')(k), Musc.('Tc')(k), Musc.('gnorm'), ...
                tincr, t_SMUforce_frames );
            tAP_prior(k) = tAP_last;
            F_indexrange = istartframe(j-N_prior_frames):(iendframe(min(n_frames,j+N_SMUforce_frames-1)));
            smuforce(k,F_indexrange) = smuforce(k,F_indexrange) + SMUforce_frames(1:length(F_indexrange));
        end
    end  % of stepping through motor units
    % add noise to EMG frame, also contains Mwaves
    EMG_frame(j,:) = EMG_frame(j,:) + EMG_noise(istartframe(j):iendframe(j));
    % estimate voluntary MAV from total EMG for each frame
    EMG_filter_frames(1+GSF_order, :) = []; % delete oldest data frame
    EMG_filter_frames = [EMG_frame(j,:); EMG_filter_frames]; % insert newest data frame
    if j > GSF_order
        [EMG_frame_est(j,:), m_frame_est(j,:)] = gsfm(EMG_filter_frames, GSF_order);
% % ********* Uncomment for debugging
%         figure;hold on;
%         for i = 1:GSF_order+1
%             subplot(GSF_order+1,1,i);
%             plot(EMG_filter_frames(i,:));
%         end
%         subplot(GSF_order+1,1,1);hold on;
%         plot(EMG_frame_est(j,:));
% % *********
    end
    EMG_est(istartframe(j):iendframe(j)) = EMG_frame_est(j,:);
    Emav_est(j) = mean(abs(EMG_est(istartframe(j):iendframe(j))));
end % of stepping through frames
%% calculate total force and emg produced by the whole motor neuron pool
% if A is a matrix, then sum(A) returns a row vector containing the sum of each column
F = sum(smuforce);
EMG = sum(smuemg) + EMG_noise;
Mwaves = sum(smumwave);
EMG_vol = EMG - Mwaves;
rectvolemg = abs(EMG_vol);
rectestemg = abs(EMG_est);
for j = 1:n_frames - 1
    Emav(j) = mean(rectvolemg(istartframe(j):iendframe(j)));
end
%% plots for debugging purposes
tvector = tstartframe(1):tincr:tendframe(n_frames-1);
indexvector = 1:length(tvector);
if plotlevel >= 3
    for k = 10:10:120
        figure; 
        subplot(2,1,1); 
        hold on;
        stimAPs = logical(APbystim{1,k});
        instfrq = [NaN 1./diff(tAP{1,k})];
        p1 = plot(tAP{1,k}(:),instfrq,'bo');
        p2 = plot(tvector, R(k, istartframe(1):iendframe(n_frames-1)), 'b-');
        p3 = plot(tAP{1,k}(stimAPs),instfrq(stimAPs),'ro');
        title(sprintf('SMU #%d, Rs=%d',k,freqstim));
        ylabel('instantaneous frequency, s^{-1}');
        for i = 1:n_frames-1
            line([tstartframe(i), tstartframe(i)], [0 5]);
        end
        legend([p1 p2 p3],{'neuralAP','neural rate','stimulatedAP'},'Location', 'NorthWest');
        legend('boxoff');
        subplot(2,1,2);
        hold on;
        plot(tvector, smuforce(k,indexvector),'b');
        ylabel('F_{SMU} [au]');
        xlabel('time, s');
    end
end
if plotlevel >= 2
    figure;
    subplot(3,1,1);hold on;
    plot(tstartframe(1):tincr:tendframe(n_frames-1), EMG_vol(istartframe(1):iendframe(n_frames-1)), 'b');
    ylabel('EMG_{vol}');
    subplot(3,1,2);hold on;
    plot(tstartframe(1):tincr:tendframe(n_frames-1), Mwaves(istartframe(1):iendframe(n_frames-1)), 'b');
    ylabel('Mwaves');
    subplot(3,1,3);hold on;
    plot(tstartframe(1):tincr:tendframe(n_frames-1), EMG_est(istartframe(1):iendframe(n_frames-1)), 'r');
    ylabel('EMG_{est}');
    
    figure;
    subplot(4,1,1);
    hold on;
    p4 = plot(tvector, F(1:length(tvector)), 'b-');
    ylabel('F_{total} [au]');

    subplot(4,1,2);
    hold on;
    p5 = plot(tvector, U_vol(1:length(tvector)),'b');
    p6 = plot(tstartframe, U_stim, 'r');
    legend([p5 p6], {sprintf('U_{vol}'), sprintf('U_{stim}')},'Location', 'NorthWest');
    legend('boxoff');
    ylabel(sprintf('U_{vol}, U_{stim}'));

    subplot(4,1,3);
    hold on;
    p7 = plot(tvector, rectvolemg(indexvector),'b'); 
    title(sprintf('M=%d',GSF_order));
    p8 = line([tstartframe(1) tendframe(1)], [Emav(1) Emav(1)], 'LineWidth', 2, 'Color', 'r');
    p9 = line([tstartframe(1) tendframe(1)], [Emav_est(1) Emav_est(1)], 'LineWidth', 2, 'Color', 'g');
    for j = 2:n_frames-1
        line([tstartframe(j) tendframe(j)], [Emav(j) Emav(j)], 'LineWidth', 2, 'Color', 'r');
        line([tstartframe(j) tendframe(j)], [Emav_est(j) Emav_est(j)], 'LineWidth', 2, 'Color', 'g');
    end
    legend([p7 p8 p9], {'|EMG_{vol}|', 'MAV', 'MAV_{est}'},'Location', 'NorthWest');
    legend('boxoff');
    ylabel('|EMG_{vol}| , MAV [au]');

    subplot(4,1,4);hold on;
    plot(tvector, Mwaves(indexvector),'k');
    ylabel('Mwaves [au]');
    xlabel('time, s');

    figure;
    hold on;
    plot(Emav, Emav_est,'ko');
    axis([0 .0002 0 .0002]);
    line([0 .0002],[0 .0002]);
    xlabel('E_{MAV}');
    ylabel('E_{MAVest}');
    title(sprintf('M=%d',GSF_order));
    % fit a straight line to the data
    %************debug
    % f=fit(Emav(5:240)', Emav_est(5:240)', 'poly1')
    %     % plot muscle force and emg responses
    %     s = sprintf('%d stimulated, %d voluntarily activated', N_stim, Nneur); %
    %     plothandle(2) = figure;
    %     hold;
    %     subplot(5,1,1)
    %     plot(t,F);
    %     ylabel('MNP force');
    %     title(s);
    %     subplot(5,1,2)
    %     plot(t,EMG);
    %     ylabel('MNP EMG');
    %     subplot(5,1,3),plot(t,mnpmwaves);
    %     ylabel('Mwaves');
    %     subplot(5,1,4),plot(t,volemg);
    %     ylabel('voluntary EMG');
    %     subplot(5,1,5),plot(t,rectvolemg);
    %     ylabel('Rectified Voluntary EMG');
    %     xlabel('time, s');
end
%% save MAT file for later plotting
% how to ID the files?,series id, max N&S levels, ramp duration, freqstim?
% Volmax = 100*max(Uvolramp);
% Stimmax = 100*max(U_stim);
if datasave
    filename = sprintf('EMG_controlled_stim_%d', freqstim);
    save(filename);
end
