% run_EMG_control_sim.m
clear;
% Set Model and Simulation Parameters
global freqstim
freqstim = 20; % /s
datasave = true;
plotlevel = 3; % 1-none, 2-analysis plots, 3-all plots
global nsamp;
% contraction timing parameters
global tvolstart;
tvolstart = 1; % for 10 pct/sec ramps
global tvolstop;
tvolstop = 2;
global tstop
tstop = tvolstop+1;

stimperiod = 1/freqstim; % this is also the frame length, s
% calculate sampling parameters
R_sample_target = 2500; % /s
samples_per_frame = round(R_sample_target/freqstim,0);
R_sample = samples_per_frame*freqstim;
tincr = 1/R_sample;
n_frames = floor(tstop*freqstim);
nsamp = samples_per_frame*n_frames;
t_stim = 0:stimperiod:tstop; % this is a vector of the stim times
nstim = floor(tstop/stimperiod); % +1 to include t=0

% load previously saved muscle model parameters for FDI
[FDI] = load_FDI();

%% generate array (smuap) of scaled single motor unit action potentials
muapduration = 0.012; % s, duration of muaps
lambda = muapduration/5; % value 5 determined by trial and error to get correct duration
tsteps = -muapduration/2:tincr:muapduration/2; % s, time steps for calculating the normalized muap
[~, aplength] = size(tsteps);
smuap = zeros(FDI.('N'),aplength);
smuapnorm = tsteps.*exp(-(tsteps.^2)/lambda^2);
for i=1:FDI.('N')
    smuap(i,:) = FDI.('P')(i)*smuapnorm(:);
end

%% combined vol and modulated stim

% voluntary ramp
volmin = 0;
excit_slope = 1/(tvolstop-tvolstart) ; % 0.1  for 10 pct/sec
u_vol = [volmin, excit_slope, tvolstart, tvolstop];
tic
[F, EMG, Mwave, Emav, Mest, tAP, APbystim] =...
    rampVol_arbStim(tincr, u_vol, n_frames, samples_per_frame, t_stim, FDI, smuap, plotlevel, datasave);
% toc
