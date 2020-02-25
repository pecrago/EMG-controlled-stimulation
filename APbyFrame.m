function [ tAP_frame, APbystim_frame, integral_next_APn, threshold_next_APn, tn_open_frame, ts_open_frame ]...
    = APbyFrame(tn_open, R, integral, threshold, ts_open, stimactive, tincr, t_frame, tic, tep, tr, ipivmin, ipivcv)

% ******* assigns SMU IPI variance, but the assignment needs to be moved elsewhere
% these two variables need to be set up front in the muscle model, not here
ipivmin = .005;
% ipivcv = 0.2;
ipivcv = 0.2;
% ******* 

% This function calculates the firing times of action potentials arriving
% at the endpoint of a single axon that result from interactions between a neural
% integrator and a stimulus occuring at an intermediate location along an axon
%
% inputs to function
    % tn_open - times of APn were unresolved in prior frame
    % R - vector of analog values of firing rate of neuron to be encoded
    % integral - the current value of the threshold integral for this frame
    % threshold - calculated firing threshold for the integrator
    % ts_open - times of APs that were unresolved in prior frame
    % stimactive - true if the axon is stimulated in this frame
    % tincr - sampling interval
    % tframe - vector of sample times for this frame
    % tic - intersite conduction time between neural pulse generator and stimulation site
    % tep - conduction time from stimulation site to endpoint
    % tr - refractory period of axon at stimulation site
    %
% ***********  outputs of function
% These two vectors record the times and sources of endpoint action potentials
tAP_frame = []; % vector of times of all new endpoint APs in this frame
APbystim_frame = []; % logical vector of AP sources, T if stimulation, F if neural
% integral_next_APn = 0; % neural integral and threshold at the end of the current frame
% threshold_next_APn = 1;
tn_open_frame = []; % neural APs that are unresolved within this frame
ts_open_frame = []; % unresolved stimulated APs 
% ***********

tn = tn_open; % times of unresolved neural AP
ts = ts_open; % times of unresolved stimulated AP
t_next_frame = max(t_frame) + tincr; % this is the  next possible stimulation time
nsamp = length(t_frame); 

if stimactive % stimulation can only occur at beginning of frame
    ts = [ts t_frame(1)]; % append time of this frame's stim pulse
end

% At this point, we start with a list of existing APs that need to be
% processed
i = 1; % index to tframe
%% version 3A
increment_ts = false;
increment_tn = false;
reset = false;
while i <= nsamp
    integral = integral + R(i)*tincr; % this is a basic integrate and fire model
    if integral >= threshold
        delta_int = integral - threshold; % delta_int is how much the sum exceeds the threshold
        delta_t = delta_int/R(i) ; % interpolate to get better resolution
        integral = delta_int; % reset integrator, but include the excess
        tn = [tn, t_frame(i) - delta_t]; % this is the neural AP time
        threshold = max(ipivmin,1+ipivcv*randn); 
    end
    i = i + 1;
    % start with stimulated APs case
    if ~isempty(ts)
        if ~isempty(tn)
            if ts(1) >= tn(1)+tic 
                if ts(1) < tn(1)+tic+tr % refractory block
                    increment_ts = true;
                end
                % refractory block or not
                increment_tn = true;
                if tn(1)+tic+tr < t_next_frame
                    tAP_frame = [tAP_frame tn(1)+tic+tep];
                    % NOTE: tn(1)+tic+tep will be less than tn(1)+tic+tr if tep<tr
                    APbystim_frame = [APbystim_frame false];
                else
                    tn_open_frame = [tn_open_frame tn(1)];
                end
            elseif ts(1)+tic < tn(1) % reset
                increment_ts = true;
                if ts(1)+tic < t_next_frame
                    tAP_frame = [tAP_frame ts+tep];
                    APbystim_frame = [APbystim_frame true];
                    reset = true;
                    tn_reset = ts(1) + tic; % 
                else
                    ts_open_frame = [ts_open_frame ts(1)];
                end
            else  % collision
                increment_tn = true; % voluntary AP accounted for
                increment_ts = true; % stimulated AP accounted for
                if (ts(1)+tn(1))/2 < t_next_frame
                    tAP_frame = [tAP_frame ts + tep]; %record AP activation by stimulation
                    APbystim_frame = [APbystim_frame true];
                else
                    ts_open_frame = [ts_open_frame ts(1)];
                end
            end
        else % ts but no tn
            if ts(1)+tic >= t_next_frame
                ts_open_frame = [ts_open_frame ts(1)];
                increment_ts = true;
            elseif ts(1)+tic < t_frame(i) % reset
                increment_ts = true;
                tAP_frame = [tAP_frame ts+tep];
                APbystim_frame = [APbystim_frame true];
                reset = true;
                tn_reset = ts(1) + tic; % 
            else % continue to integrate
            end
        end
    else % empty ts
        if ~isempty(tn)
            increment_tn = true;
            if tn(1)+tic+tr < t_next_frame
                tAP_frame = [tAP_frame tn(1)+tic+tep];
                APbystim_frame = [APbystim_frame false];
            else
                tn_open_frame = [tn_open_frame tn(1)];
            end
        end
    end
    if increment_ts
        ts(1) = [];
        increment_ts = false;
    end
    if increment_tn
        increment_tn = false;
        tn(1) = [];
        if reset % figure out what the new index i should be, and the new value of each integration
            i = 1 + ceil((tn_reset-t_frame(1))/tincr); % this is the next value of i for the integration
            delta_t = (i-1)*tincr + t_frame(1) - tn_reset;
            integral = R(i-1)*delta_t; % this is what the integral would be
            reset = false;
        end
    end
end  % while i<=nsamp
%% carry these values over to next frame
integral_next_APn = integral; 
threshold_next_APn = threshold;

if ~isempty(tAP_frame)
      if tAP_frame(1) < t_frame(1)
        0;
      end
end


