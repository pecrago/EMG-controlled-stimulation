function [ smuemg_frames, smumwave_frames ] = EMGbyFrame( tAP_frame, APbystim_frame, smuap, tincr, t_frames )
% Calculates the force and EMG waveforms as functions of time for a single motor unit
% 
% Inputs
%   tAP_frame; the times of the action potentials in this frame
%   APbystim_frame; logical vector true if AP due to stimulation
%   smuap;  The motor unit action potential scaled for this motor unit
%   tincr; sample time
%   t_frames; The total sample time vector to allow for full EMG waveform for all APs that
%           take place in this frame, i.e. it extends beyond one frame
% Outputs
%   smuemg_frames, contains EMG and Mwaves
%   smumwave_frames, only contains Mwaves
lengthsmuap = length(smuap);
totsamp = length(t_frames);
smuemg_frames = zeros(1,totsamp);
smumwave_frames = zeros(1,totsamp);
% calculate individual motor unit emgs
if isempty(t_frames)
    0;
end
if ~isempty(tAP_frame) % calculate for all active motor units
    % convolve smu APs with firing times to construct smu EMG vector
    for j = 1:1:length(tAP_frame) % j indexes the stimulus times
        istart = 1 + floor((tAP_frame(j)-t_frames(1))/tincr);% index to time when AP reaches muscle
        % add scaled smuap
        endindex = min(totsamp, istart + lengthsmuap - 1);
        for i = istart:endindex %i indexes the time since the start of the impulse
            if ((i-istart+1) < 1) || ((i-istart+1) > lengthsmuap)
                0;
            end
            if i >= totsamp
                0;
            end
            smuemg_frames(i) = smuemg_frames(i) + smuap(i-istart+1);
            if APbystim_frame(j) % if AP is due to stimulation, add to mwave
                smumwave_frames(i) = smumwave_frames(i) + smuap(i-istart+1);
            end
        end
    end
end
    

