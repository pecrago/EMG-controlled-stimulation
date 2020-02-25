function [ smuforce_frames, tAP_last ] = ForcebyFrame( tAP_frame, tAP_prior, P, Tauc, gnorm, tincr, t_frames )
% Calculates the force function of time for a single motor unit
% Inputs
%   tAP_frame
%   P
%   Tauc
%   gnorm
%   tincr
%   t_frames
% Outputs
%   smuforce
%   tAP_last
totsamp = length(t_frames);
smuforce_frames = zeros(1,totsamp);
    % calculate individual motor unit forces and emgs
if ~isempty(tAP_frame) % calculate for all active motor units
    for j = 1:length(tAP_frame)        % convolve twitch forces with firing times to construct smu force vector
        % first calculate twitch scale
        if tAP_prior == 0
            g = 1;
        else
            Tnorm = Tauc/(tAP_frame(j)-tAP_prior);
            tAP_prior = tAP_frame(j);
            if Tnorm < 0.4
                g = 1; % linear scaling region for long periods
            else % nonlinear scaling region for short periods
                g = ((1-exp(-2*(Tnorm)^3))/Tnorm)/gnorm;
            end
        end
        % add scaled twitch
        istart = floor((tAP_frame(j)-t_frames(1))/tincr)+1;% tAP when AP reaches muscle
        scale = g*P/Tauc;
        for i = istart:totsamp % i indexes the time since the start of the impulse
            smuforce_frames(i) = smuforce_frames(i)+scale*(t_frames(i) - tAP_frame(j))*...
                exp(1 - (t_frames(i) - tAP_frame(j))/Tauc);
        end
        tAP_last = tAP_frame(j); % save for next frame
    end
end

