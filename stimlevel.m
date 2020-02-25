function [stim_next] = stimlevel(MAV_est)
MAV_max = .0001;
% MAV_max = .00015;
stim_next = min(1,max(0, MAV_est/MAV_max));