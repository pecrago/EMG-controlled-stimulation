function [FDI] = load_FDI()
load('FDI_model_1_1', 'N', 'P', 'Tc', 'RET', 'stimorder', 'tic_smu', 'tep_smu', 'tr', 'ipivmin', 'ipivcv');
FDI = struct;
FDI.('N') = N;
FDI.('P') = P;
FDI.('Tc') = Tc;
FDI.('RET') = RET;
FDI.('stimorder') = stimorder;
FDI.('tic') = tic_smu;
FDI.('tep') = tep_smu;
FDI.('tr') = tr;
FDI.('ipivmin') = ipivmin;
FDI.('ipivcv') = ipivcv;
FDI.('gnorm') = (1-exp(-2*(0.4)^3))/0.4; % for MU force gain normalization