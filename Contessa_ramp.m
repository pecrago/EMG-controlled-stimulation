function [neuractive, rate] = Contessa_ramp(volmin, excit_slope, ...
    tvolstart, tvolstop, deltat, tstop, graph)
% load the following variables from previously defined model
% N, RR, A,B, C, D, E, mr, br, mtheta, btheta, ipivmin, ipivcv, 
load('FDI_model_1_1.mat','N', 'RR', 'A','B', 'C', 'D', 'E', 'mr', 'br', ...
    'mtheta', 'btheta', 'ipivmin', 'ipivcv');

% output arguments
% neuractive, rate

% input arguments
% (volmin, excit_slope, tvolstart, tvolstop, deltat, tstop, graph) 

lambdar = zeros(1,N); %threshold firing rates
theta = zeros(1,N); %time constant
RET = zeros(1,N); %recruitment excitation threshold (%) of each unit
RETnorm = zeros(1,N); % fraction of max excitation
i = 1:N;    
a=log(RR)/N;
RET = exp(a.*i); %from Fuglevand 1993
RETnorm = RET/100;
exc = [RET ceil(RET(N)):1.0:100];
excnorm = exc/100;
% steps = length(exc); % number of excitation levels in exc

mp = (C - A*exp(-(excnorm/B)));
bp = D*excnorm + E;
theta = mtheta*RETnorm + btheta; % eqn 5 of D&C 2012
lambdar = mr*RETnorm + br; % eqn 3 of D&C 2012

timespan = 0:deltat:tstop;
nsamp = length(timespan);
rate = zeros(N,nsamp);
ivolstart = floor(tvolstart/deltat);
for k = 1:N
    neuractive(k) = false;
    j = 1;
    while ((~neuractive(k)) && (j <= nsamp))
        if j <= ivolstart
            excit(j) = volmin;
        else
            excit(j) = min(1,volmin+(j-1-ivolstart)*deltat*excit_slope);
        end
        if excit(j) - RETnorm(k)>0
            neuractive(k) = true;
            t_active = j*deltat;
            while j <= nsamp
                excit(j) = min(1,volmin+(j-1-ivolstart)*deltat*excit_slope);
                mpr(j) = (C - A*exp(-(excit(j)/B)));
                bpr(j) = D*excit(j) + E;
                lambdapr(j)= mpr(j)*RETnorm(k) + bpr(j);
                rate(k,j) = lambdar(k) + (lambdapr(j)-lambdar(k))*(1 - exp((t_active-j*deltat)/theta(k)));
                j = j + 1;
            end
        else
            rate(k,j) = 0;
            j = j+1;
        end
    end
end
%%
if graph
    figure; subplot(2,1,1);hold;
    for k = 1:N
        plot(timespan,rate(k,:));
    end
    titlestring = sprintf('%d%%/s ramp excitation',excit_slope*100);
    title(titlestring);
    ylabel('firing rate [1/s]');
    subplot(2,1,2);hold;
    plot(timespan,excit(:));
    axis([0 tstop 0  1]);
    ylabel('excitation');
    xlabel('time [s]');
end
