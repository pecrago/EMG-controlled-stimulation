function [y_out, m_out] = gsfm(x, M)
% general Gram-Schmidt Prediction Error Filter
% adapted from Nathan Makowsky's 'SimpleGSFilt.m'
% x is M+1 frames of row vectors of length n_f
% frames are in reverse temporal order - current frame is index 1
[~, n_per_frame] = size(x);
w = zeros(M);
eps = zeros(M+1, M+1, n_per_frame);
for i = 1:M+1
    eps(i, 1, :) = x(i,:);
end
for m = 1:M % step m = m_step-1 through 0:M-1
    for i = 1:1:(M-m+1) % step i = i_step-1 through 0:M-m-1
%         eps(i_step, m_step, :) = x(m_step,i_step,:);
        temp1(1, :) = eps(i, m, :);
        temp2(:, 1) = eps(M-m+2, m, :);
        w(m, i) = temp1*temp2/(mag(eps(M-m+2, m, :)))^2;
        eps(i, m+1, :) = eps(i, m, :) - w(m, i)*eps(M-m+2, m, :);
    end
end
y_out(1, :) = eps(1, M+1, :);
m_out(1, :) = x(1, :) - y_out(1, :);
end % of function gsfm

function y = mag(X)
y = sqrt(sum((X.^2)));
end
