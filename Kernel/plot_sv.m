function plot_sv(I)
%--------------------------------------------------------------------------
% This file is part of the d_RPCA repository from
%           https://github.com/roland1993/d_RPCA
% If you have questions contact me at
%           haase [at] mic.uni-luebeck [dot] de
% Source code is provided under the
%           MIT Open Source License
%--------------------------------------------------------------------------

% get number of images
numImg = size(I{1}, 3);

% get number of outer iterations
outerIter = sum(cellfun(@(c) ~isempty(c), I));

% get numImg colors
cmap = jet(numImg);

% cell-array for legend entries
names = cell(numImg + 1, 1);

% compute singular values per iterate
SV = zeros(numImg, outerIter);
for i = 1 : outerIter
    
    meanI = mean(I{i}, 3);
    TMP = reshape(I{i} - meanI, [], numImg);
    
    [~, D, ~] = svd(TMP, 'econ');
    D = diag(D);
    
    SV(:, i) = D;
    
end

% do the plotting
figure;
hold on;
for i = 1 : numImg
    plot(SV(i, :), '-x', 'Color', cmap(i, :));
    names{i} = ['\sigma_', num2str(i)];
end
plot(sum(SV, 1), 'k--x');
names{numImg + 1} = '\Sigma_i \sigma_i';
hold off;
xlim([0.5, outerIter + 0.5]);
xlabel('#outer iter');
title('development of singular values');
grid on;
legend(names);

end