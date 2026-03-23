%DEMO_BSPLINES  One-line description here, please.
%
%   output = demo_BSplines(input)
%
%   Example
%   demo_BSplines
%
%   See also
%
 
% ------
% Author: David Legland
% e-mail: david.legland@inra.fr
% Created: 2018-09-25,    using Matlab 8.6.0.267246 (R2015b)
% Copyright 2018 INRA - Cepia Software Platform.

% base de temps
t = linspace(0, 1, 50);


%% Basis functions

% new figure
figure; hold on;

% draw he basis functions
plot(t, BSplines.beta3_0(t), 'LineWidth', 2);
plot(t, BSplines.beta3_1(t), 'LineWidth', 2);
plot(t, BSplines.beta3_2(t), 'LineWidth', 2);
plot(t, BSplines.beta3_3(t), 'LineWidth', 2);

legend({'\beta_0', '\beta_1', '\beta_2', '\beta_3'});
title('Basis functions');
print(gcf, '-r100', 'bsplines_basisFunction.png', '-dpng');
% print(gcf, '-vector', 'bsplines_basisFunction.eps', '-depsc2');


%% Derivatives

% new figure
figure; hold on;

% draw he basis functions
plot(t, BSplines.beta3_0d(t), 'LineWidth', 2);
plot(t, BSplines.beta3_1d(t), 'LineWidth', 2);
plot(t, BSplines.beta3_2d(t), 'LineWidth', 2);
plot(t, BSplines.beta3_3d(t), 'LineWidth', 2);

legend({'\beta_0', '\beta_1', '\beta_2', '\beta_3'});
title('Derivatives');
print(gcf, '-r100', 'bsplines_derivatives.png', '-dpng');


%% Second Derivatives

% new figure
figure; hold on;

% draw he basis functions
plot(t, BSplines.beta3_0s(t), 'LineWidth', 2);
plot(t, BSplines.beta3_1s(t), 'LineWidth', 2);
plot(t, BSplines.beta3_2s(t), 'LineWidth', 2);
plot(t, BSplines.beta3_3s(t), 'LineWidth', 2);

legend({'\beta_0', '\beta_1', '\beta_2', '\beta_3'});
title('Second Derivatives');
print(gcf, '-r100', 'bsplines_secondDerivatives.png', '-dpng');



