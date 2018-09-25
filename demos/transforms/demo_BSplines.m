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

% new functions
figure; hold on;

% draw he basis functions
plot(t, BSplines.beta3_0(t));
plot(t, BSplines.beta3_1(t));
plot(t, BSplines.beta3_2(t));
plot(t, BSplines.beta3_3(t));

legend({'\beta_0', '\beta_1', '\beta_2', '\beta_3'});


%% draw the kernel
figure; hold on;
plot(t-2, BSplines.beta3_3(t), 'b');
plot(t-1, BSplines.beta3_2(t), 'b');
plot(t, BSplines.beta3_1(t), 'b');
plot(t+1, BSplines.beta3_0(t), 'b');


%% Draw the first derivatives

figure; hold on;
plot(t-2, BSplines.beta3_3d(t), 'b');
plot(t-1, BSplines.beta3_2d(t), 'b');
plot(t, BSplines.beta3_1d(t), 'b');
plot(t+1, BSplines.beta3_0d(t), 'b');


%% Draw the second derivatives

figure; hold on;
plot(t-2, BSplines.beta3_3s(t), 'b');
plot(t-1, BSplines.beta3_2s(t), 'b');
plot(t, BSplines.beta3_1s(t), 'b');
plot(t+1, BSplines.beta3_0s(t), 'b');
