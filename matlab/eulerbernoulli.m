% Løser Euler-Bernoullibjelkelikningen numerisk
% Input:
%       E = Bjelkens Youngmodulus [Pascal]
%       D = materialets massetetthet [kg/m^3]
%       w = bjelkens bredde [meter]
%       d = bjelkens tykkelse [meter]
%       L = bjelkens lengde [meter]
%       n = antall oppstykkinger
%       b = lastvektor (uten egenlast) [Newton]
% Output:
%       y = vektor som tilfredsstiller Euler-Bernoullibjelkeligningen

function y = eulerbernoulli(E, D, w, d, L, n, b)
g = -9.81;                              % tyngdekraftens akselerasjon
I = w*d^3/12;                           % treghetsmoment
h = L/n;                                % segmentlengde
if nargin < 7                           % hvis b ikke gitt
    b = zeros(n, 1);      
end
b = b + repmat(D*w*d*-9.81, n, 1);      % legger til bjelkens egenlast
b = b * (h^4 / (E*I));
A = lagA(n);
y = A\b;
end

