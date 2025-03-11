function S = waveSpectrum(spectrumType, Parameter, W, PlotFlag)
% Function computing state-of-the-art wave spectra.
%
% Usage:
%  S = waveSpectrum(spectrumType, Parameter, W, PlotFlag)
%
% Inputs:
%  spectrumType - Spectrum type identifier
%  Parameter    - Spectrum parameters
%  W            - Column vector of wave frequencies [rad/s]
%  PlotFlag     - 1 to plot the spectrum, 0 for no plot
%
% Output:
%  S - Column vector of wave spectrum values [m^2 s]; evaluated at W(k)
%
% Spectrum Type and Parameters:
%  1. Bretschneider (p1 = A, p2 = B)
%  2. Pierson-Moskowitz (p1 = Vwind20)
%  3. ITTC-Modified Pierson-Moskowitz (p1 = Hs, p2 = T0)
%  4. ITTC-Modified Pierson-Moskowitz (p1 = Hs, p2 = T1)
%  5. ITTC-Modified Pierson-Moskowitz (p1 = Hs, p2 = Tz)
%  6. JONSWAP (p1=Vwind10, p2 = Fetch)
%  7. JONSWAP (p1 = Hs, p2 = w0, p3 = gamma)
%  8. Torsethaugen (p1 = Hs, p2 = w0)
%
% 1. Bretschneither [6]
%   p1 = A
%   p2 = B
%   S(w) = A * w^(-5) * exp(-B / w^4)
%
% 2. Pierson-Moskowitz [5,6]
%   p1 = Vwind20 - Average wind speed 20m above sea level (m/s)
%   A = 8.1 ecp(-3 * g^2)  
%   B = 0.74 * (g / Vwind20)^4
%   S(w) = A*w^(-5) * exp(-B / w^4) 
%
% 3. ITTC-Modified Pierson-Moskowitz [1]
%   p1 = Hs - Significant wave height (Hs = 4 m0^1/2) (m)
%   p2 = T0 - Modal period      (T0 = 2 pi /w0)       (s)
%   A = 487 * Hs^2 / T0^4;
%   B = 1949 / T0^4;
%   S(w) = A * w^(-5) * exp(-B / w^4) 
%
% 4. ITTC-Modified Pierson-Moskowitz [1, 5]
%   p1 = Hs - Significant wave height (Hs = 4 m0^1/2) (m)
%   p2 = T1 - Average wave period (T1 = 2 pi * m0/m1) (s)
%   A = 173 * Hs^2 / T1^4
%   B = 691 / T1^4
%   S(w) = A * w^(-5) * exp(-B / w^4)
%
% 5. ITTC-Modified Pierson-Moskowitz [1, 5]
%   p1 = Hs - Significant wave height (Hs = 4 m0^1/2) (m)
%   p2 = Tz - Average zero-crossing period (T = 2 pi (m0/m1)^1/2) (s)
%   A = 123 * Hs^2 / Tz^4
%   B = 495 / Tz^4
%   S(w) = A * w^(-5) * exp(-B / w^4)
%
% 6. JONSWAP [2]
%   p1 = Vwind10 - Wind speed 10m over sea surface (m/s)
%   p2 = fetch   - Distance to geographical boundary, typically > 1000 (m)
%   xtilde = g * fetch / Vwind10^2
%   f0 = 3.5 * (g/Vwind10) * xtilde^(-0.33)
%   w0 = 2 * pi * f0
%   alpha = 0.076 * xtilde^(-0.22)
%   gamma = 3.3
%   sigma = 0.07  if w < w0, sigma=0.09 otherwise
%   S(w) = S1 * S2   
% where
%   S1 = alpha * g^2 * (W^-5) * exp(-(5/4) * (w0/w)^4)
%   S2 = gamma^(exp(-(w-w0)^2/(2*(sigma*w0)^2)))
%
% JONSWAP with DNV coorection [3]
%   p1 = Hs - Significant wave height (Hs = 4 sqrt(m0)) [m]
%   p2 = w0 - Modal Freq. [rad/sec] 
%             Reccomended: 1.25 < w0 * sqrt(Hs) < 1.75
%   p3 = gamma - Peakedness factor (Between 1 and 5, usually 3.3)
%               Set to zero to use DNV formula
%   alpha = (5/16) * Hs^2 * w0^4 / g^2
%   sigma = 0.07  if w < w0, sigma=0.09 otherwise
%   S(w) = S1 * S2   
% where
%    S1 = alpha * g^2 * W^(-5) * exp(-(5/4) * (w0/w)^4)
%    S2 = gamma^(exp( -(w-w0)^2 / (2 * (sigma * w0)^2)) ) 
%
% 8. Torsethaugen [4]
% The Torsethaugen spectrum is an empirical two peaked spectrum for swell 
% and developing sea based on experimental data from the North Sea. For 
% small peak frequencies, i.e.  0 < w0 <= 0.6 only one peak in the 
% spectrum appears. Returns the spectral density function S of the 
% Torsethaugen spectrum for the frequencies in the vector  W [rad/s].
%
%    p1 = Hs  - significant wave height (m)
%    p2 = w0  - Modal (peak) frequency (rad/s)
%    
%    S = waveSpectrum(8, [Hs, wo], W, PlotFlag) calls the function
%    torsetSpectrum(Hs, wo, W), which computes S.
%
% References: 
%
% [1] A.R.M.J LLoyd (1998) "Seakeeping: Ship Behaviour in Rough Wheather."
%     Published by A.R.M.J LLoyd, Gosport UK.ISBN 0-9532634-01
% [2] Ochi, M.K. (1998) "Ocean Waves, The stochastic Approach"
%     Cambridge Ocean Technology Series, Vol 6,
%     Cambridge University Press.
% [3] Sørensen, A.J. (2005) "Marine Cybernetics: Modelling and Control"
%     Lecture Notes for TMR4241 Marine Control Systems, NTNU, 2005.
% [4] K. Torsethaugen (1996): "Model for a Doubly Peaked Wave Spectra"
%      Sintef report no.: STF22 A96204 prepared for Norsk Hydro.
% [5] T.I. Fossen (2021) "Handbook of Marine Craft Hydrodynamics and Motion
%     Control. John Wiley & Sons Ltd., Chichester, UK.  
% [6] Lewis E.V. "Principles of Naval Architecture volume  III
%     Motions in Waves and Controllability." SNAME, 1989.
%
% Author: Thor I. Fossen
% Date: 2024-07-07

S = [];
W = W(:); % Make sure that W is a column vector

g = 9.81;       % Acceleration of gravity
epsilon = 1e-4; % Small positive value to avoid the W(1) = 0 singularity

% Replace the first value in W with eps if it is zero
if W(1) == 0
    W(1) = epsilon;
end

switch spectrumType
    case 1  % Bretschneider [6]
        A = Parameter(1);
        B = Parameter(2);
        S = arrayfun(@(w) A * w^(-5) * exp(-B / w^4), W);
        TitleStr = 'Bretschneider Spectrum';
        L1Str = sprintf('A = %.2f [m^2 s^{-4}], B = %.2f [s^{-4}]', A, B);

    case 2  % Pierson-Moskowitz [5, 6]
        Vwind20 = Parameter(1);
        A = 8.1e-3 * g^2;
        B = 0.74 * (g / Vwind20)^4;
        S = arrayfun(@(w) A * w^(-5) * exp(-B / w^4), W);
        TitleStr = 'Pierson-Moskowitz Spectrum';
        L1Str = sprintf('Vwind @20m ASL = %.2f [m/s]', Vwind20);

    case 3  % ITTC-Modified Pierson-Moskowitz (Hs, T0)
        Hs = Parameter(1);
        T0 = Parameter(2);
        A = 487 * Hs^2 / T0^4;
        B = 1949 / T0^4;
        S = arrayfun(@(w) A * w^(-5) * exp(-B / w^4), W);
        TitleStr = 'ITTC-Modified Pierson-Moskowitz Spectrum';
        L1Str = sprintf('Hs = %.2f [m], T0 = %.2f [s]', Hs, T0);

    case 4  % ITTC-Modified Pierson-Moskowitz (Hs, T1)
        Hs = Parameter(1);
        T1 = Parameter(2);
        A = 173 * Hs^2 / T1^4;
        B = 691 / T1^4;
        S = arrayfun(@(w) A * w^(-5) * exp(-B / w^4), W);
        TitleStr = 'ITTC-Modified Pierson-Moskowitz Spectrum';
        L1Str = sprintf('Hs = %.2f [m], T1 = %.2f [s]', Hs, T1);

    case 5  % ITTC-Modified Pierson-Moskowitz (Hs, Tz)
        Hs = Parameter(1);
        Tz = Parameter(2);
        A = 123 * Hs^2 / Tz^4;
        B = 495 / Tz^4;
        S = arrayfun(@(w) A * w^(-5) * exp(-B / w^4), W);
        TitleStr = 'ITTC-Modified Pierson-Moskowitz Spectrum';
        L1Str = sprintf('Hs = %.2f [m], Tz = %.2f [s]', Hs, Tz);

    case 6  % JONSWAP (Vwind10, Fetch)
        Vwind10 = Parameter(1);
        fetch = Parameter(2);
        xtilde = g * fetch / Vwind10^2;
        f0 = 3.5 * (g / Vwind10) * xtilde^-0.33;
        w0 = 2 * pi * f0;
        alpha = 0.076 * xtilde^(-0.22);

        sigma = 0.09 * ones(size(W));
        sigma(W < w0) = 0.07;

        S1 = alpha * g^2 * (W .^ -5) .* exp(-(5/4) * (w0 ./ W) .^ 4);
        S2 = 3.3 .^ exp(-(W - w0) .^ 2 ./ (2 * (sigma .* w0) .^ 2));
        S = S1 .* S2;

        TitleStr = 'JONSWAP Spectrum';
        L1Str = sprintf('Vwind @10m ASL = %.2f [m/s], Fetch = %.2f [km]', ...
            Vwind10, fetch / 1000);

    case 7  % JONSWAP (Hs, w0, gamma)
        Hs	  = Parameter(1);
        w0    = Parameter(2);
        gamma = Parameter(3);
        alpha = (5/16) * Hs^2 * w0^4 / g^2;

        % Check gamma value and adjust if necessary
        if gamma < 1 || gamma > 7
            if gamma ~= 0
                disp(['Warning: gamma value is outside validity range, ' ...
                    'using DNV formula'])
            end
            k = 2 * pi / (w0 * sqrt(Hs));
            if k <= 3.6
                gamma = 5;
            elseif k <= 5
                gamma = exp(5.75 - 1.15 * k);
            else
                gamma = 1;
            end
        end

        sigma = 0.09 * ones(size(W));
        sigma(W < w0) = 0.07;

        S1 = alpha * g^2 * (W .^ -5) .* exp(-(5/4) * (w0 ./ W) .^ 4);
        S2 = gamma .^ exp(-(W - w0) .^ 2 ./ (2 * (sigma .* w0) .^ 2));

        % DNV conversion factor
        ConversionFactor = 1 - 0.287 * log(gamma);

        S = ConversionFactor .* S1 .* S2;

        TitleStr = 'JONSWAP Spectrum';
        L1Str = sprintf('gamma = %.2f, Hs = %.2f [m], w0 = %.2f [rad/s]', ...
            gamma, Hs, w0);

    case 8  % Torsethaugen (Hs, w0)
        Hs = Parameter(1);
        w0 = Parameter(2);
        S = torsetSpectrum(Hs, w0, W);
        TitleStr = 'Torsethaugen Spectrum';
        L1Str = sprintf('Hs = %.2f (m), w0 = %.2f (rad/s)', Hs, w0);

    otherwise
        error('Wrong spectrum type identifier, spectrumType = 1, 2,..,8');

end % END switch

% Plot Spectrum if PlotFlag is set to 1, else use 0
if PlotFlag == 1
    figure(gcf); clf;
    plot(W, S, 'linewidth', 2);
    title(TitleStr);
    legend(L1Str);
    xlabel('\omega (rad/s)');
    ylabel('S(\omega) (m^2 s)');
    grid
    set(findall(gcf,'type','line'),'linewidth',2);
    set(findall(gcf,'type','text'),'FontSize',14);
    set(findall(gcf,'type','legend'),'FontSize',14);
end

end % END waveSpectrum
