function [m_ref, l, mu, cityName] = magneticField(index)
% magneticField s compatible with MATLAB and GNU Octave (www.octave.org).
% This function returns the magnetic field vector, longitude, latitude,
% and city name based on the input index (1 to 30).
%
%  [m_ref, l, mu, cityName] = magneticField(index)
%
% Inputs:
%   index - An integer between 1 and 30 corresponding to a specific city
%
% Outputs:
%   m_ref - 3x1 column vector of magnetic field components [North; East; Down] in nT
%   l     - Longitude of the city in radians
%   mu    - Latitude of the city in radians
%   cityName - Name of the city corresponding to the index
%
% The magnetic field components are hypothetical values in nanoteslas (nT)
% for demonstration purposes. For actual values, use NOAA's online tools.
%
% To get actual magnetic field data for a specific location:
%   1. Visit the NOAA Geomagnetic Calculator:
%      https://www.ngdc.noaa.gov/geomag/calculators/magcalc.shtml
%   2. Enter the latitude, longitude, altitude, and date.
%   3. Retrieve the magnetic field components (North, East, Down) in nT.
%
% Author:
%   Thor I. Fossen
% Date:
%   2024-08-20
% Revisions:

% Data: [Index, City, North (nT), East (nT), Down (nT), Longitude (deg), Latitude (deg)]
data = {
    1,  'Trondheim',    13559,  921,   50209,  10.3951, 63.4305;
    2,  'Oslo',         13500,  900,   50000,  10.7522, 59.9139;
    3,  'New York',     22000, -5000,  45000, -74.0060, 40.7128;
    4,  'Tokyo',        30000,  2000,  35000, 139.6917, 35.6895;
    5,  'Sydney',       25000, 12000, -60000, 151.2093, -33.8688;
    6,  'Cape Town',    15000,  6000, -40000,  18.4241, -33.9249;
    7,  'Moscow',       12000,   500,  55000,  37.6173, 55.7558;
    8,  'São Paulo',    25000, -10000,-20000, -46.6333, -23.5505;
    9,  'Cairo',        35000,  1000,  30000,  31.2357,  30.0444;
    10, 'Beijing',      40000,  7000,  50000, 116.4074, 39.9042;
    11, 'London',       20000,  -500,  48000,  -0.1278,  51.5074;
    12, 'Paris',        21000,   600,  47000,   2.3522,  48.8566;
    13, 'Berlin',       19000,   400,  49000,  13.4050,  52.5200;
    14, 'Rome',         24000,   800,  46000,  12.4964,  41.9028;
    15, 'Madrid',       26000, -1000,  44000,  -3.7038,  40.4168;
    16, 'Toronto',      20000, -4000,  43000, -79.3832,  43.6532;
    17, 'Los Angeles',  22000,  2000,  42000, -118.2437, 34.0522;
    18, 'Mexico City',  28000,  3000,  41000,  -99.1332, 19.4326;
    19, 'Buenos Aires', 24000, -2000, -25000, -58.3816, -34.6037;
    20, 'Rio de Janeiro',26000, -3000,-23000, -43.1729, -22.9068;
    21, 'Mumbai',       34000,  4000,  25000,  72.8777,  19.0760;
    22, 'Delhi',        33000,  3500,  26000,  77.1025,  28.7041;
    23, 'Shanghai',     38000,  5000,  30000, 121.4737,  31.2304;
    24, 'Hong Kong',    36000,  5500,  28000, 114.1694,  22.3193;
    25, 'Singapore',    32000,  4000,  10000, 103.8198,   1.3521;
    26, 'Bangkok',      31000,  3000,  29000, 100.5018,  13.7563;
    27, 'Istanbul',     23000,   700,  46000,  28.9784,  41.0082;
    28, 'Johannesburg', 16000,  6500, -42000,  28.0473, -26.2041;
    29, 'Nairobi',      18000,  5000, -30000,  36.8219,  -1.2921;
    30, 'Moscow',       12000,   500,  55000,  37.6173,  55.7558;
    };

m_ref = [data{index,3}, data{index,4}, data{index,5}]';  % NED components
l = deg2rad(data{index,6}); % Longitude in radians
mu = deg2rad(data{index,7}); % Latitude in radians
cityName = data{index,2}; % City name

end