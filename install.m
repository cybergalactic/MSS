clear all
close all

% Configure build folder structure
Simulink.fileGenControl('set', ...
    'CacheFolder', 'cache', ...
    'CodeGenFolder', 'build',...
    'createDir',true...
);

% Install by adding to path
addpath(genpath(pwd()));

% Load supply ship
load('supply.mat');
