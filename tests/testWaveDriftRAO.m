function testWaveDriftRAO()
% Regression tests for signed second-order mean wave-drift loads.

testSignedPeriodicInterpolation();
testDepthAttenuation();
testInitializationAndDirectionalSum();
testMissingDriftData();

fprintf('testWaveDriftRAO: all tests passed\n');
end

function testDepthAttenuation()
waveModel.hasDriftRAO = true;
waveModel.Amp = 2;
waveModel.Omega = 2;
waveModel.g = 10;
waveModel.raoReferenceDepth = 1;
waveModel.mu = 0;
waveModel.raoAngles = [0; 2*pi];
waveModel.driftCoeff = {[3 3], [-2 -2], [5 5]};

tauReference = waveDriftRAO(waveModel, 0, 0, 1);
tauDeep = waveDriftRAO(waveModel, 0, 0, 3);
factor = exp(-2 * (waveModel.Omega^2 / waveModel.g) * (3 - 1));
assert(norm(tauDeep - factor * tauReference, inf) < 1e-12);

% A depth above the data reference must not amplify the supplied loads.
tauShallow = waveDriftRAO(waveModel, 0, 0, 0);
assert(norm(tauShallow - tauReference, inf) < 1e-12);
end

function testSignedPeriodicInterpolation()
waveModel.hasDriftRAO = true;
waveModel.Amp = 2;
waveModel.mu = 0;
waveModel.raoAngles = [0; pi; 2*pi];
waveModel.driftCoeff = { ...
    [1 3 1], ...
    [-1 -3 -1], ...
    [4 6 4]};

% beta_wave = -pi/2 wraps to 3*pi/2, halfway between pi and 2*pi.
[tau, returnedModel] = waveDriftRAO(waveModel, 0, -pi/2);
expected = [8; -8; 0; 0; 0; 20];
assert(norm(tau - expected, inf) < 1e-12);
assert(isequaln(returnedModel, waveModel));

% The base drift data have no heave/roll/pitch coefficients. An explicit
% reference-point translation can nevertheless create roll/pitch moments.
r_CO_to_data = [0.5; -0.25; 1.2];
tau_CO = Hmtrx(r_CO_to_data)' * tau;
expectedMoment = tau(4:6) + Smtrx(r_CO_to_data) * tau(1:3);
assert(norm(tau_CO(4:6) - expectedMoment, inf) < 1e-12);
assert(tau_CO(3) == 0);
end

function testInitializationAndDirectionalSum()
vessel.main.g = 9.81;
vessel.headings = [0 pi/2 pi 3*pi/2];
vessel.forceRAO.w = [0.05 0.5 2 5];
vessel.forceRAO.amp = cell(1, 6);
vessel.forceRAO.phase = cell(1, 6);
for dof = 1:6
    vessel.forceRAO.amp{dof} = zeros(4, 4, 1);
    vessel.forceRAO.phase{dof} = zeros(4, 4, 1);
end

vessel.drift.w = vessel.forceRAO.w;
coefficient = [2 -3 5];
vessel.drift.coeff = cell(1, 3);
for j = 1:3
    vessel.drift.coeff{j} = coefficient(j) * ones(4, 4, 1);
end

environment.Hs = 1;
environment.w0 = 1;
environment.spectrumType = 'Modified PM';
environment.spreadingFlag = true;
environment.numFreqIntervals = 64;
environment.numDirections = 7;

rng(4);
[environment, waveModel] = waveInitialization(vessel, environment);
[tau, returnedModel] = waveDriftRAO(waveModel, 0.4, 1.1);
amplitudeSquaredSum = sum(environment.Amp(:).^2);
expected = amplitudeSquaredSum * [2; -3; 0; 0; 0; 5];

assert(waveModel.hasDriftRAO);
assert(norm(tau - expected, inf) < 1e-10 * max(1, norm(expected, inf)));
assert(isequaln(returnedModel, waveModel));
end

function testMissingDriftData()
waveModel.hasDriftRAO = false;
try
    waveDriftRAO(waveModel, 0, 0);
    error('testWaveDriftRAO:ExpectedError', ...
        'waveDriftRAO accepted a model without mean-drift data.');
catch ME
    assert(strcmp(ME.identifier, 'waveDriftRAO:InvalidWaveModel'));
end
end
