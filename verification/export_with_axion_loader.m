function export_with_axion_loader(loaderDir, spkFile, outFile)
% EXPORT_WITH_AXION_LOADER  Read a .spk file with Axion's own AxionFileLoader
% (MATLAB) and save every spike to a .mat file, for compare_with_axion_loader.py.
%
%   export_with_axion_loader('/path/to/AxionFileLoader', 'recording.spk', 'axion_ref.mat')
%
% Saved fields (one entry per spike, in the loader's order):
%   wellRow, wellCol, elecCol, elecRow   electrode position (ChannelMapping)
%   startS                                waveform start (s)        = Spike_v1.Start
%   triggerSample                         trigger offset (samples)  = Spike_v1.TriggerSampleOffset
%   stdDev, thrMult                       detector noise SD (V), threshold multiplier
%   data                                  raw samples, nSamples x nSpikes (Spike_v1.Data)
%   volts                                 data .* VoltageScale (V)  = GetVoltageVector
%   allTimes, allAchk, allChannel         SpikeDataSet.LoadAllSpikes (spike times, s)
%   fs, voltageScale, duration            header values

addpath(loaderDir);
A = AxisFile(spkFile);
ds = [];
for k = 1:numel(A.DataSets)
    if isa(A.DataSets(k), 'SpikeDataSet'), ds = A.DataSets(k); break; end
end
if isempty(ds), error('No spike data set in %s', spkFile); end

spikes = ds.LoadData(LoadArgs.ByPlateDimensions);     % flat array of Spike_v1, all electrodes
n = numel(spikes);
fprintf('Axion loader: %d spikes\n', n);

ch = [spikes.Channel];
wellRow = double([ch.WellRow]);  wellCol = double([ch.WellColumn]);
elecCol = double([ch.ElectrodeColumn]); elecRow = double([ch.ElectrodeRow]);
startS = [spikes.Start];
triggerSample = double([spikes.TriggerSampleOffset]);
stdDev = double([spikes.StandardDeviation]);
thrMult = double([spikes.ThresholdMultiplier]);
data = [spikes.Data];                                 % nSamples x nSpikes
volts = double(data) .* ds.VoltageScale;

[el, allTimes] = ds.LoadAllSpikes();
allAchk = double(el.Achk); allChannel = double(el.Channel);

fs = ds.SamplingFrequency; voltageScale = ds.VoltageScale;
duration = []; try, duration = ds.Duration; catch, end

save(outFile, 'wellRow', 'wellCol', 'elecCol', 'elecRow', 'startS', 'triggerSample', 'stdDev', ...
     'thrMult', 'data', 'volts', 'allTimes', 'allAchk', 'allChannel', 'fs', 'voltageScale', 'duration', '-v7');
fprintf('saved %s\n', outFile);
end
