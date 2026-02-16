%% ====================  OPTOTAGGING — Unified (2ms vs 6ms)  ====================
% Switch between two criteria with 'criterion':
%   '2ms_peak_anywhere'  -> Cerniauskas 2019
%   '6ms_fixed'          -> De Jong 2024 (orthodromic)
%
% Both modes:
%   - Build null by circularly shuffling trials 10,000× in 0–100 ms window
%   - Threshold = 99.9th percentile of shuffled maxima (matching window width)
%   - Fidelity exclusion (>= 0.1) uses the same window as the criterion
%
% I/O:
%   - Input spikes CSV: first row = neuron IDs, subsequent rows = spike times (samples)
%   - Input stim MAT: contains 'selectedPoints' (stim times in samples)
%   - Output: CSV summarizing results + quick plots for passing units
% ------------------------------------------------------------------------

%clc; clear; close all;

%% ---------- USER SETTINGS ----------
criterion            = '6ms_fixed';  % '2ms_peak_anywhere' | '6ms_fixed'
spikeCsv             = 'Curated_clusters_spike_times_columns_sorted1.csv';
stimMat              = 'Stim_R1.mat';   % must load variable 'selectedPoints' (samples)

Fs_Hz                = 30000;            % sampling rate of spike times
binSize_ms           = 1;                % PSTH bin size (ms)
shufflingWindow_ms   = [0 100];          % analysis window after stim, in ms
numShuffles          = 10000;            % circular shuffle iterations
fidelityThreshold    = 0.1;              % must spike in >=10% of trials within criterion window

%% ---------- LOAD DATA ----------
rawData = readmatrix(spikeCsv);
if isempty(rawData)
    error('Could not read spike CSV: %s', spikeCsv);
end

idRow        = rawData(1, :);              % neuron IDs
spikeSamples = rawData(2:end, :);          % spike times (samples)
spikeData_ms = spikeSamples / Fs_Hz * 1000;

S = load(stimMat);                         % must contain selectedPoints (samples)
if ~isfield(S,'selectedPoints')
    error('Expected ''selectedPoints'' in %s', stimMat);
end
selectedPoints_ms = S.selectedPoints / Fs_Hz * 1000;

%% ---------- CHOOSE NEURONS ----------
% neuronInput = input('Enter neuron IDs (e.g. 0,3,27) or "all": ', 's');  % interactive
neuronInput = 'all';  % Set to 'all' or e.g. '0,3,27'
if strcmpi(neuronInput,'all')
    neuronIDs  = idRow;                   
    neuronCols = 1:numel(idRow);
else
    neuronIDs = str2double(strsplit(neuronInput,{',',' ',';'}));
    neuronIDs = neuronIDs(~isnan(neuronIDs));
    [found, neuronCols] = ismember(neuronIDs, idRow);
    if any(~found)
        warning('IDs not found, skipped: %s', mat2str(neuronIDs(~found)));
        neuronIDs  = neuronIDs(found);
        neuronCols = neuronCols(found);
    end
end
if isempty(neuronCols)
    error('No valid neuron IDs selected.');
end

%% ---------- PREP ----------
% Results table (consistent types no matter the mode)
resultsVarNames = {'NeuronID','Criterion','ActualCount','Thresh99_9','PassSignif', ...
                   'Fidelity','PassFidelity','NumTrials', ...
                   'WinningWindowStart_ms','WinningWindowEnd_ms','LatencyInWindow_ms'};
resultsVarTypes = {'double','string','double','double','logical', ...
                   'double','logical','double', ...
                   'double','double','double'};

results = table('Size',[0 numel(resultsVarNames)], ...
                'VariableTypes',resultsVarTypes, ...
                'VariableNames',resultsVarNames);

% Binning
binEdges = shufflingWindow_ms(1) : binSize_ms : shufflingWindow_ms(2);
nBins    = numel(binEdges) - 1;

% Window width per mode
switch lower(criterion)
    case '2ms_peak_anywhere'
        realWindow_ms  = 2;
        realWindowBins = max(1, ceil(realWindow_ms / binSize_ms));
        modeLabel      = "2ms_peak_anywhere";
    case '6ms_fixed'
        realWindow_ms  = 6;
        realWindowBins = max(1, ceil(realWindow_ms / binSize_ms));
        modeLabel      = "6ms_fixed";
    otherwise
        error('Unknown criterion: %s', criterion);
end

significantNeurons = [];
allPSTH = [];        % mean PSTH per passing unit (for plotting)
allBinnedData = [];  % concatenated trial x time heatmaps for passing units

%% ---------- MAIN LOOP ----------
nTrials = numel(selectedPoints_ms);

for k = 1:numel(neuronCols)
    col    = neuronCols(k);
    thisID = neuronIDs(k);

    neuronSpikes = spikeData_ms(:, col);
    neuronSpikes = neuronSpikes(neuronSpikes > 0);  % drop zero padding
    if isempty(neuronSpikes)
        fprintf('Neuron %d: no spikes, skipped\n', thisID);
        continue
    end

    % ---- PSTH (trial x time bins) ----
    binnedData = zeros(nTrials, nBins);
    for tr = 1:nTrials
        stimTime = selectedPoints_ms(tr);
        binnedData(tr, :) = histcounts(neuronSpikes, stimTime + binEdges);
    end
    counts = sum(binnedData, 1);  % 1 x nBins

    % ---- REAL STAT & WIN BOUNDS ----
    switch lower(criterion)
        case '2ms_peak_anywhere'
            % Real: max 2-ms window ANYWHERE in 0–100 ms
            if realWindowBins == 1
                winReal = counts;
            else
                winReal = movsum(counts, realWindowBins, 2, 'Endpoints','discard');
            end
            [actualCount, idxMaxReal] = max(winReal);
            startIdx  = idxMaxReal;                      % sums counts(startIdx : startIdx+W-1)
            endIdx    = idxMaxReal + realWindowBins - 1;
            winStart  = shufflingWindow_ms(1) + (startIdx-1) * binSize_ms;
            winEnd    = shufflingWindow_ms(1) + endIdx     * binSize_ms;

            % Latency (paper-style): mean spike time within winning 2-ms window
            latencies = [];
            for tr = 1:nTrials
                stimTime = selectedPoints_ms(tr);
                sp = neuronSpikes(neuronSpikes >= stimTime + winStart & ...
                                  neuronSpikes <  stimTime + winEnd);
                if ~isempty(sp)
                    latencies = [latencies; (sp(:) - stimTime)]; %#ok<AGROW>
                end
            end
            latencyInWindow = mean(latencies);  % NaN if none

        case '6ms_fixed'
            % Real: fixed 0–6 ms sum
            if realWindowBins > nBins
                error('Real window bins exceed available bins. Check bin size/window.');
            end
            actualCount = sum(counts(1:realWindowBins));
            startIdx  = 1;
            endIdx    = realWindowBins;
            winStart  = shufflingWindow_ms(1);
            winEnd    = shufflingWindow_ms(1) + realWindow_ms;
            % Latency not defined in this paper; leave NaN
            latencyInWindow = NaN;
    end

    % ---- NULL: shuffled max over SAME-WIDTH window ----
    shuffledMax = zeros(numShuffles,1);
    for s = 1:numShuffles
        shuffled = zeros(size(binnedData));
        shifts   = randi([1, nBins], nTrials, 1);
        for tr = 1:nTrials
            shuffled(tr,:) = circshift(binnedData(tr,:), shifts(tr));
        end
        scounts = sum(shuffled,1);
        if realWindowBins == 1
            swin = scounts;
        else
            swin = movsum(scounts, realWindowBins, 2, 'Endpoints','discard');
        end
        shuffledMax(s) = max(swin);
    end
    thr = prctile(shuffledMax, 99.9);
    passSignif = (actualCount > thr);   % strict >

    % ---- FIDELITY (same window as criterion) ----
    trialHasSpike = any(binnedData(:, startIdx:endIdx), 2);
    fidelity = mean(trialHasSpike);          % proportion (0..1), no rounding in test
    passFid  = (fidelity >= fidelityThreshold);

    % ---- RECORD ----
    results = [results; { ...
        thisID, modeLabel, actualCount, thr, passSignif, ...
        round(fidelity,3), passFid, nTrials, ...
        winStart, winEnd, latencyInWindow}]; %#ok<AGROW>

    % ---- LOG & BOOKKEEP ----
    if passSignif && passFid
        fprintf('Neuron %d: PASSED  | %s  Actual=%d  Thr=%.1f  Fid=%.3f  Win=[%.1f, %.1f] ms\n', ...
            thisID, criterion, actualCount, thr, fidelity, winStart, winEnd);
        significantNeurons = [significantNeurons, thisID]; %#ok<AGROW>
        allPSTH       = [allPSTH; mean(binnedData,1)];     %#ok<AGROW>
        allBinnedData = [allBinnedData; binnedData];       %#ok<AGROW>
    else
        fprintf('Neuron %d: FAILED | %s  PassSig=%d  PassFid=%d  Actual=%d  Thr=%.1f  Fid=%.3f\n', ...
            thisID, criterion, passSignif, passFid, actualCount, thr, fidelity);
    end
end

%% ---------- SAVE & REPORT ----------
outCsv = sprintf('optotag_unified_%s.csv', criterion);
writetable(results, outCsv);

fprintf('\n--- RESULTS (%s) saved to %s ---\n', criterion, outCsv);
disp(results);

if isempty(significantNeurons)
    disp('No neurons met all tagging criteria.');
else
    fprintf('Passing neuron IDs: %s\n', mat2str(significantNeurons));
end

%% ---------- PLOTS (only if at least one passing unit) ----------
if ~isempty(significantNeurons)
    figure('Name','Aggregate PSTH (passing units)');
    plot(mean(allPSTH,1), 'LineWidth',2);
    xlabel(sprintf('%d ms bins (%d–%d ms)', binSize_ms, shufflingWindow_ms(1), shufflingWindow_ms(2)));
    ylabel('Average spike count (per bin across trials)');
    title(sprintf('Aggregate PSTH — criterion: %s', criterion));
    grid on;

    figure('Name','Heatmap — Trials x Time (passing units)');
    imagesc(allBinnedData);
    colormap('hot'); colorbar;
    xlabel(sprintf('Time bins (bin=%d ms; %d–%d ms)', binSize_ms, shufflingWindow_ms(1), shufflingWindow_ms(2)));
    ylabel('Trial (concatenated across passing units)');
    title('Heatmap — Passing Units (Trials x Time)');
end