function [ riseTime, badTraces] = findPCNArise( filtTrace, traceStats, settings )
%
[sampleSize,traceLength] = size(filtTrace);
riseTime = ones(sampleSize,1)*NaN;
riseTimeRel = ones(sampleSize,1)*NaN;
firstFrame = ones(sampleSize,1)*NaN;
lastFrame = ones(sampleSize,1)*NaN;
fullTraceStore = ones(sampleSize,traceLength)*NaN;
badTraces = false(sampleSize,1);

for i=1:size(filtTrace, 1)
    trace = filtTrace(i,:);
    trace = trace - min(trace);
    fullTraceStore(i,:) = trace;
    if settings.cycling == 1
        mitosisFrame = traceStats(i, 1);
        %startFrame = mitosisFrame + settings.buff;
        startFrame = find(trace < settings.trunc & 1:length(trace) >= mitosisFrame+settings.buff,1,'first');
        if isempty(startFrame)
            continue;
        end
        endFrame = traceStats(i,2);
    else
        %startFrame = traceStats(i,5);
        startFrame = find(trace < settings.trunc,1,'first');
        endFrame = traceStats(i,2);
    end
    firstFrame(i) = startFrame;
    lastFrame(i) = endFrame;
    
    trace = trace(startFrame:endFrame);
    trunc = find(trace > settings.trunc,1,'first');
    if ~isempty(trunc)
        trace = trace(1:trunc);
    end
    
    if length(trace) > settings.postBuffer + settings.preBuffer
        %trace = trace - min(trace);
        if settings.smooth
            signal = nansmoothm(trace, settings.smooth,'sgolay');
        else
            signal = trace;
        end
        %%% build risetime filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% general requirements
        numFrames = size(trace,2);
        presentHeight = zeros(1,numFrames);
         prevHeight = zeros(1,numFrames);
        minFutureHeight = zeros(1,numFrames);
        futureHeightInc = zeros(1,numFrames);
        allSlope = zeros(1,numFrames);
        early = zeros(1,numFrames);
        sigSlope = gradient(signal);
        for j = 1+settings.preBuffer:numFrames - settings.postBuffer
            presentHeight(j) = signal(j) < settings.lowThresh; %default 0.05
            prevHeight(j) = max(signal(j-settings.preBuffer:j)) < settings.lowThresh;            
            minFutureHeight(j) = min(signal(j+1:j+settings.postBuffer)) > signal(j);
            futureHeightInc(j) = signal(j+settings.postBuffer) > signal(j) + settings.increase; %default 0.1
            if settings.cycling
                early(j) = j + settings.buff > settings.early;
            else
                early(j) = j + startFrame -1 > settings.early;
            end
        end
        gate = presentHeight & minFutureHeight & futureHeightInc & early & prevHeight;
        
        filterScore = sigSlope.*gate;
        tempSearch=find(gate,1,'first');
        
        if ~isempty(tempSearch)
            %             while tempSearch > 2 && filterScore(tempSearch-1) > settings.thresh
            %                 tempSearch = tempSearch-1;
            %             end
            riseTimeRel(i) = tempSearch;
            riseTime(i) = tempSearch + startFrame-1; %return absolute POI rather than relative to mitosis
        elseif max(signal(1:end - settings.postBuffer))<settings.lowThresh
            riseTime(i)=NaN;
        else
            badTraces(i) = 1;
        end
    else
        trace = filtTrace(i,:);
        signal = trace;
        filterScore = [];
        badTraces(i) = 1;
    end
    
    if any(i == []) & settings.debug
        clf
        figure(1),
        subplot(4,1,1)
        plot(filtTrace(i,:)), hold on
        %ylim([-0.1 1]);
        subplot(4,1,2)
        plot(trace), hold on
        plot(signal,'--')
        
        %ylim([0 1]);
        if ~isnan(riseTimeRel(i))
            scatter(riseTimeRel(i), trace(riseTimeRel(i)),100,'r.');
            hold off
            subplot(4,1,1)
            scatter(riseTime(i),filtTrace(i,riseTime(i)))
            hold off
        end
        subplot(4,1,3)
        plot(sigSlope)
        hline(settings.thresh);
        if  length(trace) > settings.postBuffer
            subplot(4,1,4)
            plot(filterScore);
            hline(settings.thresh);
        end
        keyboard;
    end
end

%%% debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
if settings.debug
    traceIDs=[1:96*2];
    ylims=[0 200];
    POIs = {riseTime};
    POIdisplay(traceIDs,fullTraceStore, firstFrame, lastFrame,badTraces,ylims, POIs);
    keyboard;
end