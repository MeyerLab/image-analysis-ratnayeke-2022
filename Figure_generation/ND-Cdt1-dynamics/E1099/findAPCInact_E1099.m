function [ riseTime, badTraces] = findAPCInact_D103( apcTrace, traceStats, settings )
%
[sampleSize,traceLength] = size(apcTrace);
riseTime = ones(sampleSize,1)*NaN;
riseTimeRel = ones(sampleSize,1)*NaN;

badTraces = false(sampleSize,1);

for i=1:size(apcTrace, 1)
    trace = apcTrace(i,:);
    trace = trace - min(trace);
    if settings.cycling == 1
        mitosisFrame = traceStats(i, 1);
        %startFrame = mitosisFrame + settings.buff;
        startFrame = find(trace < settings.trunc & 1:length(trace) >= mitosisFrame+settings.buff,1,'first');
        endFrame = traceStats(i,2);
    else
        %startFrame = traceStats(i,5);
        startFrame = find(trace < settings.trunc +settings.buff,1,'first');
        endFrame = traceStats(i,2);
    end
    trace = trace(startFrame:endFrame);
    trunc = find(trace > settings.trunc,1,'first');
    if ~isempty(trunc)
        trace = trace(1:trunc);
    end
    
    if length(trace) > settings.postBuffer
        trace = trace - min(trace);
        signal = nansmooth(trace, settings.smooth);
        %%% build risetime filter %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        %%%%%% general requirements
        numFrames = size(trace,2);
        gate = zeros(1,numFrames);
        for j = 1:numFrames - settings.postBuffer
            presentHeight = signal(j) < settings.lowThresh; %default 0.05
            minFutureHeight = min(signal(j+1:end)) > signal(j);
            bufferFutureHeight = signal(j+settings.postBuffer) > signal(j) + settings.increase; %default 0.1
            gate(j) = presentHeight & minFutureHeight & bufferFutureHeight;
        end
        
        sigFwdSlope = 10*getslope_forward_avg(signal, 1:settings.postBuffer);
        sigTime = (1:length(signal))/traceLength;
        filterScore = 4*sigFwdSlope - 2*signal + sigTime + 1;
        filterScore = filterScore.*gate;
        if settings.medfilt
            filterScore = medfilt1(filterScore,settings.medfilt);
        end
        %tempsearch=find(sig_fwdslope>0.05 & abs(signal)<0.03,1,'last');
        tempSearch=find(filterScore > settings.thresh,1,'last');

        
        %     filtermax=max(filterScore);
        %     tempSearch=find(filterScore==filtermax,1,'first');
        %if isempty(tempsearch) || signal(end)<0.05 %0.05
        
        if ~isempty(tempSearch) %&& filtermax>0
            while tempSearch > 2 && filterScore(tempSearch-1) > settings.thresh
                tempSearch = tempSearch-1;
            end
            riseTimeRel(i) = tempSearch;
            riseTime(i) = tempSearch + startFrame-1; %return absolute POI rather than relative to mitosis
            %elseif signal(end)<0.1
        elseif max(signal(settings.preBuffer:end))<settings.lowThresh
            riseTime(i)=NaN;
        else
            badTraces(i) = 1;
        end
    else
        trace = apcTrace(i,:);
        signal = trace;
        filterScore = [];
        badTraces(i) = 1;
    end
    
    
%         clf
%         figure(1),
%         subplot(3,1,1)
%         plot(apcTrace(i,:)), hold on
%         %ylim([-0.1 1]);
%         subplot(3,1,2)
%         plot(trace), hold on
%         plot(signal,'--')
%         %ylim([0 1]);
%         if ~isnan(riseTimeRel(i))
%             scatter(riseTimeRel(i), trace(riseTimeRel(i)),200,'r.');
%             hold off
%             subplot(3,1,1)
%             scatter(riseTime(i),apcTrace(i,riseTime(i)),200,'r.')
%             hold off
%         end
%         if  length(trace) > settings.postBuffer
%             subplot(3,1,3)
%             plot(filterScore);
%             hline(settings.thresh);
%         end
%         keyboard;
    
end
