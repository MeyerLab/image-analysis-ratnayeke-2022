function [Sstart,badtraces]=getPCNAstart(sampletraces,samplestats,settings) 
initoffset = settings.initoffset;
punctathresh =settings.thresh;
punctalowthresh = settings.lowthresh;

%This function will return degstart assuming trace begins with mitosis
[samplesize,tracelength]=size(sampletraces);
firstframe=ones(samplesize,1)*NaN;
lastframe=ones(samplesize,1)*NaN;
Sstart=ones(samplesize,1)*NaN;
Sstartdb=Sstart;
badtraces=zeros(samplesize,1);
sigstore=ones(samplesize,tracelength)*NaN;

for i=1:samplesize
    signal=sampletraces(i,:);
    firstframe(i)=samplestats(i,1);
    lastframe(i)=samplestats(i,2);
    numframes=lastframe(i)-firstframe(i)+1;
%     signal=smooth(signal(firstframe(i):lastframe(i)),3)' ; %incoming signal always raw (unsmoothened)
%     signal = signal - min(signal);
    signal=signal(firstframe(i):lastframe(i)) - min(signal(firstframe(i):lastframe(i)));
    sigstore(i,firstframe(i):lastframe(i))=signal;
    
    %%%%%% general requirements %%%%%%%
    prebuffer=initoffset-1;
    gate=zeros(1,numframes);
    for j=initoffset:numframes
        pastmax=max(signal(j-prebuffer:j-1));
        presentgate=signal(j)>=punctathresh;
        gate(j)=signal(j)>1*pastmax & presentgate;
    end
    punctastart=find(gate,1,'first');
    lowerthresh=punctathresh*0.75;
    if isempty(punctastart) && signal(end)>lowerthresh
        badtraces(i)=1;
    elseif ~isempty(punctastart)
        %%% Find last zero before punctastart call %%%%%%%%%%%%%%%%%%%%%%%%
        earlypuncta=find(signal(1:punctastart-1)<=punctalowthresh & 1:(punctastart-1) > initoffset,1,'last');
        if ~isempty(earlypuncta)
            punctastart_new=earlypuncta;
            Sstartdb(i)=punctastart_new;
            Sstart(i)=firstframe(i)+punctastart_new-1; %account for index
        end
    end
    if any(i == [0]) & settings.debug
        %         clf
        %         figure(1),
        %         subplot(4,1,1)
        %         plot(apcTrace(i,:)), hold on
        %         %ylim([-0.1 1]);
        %         subplot(4,1,2)
        %         plot(trace), hold on
        %         plot(signal,'--')
        %
        %         %ylim([0 1]);
        %         if ~isnan(riseTimeRel(i))
        %             scatter(riseTimeRel(i), trace(riseTimeRel(i)),100,'r.');
        %             hold off
        %             subplot(4,1,1)
        %             scatter(riseTime(i),apcTrace(i,riseTime(i)))
        %             hold off
        %         end
        %         subplot(4,1,3)
        %         plot(sigSlope)
        %         hline(settings.thresh);
        %         if  length(trace) > settings.postBuffer
        %             subplot(4,1,4)
        %             plot(filterScore);
        %             hline(settings.thresh);
        %         end
        keyboard;
    end
end

if settings.debug
traceIDs=1:96;
ylims=[0 200];
POIs = {Sstart};
POIdisplay(traceIDs,sigstore, firstframe, lastframe,badtraces,ylims, POIs);

keyboard;
end
