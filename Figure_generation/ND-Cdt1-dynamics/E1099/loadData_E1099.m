function [S] = loadData_E1099(conditions, dataDir)

%%% Analysis Options %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
setting.motherOption = 0;    %0:no gating 1:mothers 2:no mothers
setting.daughterOption = 0;  %0:no gating 1:daughters 2:no daughters
setting.quiescentAnalysis = 1;   %0:cycling cells 1:serum starved cells
setting.minTraceFrac = .5;
setting.nuc = 'CFP';       %0:diff sensor 1:CDk2 sensor
setting.cdk = '';       %0:diff sensor 1:CDk2 sensor
setting.apc = 'YFP';       %0:diff sensor 1:CDk2 sensor
setting.crl = 'RFP';       %0:diff sensor 1:CDk2 sensor
setting.IFoption = 1;        %0:No IF 1:IF
setting.IFlabel = 'IF_';
setting.poiCdk = 0;
setting.poiApc = 1;
setting.poiCrl = 1;
startFrame = 2;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

daughterOption = setting.daughterOption;
motherOption = setting.motherOption;
quiescentAnalysis = setting.quiescentAnalysis;
if quiescentAnalysis
    motherOption = 2; daughterOption = 2;
end

load([dataDir,'settings.mat'],'names');
names = names(2,:);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allNames = conditions(:,1);
[~,uidx] = unique(allNames,'first');
uniqueNames = allNames(sort(uidx));
uniqueCondnum = numel(uniqueNames);
condNum = size(conditions,1);
for i = 1:uniqueCondnum
    condRow = find(ismember(conditions(:,1),uniqueNames{i}));
    S(i).traceData = [];
    S(i).traceStats = [];
    S(i).motherStats = [];
    S(i).IFdata = [];
    %Traces(i).IFjitter = [];
    S(i).wellindex = [];
    S(i).cellID = [];
    S(i).shot = [];
    for c = condRow'
        rowMat = cell2mat(conditions(c,2));
        colMat = cell2mat(conditions(c,3));
        siteMat = cell2mat(conditions(c,4));
        for row = rowMat
            for col = colMat
                for site = siteMat
                    shot = [num2str(row),'_',num2str(col),'_',num2str(site)];
                    if exist([dataDir,'traceData_',shot,'.mat']) && (~setting.IFoption || exist([dataDir,setting.IFlabel,shot,'.mat']))
                        [traceDatatemp,traceStatstemp,motherStatstemp,IFdatatemp,jitters,samplecellsID] = ...
                            gathertracedata_1_rev05_NR(dataDir,shot,setting.motherOption,setting.daughterOption,setting.IFoption,setting.IFlabel);
                        S(i).traceData = [S(i).traceData; traceDatatemp];
                        S(i).traceStats = [S(i).traceStats; traceStatstemp];
                        S(i).motherStats = [S(i).motherStats; motherStatstemp];
                        S(i).IFdata = [S(i).IFdata; IFdatatemp];
                        S(i).cellID = [S(i).cellID; samplecellsID];
                        wellindexTemp = ones(size(traceDatatemp,1),3);
                        wellindexTemp(:,1) = wellindexTemp(:,1)*row;wellindexTemp(:,2) = wellindexTemp(:,2)*col;wellindexTemp(:,3) = wellindexTemp(:,3)*site;
                        S(i).wellindex = [S(i).wellindex; wellindexTemp];
                        S(i).shot = [S(i).shot; repmat({shot},size(traceDatatemp,1),1)];
                    end
                end
            end
        end
    end
    
    %% Extract Nuclear Channels
%     S(i).traceData = S(i).traceData(:,1:end-1,:);
%     S(i).traceStats(S(i).traceStats == 156) = 155;
%     S(i).traceStats(:,3) = S(i).traceStats(:,2) - S(i).traceStats(:,1) + 1;
    
    S(i).area = S(i).traceData(:,:,ismember(names,'nuclear area'));
    S(i).nucMean = S(i).traceData(:,:,ismember(names,[setting.nuc '_mean']));
    S(i).mass = S(i).area.*S(i).nucMean;
    S(i).massNorm = S(i).mass./repmat(max(S(i).mass,[],2),1,size(S(i).area,2));
    S(i).POI(:, 1) = S(i).traceStats(:,1);
    
    %% Gate on length
    numFrames = size(S(i).traceData,2);
    minLengthTrace = ceil(numFrames*setting.minTraceFrac);
    if setting.quiescentAnalysis
        minLengthTrace =  (numFrames-startFrame)*.75;
    end
    S(i).traceStats(:,5) = findInMat(~isnan(S(i).area));
    badlengths = S(i).traceStats(:,2) - S(i).traceStats(:,5) < minLengthTrace; %| sensor(i).motherStats(:,3)<5;
    S(i)=gateout_all(S(i),~badlengths);
    
    %% Gate nuclear
    noiseThresh = .05; %(.07)
    badNoise = gatenoisy(S(i).massNorm, S(i).traceStats, daughterOption, quiescentAnalysis, noiseThresh, 4,4);
    S(i) = gateout_all(S(i),~badNoise);
    
    %% Extract and gate cdk
    if ~isempty(setting.cdk)
        S(i).cdkNuc = S(i).traceData(:,:,ismember(names,[setting.cdk '_mean']));
        S(i).cdkCyt = S(i).traceData(:,:,ismember(names,[setting.cdk '_cyto ring']));
        S(i).cdkLocalBg = S(i).traceData(:,:,ismember(names, [setting.cdk '_block bg']));
        
        maxThresh = 200; %threshold above which max of each trace must be  %150
        noiseThresh = .5;%0.20; %threshold rate of DHBratio change (not absolute) above which trace is considered noisy
        smoothWindow = 5;
        [S(i).cdk,badTracesCdk] = gate_Cdk2_NR(S(i).cdkNuc,S(i).cdkCyt,maxThresh,noiseThresh,smoothWindow);
        %sensor(i).cdk = sensor(i).cdkCyt./sensor(i).cdkNuc;
        S(i) = gateout_all(S(i),~badTracesCdk);
    end
    
    %% Extract and gate apc
    if ~isempty(setting.apc)
        S(i).apcLocalBg = S(i).traceData(:,:,ismember(names, [setting.apc '_block bg']));
        %[S(i).apcLocalBg, badBg] = fillTraceVals(S(i).apcLocalBg, S(i).traceStats, 5);
        S(i).apcNuc = S(i).traceData(:,:,ismember(names,[setting.apc '_mean']));% - S(i).apcLocalBg;
        [S(i).apcNuc, badShot] = correctBlankFrame(S(i).apcNuc, S(i).shot, 10:95, -10);
        S(i).apcNuc = fillTraceVals(S(i).apcNuc, S(i).traceStats, 5);
        S(i).apcNuc(S(i).apcNuc <=0) = .01;

        
        % Gate out traces on noise and misexpression
        noiseThresh = 200;
        maskTrace = [2:95];
        meanWindow = 2:30;
        lowThresh = 15;
        highStart = S(i).apcNuc(:,1) > lowThresh;
        noiseMask = ones(1,numFrames-1);
        noiseMask(maskTrace)=0;
        clear highTracestart highNoise;
        for n = 1:size(S(i).apcNuc,1)
            highTracestart(n,1) = S(i).apcNuc(n,S(i).traceStats(n,1)) > lowThresh;
            highNoise(n,1) = any(abs(diff(S(i).apcNuc(n,:))) > noiseThresh & noiseMask);
        end
        highMean = nanmean(S(i).apcNuc(:,meanWindow),2) > lowThresh;
        S(i) = gateout_all(S(i),~(highStart | highTracestart | highMean | highNoise));% | badBg));
        
        % Transform data
        %S(i).apcNuc = S(i).apcNuc - min(S(i).apcNuc(:,5:end),[],2);
        S(i).apcArea = (S(i).apcNuc).*S(i).area;
        for n = 1:size(S(i).apcArea,1)
            S(i).apcSmooth(n,:) = nansmooth(S(i).apcArea(n,:),5);
        end
        
        % Find apc inactivation
        if setting.poiApc
            lowThresh = 500;
            highThresh = 300;
            sigThresh = [6000 3000];
            blockSize = 15;
            blockThresh = 2;
            percentile = 95;
             if setting.poiApc
            if setting.quiescentAnalysis
                minstartframe = 30;
                [S(i).POI(:,2), badTracesAPC] = findAPCInact_E1099(S(i).apcNuc, S(i).traceStats, ...
                    struct('postBuffer',6,'smooth',5,'cycling',0,...
                    'buff',3,'thresh',.9,'preBuffer',1,'lowThresh',15, 'increase',5,'trunc',50, 'medfilt',3));
            else
                [S(i).POI(:,2), badTracesAPC] = findAPCInact_E1099(S(i).apcNormM, S(i).traceStats, ...
                    struct('postBuffer',5,'smooth',3,'cycling',1,...
                    'buff',1,'thresh',.9,'preBuffer',1,'lowThresh',.1, 'increase',.025,'trunc',.5, 'medfilt',5));
            end
        end
        end
    end
    
    %% Extract and gate crl
    if ~isempty(setting.crl)
%         S(i).crlLocalBg = S(i).traceData(:,:,ismember(names, [setting.crl '_block bg']));
%         [S(i).crlLocalBg, badBg] = fillTraceVals(S(i).crlLocalBg, S(i).traceStats, 5);
        S(i).crlNuc = S(i).traceData(:,:,ismember(names,[setting.crl '_mean']));% - S(i).crlLocalBg;
        S(i).crlNuc = correctBlankFrame(S(i).crlNuc, S(i).shot, 10:40, -10);
        S(i).crlNuc = fillTraceVals(S(i).crlNuc, S(i).traceStats, 5);
        S(i).crlNuc(S(i).crlNuc <=0) = .01;
        S(i).crlArea = (S(i).crlNuc).*S(i).area;
        
        % Gate out traces on noise and misexpression
        noiseThresh = 20000;
        maskTrace = [];
        meanWindow = 1:30;
        lowThresh = 50;
        noiseMask = ones(1,numFrames-1);
        noiseMask(maskTrace)=0;
        clear lowTracestart highNoise;
        for n = 1:size(S(i).crlNuc,1)
            highNoise(n,1) = any(abs(diff(S(i).crlNuc(n,:)))>noiseThresh & noiseMask);
        end
        lowMean = nanmean(S(i).crlNuc(:,meanWindow),2) < lowThresh;
        %S(i) = gateout_all(S(i),~( lowMean | highNoise | badBg));
        S(i) = gateout_all(S(i),~( highNoise));
        S(i).crlExpress = ~lowMean;
        % Transform data
        for n = 1:size(S(i).crlArea,1)
            ignore_norm = 5;
            S(i).crlNorm(n,:) =  S(i).crlArea(n,:)./max(S(i).crlArea(n,ignore_norm:end-1));
            S(i).crlSmooth(n,:) = nansmooth(S(i).crlArea(n,:), 5);
            %S(i).crlDiff(n,:) = gradient(S(i).crlSmooth(n,:));
        end
        %S(i).crl4Activity = calcCRL4Activity(S(i).crlArea,1:10, 0, 1, 3,2e3);

        % Find crl inactivation
        if setting.poiCrl
            S(i).POI(:,3) = findCRL4Act_E1099(S(i).crlNorm,S(i).traceStats, ...
                struct('cycling',0,'low',.1,'smooth',5,'firstD',-.05,'secD',-.005,'early', 0,'falseCall', .5,'buff',0));
            %S(i).POI(:,4) = aprioriCRL4act(S(i).crl4Activity, S(i).POI(:,3), .008, 5, S(i).crlSmooth);

        end
        
        S(i).crlNormAct = NaN*ones(size(S(i).crlNuc));
        for n = 1:size(S(i).crlArea,1)
            if ~isnan(S(i).POI(n,3))
                S(i).crlNormAct(n,:) =  S(i).crlNuc(n,:)/S(i).crlNuc(n,S(i).POI(n,3));
            else
                S(i).crlNormAct(n,:) = S(i).crlNorm(n,:);
            end
        end
        
        
        
    end
    
    %% transform and gate IF data %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if setting.IFoption
        load([dataDir, 'settings_IF.mat'],'header');
        IFnames = header(2,:);
        
        S(i).IFarea = S(i).IFdata(:, find(ismember(IFnames,'nuclear area')));
        S(i).DAPI = S(i).IFdata(:, find(ismember(IFnames,'DAPI_mean')));
        S(i).RFP = S(i).IFdata(:, find(ismember(IFnames,'RFP_mean')));
        S(i).FarRed = S(i).IFdata(:, find(ismember(IFnames,'FarRed_mean')));
        
        S(i).DAPIbg = S(i).IFdata(:, find(ismember(IFnames,'DAPI_block bg')));
        S(i).RFPbg = S(i).IFdata(:, find(ismember(IFnames,'RFP_block bg')));
        S(i).FarRedbg = S(i).IFdata(:, find(ismember(IFnames,'FarRed_block bg')));
        
        S(i).dna = S(i).DAPI .* S(i).IFarea;
        S(i).x = S(i).IFdata(:, find(ismember(IFnames,'x')));
        S(i).y = S(i).IFdata(:, find(ismember(IFnames,'y')));
        
    end
    
    
end

%% Extra transformations
for i = 1:length(S)
end


save([dataDir 'sensordata.mat'],'S','-v7.3');


end

%% Extra code
% Median filter to get rid of single frame noise
%     for numcell = 1:size(S(i).apcNuc,1)
%         S(i).apcNuc_filt(numcell,:) = S(i).apcNuc(numcell,:);
%         ind = ~isnan(S(i).apcNuc_filt(numcell,:));
%         S(i).apcNuc_filt(numcell,ind) = medfilt1(S(i).apcNuc_filt(numcell,ind),3);
%     end
