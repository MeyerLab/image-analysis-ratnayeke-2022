function sensor = loadData_C183(conditions, pth)

%% Combine wells%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allnames = conditions(:,1);
[~,uidx] = unique(allnames,'first');
uniquenames = allnames(sort(uidx));
uniquecondnum = numel(uniquenames);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
condnum = size(conditions,1);
load([pth, 'settings_live.mat'],'names');
liveheader = names;
namesLive = liveheader(2,:);

load([pth, 'settings_live_IF.mat'],'header');
IFheader = header;
namesIF = IFheader(2,:);
for ii = 1:uniquecondnum
    condrow = find(ismember(conditions(:,1), uniquenames{ii}));
    tracedata = [];
    tracestats = [];
    motherstats = [];
    IFdata_combined = [];
    tracedata_combined = [];
    wellindex_combined = [];
    shotmat_combined = [];
    cc = 0;
    for c = condrow'
        rowmat = cell2mat(conditions(c,2));
        colmat = cell2mat(conditions(c,3));
        sitemat = cell2mat(conditions(c,4));
        for row = rowmat
            for col = colmat
                for site = sitemat
                    cc = cc+1;
                    shot = [num2str(row),'_',num2str(col), '_', num2str(site)];
                    if exist([pth,'IF_', shot, '.mat']) & exist([pth,'tracedata_', shot, '.mat']);
                        load([pth,'IF_', shot, '.mat'], 'IFdata');
                        shotmat = repmat({shot}, size(IFdata,1), 1);
                        %IFdata(:, find(ismember(names,'1_DAPI_mean'))) = IFdata(:, find(ismember(names,'1_DAPI_mean')))./median(IFdata(:, find(ismember(names,'1_DAPI_mean'))));
                        IFdata_combined = [IFdata_combined; IFdata];
                        shotmat_combined = [shotmat_combined; shotmat];
                        wellindextemp = ones(size(IFdata,1),3);
                        wellindextemp(:,1) = wellindextemp(:,1)*row;
                        wellindextemp(:,2) = wellindextemp(:,2)*col;
                        wellindextemp(:,3) = wellindextemp(:,3)*site;
                        wellindex_combined = [wellindex_combined; wellindextemp];
                        load([pth,'tracedata_', shot, '.mat'], 'tracedata');
                        tracedata_combined = [tracedata_combined; squeeze(tracedata)];
                        
                    end
                end
            end
        end
    end
    IFdata_combined(IFdata_combined(:)<0) = 1;
    
    %%% Extract Channels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sensor(ii).area = IFdata_combined(:, find(ismember(namesIF,'nuclear area')));
    sensor(ii).DAPI1 = IFdata_combined(:, find(ismember(namesIF,'1_DAPI_mean')));
    sensor(ii).FarRed1 = IFdata_combined(:, find(ismember(namesIF,'1_FarRed_mean')));
    sensor(ii).YFP2 = IFdata_combined(:, find(ismember(namesIF,'2_YFP_median')));
    sensor(ii).FarRed2 = IFdata_combined(:, find(ismember(namesIF,'2_FarRed_median')));
    
    sensor(ii).YFPPunctaArea = IFdata_combined(:,contains(namesIF,'2_YFP_puncta area'));
    sensor(ii).YFPIntensity = IFdata_combined(:,contains(namesIF,'2_YFP_puncta intensity'));
    sensor(ii).YFPPunctaNum = IFdata_combined(:,contains(namesIF,'2_YFP_puncta number'));
        
    sensor(ii).YFPnuclive = tracedata_combined(:, find(ismember(namesLive,'YFP_mean')));
    sensor(ii).YFPcytlive = tracedata_combined(:, find(ismember(namesLive,'YFP_cyto ring')));
    sensor(ii).cdk2 = sensor(ii).YFPcytlive./sensor(ii).YFPnuclive;
    
    sensor(ii).shot = shotmat_combined(:);
    sensor(ii).pos = IFdata_combined(:, 1:2);
    sensor(ii).wellindex = wellindex_combined;
    
    
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
condnum = length(sensor);
for i = 1:condnum
    
    sensor(i).dna = (sensor(i).area.*sensor(i).DAPI1);
    ind = sensor(i).area > 1000 &  sensor(i).area < 3000 ...
    &(sensor(i).YFPnuclive + sensor(i).YFPcytlive) > 200 & sensor(i).cdk2 < 2.5 & sensor(i).cdk2 > .2;
    sensor(i)=gateout_all(sensor(i), ind);
    sensor(i).FarRed2corr = sensor(i).FarRed2 - (sensor(i).FarRed1*.0131);
end
end