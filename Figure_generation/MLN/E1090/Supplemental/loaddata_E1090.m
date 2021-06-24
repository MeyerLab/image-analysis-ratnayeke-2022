function sensor = loaddata_E1090(conditions, pth)
 
%% Combine wells%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
allnames = conditions(:,1);
[~,uidx] = unique(allnames,'first');
uniquenames = allnames(sort(uidx));
uniquecondnum = numel(uniquenames);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
condnum = size(conditions,1);
load([pth, 'settings.mat'],'header');

names = header(2,:);
for ii = 1:uniquecondnum
    condrow = find(ismember(conditions(:,1), uniquenames{ii}));
    tracedata = [];
    tracestats = [];
    motherstats = [];
    IFdata_combined = [];
    wellindex_combined = [];
    well_ID_combined = {};
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
                    if exist([pth,'IFdata_', shot, '.mat']);
                        load([pth,'IFdata_', shot, '.mat'], 'IFdata');
                        shotmat = repmat({shot}, size(IFdata,1), 1);
                        IFdata(:, find(ismember(names,'1_DAPI_mean'))) = IFdata(:, find(ismember(names,'1_DAPI_mean')))./median(IFdata(:, find(ismember(names,'1_DAPI_mean'))));
                        IFdata_combined = [IFdata_combined; IFdata];
                        shotmat_combined = [shotmat_combined; shotmat];
                        wellindextemp = ones(size(IFdata,1),3);
                        wellindextemp(:,1) = wellindextemp(:,1)*row;
                        wellindextemp(:,2) = wellindextemp(:,2)*col;
                        wellindextemp(:,3) = wellindextemp(:,3)*site;
                        wellindex_combined = [wellindex_combined; wellindextemp];
                        well_ID = repmat({[num2str(row) '_' num2str(col)]},size(IFdata,1),1);
                        well_ID_combined = [well_ID_combined; well_ID];

                    end
                end
            end
        end
    end
    IFdata_combined(IFdata_combined(:)<0) = 1;
    
    %%% Extract Channels %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    sensor(ii).area = IFdata_combined(:, find(ismember(names,'nuclear area')));
    sensor(ii).DAPI1 = IFdata_combined(:, find(ismember(names,'1_DAPI_mean')));
    sensor(ii).Cy31 = IFdata_combined(:, find(ismember(names,'1_Cy3_mean')));
    sensor(ii).Cy31bg = IFdata_combined(:, find(ismember(names,'1_Cy3_block bg')));
    sensor(ii).FarRed1 = IFdata_combined(:, find(ismember(names,'1_FarRed_mean')));
    sensor(ii).CDK2nuc = IFdata_combined(:, find(ismember(names,'2_YFP_median')));
    sensor(ii).CDK2cyt = IFdata_combined(:, find(ismember(names,'2_YFP_cyto ring')));
    sensor(ii).CDK4nuc = IFdata_combined(:, find(ismember(names,'2_RFP_median')));
    sensor(ii).CDK4cyt = IFdata_combined(:, find(ismember(names,'2_RFP_cyto ring')));
    sensor(ii).FarRed2 = IFdata_combined(:, find(ismember(names,'2_FarRed_mean')));
    sensor(ii).CFP2 = IFdata_combined(:, find(ismember(names,'2_CFP_mean')));
    sensor(ii).FarRed2bg = IFdata_combined(:, find(ismember(names,'2_FarRed_block bg')));
    sensor(ii).CFP2bg = IFdata_combined(:, find(ismember(names,'2_CFP_block bg')));


    sensor(ii).shot = shotmat_combined(:);
    sensor(ii).pos = IFdata_combined(:, 1:2);
    sensor(ii).wellindex = wellindex_combined;
    sensor(ii).well_ID = well_ID_combined;
    
    
end
%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
condnum = length(sensor);
for i = 1:condnum
    sensor(i).dna = (sensor(i).area.*sensor(i).DAPI1);
    ind = sensor(i).area > 200 & sensor(i).area <1000;
    sensor(i)=gateout_all(sensor(i), ind);

    sensor(i).CDK2 = sensor(i).CDK2cyt./sensor(i).CDK2nuc;
    ind = sensor(i).CDK2 < 4 & sensor(i).CDK2 >.1;
    sensor(i)=gateout_all(sensor(i), ind);
    
    sensor(i).CDK4 = sensor(i).CDK4cyt./sensor(i).CDK4nuc;
    ind = sensor(i).CDK2 < 4 & sensor(i).CDK2 >.1;
    sensor(i)=gateout_all(sensor(i), ind);
    
    ind = log2(sensor(i).FarRed2) > 9;
    sensor(i)=gateout_all(sensor(i), ind);
    
    sensor(i).normGem = sensor(i).CFP2 ./ sensor(i).FarRed2;
    sensor(i).CFP2bgsub = sensor(i).CFP2-sensor(i).CFP2bg;
    sensor(i).normGembgsub = sensor(i).CFP2bgsub./sensor(i).FarRed2;
    sensor(i).Cy31bgsub = sensor(i).Cy31-sensor(i).Cy31bg;
end
