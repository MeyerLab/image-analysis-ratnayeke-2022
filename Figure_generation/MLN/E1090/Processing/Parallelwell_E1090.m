clear all;close all;clc;
%%% designate wells to analyze %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
rows  = [1:6];     %1
cols  =  [1:12];    %4:8
sites =  [1:25];     

manualwells  =  [
    1 10 1;

    ];

manualcontrol  =  0;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
numrows = length(rows);
numcols = length(cols); 
numsites = length(sites);
shots = numrows*numcols*numsites;
if manualcontrol == 1
    shots = size(manualwells,1);
end
time1 = tic;
parfor shot = 1:shots
    if manualcontrol == 1
        row = manualwells(shot,1);
        col = manualwells(shot,2);
        site = manualwells(shot,3);
    else
        siteidx = mod(shot,numsites);
        if siteidx == 0
            siteidx = numsites;
        end
        site = sites(siteidx);
        colidx = mod(ceil(shot/numsites),numcols);
        if colidx == 0
            colidx = numcols;
        end
        col = cols(colidx);
        rowidx = ceil(shot/(numcols*numsites));
        row = rows(rowidx);
    end
    fprintf([num2str(row),'_',num2str(col),'_',num2str(site),'\n']);
    try 
    %%% Timelapse %%%%%%%%%%%%%%%%%%%%%%%%%
    IF_multi_E1090(row,col,site);
    catch
        disp(['Error: ',num2str(row),'_',num2str(col),'_',num2str(site),'------------------']);
    end

end
toc(time1)
