function reformat3i
paths = {'F:\Data\3i Data\PCNA 45\C143-PCNA3i\1_1_1\'};
rows = [1];
cols = [1];
channels = {'C0_','C1_','C2_'};

outchannel = {'CFP_','YFP_','RFP_'};
outpath = 'F:\Data\3i Data\PCNA 45\C143-PCNA3i\Raw\';
copy = true;
            

if ~exist(outpath)
    mkdir(outpath)
end

for p = 1:length(paths)
    path = paths{p};
    files = dir([path '*.tiff']);
    filenames = extractfield(files,'name');
    for r = 1:length(rows)
        row = rows(r);
        for c = 1:length(cols)
            col = cols(c); 
            dest_folder = [outpath num2str(p) '_' num2str(row) '_' num2str(col) '\'];
            if ~exist(dest_folder)
                mkdir(dest_folder)
            end
            timepoints = 0:90;
            for chan = 1:length(channels)
                channel = channels{chan};
                for t = 1:length(timepoints)
                    old_file = strcat(path,channel, sprintf('%02d',timepoints(t)), '.tiff');
                    new_file = strcat(dest_folder,num2str(p),'_',num2str(row),'_',num2str(col),'_', outchannel{chan}, num2str(t), '.tif');
                    if copy
                        copyfile(old_file,new_file,'f');
                    else
                        movefile(old_file,new_file,'f');
                    end
                end
            end
        end
    end
end
end
