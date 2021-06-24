%%
%file_path = 'E:\Nikon\HCA_Live2\20190523_131630_527\WellB02_ChannelCFP,YFP,mCherry_Seq0000.nd2';
%file_path = 'E:\Nikon\HCA_Live1\20190523_111647_383\WellB02_ChannelCFP,YFP,mCherry_Seq0000.nd2';
file_path = 'E:\Nikon\HCA_Live1\20190522_193505_888\WellB02_ChannelCFP,YFP,mCherry_Seq0000.nd2';
%file_path = 'E:\Nikon\HCA_Live_fixed\20190524_215451_884\WellB02_ChannelDAPI,Cy5,YFP_Seq0000.nd2';

r = loci.formats.Memoizer(bfGetReader(), 0);
r.setId(file_path);
r.close();
% s = javaObject('loci.formats.in.DynamicMetadataOptions');
% s.set('nativend2.chunkmap', 'false');
% reader.setMetadataOptions(s);

%%
%data = bfopen(file_path);
%%
% %%
% clear I;
% iSite = 1;
% iTime = 1;
% iSeries = 1;
% 
% for i = 1:3
% iZ = 1;
% iC = i;
% iT = 5;
% reader.setSeries(iSeries - 1);
% iPlane = reader.getIndex(iZ - 1, iC -1, iT - 1) + 1;
% I(:,:,i) = imadjust(mat2gray(bfGetPlane(reader, iPlane)));
% 
% end
% figure, imshow(I);
% 

%%
bfInitLogging('INFO');
reader = javaObject('loci.formats.Memoizer', bfGetReader(),0);
reader.setId(file_path);

reader.setSeries(0);
sites = reader.getSeriesCount();
channels = reader.getSizeC();
time = reader.getSizeT();

STmat = nd2Map(time,sites);

% %%
clear I;
iSite = 1;
iTime = 3;

for i = 1:channels
iZ = 1;
iC = i;
iT = STmat{iSite,iTime}(2);
iSeries = STmat{iSite,iTime}(1);
reader.setSeries(iSeries - 1);
iPlane = reader.getIndex(iZ - 1, iC -1, iT - 1) + 1;
I(:,:,i) = imadjust(mat2gray(bfGetPlane(reader, iPlane)));

end
figure, imshow(I);


%% Generate movie
clear I;
iSite = 1;
iC = 1;
iZ = 1;
cmap = gray(256);

for t = 1:90

iTime = t;
iT = STmat{iSite,iTime}(2);
iSeries = STmat{iSite,iTime}(1);
reader.setSeries(iSeries - 1);
iPlane = reader.getIndex(iZ - 1, iC -1, iT - 1) + 1;
image = bfGetPlane(reader, iPlane);
image = imresize(image,.25);
%figure, imagesc(image);
M(t) = im2frame(uint8(image), cmap);
end
movie(M);