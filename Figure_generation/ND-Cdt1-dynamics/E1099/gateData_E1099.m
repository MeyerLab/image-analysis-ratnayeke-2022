function [ S_gated ] = gateData_E1099( S, dataDir )

condNum = size(S,2);
numFrames = size(S(1).area,2);

S_gated = S;

% % GATE bad Mitosis
% for i=1:condNum
%     POI_align = 1;
%     inds=false(size(S_gated(i).apcNormM,1),1);
%     for numcell=1:size(S_gated(i).apcNormM,1)
%         if(S_gated(i).POI(numcell,POI_align)+2<=numFrames & S_gated(i).POI(numcell,POI_align)-5>0)
%             cond1=S_gated(i).apcNormM(numcell,S_gated(i).POI(numcell,POI_align)+2)<.1;
%             cond2=S_gated(i).apcNormM(numcell,S_gated(i).POI(numcell,POI_align)-5)>.6;
%         else
%             cond1=0;
%             cond2=0;
%         end
%         inds(numcell)=cond1 & cond2;
%     end
%     S_gated(i)=gateout_all(S_gated(i),inds);
%     
% end

% GATE bad APC
for i=1:condNum
    POI_align = 2;
    inds1=false(size(S_gated(i).apcNuc,1),1);
    inds2=false(size(S_gated(i).apcNuc,1),1);
    for numcell=1:size(S_gated(i).apcNuc,1)
        if(S_gated(i).POI(numcell,POI_align)+15<=numFrames & S_gated(i).POI(numcell,POI_align)-0>0)
            cond1=S_gated(i).apcNuc(numcell,S_gated(i).POI(numcell,POI_align)+15)> 15 ;
           
            cond2=S_gated(i).apcNuc(numcell,S_gated(i).POI(numcell,POI_align)-0)<15;
        elseif isnan(S_gated(i).POI(numcell,POI_align))
            cond1=1;
            cond2=1;
        else
            %S_gated(i).POI(numcell,POI_align)=NaN;
            cond1=1;
            cond2=1;
        end
        inds1(numcell)=cond1 & cond2;
%         
%         if(S_gated(i).POI(numcell,POI_align)+10<=numFrames & S_gated(i).POI(numcell,POI_align)-0>0)
%             cond1=S_gated(i).apcNormM(numcell,S_gated(i).POI(numcell,POI_align)+10)>.07;
%             cond2=S_gated(i).apcNormM(numcell,S_gated(i).POI(numcell,POI_align)-0)<.03;
%         elseif isnan(S_gated(i).POI(numcell,POI_align))
%             cond1=1;
%             cond2=1;
%         else
%             %S_gated(i).POI(numcell,POI_align)=NaN;
%             cond1=1;
%             cond2=1;
%         end
%         inds2(numcell)=cond1 & cond2;
    end
    S_gated(i)=gateout_all(S_gated(i),inds1);
    
end

% GATE bad CRL4 on
for i=1:condNum
    POI_align = 3;
    inds=false(size(S_gated(i).crlNormAct,1),1);
    for numcell=1:size(S_gated(i).crlNormAct,1)
        if(S_gated(i).POI(numcell,POI_align)+7<=numFrames & S_gated(i).POI(numcell,POI_align)-5>0)
            cond1=S_gated(i).crlNormAct(numcell,S_gated(i).POI(numcell,POI_align)+7) < .2;
            cond2=S_gated(i).crlNormAct(numcell,S_gated(i).POI(numcell,POI_align)-5)<1;
        elseif isnan(S_gated(i).POI(numcell,POI_align))
            cond1=1;
            cond2=1;
        else
            S_gated(i).POI(numcell,POI_align)=NaN;
            cond1=1;
            cond2=1;
        end
        inds(numcell)=cond1 & cond2;
    end
    S_gated(i)=gateout_all(S_gated(i),inds);
    
end

% % GATE bad CRL4 off
% for i=1:condNum
%     POI_align = 4;
%     inds=false(size(S_gated(i).crlNormPostM,1),1);
%     for numcell=1:size(S_gated(i).crlNormPostM,1)
%         if(S_gated(i).POI(numcell,POI_align)+3<=numFrames & S_gated(i).POI(numcell,POI_align)-0>0)
%             cond1=S_gated(i).crlNormPostM(numcell,S_gated(i).POI(numcell,POI_align)+3) > ...
%                 S_gated(i).crlNormPostM(numcell,S_gated(i).POI(numcell,POI_align)) + .025;
%             cond2=S_gated(i).crlNormPostM(numcell,S_gated(i).POI(numcell,POI_align)-0)<.1;
%         elseif isnan(S_gated(i).POI(numcell,POI_align))
%             cond1=1;
%             cond2=1;
%         else
%             S_gated(i).POI(numcell,POI_align)=NaN;
%             cond1=1;
%             cond2=1;
%         end
%         inds(numcell)=cond1 & cond2;
%     end
%     S_gated(i)=gateout_all(S_gated(i),inds);
%     
% end
save([dataDir 'sensordata_gated.mat'],'S_gated','-v7.3');


end

