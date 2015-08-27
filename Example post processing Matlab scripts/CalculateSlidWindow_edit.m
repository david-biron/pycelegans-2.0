SegNum = 18; 
% frameRate = 10; 
scrsz = get(0,'screensize'); 

TotSeg = 20;                    %how many segments to divide worm by
SegNum = TotSeg - 2;            %number of segments in the analysis
SegLngth = 100/TotSeg;          %the point length of each segment along worm spline.
SEGMENTS = 1:SegNum;            %set the range of segments for the analysis
FRAMES = 1:NumFrames;           %set the range to be analyzed

HalfWidth = 5*60*frameRate;
%% Compute Sliding Window of Quiescence

QBoutDurSW = zeros(2,length(RT_L)) * NaN; %Orig,FullQuiesc,BodyQuiesc
QBoutNumSW = zeros(2,length(RT_L));
QFSW = zeros(2,length(RT_L))*NaN;
BehaviorPrc10min = zeros(4,length(RT_L))*NaN; %FRW,REV,DWELL,QUIESC
WL_Slid = zeros(1,length(RT_L))*NaN;
for i = HalfWidth+1:length(RT_L)-HalfWidth
    for j = 1:2
        NaNfrac = sum(RT_L(i-HalfWidth:i+HalfWidth)==0)/length(RT_L(i-HalfWidth:i+HalfWidth));
        % If 33% of the data is missing, don't compute
        if NaNfrac < (1/3)
            temp = QBoutDur(j,i-HalfWidth:i+HalfWidth);
            QBoutDurSW(j,i) = sum(temp(isfinite(temp)))/length(temp(isfinite(temp)));
            QBoutNumSW(j,i) = sum(QBoutNum(j,i-HalfWidth:i+HalfWidth));
            temp = QboutAll_LRT(j,i-HalfWidth:i+HalfWidth);
            tt = (RT_L(i-HalfWidth:i+HalfWidth)==1); temp = temp(tt); %take out points where no data
            
            QFSW(j,i) = sum(temp(isfinite(temp)))/length(temp(isfinite(temp)));
        end
        
    end
    if NaNfrac < (1/3)
        for j = 1:4
            temp = Behavior(j,i-HalfWidth:i+HalfWidth);
            BehaviorPrc10min(j,i) = sum(temp(isfinite(temp)))/length(temp(isfinite(temp)));
        end
        temp = WormLengthRT(i-HalfWidth:i+HalfWidth);
        WL_Slid(i) = nanmean(temp);
    end
end
%%
save([Worm ' BehaviorSlidWindow'],'frameRate','scrsz','Worm','SegNum','FRAMES','SegLngth','NumFrames','FRAMES','SEGMENTS','frames','RT','RT_L',...
    'QBoutDurSW','QBoutNumSW','QFSW','BehaviorPrc10min','WL_Slid');


