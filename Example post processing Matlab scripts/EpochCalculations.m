scrsz = get(0,'screensize');
getsize = @(x)x(1);
HourXaxis = RT/10/60/60;

%%
% load('F:\PyCelegansW\Analysis\N2\W20009_N2\N2 20009 Behavior.mat');
% load('F:\PyCelegansW\Analysis\N2\W20009_N2\N2 20009 BehaviorSlidWindow.mat');
%%
BehaviorPrc1min = zeros(4,length(RT_L))*NaN; %FRW,REV,DWELL,QUIESC
EdgePrcSW1min = NaN(1,length(RT));
for i = 301:length(RT_L)-300
    NaNfrac = sum(RT_L(i-300:i+300)==0)/length(RT_L(i-300:i+300));
    % If 33% of the data is missing, don't compute
    if NaNfrac < (1/3)
        for j = 1:4
            temp = Behavior(j,i-300:i+300);
            BehaviorPrc1min(j,i) = sum(temp(isfinite(temp)))/length(temp(isfinite(temp)));
            temp = EdgePrc(i-300:i+300);
            EdgePrcSW1min(i) = mean(temp(isfinite(temp)));
        end
    end
end
BehaviorPrc2min = zeros(4,length(RT_L))*NaN; %FRW,REV,DWELL,QUIESC
EdgePrcSW2min = NaN(1,length(RT));
for i = 601:length(RT_L)-600
    NaNfrac = sum(RT_L(i-600:i+600)==0)/length(RT_L(i-600:i+600));
    % If 33% of the data is missing, don't compute
    if NaNfrac < (1/3)
        for j = 1:4
            temp = Behavior(j,i-600:i+600);
            BehaviorPrc2min(j,i) = sum(temp(isfinite(temp)))/length(temp(isfinite(temp)));
            temp = EdgePrc(i-600:i+600);
            EdgePrcSW2min(i) = mean(temp(isfinite(temp)));
        end
    end
end
BehaviorPrc5min = zeros(4,length(RT_L))*NaN; %FRW,REV,DWELL,QUIESC
EdgePrcSW5min = NaN(1,length(RT));
for i = 1501:length(RT_L)-1500
    NaNfrac = sum(RT_L(i-1500:i+1500)==0)/length(RT_L(i-1500:i+1500));
    % If 33% of the data is missing, don't compute
    if NaNfrac < (1/3)
        for j = 1:4
            temp = Behavior(j,i-1500:i+1500);
            BehaviorPrc5min(j,i) = sum(temp(isfinite(temp)))/length(temp(isfinite(temp)));
            temp = EdgePrc(i-1500:i+1500);
            EdgePrcSW5min(i) = mean(temp(isfinite(temp)));
        end
    end
end
BehaviorPrc20min = zeros(4,length(RT_L))*NaN; %FRW,REV,DWELL,QUIESC
EdgePrcSW20min = NaN(1,length(RT));
for i = 6001:length(RT_L)-6000
    NaNfrac = sum(RT_L(i-6000:i+6000)==0)/length(RT_L(i-6000:i+6000));
    % If 33% of the data is missing, don't compute
    if NaNfrac < (1/3)
        for j = 1:4
            temp = Behavior(j,i-6000:i+6000);
            BehaviorPrc20min(j,i) = sum(temp(isfinite(temp)))/length(temp(isfinite(temp)));
            temp = EdgePrc(i-6000:i+6000);
            EdgePrcSW20min(i) = mean(temp(isfinite(temp)));
        end
    end
end
%%

PrimaryBehavior1min = zeros(1,length(RT))*NaN;
for i = 1:length(RT)
    if isfinite(BehaviorPrc1min(:,i))
        PrimaryBehavior1min(i) = find(BehaviorPrc1min(:,i) == max(BehaviorPrc1min(:,i)),1,'first');
    end
end
Edges1min = (find(diff(PrimaryBehavior1min)));
Epochs1min = {};
kF = 1; kR = 1; kD = 1; kQ = 1;
for i = 1:length(Edges1min)-1
    if PrimaryBehavior1min(Edges1min(i)+1) == 1
        Epochs1min{1}(kF,:) = [Edges1min(i)+1 Edges1min(i+1)+1];
        kF = kF+1;
    end
    if PrimaryBehavior1min(Edges1min(i)+1) == 2
        Epochs1min{2}(kR,:) = [Edges1min(i)+1 Edges1min(i+1)+1];
        kR = kR+1;
    end
    if PrimaryBehavior1min(Edges1min(i)+1) == 3
        Epochs1min{3}(kD,:) = [Edges1min(i)+1 Edges1min(i+1)+1];
        kD = kD+1;
    end
    if PrimaryBehavior1min(Edges1min(i)+1) == 4
        Epochs1min{4}(kQ,:) = [Edges1min(i)+1 Edges1min(i+1)+1];
        kQ = kQ+1;
    end
end
if kF == 1; Epochs1min{1} = []; end
if kR == 1; Epochs1min{2} = []; end
if kD == 1; Epochs1min{3} = []; end
if kQ == 1; Epochs1min{4} = []; end
%%%%%%
PrimaryBehavior2min = zeros(1,length(RT))*NaN;
for i = 1:length(RT)
    if isfinite(BehaviorPrc2min(:,i))
        PrimaryBehavior2min(i) = find(BehaviorPrc2min(:,i) == max(BehaviorPrc2min(:,i)),1,'first');
    end
end
Edges2min = (find(diff(PrimaryBehavior2min)));
Epochs2min = {};
kF = 1; kR = 1; kD = 1; kQ = 1;
for i = 1:length(Edges2min)-1
    if PrimaryBehavior2min(Edges2min(i)+1) == 1
        Epochs2min{1}(kF,:) = [Edges2min(i)+1 Edges2min(i+1)+1];
        kF = kF+1;
    end
    if PrimaryBehavior2min(Edges2min(i)+1) == 2
        Epochs2min{2}(kR,:) = [Edges2min(i)+1 Edges2min(i+1)+1];
        kR = kR+1;
    end
    if PrimaryBehavior2min(Edges2min(i)+1) == 3
        Epochs2min{3}(kD,:) = [Edges2min(i)+1 Edges2min(i+1)+1];
        kD = kD+1;
    end
    if PrimaryBehavior2min(Edges2min(i)+1) == 4
        Epochs2min{4}(kQ,:) = [Edges2min(i)+1 Edges2min(i+1)+1];
        kQ = kQ+1;
    end
end
if kF == 1; Epochs2min{1} = []; end
if kR == 1; Epochs2min{2} = []; end
if kD == 1; Epochs2min{3} = []; end
if kQ == 1; Epochs2min{4} = []; end
%%%%
PrimaryBehavior5min = zeros(1,length(RT))*NaN;
for i = 1:length(RT)
    if isfinite(BehaviorPrc5min(:,i))
        PrimaryBehavior5min(i) = find(BehaviorPrc5min(:,i) == max(BehaviorPrc5min(:,i)),1,'first');
    end
end
Edges5min = (find(diff(PrimaryBehavior5min)));
Epochs5min = {};
kF = 1; kR = 1; kD = 1; kQ = 1;
for i = 1:length(Edges5min)-1
    if PrimaryBehavior5min(Edges5min(i)+1) == 1
        Epochs5min{1}(kF,:) = [Edges5min(i)+1 Edges5min(i+1)+1];
        kF = kF+1;
    end
    if PrimaryBehavior5min(Edges5min(i)+1) == 2
        Epochs5min{2}(kR,:) = [Edges5min(i)+1 Edges5min(i+1)+1];
        kR = kR+1;
    end
    if PrimaryBehavior5min(Edges5min(i)+1) == 3
        Epochs5min{3}(kD,:) = [Edges5min(i)+1 Edges5min(i+1)+1];
        kD = kD+1;
    end
    if PrimaryBehavior5min(Edges5min(i)+1) == 4
        Epochs5min{4}(kQ,:) = [Edges5min(i)+1 Edges5min(i+1)+1];
        kQ = kQ+1;
    end
end
if kF == 1; Epochs5min{1} = []; end
if kR == 1; Epochs5min{2} = []; end
if kD == 1; Epochs5min{3} = []; end
if kQ == 1; Epochs5min{4} = []; end
%%%%
PrimaryBehavior10min = zeros(1,length(RT))*NaN;
for i = 1:length(RT)
    if isfinite(BehaviorPrc10min(:,i))
        PrimaryBehavior10min(i) = find(BehaviorPrc10min(:,i) == max(BehaviorPrc10min(:,i)),1,'first');
    end
end
Edges10min = (find(diff(PrimaryBehavior10min)));
Epochs10min = {};
kF = 1; kR = 1; kD = 1; kQ = 1;
for i = 1:length(Edges10min)-1
    if PrimaryBehavior10min(Edges10min(i)+1) == 1
        Epochs10min{1}(kF,:) = [Edges10min(i)+1 Edges10min(i+1)+1];
        kF = kF+1;
    end
    if PrimaryBehavior10min(Edges10min(i)+1) == 2
        Epochs10min{2}(kR,:) = [Edges10min(i)+1 Edges10min(i+1)+1];
        kR = kR+1;
    end
    if PrimaryBehavior10min(Edges10min(i)+1) == 3
        Epochs10min{3}(kD,:) = [Edges10min(i)+1 Edges10min(i+1)+1];
        kD = kD+1;
    end
    if PrimaryBehavior10min(Edges10min(i)+1) == 4
        Epochs10min{4}(kQ,:) = [Edges10min(i)+1 Edges10min(i+1)+1];
        kQ = kQ+1;
    end
end
if kF == 1; Epochs10min{1} = []; end
if kR == 1; Epochs10min{2} = []; end
if kD == 1; Epochs10min{3} = []; end
if kQ == 1; Epochs10min{4} = []; end
%%%%%%
PrimaryBehavior20min = zeros(1,length(RT))*NaN;
for i = 1:length(RT)
    if isfinite(BehaviorPrc20min(:,i))
        PrimaryBehavior20min(i) = find(BehaviorPrc20min(:,i) == max(BehaviorPrc20min(:,i)),1,'first');
    end
end
Edges20min = (find(diff(PrimaryBehavior20min)));
Epochs20min = {};
kF = 1; kR = 1; kD = 1; kQ = 1;
for i = 1:length(Edges20min)-1
    if PrimaryBehavior20min(Edges20min(i)+1) == 1
        Epochs20min{1}(kF,:) = [Edges20min(i)+1 Edges20min(i+1)+1];
        kF = kF+1;
    end
    if PrimaryBehavior20min(Edges20min(i)+1) == 2
        Epochs20min{2}(kR,:) = [Edges20min(i)+1 Edges20min(i+1)+1];
        kR = kR+1;
    end
    if PrimaryBehavior20min(Edges20min(i)+1) == 3
        Epochs20min{3}(kD,:) = [Edges20min(i)+1 Edges20min(i+1)+1];
        kD = kD+1;
    end
    if PrimaryBehavior20min(Edges20min(i)+1) == 4
        Epochs20min{4}(kQ,:) = [Edges20min(i)+1 Edges20min(i+1)+1];
        kQ = kQ+1;
    end
end
if kF == 1; Epochs20min{1} = []; end
if kR == 1; Epochs20min{2} = []; end
if kD == 1; Epochs20min{3} = []; end
if kQ == 1; Epochs20min{4} = []; end

%%

figure;hold on;
for i = 1:getsize(size(Epochs1min{1}))
    rectangle('Position',[(Epochs1min{1}(i,1))/10/60/60,0,(Epochs1min{1}(i,2)-Epochs1min{1}(i,1))/10/60/60,1],'FaceColor','g','EdgeColor','g')
end
for i = 1:getsize(size(Epochs1min{2}))
    rectangle('Position',[(Epochs1min{2}(i,1))/10/60/60,0,(Epochs1min{2}(i,2)-Epochs1min{2}(i,1))/10/60/60,1],'FaceColor','r','EdgeColor','r')
end
for i = 1:getsize(size(Epochs1min{3}))
    rectangle('Position',[(Epochs1min{3}(i,1))/10/60/60,0,(Epochs1min{3}(i,2)-Epochs1min{3}(i,1))/10/60/60,1],'FaceColor','k','EdgeColor','k')
end
for i = 1:getsize(size(Epochs1min{4}))
    rectangle('Position',[(Epochs1min{4}(i,1))/10/60/60,0,(Epochs1min{4}(i,2)-Epochs1min{4}(i,1))/10/60/60,1],'FaceColor','b','EdgeColor','b')
end
ylim([-0.5 1.5]);
plot((RT)/10/60/60,QFSW(1,:),'r');
plot((RT)/10/60/60,QFSW(2,:),'m');
xlabel('Time (hr)')
title([Worm ' Behavioral Epochs (1min)']);
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' BehaviorEpochs1min.fig']);
saveas(gcf,[Worm ' BehaviorEpochs1min.pdf']);
close;
%%%%%%%%%%%%%%%%%
figure;hold on;
for i = 1:getsize(size(Epochs2min{1}))
    rectangle('Position',[(Epochs2min{1}(i,1))/10/60/60,0,(Epochs2min{1}(i,2)-Epochs2min{1}(i,1))/10/60/60,1],'FaceColor','g','EdgeColor','g')
end
for i = 1:getsize(size(Epochs2min{2}))
    rectangle('Position',[(Epochs2min{2}(i,1))/10/60/60,0,(Epochs2min{2}(i,2)-Epochs2min{2}(i,1))/10/60/60,1],'FaceColor','r','EdgeColor','r')
end
for i = 1:getsize(size(Epochs2min{3}))
    rectangle('Position',[(Epochs2min{3}(i,1))/10/60/60,0,(Epochs2min{3}(i,2)-Epochs2min{3}(i,1))/10/60/60,1],'FaceColor','k','EdgeColor','k')
end
for i = 1:getsize(size(Epochs2min{4}))
    rectangle('Position',[(Epochs2min{4}(i,1))/10/60/60,0,(Epochs2min{4}(i,2)-Epochs2min{4}(i,1))/10/60/60,1],'FaceColor','b','EdgeColor','b')
end
ylim([-0.5 1.5]);
plot((RT)/10/60/60,QFSW(1,:),'r');
plot((RT)/10/60/60,QFSW(2,:),'m');
xlabel('Time (hr)')
title([Worm ' Behavioral Epochs (2min)']);
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' BehaviorEpochs2min.fig']);
saveas(gcf,[Worm ' BehaviorEpochs2min.pdf']);
close;
%%%%%%%%%%%%%
figure;hold on;
for i = 1:getsize(size(Epochs5min{1}))
    rectangle('Position',[(Epochs5min{1}(i,1))/10/60/60,0,(Epochs5min{1}(i,2)-Epochs5min{1}(i,1))/10/60/60,1],'FaceColor','g','EdgeColor','g')
end
for i = 1:getsize(size(Epochs5min{2}))
    rectangle('Position',[(Epochs5min{2}(i,1))/10/60/60,0,(Epochs5min{2}(i,2)-Epochs5min{2}(i,1))/10/60/60,1],'FaceColor','r','EdgeColor','r')
end
for i = 1:getsize(size(Epochs5min{3}))
    rectangle('Position',[(Epochs5min{3}(i,1))/10/60/60,0,(Epochs5min{3}(i,2)-Epochs5min{3}(i,1))/10/60/60,1],'FaceColor','k','EdgeColor','k')
end
for i = 1:getsize(size(Epochs5min{4}))
    rectangle('Position',[(Epochs5min{4}(i,1))/10/60/60,0,(Epochs5min{4}(i,2)-Epochs5min{4}(i,1))/10/60/60,1],'FaceColor','b','EdgeColor','b')
end
ylim([-0.5 1.5]);
plot((RT)/10/60/60,QFSW(1,:),'r');
plot((RT)/10/60/60,QFSW(2,:),'m');
xlabel('Time (hr)')
title([Worm ' Behavioral Epochs (5min)']);
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' BehaviorEpochs5min.fig']);
saveas(gcf,[Worm ' BehaviorEpochs5min.pdf']);
close;
%%%%%%%%%%%%%%
figure;hold on;
for i = 1:getsize(size(Epochs10min{1}))
    rectangle('Position',[(Epochs10min{1}(i,1))/10/60/60,0,(Epochs10min{1}(i,2)-Epochs10min{1}(i,1))/10/60/60,1],'FaceColor','g','EdgeColor','g')
end
for i = 1:getsize(size(Epochs10min{2}))
    rectangle('Position',[(Epochs10min{2}(i,1))/10/60/60,0,(Epochs10min{2}(i,2)-Epochs10min{2}(i,1))/10/60/60,1],'FaceColor','r','EdgeColor','r')
end
for i = 1:getsize(size(Epochs10min{3}))
    rectangle('Position',[(Epochs10min{3}(i,1))/10/60/60,0,(Epochs10min{3}(i,2)-Epochs10min{3}(i,1))/10/60/60,1],'FaceColor','k','EdgeColor','k')
end
for i = 1:getsize(size(Epochs10min{4}))
    rectangle('Position',[(Epochs10min{4}(i,1))/10/60/60,0,(Epochs10min{4}(i,2)-Epochs10min{4}(i,1))/10/60/60,1],'FaceColor','b','EdgeColor','b')
end
ylim([-0.5 1.5]);
plot((RT)/10/60/60,QFSW(1,:),'r');
plot((RT)/10/60/60,QFSW(2,:),'m');
xlabel('Time (hr)')
title([Worm ' Behavioral Epochs (10min)']);
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' BehaviorEpochs10min.fig']);
saveas(gcf,[Worm ' BehaviorEpochs10min.pdf']);
close;
%%%%%%%%%%%%%
figure;hold on;
for i = 1:getsize(size(Epochs20min{1}))
    rectangle('Position',[(Epochs20min{1}(i,1))/10/60/60,0,(Epochs20min{1}(i,2)-Epochs20min{1}(i,1))/10/60/60,1],'FaceColor','g','EdgeColor','g')
end
for i = 1:getsize(size(Epochs20min{2}))
    rectangle('Position',[(Epochs20min{2}(i,1))/10/60/60,0,(Epochs20min{2}(i,2)-Epochs20min{2}(i,1))/10/60/60,1],'FaceColor','r','EdgeColor','r')
end
for i = 1:getsize(size(Epochs20min{3}))
    rectangle('Position',[(Epochs20min{3}(i,1))/10/60/60,0,(Epochs20min{3}(i,2)-Epochs20min{3}(i,1))/10/60/60,1],'FaceColor','k','EdgeColor','k')
end
for i = 1:getsize(size(Epochs20min{4}))
    rectangle('Position',[(Epochs20min{4}(i,1))/10/60/60,0,(Epochs20min{4}(i,2)-Epochs20min{4}(i,1))/10/60/60,1],'FaceColor','b','EdgeColor','b')
end
ylim([-0.5 1.5]);
plot((RT)/10/60/60,QFSW(1,:),'r');
plot((RT)/10/60/60,QFSW(2,:),'m');
xlabel('Time (hr)')
title([Worm ' Behavioral Epochs (20min)']);
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' BehaviorEpochs20min.fig']);
saveas(gcf,[Worm ' BehaviorEpochs20min.pdf']);
close;
%%
save([Worm ' Epochs'], 'frameRate','scrsz','Worm','SegNum','SegLngth','NumFrames','FRAMES','SEGMENTS','frames','RT','RT_L','IN','OUT',...
                                    'Epochs1min','Epochs2min','Epochs5min','Epochs10min','Epochs20min',...
                                    'BehaviorPrc1min','BehaviorPrc2min','BehaviorPrc5min','BehaviorPrc10min','BehaviorPrc20min',...
                                    'EdgePrcSW1min','EdgePrcSW2min','EdgePrcSW5min','EdgePrcSW10min','EdgePrcSW20min');
