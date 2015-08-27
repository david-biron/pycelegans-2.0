% frameRate = 10;
%%
cJet = colormap(jet(SegNum));
c1b = colormap(winter(SegNum));
c1c = colormap(autumn(SegNum));
close all;
scrsz = get(0,'screensize');
HourXaxis = RT/frameRate/60/60;
getsize = @(x)x(1);
%% Plot Angles
figure;hold on;
for i = SEGMENTS
    plot((frames(FRAMES)),AHA_S05{5}(i,FRAMES),'Color',cJet(i,:))
    hold on;
end
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' Angles.fig']);
close;

%% Bends
figure;hold on;
for i = 1:length(Waves1)
   plot(Waves1{i}.T,Waves1{i}.S,'k.','MarkerSize',10); plot(Waves1{i}.T(1),Waves1{i}.S(1),'k.','MarkerSize',30);
   plot(Waves1{i}.FRWT,Waves1{i}.FRWS,'g.','MarkerSize',40);
   plot(Waves1{i}.REVT,Waves1{i}.REVS,'r.','MarkerSize',30);
   plot(Waves1{i}.QuiescT,Waves1{i}.QuiescS,'c.','MarkerSize',20);
end
for i = 1:length(Waves2)
    plot(Waves2{i}.T,Waves2{i}.S,'k.','MarkerSize',10); plot(Waves2{i}.T(1),Waves2{i}.S(1),'k.','MarkerSize',30);
    plot(Waves2{i}.FRWT,Waves2{i}.FRWS,'g.','MarkerSize',40);
    plot(Waves2{i}.REVT,Waves2{i}.REVS,'r.','MarkerSize',30);
    plot(Waves2{i}.QuiescT,Waves2{i}.QuiescS,'c.','MarkerSize',20);
end
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' Waves.fig']);
close

%% Quiescence Analysis
figure;hold on;
plot(HourXaxis,QFSW(1,:),'k')
plot(HourXaxis,QFSW(2,:),'b')
ylim([0 1]);
xlabel('Time (hr)')
ylabel('Quiescence Fraction')
title([Worm ' Quiescence Fraction']);
legend('Standard','Body')
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' QF.fig']);
saveas(gcf,[Worm ' QF.pdf']);
close;

figure;hold on;
plot(HourXaxis,QBoutDurSW(1,:),'k')
plot(HourXaxis,QBoutDurSW(2,:),'b')
xlabel('Time (hr)')
ylabel('Qbout Duration (sec)')
title([Worm ' Qbout Duration']);
legend('Standard','Body')
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' QboutDur.fig']);
saveas(gcf,[Worm ' QboutDur.pdf']);
close;

%%
figure;hold on;
plot(HourXaxis,BehaviorPrc10min(1,:),'g');
plot(HourXaxis,BehaviorPrc10min(2,:),'r');
plot(HourXaxis,BehaviorPrc10min(3,:),'k');
plot(HourXaxis,BehaviorPrc10min(4,:),'b');
ylim([0 1]);
xlabel('Time (hr)')
ylabel('Fraction')
title([Worm ' Behavior Percentage']);
legend('FRW','REV','DWELL','QUIESC')
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' BehaviorPrc10min.fig']);
saveas(gcf,[Worm ' BehaviorPrc10min.pdf']);
close;
%%

figure;hold on;
plot(frames(FRAMES)/frameRate/60/60,WormLengths(FRAMES),'k')
plot(HourXaxis,WL_Slid,'r');
xlabel('Time (hr)')
ylabel('Length (pixels)')
title([Worm ' Worm Length']);
set(gcf,'position',scrsz);
saveas(gcf,[Worm ' WormLength.fig']);
saveas(gcf,[Worm ' WormLength.pdf']);
close;
















