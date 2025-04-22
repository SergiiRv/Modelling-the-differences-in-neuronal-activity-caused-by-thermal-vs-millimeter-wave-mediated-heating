% REMEMBER to update the Ksrink_coefCa range of processing in the SliseVec

clear Probabilities ProbClass ProbSort 
Time_preserved = T;
T = [1000, 1500, 2000, 3000];
Probabilities = zeros(indexx, length(T));
ProbSort = zeros(indexx, length(T)+1);
ProbClass = cell(indexx,1);
% calculating the probabilities along the time bids with relation to the
% first non affected region of 1 sec (or 2 sec)
for j =1:indexx
    Probabilities(j,:) = BineriesProb(VlocsM(:,j), VlocsM(:,j), T);%BineriesProb(VlocsM(:,j), V_ISI_M(:,j), T);
end;

ProbSort(:,1) = Probabilities(:,1);
ProbSort(:,2:end) = Probabilities;
% formation of coloring
for j=1:indexx
    if Probabilities(j,1)<Probabilities(j,length(T))
        ProbClass(j) = {'r'};
    elseif Probabilities(j,1)==Probabilities(j,length(T))
        ProbClass(j) = {'b'};%{'k'};
    else
        ProbClass(j) = {'b'};
    end
end;
for j=1:indexx
    if sum(Probabilities(j,(end-1):end))==0
        ProbClass(j) = {'k'};
    end
end;

T_=zeros(length(T)+1,1); T_(1) = 0; T_(2:end) = T;

figure('Name','Probabilities over all data ','NumberTitle','off')
for j=1:indexx
    plot(T_,ProbSort(j,:), ProbClass{j});
    hold on
end
hold off
title('AP probabilities at different transient parameters changes')
xlabel('time, msec'); set(gcf,'color','w'); set(gca,'fontsize',16); ylabel('Probability, a.u.');
saveas(gcf,'C:\the_model_2\Probabilities over all data'); 
saveas(gcf,'C:\the_model_2\Probabilities over all data.jpg');

figure('Name','Probabilities over all set of data (Normalized)','NumberTitle','off')
for j=1:indexx
    plot(T_,ProbSort(j,:)./ProbSort(j,1), ProbClass{j});
    hold on
end
hold off
title('AP probabilities at different transient parameters changes')
xlabel('time, msec'); set(gcf,'color','w'); set(gca,'fontsize',16); ylabel('Probability, a.u.');
saveas(gcf,'C:\the_model_2\Probabilities over all set of data (Normalized)'); 
saveas(gcf,'C:\the_model_2\Probabilities over all set of data (Normalized).jpg');

%formation of the indexing vectors
ProbClassD = zeros(length(ProbClass),1);
ProbClassSortColTriad = zeros(length(ProbClass),3);
for j=1:indexx
    if ProbClass{j}=='b'
        ProbClassD(j) = 1;
        ProbClassSortColTriad (j,:) = [0 0 1];
    elseif ProbClass{j}=='r'
        ProbClassD(j) = 2;
        ProbClassSortColTriad (j,:) = [1 0 0];
    else
        ProbClassD(j) = 0;
        ProbClassSortColTriad (j,:) = [0 1 0];
    end
end;
% indexing vector for sorting
ProbClassSortD = ProbClassD;
[ProbClassSortOut,ProbClassSortIndex] =  sort(ProbClassSortD,'descend');
ProbClassSortColTriadSorted = ProbClassSortColTriad(ProbClassSortIndex,:);
ProbClassSort = ProbClass(ProbClassSortIndex);
PCSO2min = min(find(ProbClassSortOut==2));
PCSO2max = max(find(ProbClassSortOut==2));
PCSO1min = min(find(ProbClassSortOut==1));
PCSO1max = max(find(ProbClassSortOut==1));
PCSO0min = min(find(ProbClassSortOut==0));
PCSO0max = max(find(ProbClassSortOut==0));

ParamMxSort = ParamMx(ProbClassSortIndex,:);
% here is the devision *100% and /2 sec to get the persentage rate %/sec of
% change
figure('Name','AP probabilities at different transient parameters changes throughout the paramrun','NumberTitle','off')
scatter3(ParamMxSort(:,1).*100, ParamMxSort(:,2).*100, ParamMxSort(:,3).*100,15,ProbClassSortColTriadSorted);
title('AP probabilities at different transient parameters changes throughout the paramrun')
xlabel('IK(Ca), %'); set(gcf,'color','w'); set(gca,'fontsize',16); ylabel('Na/K ATP, %'); zlabel('Ca-ATP, %');
saveas(gcf,'C:\the_model_2\Probabilities over all set of data 3D'); 
saveas(gcf,'C:\the_model_2\Probabilities over all set of data 3D.jpg');

SliseVec = [0:0.2:2];%[-0.1:-0.1:-1];%[2:1:12];   %Ksrink_coefCa;
figure('Name','AP probabilities at different transient parameters changes throughout the paramrun','NumberTitle','off')
for j=1:2:length(SliseVec)
    matrix = ParamMx(find(ParamMx(:,3)==SliseVec(j)),:);
    x = matrix(:,1).*100;
    y = matrix(:,2).*100;
    z = ProbClassD(find(ParamMx(:,3)==SliseVec(j))); %.*10;
    xi=linspace(min(x),max(x),100);
    yi=linspace(min(y),max(y),100);
    [XI YI]=meshgrid(xi,yi);
    ZI = griddata(x,y,z,XI,YI,'cubic'); %'cubic'
    ZI(1,1) = 2; ZI(1,100) = 1; ZI(100,1) = 0; 
    subplot(length(SliseVec),1,j);
    contourf(XI,YI,ZI);
    colorbar('Ticks',[0,1,2], 'TickLabels',{'AP.gen.fail','P.Dcrease','P.Increase'})
    colormap(jet); %colormap(prism(4))map; 
    title(['Ca-ATP = ', num2str(SliseVec(j)*100), '%']);
    %colorbar('Ticks',[0,1,2], 'TickLabels',{'AP.gen.fail','P.Dcrease','P.Increase'})
    xlabel('IK(Ca), %'); ylabel('Na/K ATP, %');
end;
set(gcf,'color','w');
saveas(gcf,'C:\the_model_2\AP probabilities at different transient parameters changes throughout the paramrun'); 
saveas(gcf,'C:\the_model_2\AP probabilities at different transient parameters changes throughout the paramrun.jpg');
% ParamMxSort2 = ParamMx(ProbClassSortIndex(PCSO2min:PCSO2max),:);
% ProbClassSort2 
% reprocessing of the ISI (because it got negative)

V_ISI_M_mod = (abs(V_ISI_M)).^-1;
BinneryOutF_mod = BinneryOutF;
BinneryOutI_mod = BinneryOutI;
BinneryOutHW_mod = BinneryOutHW;
BinneryOutAm_mod = BinneryOutAm;
JujeOut_mod = JujeOut;
deltaD = 50; %250;
for j = 1:indexx
    BinneryOutF_mod(j,:) = Bineries_mod(VlocsM(:,j), V_ISI_M_mod(:,j), T+deltaD);
    BinneryOutI_mod(j,:) = Bineries_mod(IlocsM(:,j), IpromsM(:,j), T+deltaD);
    BinneryOutHW_mod(j,:) = Bineries_mod(VlocsM(:,j), VwidthsM(:,j), T+deltaD);
    BinneryOutAm_mod(j,:) = Bineries_mod(VlocsM(:,j), VpromsM(:,j), T+deltaD);
    JujeOut_mod(j,:) = [Clesifyier_mod(BinneryOutF_mod(j,:)), Clesifyier_mod(BinneryOutI_mod(j,:)), Clesifyier_mod(BinneryOutHW_mod(j,:)), Clesifyier_mod(BinneryOutAm_mod(j,:))]; 
end;
 
% assingn the slising vector 
% here it is the temp rate

%SliseVec = HeatRate;

map = [0,0,1;
       0,1,1;
       0,1,0;
       1,1,0;
       1,0,0];
% BinniesF
figure('Name','Inter-Spike Interval model parameter sensitivity','NumberTitle','off')
for j=1:2:length(SliseVec)
    matrix = ParamMx(find(ParamMx(:,3)==SliseVec(j)),:);
    x = matrix(:,1).*100;
    y = matrix(:,2).*100;
    z = JujeOut_mod(find(ParamMx(:,3)==SliseVec(j)),1); %.*10;
    xi=linspace(min(x),max(x),100);
    yi=linspace(min(y),max(y),100);
    [XI YI]=meshgrid(xi,yi);
    ZI = griddata(x,y,z,XI,YI,'cubic'); %'cubic'
    %ZI(1,1) = -2; ZI(1,100) = -1; ZI(100,1) = 0; ZI(100,100) = 1; ZI(100,99) = 2;
    subplot(length(SliseVec),1,j);
    contourf(XI,YI,ZI);
    colormap(jet); %colormap(prism(4))map; 
    title(['Ca-ATP = ', num2str(SliseVec(j)*100), '%']);
    colorbar('Ticks',[-2,-0.5,1,2], 'TickLabels',{'Mon.Decrease','NoChange','Bi-Phasic','Mon.Increase'})
    caxis([-2 2]); xlabel('IK(Ca), %'); ylabel('Na/K ATP, %');
end;
set(gcf,'color','w');
saveas(gcf,'C:\the_model_2\Inter-Spike Interval model parameter sensitivity'); 
saveas(gcf,'C:\the_model_2\Inter-Spike Interval model parameter sensitivity.jpg');
% BinniesI
figure('Name','AP current Am model parameter sensitivity','NumberTitle','off')
for j=1:2:length(SliseVec)
    matrix = ParamMx(find(ParamMx(:,3)==SliseVec(j)),:);
    x = matrix(:,1).*100;
    y = matrix(:,2).*100;
    z = JujeOut_mod(find(ParamMx(:,3)==SliseVec(j)),2); %.*10;
    xi=linspace(min(x),max(x),100);
    yi=linspace(min(y),max(y),100);
    [XI YI]=meshgrid(xi,yi);
    ZI = griddata(x,y,z,XI,YI,'cubic'); %'cubic'
    %ZI(1,1) = -2; ZI(1,100) = -1; ZI(100,1) = 0; ZI(100,100) = 1; ZI(100,99) = 2;
    subplot(length(SliseVec),1,j);
    contourf(XI,YI,ZI);
    colormap(jet); %colormap(prism(4))map; 
    title(['Ca-ATP = ', num2str(SliseVec(j)*100), '%']);
    colorbar('Ticks',[-2,-0.5,1,2], 'TickLabels',{'Mon.Decrease','NoChange','Bi-Phasic','Mon.Increase'})
    caxis([-2 2]); xlabel('IK(Ca), %'); ylabel('Na/K ATP, %');
end;
set(gcf,'color','w');
saveas(gcf,'C:\the_model_2\AP current Am model parameter sensitivity'); 
saveas(gcf,'C:\the_model_2\AP current Am model parameter sensitivity.jpg');
% BinniesHW
figure('Name','AP Half-Width model parameter sensitivity','NumberTitle','off')
for j=1:2:length(SliseVec)
    matrix = ParamMx(find(ParamMx(:,3)==SliseVec(j)),:);
    x = matrix(:,1).*100;
    y = matrix(:,2).*100;
    z = JujeOut_mod(find(ParamMx(:,3)==SliseVec(j)),3); %.*10;
    xi=linspace(min(x),max(x),100);
    yi=linspace(min(y),max(y),100);
    [XI YI]=meshgrid(xi,yi);
    ZI = griddata(x,y,z,XI,YI,'cubic'); %'cubic'
    %ZI(1,1) = -2; ZI(1,100) = -1; ZI(100,1) = 0; ZI(100,100) = 1; ZI(100,99) = 2;
    subplot(length(SliseVec),1,j);
    contourf(XI,YI,ZI);
    colormap(jet); %colormap(prism(4))map; 
    title(['Ca-ATP = ', num2str(SliseVec(j)*100), '%']);
    colorbar('Ticks',[-2,-0.5,1,2], 'TickLabels',{'Mon.Decrease','NoChange','Bi-Phasic','Mon.Increase'})
    caxis([-2 2]); xlabel('IK(Ca), %'); ylabel('Na/K ATP, %');
end;
set(gcf,'color','w');
saveas(gcf,'C:\the_model_2\AP Half-Width model parameter sensitivity'); 
saveas(gcf,'C:\the_model_2\AP Half-Width model parameter sensitivity.jpg');
% BinniesAm
figure('Name','AP Amplitude model parameter sensitivity','NumberTitle','off')
for j=1:2:length(SliseVec)
    matrix = ParamMx(find(ParamMx(:,3)==SliseVec(j)),:);
    x = matrix(:,1).*100;
    y = matrix(:,2).*100;
    z = JujeOut_mod(find(ParamMx(:,3)==SliseVec(j)),4); %.*10;
    xi=linspace(min(x),max(x),100);
    yi=linspace(min(y),max(y),100);
    [XI YI]=meshgrid(xi,yi);
    ZI = griddata(x,y,z,XI,YI,'cubic'); %'cubic'
    %ZI(1,1) = -2; ZI(1,100) = -1; ZI(100,1) = 0; ZI(100,100) = 1; ZI(100,99) = 2;
    subplot(length(SliseVec),1,j);
    contourf(XI,YI,ZI);
    colormap(jet); %colormap(prism(4))map; 
    title(['Ca-ATP = ', num2str(SliseVec(j)*100), '%']);
    colorbar('Ticks',[-2,-0.5,1,2], 'TickLabels',{'Mon.Decrease','NoChange','Bi-Phasic','Mon.Increase'})
    caxis([-2 2]); xlabel('IK(Ca), %'); ylabel('Na/K ATP, %');
end;
set(gcf,'color','w');
saveas(gcf,'C:\the_model_2\AP Amplitude model parameter sensitivity'); 
saveas(gcf,'C:\the_model_2\AP Amplitude model parameter sensitivity.jpg');

%__________________________________________________________________________
%__________________________________________________________________________
JujeOut_Aver = JujeOut;
Zup = 2;
deltaD = 50;
for j = 1:indexx
    MeanOutF = AverBineriesProb(VlocsM(:,j), V_ISI_M_mod(:,j), T+deltaD);
    MeanOutI = AverBineriesProb(IlocsM(:,j), IpromsM(:,j), T+deltaD);
    MeanOutHW = AverBineriesProb(VlocsM(:,j), VwidthsM(:,j), T+deltaD);
    MeanOutAm = AverBineriesProb(VlocsM(:,j), VpromsM(:,j), T+deltaD);
    JujeOut_Aver(j,:) = [MeanOutF, MeanOutI, MeanOutHW, MeanOutAm]; 
end;

% BinniesF
figure('Name','Inter-Spike Interval Aver. Value model parameter sensitivity','NumberTitle','off')
for j=1:2:length(SliseVec)
    matrix = ParamMx(find(ParamMx(:,3)==SliseVec(j)),:);
    x = matrix(:,1).*100;
    y = matrix(:,2).*100;
    z = JujeOut_Aver(find(ParamMx(:,3)==SliseVec(j)),1).*100; %.*10;
    xi=linspace(min(x),max(x),100);
    yi=linspace(min(y),max(y),100);
    [XI YI]=meshgrid(xi,yi);
    ZI = griddata(x,y,z,XI,YI,'cubic'); %'cubic'
     
    subplot(length(SliseVec),1,j);
    contourf(XI,YI,ZI);
    colormap(jet); %colormap(prism(4))map; 
    title(['Ca-ATP = ', num2str(SliseVec(j)*100), '%']);
    colorbar%('Ticks',[-2,-0.5,1,2], 'TickLabels',{'Mon.Decrease','NoChange','Bi-Phasic','Mon.Increase'})
    caxis([300 2700]); xlabel('IK(Ca), %'); ylabel('Na/K ATP, %');
end;
set(gcf,'color','w');
saveas(gcf,'C:\the_model_2\Inter-Spike Interval Aver Value model parameter sensitivity'); 
saveas(gcf,'C:\the_model_2\Inter-Spike Interval Aver Value model parameter sensitivity.jpg');
% BinniesI
figure('Name','AP current Am Aver. Value model parameter sensitivity','NumberTitle','off')
for j=1:2:length(SliseVec)
    matrix = ParamMx(find(ParamMx(:,3)==SliseVec(j)),:);
    x = matrix(:,1).*100;
    y = matrix(:,2).*100;
    z = JujeOut_Aver(find(ParamMx(:,3)==SliseVec(j)),2).*100; %.*10;
    xi=linspace(min(x),max(x),100);
    yi=linspace(min(y),max(y),100);
    [XI YI]=meshgrid(xi,yi);
    ZI = griddata(x,y,z,XI,YI,'cubic'); %'cubic'
    
    subplot(length(SliseVec),1,j);
    contourf(XI,YI,ZI);
    colormap(jet); %colormap(prism(4))map; 
    title(['Ca-ATP = ', num2str(SliseVec(j)*100), '%']);
    colorbar%('Ticks',[-2,-0.5,1,2], 'TickLabels',{'Mon.Decrease','NoChange','Bi-Phasic','Mon.Increase'})
    caxis([160 210]); xlabel('IK(Ca), %'); ylabel('Na/K ATP, %');
end;
set(gcf,'color','w');
saveas(gcf,'C:\the_model_2\AP current Am Aver Value model parameter sensitivity'); 
saveas(gcf,'C:\the_model_2\AP current Am Aver Value model parameter sensitivity.jpg');
% BinniesHW
figure('Name','AP Half-Width Aver. Value model parameter sensitivity','NumberTitle','off')
for j=1:2:length(SliseVec)
    matrix = ParamMx(find(ParamMx(:,3)==SliseVec(j)),:);
    x = matrix(:,1).*100;
    y = matrix(:,2).*100;
    z = JujeOut_Aver(find(ParamMx(:,3)==SliseVec(j)),3).*100; %.*10;
    xi=linspace(min(x),max(x),100);
    yi=linspace(min(y),max(y),100);
    [XI YI]=meshgrid(xi,yi);
    ZI = griddata(x,y,z,XI,YI,'cubic'); %'cubic'
    
    subplot(length(SliseVec),1,j);
    contourf(XI,YI,ZI);
    colormap(jet); %colormap(prism(4))map; 
    title(['Ca-ATP = ', num2str(SliseVec(j)*100), '%']);
    colorbar%('Ticks',[-2,-0.5,1,2], 'TickLabels',{'Mon.Decrease','NoChange','Bi-Phasic','Mon.Increase'})
    caxis([35 40]);xlabel('IK(Ca), %'); ylabel('Na/K ATP, %');
end;
set(gcf,'color','w');
saveas(gcf,'C:\the_model_2\AP Half-Width Aver Value model parameter sensitivity'); 
saveas(gcf,'C:\the_model_2\AP Half-Width Aver Value model parameter sensitivity.jpg');
% BinniesAm
figure('Name','AP Amplitude Aver. Value model parameter sensitivity','NumberTitle','off')
for j=1:2:length(SliseVec)
    matrix = ParamMx(find(ParamMx(:,3)==SliseVec(j)),:);
    x = matrix(:,1).*100;
    y = matrix(:,2).*100;
    z = JujeOut_Aver(find(ParamMx(:,3)==SliseVec(j)),4).*100; %.*10;
    xi=linspace(min(x),max(x),100);
    yi=linspace(min(y),max(y),100);
    [XI YI]=meshgrid(xi,yi);
    ZI = griddata(x,y,z,XI,YI,'cubic'); %'cubic'
    
    subplot(length(SliseVec),1,j);
    contourf(XI,YI,ZI);
    colormap(jet); %colormap(prism(4))map; 
    title(['Ca-ATP = ', num2str(SliseVec(j)*100), '%']);
    colorbar%('Ticks',[-2,-0.5,1,2], 'TickLabels',{'Mon.Decrease','NoChange','Bi-Phasic','Mon.Increase'})
    caxis([85 125]); xlabel('IK(Ca), %'); ylabel('Na/K ATP, %');
end;
set(gcf,'color','w');
saveas(gcf,'C:\the_model_2\AP Amplitude Aver Value model parameter sensitivity'); 
saveas(gcf,'C:\the_model_2\AP Amplitude Aver Value model parameter sensitivity.jpg');