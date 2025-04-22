%perameters of the run
% for I-K-Ca
gkcabar_coef = [0.1:0.3:4.6]; %
% Na/K ATP
Ksrink_coefNaK = [0:0.2:2]; %
% for I-Ca
gcabar_coef = [0.5:0.5:2.5];
% for Ca-ATP
Ksrink_coefCa = [0:0.2:2]; %
% Heating Rate
HeatRate = [0.1:0.1:1];

% making the matrix of combinations PeramMX
ParamMx = ones(length(gkcabar_coef)*length(Ksrink_coefNaK)*length(Ksrink_coefCa), 3);
indexx = 0;
for i =1:length(gkcabar_coef)
    for j=1:length(Ksrink_coefNaK)
        for k =1:length(Ksrink_coefCa)
               indexx = indexx+1;
               ParamMx(indexx,1) = gkcabar_coef(i);
               ParamMx(indexx,2) = Ksrink_coefNaK(j);
               ParamMx(indexx,3) = Ksrink_coefCa(k);    
        end
    end
end;

BinneryOutF = zeros(length(ParamMx), 7);
BinneryOutI = zeros(length(ParamMx), 7);
BinneryOutHW = zeros(length(ParamMx), 7);
BinneryOutAm = zeros(length(ParamMx), 7);
JujeOut = zeros(length(ParamMx), 4);

VpksM = zeros(1000, indexx);
VlocsM = zeros(1000, indexx);
VwidthsM = zeros(1000, indexx);
VpromsM = zeros(1000, indexx);
IpksM = zeros(1000, indexx);
IlocsM = zeros(1000, indexx);
IwidthsM = zeros(1000, indexx);
IpromsM = zeros(1000, indexx);
V_ISI_M = zeros(1000, indexx);
%indexx = 5;
for i = 1:indexx
    %setting up the Flag file
    Flag = 0;
    'zero flag'
    fileID = fopen('C:\the_model_2\FlagFile.txt','w+');
    fprintf(fileID,'%e',Flag); %fprintf(fileID,'%d\r\n',0);
    fclose(fileID);
    type 'FlagFile.txt'
    % writing down the parameters 
    fileID = fopen('C:\the_model_2\ParameterFile.txt','w+');
    fprintf(fileID,'%1.2f\r\n',ParamMx(i,1:end));
    fclose(fileID);
    type 'C:\the_model_2\ParameterFile.txt'
    %!C:\nrn\bin\nrniv.exe -nogui C:\the_model\test-7-Transient_Temp_just_soma.hoc &
    %!C:\nrn\bin\nrniv.exe -nogui C:\the_model\FlagNRNoutput.hoc &
    winopen test-7-Transient_Temp_just_soma.hoc
    % control of execution and prevention of the hoc jamming in execution
    pause(1)
    
    'new flag'
    type 'FlagFile.txt'
    FListOpen=fopen('all'); %get the names of the open files
        while sum(ismember('C:\the_model_2\FlagFile.txt',FListOpen)) > 0
            pause(0.5)
        end;
    fileID = fopen('C:\the_model_2\FlagFile.txt','r');
    Flag = fscanf(fileID,'%e');
    fclose(fileID);
    while Flag == 0
        fileID = fopen('C:\the_model_2\FlagFile.txt','r');
        Flag = fscanf(fileID,'%e');
        fclose(fileID);
        if Flag == 0
            pause(1)
        end
    end; 
    pause(1)
    system('"C:\Windows\System32\taskkill.exe" /im cmd.exe &');
    % reading off the result
    'reading the file of data'
    % 1-time; 2-Volt; 3-temp; 4-I; 5-Leak
    DATASomaTemp = importfile('DATA_Soma_Temp.txt', 2, 400000);
    indext = find(DATASomaTemp(:,1)>0);
    [Vpks,Vlocs,Vwidths,Vproms] = findpeaks(DATASomaTemp(indext,2),DATASomaTemp(indext,1),'MinPeakProminence',5);
    [Ipks,Ilocs,Iwidths,Iproms] = findpeaks(DATASomaTemp(indext,4),DATASomaTemp(indext,1),'MinPeakProminence',0.05);
    'processing'
    % making the regions of calculus (overall one control region and 6 calc,=> 7 piints)
    dt = DATASomaTemp(10,1) - DATASomaTemp(9,1);
    T1= 1000;
    T2= 1200;
    T3= 1400;
    T4= 1600;
    T5= 1800;
    T6= 2000;
    T7= 3000;
    T = [T1, T2, T3, T4, T5, T6, T7];
    if isempty(Vlocs)==0
        %interspike interval
        V_ISI = zeros(1,length(Vlocs));
        V_ISI(1) = Vlocs(1);
        V_ISI(2:end) = Vlocs(1:end-1)-Vlocs(2:end);
        V_ISI(1) = V_ISI(2);
        % making the regions of calculus
        BinniesF = Bineries(Vlocs, V_ISI, T);
        BinniesI = Bineries(Vlocs, Iproms, T);
        BinniesHW = Bineries(Vlocs, Vwidths, T);
        BinniesAm = Bineries(Vlocs, Vproms, T);
    else
        BinniesF = zeros(1,7);
        BinniesI = zeros(1,7);
        BinniesHW = zeros(1,7);
        BinniesAm = zeros(1,7);
    end;
    'saving the data'
    VpksM(1:length(Vpks),i) = Vpks;
    VlocsM(1:length(Vlocs),i) = Vlocs;
    VwidthsM(1:length(Vwidths),i) = Vwidths;
    VpromsM(1:length(Vproms),i) = Vproms;
    IpksM(1:length(Ipks),i) = Ipks;
    IlocsM(1:length(Ilocs),i) = Ilocs;
    IwidthsM(1:length(Iwidths),i) = Iwidths;
    IpromsM(1:length(Iproms),i) = Iproms;
    V_ISI_M(1:length(V_ISI),i) = V_ISI';
    'classification done'
    BinneryOutF(i,:) = BinniesF;
    BinneryOutI(i,:) = BinniesI;
    BinneryOutHW(i,:) = BinniesHW;
    BinneryOutAm(i,:) = BinniesAm;
    JujeOut(i,:) = [Clesifyier(BinniesF), Clesifyier(BinniesI), Clesifyier(BinniesHW), Clesifyier(BinniesAm)];  
    
    % progress line
    Percent = i/indexx*100;
    barh(Percent, 'FaceColor','r');
    xlim([0 100]);
    
end;




