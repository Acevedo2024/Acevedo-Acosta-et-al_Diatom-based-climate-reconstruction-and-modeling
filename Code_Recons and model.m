%%   Diatom-based climate reconstruction and modeling: modern calibration set and paleodata
%    J.D. Acevedo-Acosta(1), A. Martínez López(1*), and L. Carro-Calvo(2) 
% 1. Instituto Politécnico Nacional, Centro Interdisciplinario de Ciencias Marinas–IPN (CICIMAR), Av. Instituto Politécnico Nacional s/n, Col. Playa Palo de Santa Rita, La Paz, B.C.S. 23096, México. 
% 2. Departamento de Teoría de la Señal y las Comunicaciones y Sistemas Telemáticos y Computación, Universidad Rey Juan Carlos, Madrid – España.  
%    Corresponding Author *: Aída Martínez López (amartin@ipn.mx) 

%% Load data 

Years = xlsread("Inst_Biol_DataReconst.xlsx","Years")'
x =     xlsread("Inst_Biol_DataReconst.xlsx","Input")'
t =     xlsread("Inst_Biol_DataReconst.xlsx","Target")'
x1 =    xlsread("Inst_Biol_DataReconst.xlsx","Sedcore")'
Ins =   xlsread("Inst_Biol_DataReconst.xlsx","Ints_Data")'

Model = struct ('Layers',{},'Arch',{},'Regr',{},'RMSEsst',{},'RMSEtemp',{},...
    'RMSEprec',{},'Rsst',{},'Rsstm',{},'Rtemp',{},'Rtempm',{},'RATS',{},'RATSm',{},...
    'RMiTS',{},'RMiTSm',{}, 'RMaTS',{},'RMaTSm',{},'RATW',{},'RATWm',{}, 'RMiTW',{},'RMiTWm',{},...
    'RMaTW',{},'RMaTWm',{});
RecoModel = struct ('Arch',{},'SST',{},'Temp',{},'Preci',{},'SSTm',{},'Temm',{},'Precm',{},'SSTm3',{},'Tempm3',{},'Precim3',{});


%% run neuronal network 

for k= 3:5

    Yrec= [];
    Ymod = [];
    B= [];
    C= [];
    D= [];
    RMSEMAT = []
    regre = []

    for i= 1:10
        trainFcn = 'trainbr';
        %net=feedforwardnet([k k],trainFcn);
        net = feedforwardnet(k,trainFcn);
        net.divideFcn = 'dividerand';
        net.divideMode = 'time';
        net.divideParam.trainRatio = 70/100;
        net.divideParam.valRatio = 15/100;
        net.divideParam.testRatio = 15/100;
        net.performFcn = 'mse';
        net = train(net,x,t);
        y = net(x);
        e = gsubtract(t,y);
        performance = perform(net,t,y)
        [r,m,b] = regression(t,y,'one');
       
        if r < 0.9
            trainFcn = 'trainbr';
            %net=feedforwardnet([k k],trainFcn);
            net = feedforwardnet(k,trainFcn);
            net.divideFcn = 'dividerand';
            net.divideMode = 'time';
            net.divideParam.trainRatio = 70/100;
            net.divideParam.valRatio = 15/100;
            net.divideParam.testRatio = 15/100;
            net.performFcn = 'mse';
            net = train(net,x,t);
            y = net(x);
            [r,m,b] = regression(t,y,'one');
        end
        regre = [regre; r];
        proxy= x1;
        y1=sim(net,proxy);
        Yrec= [Yrec ; y1];
        Ymod= [Ymod ; y];
        D= [D; r];
        B= [B; performance];
        C= [C; y];
        %% Calculate the RMSE of modeling
        for j = 1:3
            VC1=[]; 
            VC2=[] ;
            xi= t(j,:);
            xo= y(j,:);
            n= length(xi);
            for i = 1:87
                a=(xi(i)-xo(i))^2;
                VC1= [VC1; a];
            end
            VC1 = sum(VC1);
            d=(1/n)*VC1;
            ee= sum(sqrt(d));
            RMSEMAT = [RMSEMAT; ee];
        end
        RMSEtsm = mean(RMSEMAT([1:3:end],:));
        RMSEtemp = mean(RMSEMAT([2:3:end],:));
        RMSEpreci = mean(RMSEMAT([3:3:end],:));

    end

    %% Averages and Mobile averages of the variables
    SST = Yrec([1:3:end],:);
    Temp = Yrec([2:3:end],:);
    Preci = Yrec([3:3:end],:);
    SSTmean = mean(SST);
    Tempmean = mean(Temp);
    Precimean = mean(Preci);
    SST_movil3= movmean(SSTmean,3);
    Temp_movil3= movmean(Tempmean,3);
    Preci_movil3= movmean(Precimean,3);

    %% correlation between observed values ​​and expected
    %SST
    [RhoSST,PsSST] = corrcoef(SSTmean(1:1:21),Ins(2,1:21));
    [RhoSSTm,PsSSTm] = corrcoef(SST_movil3(1:1:21),Ins(2,1:21));
    %Temp
    [RhoTem,PsTem] = corrcoef(Tempmean(1:1:24),Ins(3,1:24));
    [RhoTemm,PsTemm] = corrcoef(Temp_movil3(1:1:24),Ins(3,1:24));
    %Temp_summer
    [RhoATS,PsATS] = corrcoef(Tempmean(1:1:24),Ins(4,1:24));
    [RhoATSm,PsATSm] = corrcoef(Temp_movil3(1:1:24),Ins(4,1:24));
    [RhoMiTS,PsMiTS] = corrcoef(Tempmean(1:1:24),Ins(5,1:24));
    [RhoMiTSm,PsMiTSm] = corrcoef(Temp_movil3(1:1:24),Ins(5,1:24));
    [RhoMaTS,PsMaTS] = corrcoef(Tempmean(1:1:24),Ins(6,1:24));
    [RhoMaTSm,PsMaTSm] = corrcoef(Temp_movil3(1:1:24),Ins(6,1:24));
    %Temp_Winter
    [RhoATW,PsATW] = corrcoef(Tempmean(1:1:24),Ins(7,1:24));
    [RhoATWm,PsATWm] = corrcoef(Temp_movil3(1:1:24),Ins(7,1:24));
    [RhoMiTW,PsMiTW] = corrcoef(Tempmean(1:1:24),Ins(8,1:24));
    [RhoMiTWm,PsMiTWm] = corrcoef(Temp_movil3(1:1:24),Ins(8,1:24));
    [RhoMaTW,PsMaTW] = corrcoef(Tempmean(1:1:24),Ins(9,1:24));
    [RhoMaTWm,PsMaTWm] = corrcoef(Temp_movil3(1:1:24),Ins(9,1:24));

    for i = k-2
        Model(i).Layers= 1
        %Modelos(i).Arch = k
        Model(i).Arch = [k k]
        Model(i).regr = mean(regre)
        Model(i).RMSEsst = mean(RMSEMAT([1:3:end],:));
        Model(i).RMSEtemp = mean(RMSEMAT([2:3:end],:));
        Model(i).RMSEprec = mean(RMSEMAT([3:3:end],:));
        %SST
        Model(i).Rsst = [RhoSST(1,2) PsSST(1,2)]
        Model(i).Rsstm = [RhoSSTm(1,2) PsSSTm(1,2)]
        %Temp
        Model(i).Rtemp = [RhoTem(1,2) PsTem(1,2)]
        Model(i).Rtempm = [RhoTemm(1,2) PsTemm(1,2)]
       %Temp_Summer
        Model(i).RATS = [RhoATS(1,2) PsATS(1,2)]
        Model(i).RATSm = [RhoATSm(1,2) PsATSm(1,2)]
        Model(i).RMiTS = [RhoMiTS(1,2) PsMiTS(1,2)]
        Model(i).RMiTSm = [RhoMiTSm(1,2) PsMiTSm(1,2)]
        Model(i).RMaTS = [RhoMaTS(1,2) PsMaTS(1,2)]
        Model(i).RMaTSm = [RhoMaTSm(1,2) PsMaTSm(1,2)] 
        %Temp_Winter
        Model(i).RATW = [RhoATW(1,2) PsATW(1,2)]
        Model(i).RATWm = [RhoATWm(1,2) PsATWm(1,2)]
        Model(i).RMiTW = [RhoMiTW(1,2) PsMiTW(1,2)]
        Model(i).RMiTWm = [RhoMiTWm(1,2) PsMiTWm(1,2)]
        Model(i).RMaTW = [RhoMaTW(1,2) PsMaTW(1,2)]
        Model(i).RMaTWm = [RhoMaTWm(1,2) PsMaTWm(1,2)] 
        
        RecoModel(i).Arch = k
        RecoModel(i).Arch = [k k]
        RecoModel(i).SST = SST
        RecoModel(i).SSTm = SSTmean
        RecoModel(i).SSTm3 = SST_movil3
        RecoModel(i).Temp = Temp
        RecoModel(i).Tempm = Tempmean
        RecoModel(i).Tempm3 = Temp_movil3
        RecoModel(i).Preci = Preci
        RecoModel(i).Preciprom = Precimean
        RecoModel(i).Precim3 = Preci_movil3

    end
end
save('filename.mat', 'Model')
save('RecModel.mat', 'RecoModel')

%% Graph the reconstruction results
amatrix = RecoModel(3).Temp;
figure
[aline(1), aFill(1)] = stdshade(amatrix,0.3,[0.4 0.4 0.4], 1:size(amatrix,2), 3);
hold on
ax = gca;
set(ax,'fontname','Times New Roman','fontsize',12,'TickDir','out')
ylim([15 37])
plot(Ins(10,:),'blue')
hold on
ylim([15 37])
plot(RecoModel(3).Tempm3(1:199),'red')
ax.XDir = 'reverse';
xticks(1:18:200)
xticklabels(Years)
xlabel('Years')
ax = gca;
set(ax,'fontname','Times New Roman','fontsize',12,'TickDir','out')
ylabel('Atm. Temperature (°C)')
print -djpeg -r600 -cmyk Pinstru5n

%% Graph the reconstruction results (1940-2011) 
Years_01= [2011 2005 1999 1993 1987 1981 1975 1969 1963 1957 1951 1945]
amatrix = RecoModel(3).SST;
amatrix = amatrix(1:end,1:24);
figure
[aline(1), aFill(1)] = stdshade(amatrix,0.3,[0.4 0.4 0.4], 1:size(amatrix,2), 3);
hold on
ax = gca;
set(ax,'fontname','Times New Roman','fontsize',12,'TickDir','out')
ylim([15 37])
plot(Ins(2,:),'blue')
hold on
ylim([15 37])
plot(RecoModel(3).SSTm3(1:24),'black')
ax.XDir = 'reverse';
xticks(1:2:24)
xticklabels(Years_01)
xlabel('Years')
ax = gca;
set(ax,'fontname','Times New Roman','fontsize',12,'TickDir','out')
ylabel('SST (°C)')
print -djpeg -r600 -cmyk Pinstru5n

%% Second part [Climate Modelling]

close all
clear all
clc 
% load the second file ("Data_Model.xlsx")
Years = xlsread("Data_Model.xlsx","Years")';
X = xlsread("Data_Model.xlsx","Input")';
T = xlsread("Data_Model.xlsx","Target")';
Ins = xlsread("Data_Model.xlsx","InstData")';
Ins = Ins(2,:);
Years = num2cell(Years);
Xo = num2cell(X);
To = num2cell(T);

Modpredic= struct('Layers',{},'Arch',{},'Regr',{},'perforModel',{},'PerforPredic',{},...
    'Predic',{},'Mean',{},'Meanmov5',{},'Corr',{},'Corrm',{});

for k = 3:10
    Tpredict = [];
    regre = [];
    Perfor = [];
    PerforP = [];
    for i= 1:6
        X=Xo;
        T=To;  
        trainFcn = 'trainbr'; 
        inputDelays = 1:6;    
        feedbackDelays = 1:6; 
        hiddenLayerSize =[k];
        %hiddenLayerSize =[k k];

        %% Model training
        net = narxnet(inputDelays,feedbackDelays,hiddenLayerSize,'open',trainFcn);
        [x,xi,ai,t] = preparets(net,X,{},T);
        net.divideFcn = 'dividerand';
        net.divideMode = 'time';
        net.divideParam.trainRatio = 70/100;
        net.divideParam.valRatio = 15/100;
        net.divideParam.testRatio = 15/100;
        net.performFcn = 'mse';
        [net,tr] = train(net,x,t,xi,ai);
        y = net(x,xi,ai);
        e = gsubtract(t,y);
        performance = perform(net,t,y);
        Perfor = [Perfor; performance];
        [r,m,b] = regression(t,y,'one');
        regre = [regre; r];
        for i= 1:8 % Far steps in prediction
            nets = removedelay(net);
            nets.name = [net.name ' - Predict One Step Ahead'];
            [xs,xis,ais,ts] = preparets(nets,X,{},T);
            ys = nets(xs,xis,ais);
            X = [(xi(1,:)) (xs(1,:))];
            T = [(xi(2,:)) ys];
        end

        stepAheadPerformance = perform(nets,ts,ys);
        PerforP = [PerforP; stepAheadPerformance];
        Tpredict = [Tpredict; T];

    end
    Var = mean(cell2mat(Tpredict));
    Var5= movmean(Var,3);
    [Rho,Ps] = corrcoef(Var(10:34),Ins(10:34));
    [Rhom,Psm] = corrcoef(Var5(10:34),Ins(10:34));

    %% Modelized data graphics
    grid on, hold on
    ax = gca;
    plot((1:207),Var(1,:),'red',(1:207),Var5(1,:),'blue',(1:207),Ins(1,:),'black')
    xticks(0:14:207)
    xticklabels(Years)
    ax = gca;
    ax.FontSize = 12;
    ax.XDir = 'reverse';
    xlabel('Years')
    ylabel('SST')
    title(['RNAn ',num2str(k),' regre: ',num2str(mean(regre)),'r: ',num2str(Rho(1,2)),...
        'rm5: ',num2str(Rhom(1,2))])
    saveas(gcf,sprintf('2Cpredict%d%d.jpg',k));
    close

    for i = k-2
        Modpredic(i).Layers = 1;
        Modpredic(i).Arch = hiddenLayerSize;
        Modpredic(i).Regr = mean(regre);
        Modpredic(i).perforModel = mean(Perfor);
        Modpredic(i).PerforPredic= mean(PerforP);
        Modpredic(i).Predic = Tpredict;
        Modpredic(i).Mean = Var;
        Modpredic(i).Meanmov5 = Var5;
        Modpredic(i).Corr = [Rho(1,2) Ps(1,2)];
        Modpredic(i).Corrm = [Rhom(1,2) Psm(1,2)];
    end
end
save('filename.mat', 'Modpredic')

%% %% Graph the model results (1951-2041) 
Years2 = [2041 2032 2023 2014 2005 1996 1987 1978 1969 1960 1951];
for i = 1:8
    amatrix = Modpredic(i).Predic;
    amatrix =cell2mat(amatrix(:,(1:29)))
    figure
    [aline(1), aFill(1)] = stdshade(amatrix,0.3,'r', 1:size(amatrix,2), 3);
    hold on
    ax = gca;
    set(ax,'fontname','Times New Roman','fontsize',12,'TickDir','out')
    ylim([21 30])
    plot(Ins,'blue')
    hold on
    plot(mean(amatrix),'black')
    ax.XDir = 'reverse';
    xticks(0:3:32)
    xticklabels(Years2)
    xlabel('Years')
    ax = gca;
    set(ax,'fontname','Times New Roman','fontsize',12,'TickDir','out')
    ylim([21 31])
    ax.FontSize = 12;
    ylabel('SST (°C)')
    saveas(gcf,sprintf('fig_Art%d%d.jpg',i));
    close
end


