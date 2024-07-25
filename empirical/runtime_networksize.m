%% Parameter definitions
addpath('EWC/functions/')
N = [3 10 20 50 100];
Fs = 2035;        % Sampling frequency (Hz)
ts = 1/Fs;          % Timestep (s)
d=0.015;            % Conduction delay (s)
repeat = 10;
javaaddpath("path/to/JIDT/infodynamics.jar");
implementingClass = 'infodynamics.measures.continuous.gaussian.TransferEntropyCalculatorGaussian';
teCalc = javaObject(implementingClass);
implementingClass = 'infodynamics.measures.continuous.gaussian.ConditionalMutualInfoCalculatorMultiVariateGaussian';
cmiCalc = javaObject(implementingClass);
implementingClass = 'infodynamics.measures.continuous.gaussian.MutualInfoCalculatorMultiVariateGaussian';
miCalc = javaObject(implementingClass);
runtimeTE=zeros(repeat,length(N));
runtimecMIEWC=zeros(repeat,length(N));
runtimeMI=zeros(repeat,length(N));
runtimeTEEWC=zeros(repeat,length(N));
sub=20
E=load(sprintf('%d_resting.csv',sub));
E=E(:,20*Fs+1:220*Fs); % Limit to 200 seconds of data (excluding first 20 seconds)
E=E';
for n=1:length(N)
    N(n)
    C=ones(N(n))-diag(diag(ones(N(n))));
    delay=floor(C*d*Fs);
    for rep=1:repeat
        rep
        main_data = E(:,randi(100,1,N(n)));
        epoch = 10*Fs;
        win = 1*Fs;                                 % 1 second window
        numepochs=floor(size(main_data,1)/epoch);
        epoch = 10*Fs;
        win = 1*Fs;                                 % 1 second window
        numepochs=floor(size(main_data,1)/epoch);
        TE_epoch = zeros(N(n),N(n),numepochs);
        cMI_epoch = zeros(N(n),N(n),numepochs);
        P_epoch = zeros(N(n),N(n),numepochs);
        %% TE full
        tic
        for s = 1:numepochs
            data = main_data(((s-1)*epoch)+1:(s*epoch),:);
            TE_epoch(:,:,s)=TransferEntFull(data,N(n),delay,teCalc);
        end
        runtimeTE(rep,n)=toc;
        %% cMI EWC
        tic
        for s = 1:numepochs
            data = main_data(((s-1)*epoch)+1:(s*epoch),:);
            [cMI_epoch(:,:,s),~]=CondMutInfo(data,win,N(n),delay,cmiCalc);
        end
        runtimecMIEWC(rep,n)=toc;
        %% TE EWC
        tic
        for s=1:numepochs
            data = main_data(((s-1)*epoch)+1:(s*epoch),:);
            [TE_epoch(:,:,s),~]=TransferEnt(data,win,N(n),delay,teCalc);
        end
        runtimeTEEWC(rep,n)=toc;
        %% MI full
        tic
        for s=1:numepochs
            data = main_data(((s-1)*epoch)+1:(s*epoch),:);
            cMI_epoch(:,:,s)=MIFull(data,N(n),delay,miCalc);
        end
        runtimeMI(rep,n)=toc;
        %% Pearson EWC
        tic
        for s=1:numepochs
            data = main_data(((s-1)*epoch)+1:(s*epoch),:);
            [P_epoch(:,:,s),~]=PearsonEWC(data,win,N(n),delay);
        end
        runtimeP(rep,n)=toc;
    end
end
save(sprintf('runtime%d.mat',sub),'-mat','runtimeTE','runtimeMI','runtimeTEEWC','runtimecMIEWC','runtimeP')
figure;
boxchart(runtimeTE);
hold on
boxchart(runtimecMIEWC)
boxchart(runtimeTEEWC)
boxchart(runtimeMI)
boxchart(runtimeP)
hold off
xticklabels({'3','10','20','50','100'})
xlabel("N")
ylabel("Network inference time for a 200s recording (s)")
legend("TE-Full","cMI-EWC","TE-EWC","MI-Full","PC-EWC")
fontsize(gcf,18,"points")
for i=1:5
    mediantimediff(i)=(median(runtimeTE1(:,i))/median(runtimeP1(:,i)));
    lessertime(:,i)=100*ones(30,1)-(runtimeP1(:,i)./runtimeTE1(:,i))*100;
end
figure;boxchart(-lessertime)
xticklabels({'3','10','20','50','100'})
xlabel("N")
ylabel("Relative time taken to compute PC-EWC (% of time to compute TE-Full)")
fontsize(gcf,18,"points")
