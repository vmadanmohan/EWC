% ----------Code to run a Linear Stochastic Model (with delays)------------
% ------------------------------(2 sources)--------------------------------
% ---------------------------(Noise Variation)-----------------------------
% ---------------------------------(EWC)-----------------------------------
%% Parameter definitions
addpath("/path/to/EWC/functions")
N = 4;              % Number of region
Fs = 2035;        % Sampling frequency (Hz)
ts = 1/Fs;          % Timestep (s)
T = 205;             % Total time (s)
dur = T*Fs;         % Number of samples (Total time * sampling frequency)
time = (0:dur-1)*ts;    % Time vector (for plots, interpretation etc.)
d12=0.015;
d32=0.015;
d42=d12;

C=[0 1 0 0;1 0 1 0;0 1 0 0;0 0 0 0];
K = 1;                 % Coupling strength
C=K*C;
C = C./(1+max(eig(C)));     % Stability scaling           

sigma = 0.05:0.025:1;           % Noise level
source = [1,3,4];                 % Index of source region
firing_rate = 0.2;          % Firing rate in Hz
pulse_amp = 0.1;
repeat = 20;
javaaddpath("/path/to/JIDT/infodynamics.jar");
implementingClass = 'infodynamics.measures.continuous.gaussian.ConditionalMutualInfoCalculatorMultiVariateGaussian';
cmiCalc = javaObject(implementingClass);

cMI = zeros(N,N,length(sigma));
cMI_std = zeros(N,N,length(sigma));
contrastEWC=zeros(3,repeat,length(sigma));
for p=1:length(sigma)
    p
    cMI_rep=zeros(N,N,repeat);
    cMI_repstd=zeros(N,N,repeat);
    for rep=1:repeat
        %% Dynamics
        delay=floor([0 d12 0 0;d12 0 d32 0;0 d32 0 0;0 d42 0 0]*Fs);    % Delay between regions, in seconds (ensure delay is not greater than stabilising duration), converted to timesteps
        E = zeros(dur,N);
        t=2;
        c=1;
        while c<=5*Fs              % Stabilise for 5 seconds
            for k = 1:N
                E(t,k) = E(t-1,k) + (ts)*(-E(t-1,k) + E(t-1,:)*C(:,k) + sigma(p)*randn); % No delay
            end
            c=c+1;
            t=t+1;
        end
        c=0;
        while t<=dur                % Stabilise for the rest of the duration
            for k=1:N
                if k==2
                    s = 0;
                    for i=1:3 % region 4 is isolated
                        s=s+E(t-delay(i,k),i)*C(i,k);
                    end
                    E(t,k) = E(t-1,k) + (ts)*(-E(t-1,k) + s + sigma(p)*randn);
                elseif k==4
                    pulse = (rand<=firing_rate*ts);
                    if pulse
                        E(t,k)=pulse_amp;
                    else
                        E(t,k) = E(t-1,k) + (ts)*(-E(t-1,k) + sigma(p)*randn);
                    end
                else
                    pulse = (rand<=firing_rate*ts);
                    if pulse
                        E(t,k)=pulse_amp;
                    else
                        s = 0;
                        for i=1:3
                            s=s+E(t-delay(i,k),i)*C(i,k);
                        end
                        E(t,k) = E(t-1,k) + (ts)*(-E(t-1,k) + s + sigma(p)*randn);
                    end
                end
            end
            t=t+1;
        end
        %% MI analysis
        delay=floor([0 d12 nan nan;d12 0 d32 d42;nan d32 0 nan; nan d42 nan 0]*Fs); % Delay between regions, in seconds (ensure delay is not greater than stabilising duration), converted to timesteps
        epoch = 10*Fs;              % 10 second epochs
        win = 1*Fs;                  % 1 second window
        main_data = E;
        main_data = main_data(5*Fs:end,:);  % Removing 5 seconds of data from the start
        numepochs=floor(size(main_data,1)/epoch);
        cMI_epoch = zeros(N,N,numepochs);
        cMI_std_epoch = zeros(N,N,numepochs);
        for s=1:numepochs
            data = main_data(((s-1)*epoch)+1:(s*epoch),:);
            [cMI_epoch(:,:,s),~]=CondMutInfo(data,win,N,delay,cmiCalc);  % Mean cMI over all significant events in target region
        end
        cMI_rep(:,:,rep)=mean(cMI_epoch,3);
        cMI_repstd(:,:,rep)=std(cMI_epoch,0,3);
        contrastEWC(1,rep,p)=(cMI_rep(1,2,rep)-cMI_rep(4,2,rep))/sqrt((cMI_repstd(1,2,rep)^2+cMI_repstd(4,2,rep)^2)/2);
        contrastEWC(2,rep,p)=(cMI_rep(3,2,rep)-cMI_rep(4,2,rep))/sqrt((cMI_repstd(3,2,rep)^2+cMI_repstd(4,2,rep)^2)/2);
        contrastEWC(3,rep,p)=(cMI_rep(1,2,rep)-cMI_rep(2,1,rep))/sqrt((cMI_repstd(1,2,rep)^2+cMI_repstd(2,1,rep)^2)/2);
    end
    cMI(:,:,p)=mean(cMI_rep,3);
    cMI_std(:,:,p)=std(cMI_rep,0,3);
end
figure;
hold on
shadedErrorBar(sigma/pulse_amp,smooth(squeeze(cMI(1,2,:))',100),smooth(squeeze(cMI_std(1,2,:))'/sqrt(repeat),100),'lineProps','r')
shadedErrorBar(sigma/pulse_amp,smooth(squeeze(cMI(2,1,:))',100),smooth(squeeze(cMI_std(2,1,:))'/sqrt(repeat),100),'lineProps','g')
shadedErrorBar(sigma/pulse_amp,smooth(squeeze(cMI(2,3,:))',100),smooth(squeeze(cMI_std(2,3,:))'/sqrt(repeat),100),'lineProps','k')
shadedErrorBar(sigma/pulse_amp,smooth(squeeze(cMI(3,2,:))',100),smooth(squeeze(cMI_std(3,2,:))'/sqrt(repeat),100),'lineProps','b')
shadedErrorBar(sigma/pulse_amp,smooth(squeeze(cMI(4,2,:))',100),smooth(squeeze(cMI_std(4,2,:))'/sqrt(repeat),100),'lineProps','c')
legend('1 \rightarrow 2','2 \rightarrow 1','2 \rightarrow 3','3 \rightarrow 2','4 \rightarrow 2')
xlabel('Noise amplitude / Pulse amplitude')
ylabel('Conditional Mutual Information (EWC)')
xlim([0.5,10])
ylim([-0.01,0.6])
fontsize(gca,15,"points")
hold off
contrastmeanEWC=mean(contrastEWC,2);
contraststdEWC=std(contrastEWC,0,2);
figure;
hold on
shadedErrorBar(sigma/pulse_amp,smooth(squeeze(contrastmeanEWC(1,1,:))',100),smooth(squeeze(contraststdEWC(1,1,:))'/sqrt(repeat),100),'lineProps','r')
shadedErrorBar(sigma/pulse_amp,smooth(squeeze(contrastmeanEWC(2,1,:))',100),smooth(squeeze(contraststdEWC(2,1,:))'/sqrt(repeat),100),'lineProps','b')
shadedErrorBar(sigma/pulse_amp,smooth(squeeze(contrastmeanEWC(3,1,:))',100),smooth(squeeze(contraststdEWC(3,1,:))'/sqrt(repeat),100),'lineProps','g')
yline(0)
legend('(1 \rightarrow 2) - (4 \rightarrow 2)','(3 \rightarrow 2) - (4 \rightarrow 2)','(1 \rightarrow 2) - (2 \rightarrow 1)')
xlim([0.5,10])
hold off
ylim([-0.1,2.15])
xlabel('Noise amplitude / Pulse amplitude')
ylabel('Standardised contrast')
fontsize(gca,15,"points")
