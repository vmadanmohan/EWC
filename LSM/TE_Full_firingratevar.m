% ----------Code to run a Linear Stochastic Model (with delays)------------
% ------------------------------(2 sources)--------------------------------
% -------------------------(Firing Rate Variation)-------------------------
% -------------------------------(TE-Full)---------------------------------
%% Parameter definitions
% clear java
addpath("/path/to/EWC/functions")
N = 3;              % Number of region
Fs = 2035;        % Sampling frequency (Hz)
ts = 1/Fs;          % Timestep (s)
T = 205;             % Total time (s)
dur = T*Fs;         % Number of samples (Total time * sampling frequency)
time = (0:dur-1)*ts;    % Time vector (for plots, interpretation etc.)
d12=0.015;
d32=0.015;

C=[0 1 0;1 0 1;0 1 0];
K = 1;                 % Coupling strength
C = C*K;
C = C./(1+max(eig(C)));     % Stability scaling

sigma = 0.2;           % Noise level
source = [1,3];                 % Index of source region
firing_rate = 0.1:0.05:1;          % Firing rate in Hz
pulse_amp = 0.1;
repeat = 10;
javaaddpath("/Users/vmadanmohan/Documents/MATLAB/infodynamics-dist-1.6/infodynamics.jar");
implementingClass = 'infodynamics.measures.continuous.gaussian.TransferEntropyCalculatorGaussian';
teCalc = javaObject(implementingClass);

TE = zeros(N,N,length(firing_rate));
TE_std = zeros(N,N,length(firing_rate));
for p=1:length(firing_rate)
    p
    TE_rep=zeros(N,N,repeat);
    TE_repstd=zeros(N,N,repeat);
    for rep=1:repeat
        %% Dynamics
        delay=floor([0 d12 0;d12 0 d32;0 d32 0]*Fs);    % Delay between regions, in seconds (ensure delay is not greater than stabilising duration), converted to timesteps
        E = zeros(dur,N);
        t=2;
        c=1;
        while c<=5*Fs              % Stabilise for 5 seconds
            for k = 1:N
                E(t,k) = E(t-1,k) + (ts)*(-E(t-1,k) + E(t-1,:)*C(:,k) + sigma*randn); % No delay
            end
            c=c+1;
            t=t+1;
        end
        c=1;
        while t<=dur                % Stabilise for the rest of the duration
            for k=1:N
                if k==2
                    s = 0;
                    for i=1:N
                        s=s+E(t-delay(i,k),i)*C(i,k);
                    end
                    E(t,k) = E(t-1,k) + (ts)*(-E(t-1,k) + s + sigma*randn);
                else
                    pulse = (rand<=firing_rate(p)*ts);
                    if pulse
                        E(t,k)=pulse_amp;
                    else
                        s = 0;
                        for i=1:N
                            s=s+E(t-delay(i,k),i)*C(i,k);
                        end
                        E(t,k) = E(t-1,k) + (ts)*(-E(t-1,k) + s + sigma*randn);
                    end
                end
            end
            t=t+1;
        end
        %% MI analysis
        delay=floor([0 d12 nan;d12 0 d32;nan d32 0]*Fs);
        epoch = 10*Fs;              % 10 second epochs
        win = 1*Fs;                 % 1 second window
        main_data = E;
        main_data = main_data(5*Fs:end,:);  % Removing 5 seconds of data from the start
        numepochs=floor(size(main_data,1)/epoch);
        TE_epoch = zeros(N,N,numepochs);
        TE_std_epoch = zeros(N,N,numepochs);
        for s=1:numepochs
            data = main_data(((s-1)*epoch)+1:(s*epoch),:);
            TE_epoch(:,:,s)=TransferEntFull(data,N,delay,teCalc);
        end
        TE_rep(:,:,rep)=mean(TE_epoch,3);
        TE_repstd(:,:,rep)=std(TE_epoch,0,3);
    end
    TE(:,:,p)=mean(TE_rep,3);
    TE_std(:,:,p)  =std(TE_rep,0,3);
end
smoothfactor=10;
figure;
hold on
shadedErrorBar(firing_rate,smooth(squeeze(TE(1,2,:))',smoothfactor),smooth(squeeze(TE_std(1,2,:))'/sqrt(repeat),smoothfactor),'lineProps','r')
shadedErrorBar(firing_rate,smooth(squeeze(TE(2,1,:))',smoothfactor),smooth(squeeze(TE_std(2,1,:))'/sqrt(repeat),smoothfactor),'lineProps','g')
shadedErrorBar(firing_rate,smooth(squeeze(TE(2,3,:))',smoothfactor),smooth(squeeze(TE_std(2,3,:))'/sqrt(repeat),smoothfactor),'lineProps','k')
shadedErrorBar(firing_rate,smooth(squeeze(TE(3,2,:))',smoothfactor),smooth(squeeze(TE_std(3,2,:))'/sqrt(repeat),smoothfactor),'lineProps','b')
legend('1 \rightarrow 2','2 \rightarrow 1','2 \rightarrow 3','3 \rightarrow 2')
xlabel('Source firing rate (Hz)')
ylabel('Transfer Entropy (Full)')
xlim([0.1,1])
ylim([-0.0001,0.0045])
fontsize(gca,18,"points")
hold off
