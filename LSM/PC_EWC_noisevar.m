% ----------Code to run a Linear Stochastic Model (with delays)------------
% ------------------------------(3 sources)--------------------------------
% ---------------------------(Noise Variation)-----------------------------
% -------------------------------(PC-EWC)----------------------------------
%% Parameter definitions
addpath("path/to/EWC/functions")
N = 4;                % Number of region
Fs = 2035;            % Sampling frequency (Hz)
ts = 1/Fs;            % Timestep (s)
T = 205;              % Total time (s)
dur = T*Fs;           % Number of samples (Total time * sampling frequency)
time = (0:dur-1)*ts;    % Time vector (for plots, interpretation etc.)
d12=0.015;          % Delay between regions 1 and 2
d32=0.015;          % Delay between regions 3 and 2
d42=d12;            % Delay between regions 4 and 2

C=[0 1 0 0;1 0 1 0;0 1 0 0;0 0 0 0];    % Adjacency matrix
K = 1;                 % Coupling strength
C=K*C;
C = C./(1+max(eig(C)));      % Stability scaling           

sigma = 0.05:0.025:1;        % Noise level
source = [1,3,4];            % Index of source region
firing_rate = 0.2;            % Firing rate in Hz
pulse_amp = 0.1;              % Source pulse amplitude
repeat = 20;                  % Number of trials

P = zeros(N,N,length(sigma));
P_std = zeros(N,N,length(sigma));
contrastEWC=zeros(3,repeat,length(sigma));
for p=1:length(sigma)
    p
    P_rep=zeros(N,N,repeat);
    P_repstd=zeros(N,N,repeat);
    for rep=1:repeat
        %% Dynamics
        delay=floor([0 d12 0 0;d12 0 d32 0;0 d32 0 0;0 d42 0 0]*Fs);    % Delay between regions, in seconds (ensure delay is not greater than stabilising duration), converted to timesteps
        E = zeros(dur,N);
        t=2;
        c=1;
        while c<=5*Fs              % Stabilise for 5 seconds (No source firing)
            for k = 1:N
                E(t,k) = E(t-1,k) + (ts)*(-E(t-1,k) + E(t-1,:)*C(:,k) + sigma(p)*randn); % No delay
            end
            c=c+1;
            t=t+1;
        end
        c=0;
        while t<=dur                % Poisson sources active for the rest of the duration
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
        delay=floor([0 d12 nan nan;d12 0 d32 d42;nan d32 0 nan; nan d42 nan 0]*Fs); % Delay between regions, in seconds, converted to timesteps (setting delay as NaN excludes it from the EWC protocol)
        epoch = 10*Fs;                % 10 second epochs
        win = 1*Fs;                    % 1 second window
        main_data = E;
        main_data = main_data(5*Fs:end,:);  % Removing 5 seconds of data from the start (stabilising duration)
        numepochs=floor(size(main_data,1)/epoch);
        P_epoch = zeros(N,N,numepochs);
        P_std_epoch = zeros(N,N,numepochs);
        for s=1:numepochs
            data = main_data(((s-1)*epoch)+1:(s*epoch),:);
            [P_epoch(:,:,s),~]=PearsonEWC(data,win,N,delay);
        end
        P_rep(:,:,rep)=mean(P_epoch,3);
        P_repstd(:,:,rep)=std(P_epoch,0,3);
        contrastEWC(1,rep,p)=(abs(P_rep(1,2,rep))-abs(P_rep(4,2,rep)))/sqrt((P_repstd(1,2,rep)^2+P_repstd(4,2,rep)^2)/2);  % Difference between a true connection and an absent connection, normalised by their pooled STD
        contrastEWC(2,rep,p)=(abs(P_rep(3,2,rep))-abs(P_rep(4,2,rep)))/sqrt((P_repstd(3,2,rep)^2+P_repstd(4,2,rep)^2)/2);  % Difference between a true connection and an absent connection, normalised by their pooled STD
        contrastEWC(3,rep,p)=(abs(P_rep(1,2,rep))-abs(P_rep(2,1,rep)))/sqrt((P_repstd(1,2,rep)^2+P_repstd(2,1,rep)^2)/2);  % Difference between a true connection and an absent connection, normalised by their pooled STD
    end
    P(:,:,p)=mean(P_rep,3);
    P_std(:,:,p)=std(P_rep,0,3);
end
smoothfactor=20;
figure;
hold on
shadedErrorBar(sigma/pulse_amp,smooth(abs(squeeze(P(1,2,:))'),smoothfactor),smooth(squeeze(P_std(1,2,:))'/sqrt(repeat),smoothfactor),'lineProps','r')
shadedErrorBar(sigma/pulse_amp,smooth(abs(squeeze(P(2,1,:))'),smoothfactor),smooth(squeeze(P_std(2,1,:))'/sqrt(repeat),smoothfactor),'lineProps','g')
shadedErrorBar(sigma/pulse_amp,smooth(abs(squeeze(P(2,3,:))'),smoothfactor),smooth(squeeze(P_std(2,3,:))'/sqrt(repeat),smoothfactor),'lineProps','k')
shadedErrorBar(sigma/pulse_amp,smooth(abs(squeeze(P(3,2,:))'),smoothfactor),smooth(squeeze(P_std(3,2,:))'/sqrt(repeat),smoothfactor),'lineProps','b')
shadedErrorBar(sigma/pulse_amp,smooth(abs(squeeze(P(4,2,:))'),smoothfactor),smooth(squeeze(P_std(4,2,:))'/sqrt(repeat),smoothfactor),'lineProps','c')
legend('1 \rightarrow 2','2 \rightarrow 1','2 \rightarrow 3','3 \rightarrow 2','4 \rightarrow 2')
xlabel('Noise amplitude / Pulse amplitude')
ylabel('Partial correlation (EWC)')
xlim([0.5,10])
ylim([-0.01,0.6])
fontsize(gca,15,"points")
hold off

contrastmeanEWC=mean(contrastEWC,2);
contraststdEWC=std(contrastEWC,0,2);
figure;
hold on
shadedErrorBar(sigma/pulse_amp,smooth(squeeze(contrastmeanEWC(1,1,:))',smoothfactor),smooth(squeeze(contraststdEWC(1,1,:))'/sqrt(repeat),smoothfactor),'lineProps','r')
shadedErrorBar(sigma/pulse_amp,smooth(squeeze(contrastmeanEWC(2,1,:))',smoothfactor),smooth(squeeze(contraststdEWC(2,1,:))'/sqrt(repeat),smoothfactor),'lineProps','b')
shadedErrorBar(sigma/pulse_amp,smooth(squeeze(contrastmeanEWC(3,1,:))',smoothfactor),smooth(squeeze(contraststdEWC(3,1,:))'/sqrt(repeat),smoothfactor),'lineProps','g')
yline(0)
legend('(1 \rightarrow 2) - (4 \rightarrow 2)','(3 \rightarrow 2) - (4 \rightarrow 2)','(1 \rightarrow 2) - (2 \rightarrow 1)')
xlim([0.5,10])
ylim([-0.1,2.15])
hold off
xlabel('Noise amplitude / Pulse amplitude')
ylabel('Standardised contrast')
fontsize(gca,15,"points")
