% Main script that calls EWC and neural osc functions - generates the communication principles for 1 subject (NOT SURROGATE CORRECTED) - repeat over subs to generate observed/original dataset

sub=1;			% Subject index
N = 150;              % Number of regions % 100 for schaefer (HCP), 148 for Destrieux (OMEGA)
Fs = 2035;        % Sampling frequency (Hz) % 2400 for OMEGA, 2035 for HCP
thresh=3;
recording=load(sprintf("/path/to/recordings/%d_resting.mat",sub));
main_data=recording.main_data;
delay=recording.delay; % Delay between regions, in seconds, converted to timesteps
SC=load("/path/to/SC/group_SC.txt")
SC=threshold_proportional(SC,0.15); % Keep top 15% connections
SC(SC>0)=1;                         % Binarize 
SC(1:N/2,N/2+1:end)=0;              % Remove inter hemispheric connections
SC(N/2+1:end,1:N/2)=0;              % ^^^
main_bandpass(:,:,1)=bandpass(main_data,[4,6],Fs);
main_bandpass(:,:,2)=bandpass(main_data,[8,12],Fs);
main_bandpass(:,:,3)=bandpass(main_data,[14,24],Fs);
main_bandpass(:,:,4)=bandpass(main_data,[30,59],Fs);
main_bandpass(:,:,5)=bandpass(main_data,[60,80],Fs);
epoch = 10*Fs;                      % 10 second epochs
win = 1*Fs;                         % 1 second window
numepochs=floor(size(main_data,1)/epoch);
P = zeros(N,N,numepochs);
P_std = zeros(N,N,numepochs);
Tpow = zeros(N,N,numepochs);
Apow = zeros(N,N,numepochs);
Bpow = zeros(N,N,numepochs);
Glopow = zeros(N,N,numepochs);
Ghipow = zeros(N,N,numepochs);
PLV_theta = zeros(N,N,numepochs);
PLV_alpha = zeros(N,N,numepochs);
PLV_beta = zeros(N,N,numepochs);
PLV_gammalo = zeros(N,N,numepochs);
PLV_gammahi = zeros(N,N,numepochs);
for s=1:numepochs
    data = main_data(((s-1)*epoch)+1:(s*epoch),:);
    bandpassSig=hilbert(main_bandpass(((s-1)*epoch)+1:(s*epoch),:,:));
    sigevents = eventiden(data,delay,win,thresh);
    [P(:,:,s),~]=PearsonEWCNetwork(data,SC,win,N,delay,sigevents);  % Mean PC over all significant events in target region
    [Tpow(:,:,s),Apow(:,:,s),Bpow(:,:,s),Glopow(:,:,s),Ghipow(:,:,s),PLV_theta(:,:,s),PLV_alpha(:,:,s),PLV_beta(:,:,s),PLV_gammalo(:,:,s),PLV_gammahi(:,:,s)] = neural_osc_network(data,bandpassSig,SC,win,N,delay,sigevents,Fs);
end
P = abs(P);
comm_principle_mean = zeros(N,10);
for i=1:N
    comm_principle_mean(i,1)=corr(reshape(squeeze(P(i,find(SC(i,:)),:)),[],1),reshape(squeeze(Tpow(i,find(SC(i,:)),:)),[],1));
    comm_principle_mean(i,2)=corr(reshape(squeeze(P(i,find(SC(i,:)),:)),[],1),reshape(squeeze(Apow(i,find(SC(i,:)),:)),[],1));
    comm_principle_mean(i,3)=corr(reshape(squeeze(P(i,find(SC(i,:)),:)),[],1),reshape(squeeze(Bpow(i,find(SC(i,:)),:)),[],1));
    comm_principle_mean(i,4)=corr(reshape(squeeze(P(i,find(SC(i,:)),:)),[],1),reshape(squeeze(Glopow(i,find(SC(i,:)),:)),[],1));
    comm_principle_mean(i,5)=corr(reshape(squeeze(P(i,find(SC(i,:)),:)),[],1),reshape(squeeze(Ghipow(i,find(SC(i,:)),:)),[],1));
    comm_principle_mean(i,6)=corr(reshape(squeeze(P(i,find(SC(i,:)),:)),[],1),reshape(squeeze(PLV_theta(i,find(SC(i,:)),:)),[],1));
    comm_principle_mean(i,7)=corr(reshape(squeeze(P(i,find(SC(i,:)),:)),[],1),reshape(squeeze(PLV_alpha(i,find(SC(i,:)),:)),[],1));
    comm_principle_mean(i,8)=corr(reshape(squeeze(P(i,find(SC(i,:)),:)),[],1),reshape(squeeze(PLV_beta(i,find(SC(i,:)),:)),[],1));
    comm_principle_mean(i,9)=corr(reshape(squeeze(P(i,find(SC(i,:)),:)),[],1),reshape(squeeze(PLV_gammalo(i,find(SC(i,:)),:)),[],1));
    comm_principle_mean(i,10)=corr(reshape(squeeze(P(i,find(SC(i,:)),:)),[],1),reshape(squeeze(PLV_gammahi(i,find(SC(i,:)),:)),[],1));
end
comm_principle_mean(isnan(comm_principle_mean))=0;
save(sprintf("%d_results",sub),"comm_principle_mean","Tpow","Apow","Bpow","Glopow","Ghipow","PLV_theta","PLV_alpha","PLV_beta","PLV_gammalo","PLV_gammahi","P")
