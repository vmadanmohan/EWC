function [Tpow,Apow,Bpow,Glopow,Ghipow,PLV_theta,PLV_alpha,PLV_beta,PLV_gammalo,PLV_gammahi] = neural_osc_network_cyclicsurr(orig_data,orig_bandpassSig,SC,win,N,delay,sigevents,sampling_rate,r)
    Tpow=zeros(N,N);
    Apow=zeros(N,N);
    Bpow=zeros(N,N);
    Glopow=zeros(N,N);
    Ghipow=zeros(N,N);
    PLV_theta=zeros(N,N);
    PLV_alpha=zeros(N,N);
    PLV_beta=zeros(N,N);
    PLV_gammalo=zeros(N,N);
    PLV_gammahi=zeros(N,N);
    orig_bandpassSig=angle(orig_bandpassSig);
    for i=1:N          % Loop over source regions
	data=datashift_randcyc(orig_data,i,r(i));
	bandpassSig=zeros(size(orig_bandpassSig,1),size(orig_bandpassSig,2),size(orig_bandpassSig,3));
	for band=1:5
		bandpassSig(:,:,band)=datashift_randcyc(orig_bandpassSig(:,:,band),i,r(i));
	end
        reldelays = delay(i,:);  % Delays relative to the source region
        signum = sum(sigevents(:,i),1);  % Number of significant events
        if signum~=0
            sigind = find(sigevents(:,i));  % Indices of significant events
            totevents=signum;
            Tpow_event = zeros(totevents,N);
            Apow_event = zeros(totevents,N);
            Bpow_event = zeros(totevents,N);
            Glopow_event = zeros(totevents,N);
	    Ghipow_event = zeros(totevents,N);
	    PLV_event = zeros(totevents,N,5);
            while signum            % While loop starting from the last significant event and moving backwards
                for j=find(SC(i,:))   % Loop over all possible targets
                    [cxy,~] = mscohere(data(sigind(signum):sigind(signum)+win-1,i),data(sigind(signum)+reldelays(j):sigind(signum)+reldelays(j)+win-1,j),floor(sampling_rate/2),[],[],sampling_rate); % Coherence spectrum between source and target within the communication window
                    [pxx,~] = pwelch(data(sigind(signum)+reldelays(j):sigind(signum)+reldelays(j)+win-1,j),floor(sampling_rate/2),[],[],sampling_rate); % Target's power spectral density within the communication window
                    Tpow_event(signum,j)=trapz(pxx(3:4))/trapz(pxx); % RELATIVE Theta power (Theta power of target within the window/total power of the target within the window)
                    Apow_event(signum,j)=trapz(pxx(5:7))/trapz(pxx); % RELATIVE Alpha power (Alpha power of target within the window/total power of the target within the window)
                    Bpow_event(signum,j)=trapz(pxx(8:13))/trapz(pxx); % RELATIVE Beta power (Beta power of target within the window/total power of the target within the window)
                    Glopow_event(signum,j)=trapz(pxx(16:31))/trapz(pxx); % RELATIVE Gamma power (Gamma power of target within the window/total power of the target within the window)
                    Ghipow_event(signum,j)=trapz(pxx(32:41))/trapz(pxx);
		    for band=1:5
                    	PLV_event(signum,j,band)=abs(sum(exp(1i*(bandpassSig(sigind(signum):sigind(signum)+win-1,i,band)-bandpassSig(sigind(signum)+reldelays(j):sigind(signum)+reldelays(j)+win-1,j,band))))/win);
                    end
                end
                signum=signum-1;
            end
        end
        Tpow(i,:)=mean(Tpow_event,1);
        Apow(i,:)=mean(Apow_event,1);
        Bpow(i,:)=mean(Bpow_event,1);
        Glopow(i,:)=mean(Glopow_event,1);
	Ghipow(i,:)=mean(Ghipow_event,1);
	PLV_theta(i,:)=mean(PLV_event(:,:,1),1);
        PLV_alpha(i,:)=mean(PLV_event(:,:,2),1);
        PLV_beta(i,:)=mean(PLV_event(:,:,3),1);
        PLV_gammalo(i,:)=mean(PLV_event(:,:,4),1);
        PLV_gammahi(i,:)=mean(PLV_event(:,:,5),1);
    end
