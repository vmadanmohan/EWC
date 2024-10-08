addpath("/path/to/EWC/functions")
N = 100;              % Number of regions
Fs = 2035;        % Sampling frequency (Hz)
subjects = 30;
P_mean=zeros(N,N,subjects);
P_mean_subjectlevelstd=zeros(N,N,subjects);
numepochs=zeros(subjects,1);
parfor sub=1:subjects
    sub
    recording=load(sprintf("/path/to/%d_resting.mat",sub));
    main_data=recording.main_data;
    delay=recording.delay;
    epoch = 10*Fs;              % 10 second epochs
    win = 1*Fs;                  % 1 second window
    numepochs(sub)=floor(size(main_data,1)/epoch);
    P = zeros(N,N,numepochs(sub));
    P_std = zeros(N,N,numepochs(sub));
    for s=1:numepochs(sub)
        data = main_data(((s-1)*epoch)+1:(s*epoch),:);
        [P(:,:,s),P_std(:,:,s)]=PearsonEWC(data,win,N,delay);
    end
    P_mean(:,:,sub)=mean(P,3);
    P_mean_subjectlevelstd(:,:,sub)=std(P,0,3);
end
