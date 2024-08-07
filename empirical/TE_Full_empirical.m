addpath("/path/to/EWC/functions")
N = 100;              % Number of regions
Fs = 2035;        % Sampling frequency (Hz)
subjects = 30;
javaaddpath("/path/to/JIDT/infodynamics.jar");
implementingClass = 'infodynamics.measures.continuous.gaussian.TransferEntropyCalculatorGaussian';
teCalc = javaObject(implementingClass);
TE_mean=zeros(N,N,subjects);
TE_mean_subjectlevelstd=zeros(N,N,subjects);
numepochs=zeros(subjects,1);
parfor sub=1:subjects
    recording=load(sprintf("/path/to/%d_resting.mat",sub));
    main_data=recording.main_data;
    delay=recording.delay;
    epoch = 10*Fs;              % 10 second epochs
    win = 1*Fs;                  % 1 second window
    numepochs(sub)=floor(size(main_data,1)/epoch);
    TE = zeros(N,N,numepochs(sub));
    for s=1:numepochs(sub)
        data = main_data(((s-1)*epoch)+1:(s*epoch),:);
        TE(:,:,s)=TransferEntFull(data,N,delay,teCalc);
    end
    TE_mean(:,:,sub)=mean(TE,3);
    TE_mean_subjectlevelstd(:,:,sub)=std(TE,0,3);
end
