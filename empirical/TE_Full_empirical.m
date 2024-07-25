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
    mapping=load('mapping.txt');
    delay=load(sprintf("/path/to/%d_delayfile.txt",sub));
    epoch = 10*Fs;              % 10 second epochs
    win = 1*Fs;                  % 1 second window
    numepochs(sub)=floor(size(main_data,1)/epoch);
    TE = zeros(N,N,numepochs(sub));
    for s=1:numepochs(sub)
        data = main_data(((s-1)*epoch)+1:(s*epoch),:);
        TE(:,:,s)=TransferEntFull(data,N,delay,teCalc);
    end
    TE_mean(:,:,sub)=matrestruct(mean(TE,3),mapping);
    TE_mean_subjectlevelstd(:,:,sub)=matrestruct(std(TE,0,3),mapping);
end
