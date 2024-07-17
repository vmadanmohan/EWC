# save orthogonalised source localised recordings as .mat file
import numpy as np
import scipy as sci
from osl.source_recon import parcellation
from scipy.io import savemat

def delete_bad_segments(bad,data,numepochs):
    epoch_markers = [s*epoch for s in range(numepochs+1)]
    if len(bad.shape)==2:
        a=bad.shape[0]
        for i in range(a-1,-1,-1):
            j=0
            while epoch_markers[j]<int(np.floor(bad[i][0]*sampling_rate)):
                j=j+1
            start=epoch_markers[j-1]
            if bad[i][1]*sampling_rate<epoch_markers[-1]:
                while epoch_markers[j]<int(np.ceil(bad[i][1]*sampling_rate)):
                    j=j+1
                stop=epoch_markers[j]
                data=np.delete(data,np.arange(start,stop),0)
            else:
                data=np.delete(data,np.arange(start,data.shape[0]),0)
    else:
        j=0
        while epoch_markers[j]<int(np.floor(bad[0]*sampling_rate)):
            j=j+1
        start=epoch_markers[j-1]
        if bad[1]*sampling_rate<epoch_markers[-1]:
            while epoch_markers[j]<int(np.ceil(bad[1]*sampling_rate)):
                j=j+1
            stop=epoch_markers[j]
            data=np.delete(data,np.arange(start,stop),0)
        else:
            data=np.delete(data,np.arange(start,data.shape[0]),0)
    return data

N = 100
epoch_dur = 10    # in seconds
sampling_rate = 2035
T = 1         # window time in seconds
window = int(T*sampling_rate)   # window time in timesteps

sub_IDs = open("subIDs.txt", "r") 
data = sub_IDs.read() 
subs = data.split("\n") 
sub_IDs.close()

for sub in range(len(subs)):
    epoch = int(epoch_dur*sampling_rate)
    hemi = ['LLRLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLLRRRRRRRRRLRRRRRRRRRRLRRRRRRRRRRLRRRRRRRRRRLRRRRRRRRRR']
    tic = perf_counter();
    delay = np.loadtxt(f'{sub}_delayfile.txt',delimiter=',',dtype=int)
    main_data = np.loadtxt(f'{sub}_restin3.csv',delimiter=',')
    leakage_corr_ts = parcellation.symmetric_orthogonalise(main_data,maintain_magnitudes=True)
    main_data = leakage_corr_ts.T
    toc=perf_counter();
    print(f"file load time = {toc-tic} seconds")
    numepochs=int(main_data.shape[0]/epoch)
    bad = np.loadtxt(f'{sys.argv[1]}_badseg.txt',delimiter=',')
	if bad.size==0:
		bad = np.empty((1,2))
    main_data = delete_bad_segments(bad,main_data,numepochs)
    main_data = main_data[10*sampling_rate:-1-(10*sampling_rate)+1,:]  # Removing 10 seconds of data from the start and end to remove any transients from processing
    savemat(f'{sub}_restin.mat',{'main_data':main_data})
