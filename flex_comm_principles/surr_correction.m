% Correction of observed EWC-neuralOsc correlation by cyclic surrogate correlation - ONLY RUN AFTER generating observed dataset (for all subs) and the surrogate dataset for each subject.
sub=1; % subject index
surr_comp=[];
comm_principle_surrcorr=[];
result=load(sprintf("/path/to/observed_results/%d_results.mat",sub));
obs_comm=result.comm_principle_mean;
for surr=1:1000
    surrresult=load(sprintf("/path/to/surrogate/results/%d/%d_results_surr.mat",sub,surr));
    surr_comp(:,:,surr)=abs(surrresult.comm_principle_mean_surr)>=abs(obs_comm); % check instances where surrogate relationship is stronger than the observed relationship
end
surr_comp=sum(surr_comp,3)/1000;
for i=1:10 % 10 = total number of oscillatory measures we compute
    surr_comp(:,i)=mafdr(surr_comp(:,i),'BHFDR',true);
end
mask=surr_comp<0.05;
comm_principle_mean_surrcorr=obs_comm.*mask;
save(sprintf("/path/to/output/dir/%d_results_surrcorr",sub),"comm_principle_mean_surrcorr")
