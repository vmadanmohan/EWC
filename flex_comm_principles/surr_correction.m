% Correction of observed EWC-neuralOsc correlation by cyclic surrogate correlation - ONLY RUN AFTER generating observed dataset (for all subs) and the surrogate dataset for each subject.
sub=1; % subject index
surr_comp=[];
comm_principle_surrcorr=[];
result=load(sprintf("observed/%d_results.mat",sub));
obs_comm=result.comm_principle_mean;
for surr=1:Nsurr
    surrresult=load(sprintf("surrogate/surrogate_results/%d/%d_results_surr.mat",sub,surr));
    surr_comp(:,:,surr)=abs(surrresult.comm_principle_mean_surr)>=abs(obs_comm); % check instances where surrogate relationship is stronger than the observed relationship
end
surr_comp=sum(surr_comp,3)/Nsurr;
for i=1:10 % 10 = total number of oscillatory measures we compute
    surr_comp(:,i)=mafdr(surr_comp(:,i),'BHFDR',true);
end
mask=surr_comp<0.05;
comm_principle_mean_surrcorr=obs_comm.*mask;
save(sprintf("observed/%d_results_surrcorr",sub),"comm_principle_mean_surrcorr")
