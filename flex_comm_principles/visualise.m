function visualise
    addpath(genpath('funcs/'));    
    surf = load('funcs/parc_plotter/data/fsaverage/mat/fsaverage_inflated.mat') ;
    surfStruct = surf.surfStruct ;
    
    annots = load('funcs/parc_plotter/data/fsaverage/mat/fsaverage_annots.mat') ;
    annotMap = annots.allAnnots ;
    annotName = 'schaefer7_100' ; % see funcs/setup_data.m for all available 'annotName'
    % To see the region names
    % annotMap(annotName).combo_names;
    armean=[];
    for sub=1:30
        armean(:,:,sub)=load(sprintf("observed/%d_results_surrcorr",sub),'comm_principle_mean_surrcorr');
    end
    armean=tanh(mean(atanh(armean),3));
    for i = 1:5
        dataVec = armean(:,i);
        parc_plot(surfStruct,annotMap,annotName,dataVec,'viewcMap',1,'newFig',1,'cMap',flip(brewermap(1000,'RdYlBu'),1),'border',0,'valrange',[min(min(armean(:,1:5))) max(max(armean(:,1:5)))]);
    end
    for i = 6:10
        dataVec = armean(:,i);
        parc_plot(surfStruct,annotMap,annotName,dataVec,'viewcMap',1,'newFig',1,'cMap',flip(brewermap(1000,'RdYlBu'),1),'border',0,'valrange',[min(min(armean(:,6:10))) max(max(armean(:,6:10)))]);
    end
end
