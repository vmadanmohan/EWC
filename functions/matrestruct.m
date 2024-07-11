function mat=matrestruct(C,mapping)
    rowrestruct = zeros(size(C,1),size(C,2));
    colrestruct = zeros(size(C,1),size(C,2));
    for i =1:length(mapping)
        rowrestruct(mapping(i),:)=C(i,:);
    end
    for i=1:length(mapping)
        colrestruct(:,mapping(i))=rowrestruct(:,i);
    end
    mat=colrestruct;
