function subset = PathWay_subset(pathway,genepathways,genedata)

Ndim = size(genedata,2);

index = [];
for i = 1:Ndim
    if strcmp(genepathways(i),pathway);
        index = [index, i];
    end
end
subset = genedata(:,index);
end

   