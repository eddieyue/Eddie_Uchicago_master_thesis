function [genedata, genenames, genepathways, is_isoprenoid] = ...
    get_arabidopsis_gene_data(),
genedata_url=...
    'http://genomebiology.com/content/supplementary/gb-2004-5-11-r92-s';
genenames=cell(0);
genepathways=cell(0);
genedata=[];
is_isoprenoid=[];
for datasetnum=1:2,
    dataset=urlread(strcat(genedata_url,int2str(datasetnum),'.txt'));
    dataset=regexp(dataset,'\n','split');
    start_read=false;
    for i=1:length(dataset),
        temp=regexp(dataset(i),' ','split');
        temp=temp{1}; 
        temp=regexprep(temp,char(13),''); % removes trailing newline
        if(~start_read),
            if(any(strcmp(temp,'Pathwayname'))), % identifies header line
                start_read=1;
                length_line=length(temp);
                pathway_column=find(strcmp(temp,'Pathwayname'));
                name_column=find(strcmp(temp,'AGI'));
                first_column=find(strcmp(temp,'c1'));
                last_column=find(strcmp(temp,'c118'));
            end
        elseif(length(temp)==length_line) % identifies a line of data
            % first check if this gene has already appeared
            % (some genes are listed multiple times, under different
            % pathways)
            if(any(strcmp(genenames,temp{name_column}))),
                % this gene was already recorded -- just add this new
                % pathway name to its pathway string
                gene_num=find(strcmp(genenames,temp{name_column}));
                genepathways{gene_num}=strcat(...
                    genepathways{gene_num},',',temp{pathway_column});
            else % a gene that has not appeared before
                genenames{length(genenames)+1}=temp{name_column};
                genepathways{length(genepathways)+1}=temp{pathway_column};
                genedata=[genedata;...
                    str2num(char(temp(first_column:last_column)))'];
                is_isoprenoid=[is_isoprenoid;(datasetnum==1)];
            end
        end
    end
end
% genedata = genedata';
end
