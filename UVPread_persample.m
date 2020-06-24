%% Construct a matlab base structured in size for a given depth strata (UVP)
% by Fabien Lombard 2016-2018
% adapted for multiple depth bin by Guillaume Bourdin January 2020

% clear all
% close all

%% offset_depth : 
%currently ecotaxa does not implement the offset dept of 1.2 m between the 
%depth sensor and the imaged zone
%this will be corrected in future versions of ecotaxa (and thus depth
%offset would have to be fixed to 0m

offset_depth=1.2;
%%
%for i = 1:width(t), if iscell(t.(i)), t.(i) = cell2mat(t.(i)); end, end

A=dir('*tsv');
filenames={A.name};
% f = msgbox('select the particle base')
fprintf('Select the particle base\n')
[file,path] = uigetfile('*.mat')

load(file)  % load tha particles data base in which volume per depth intervals are storred
samplebase={base(:).profile};

% prompt = {'Enter max depth'};%
% dlg_title = 'Input';
% num_lines = 1;
% defaultans = {'200'};
% maxdepth = inputdlg(prompt,dlg_title,num_lines,defaultans);
% maxdepth=str2num(maxdepth{1});
prompt = {'Enter depth bin size above MLD (m)', 'Input MLD (enter "no" to use default = 150 m)', 'Enter depth bin size below MLD (m)'};%
dlg_title = 'Input';
num_lines = 1;
defaultans = {'50', 'Yes', '300'};
sz_depthbin = inputdlg(prompt,dlg_title,num_lines,defaultans);
up_sz_depthbin = str2double(sz_depthbin{1});
down_sz_depthbin = str2double(sz_depthbin{3});
if contains(sz_depthbin{2}, {'yes', 'Yes', 'YES', 'true', 'True', 'TRUE'})
    [MLDfile,~] = uigetfile('*.csv');
    MLD = readtable(MLDfile);
else
    MLD = table(fullfile({base.profile}'), repmat(150, size(base,1), 1), 'VariableNames', {'profile', 'mld'});
    fprintf('Default MLD used: 150 m')
end

%samplebase=samplebase';
%% caution this assumes that the same UVP was used for an entire project (not anymore= corrected )
expo=[base(:).exp];	
aa=[base(:).aa];
%pixelsize=aa.^(expo); %
pixelsize=[base(:).PixelSize]; %in mm
% pixelsize=pixelsize.^0.5; % in mm2


% sampleheaders=readtable('allsampleheaders.xlsx');
% scanheaders=readtable('allscanheaders.xls');
base_Zooscan=[];
Idlist=[];
%return
h = waitbar(0,'Please wait...');

%return

%filenames([425 466])=[];

[~,m]=size(filenames);
%return

for i=1:m
%     if i==430
%         continue
%     end
%     if i==425
%         continue
%     end
%     if i==466
%         continue
%     end
    
    S=readtable(char(filenames(i)),'Filetype','text','ReadVariableNames',1);
    sample=unique(S.sample_id);
    fprintf('Sorting %s ... ', cell2mat(sample))
    %for i = 1:width(S), if iscell(t.(i)), t.(i) = cell2mat(t.(i)); end, end
    %% cleaning anotation hierarchy
%     [n,no_use]=size(S);
    S.object_annotation_hierarchy = strrep(S.object_annotation_hierarchy,'-','_');
    S.object_annotation_hierarchy = strrep(S.object_annotation_hierarchy,'>','_');

%     for j=1:n
%         f=find(S.object_annotation_hierarchy{j,:}=='-');
%         S.object_annotation_hierarchy{j,1}(f)='_';
%         
%         f=find(S.object_annotation_hierarchy{j,:}=='>');
%         S.object_annotation_hierarchy{j,1}(f)='_';
%     end
    
    %% cleaning for messy (text) entrance of files
    if iscell(S.object_depth_max)==1
        S.object_depth_max=cellfun(@str2num,S.object_depth_max);
    end

    if iscell(S.object_major)==1
        S.object_major=cellfun(@str2num,S.object_major);
    end

    if iscell(S.object_minor)==1
        S.object_minor=cellfun(@str2num,S.object_minor);
    end

    if iscell(S.object_area_exc)==1
        S.object_area_exc=cellfun(@str2num,S.object_area_exc);
    end

    if iscell(S.object_area)==1
        S.object_area=cellfun(@str2num,S.object_area);
    end

    if iscell(S.object_feret)==1
        S.object_feret=cellfun(@str2num,S.object_feret);
    end
    
    base_Zooscan(i).SampleID=sample;
    base_Zooscan(i).DN=unique(S.sample_dn);
    %base_Zooscan(i).Idstatus=
    base_Zooscan(i).Ship=unique(S.sample_ship);
    base_Zooscan(i).Scientificprog=unique(S.sample_cruise);
    base_Zooscan(i).StationId=sample;
    temp=cell2mat(sample);
    %base_Zooscan(i).StationIdnum=num2str(temp(6:8));  % to extract number in tara
    base_Zooscan(i).StationIdnum=(temp);
    base_Zooscan(i).Date=unique(S.object_date);
    base_Zooscan(i).time=unique(S.object_time);
    %base_Zooscan(i).Datenum=
    base_Zooscan(i).Latitude=unique(S.object_lat);
    base_Zooscan(i).Longitude=unique(S.object_lon);
    base_Zooscan(i).Depth=str2num(cell2mat(unique(S.sample_bottomdepth)));
    base_Zooscan(i).CTDref=unique(S.sample_ctdrosettefilename);
    base_Zooscan(i).profileid=unique(S.sample_profileid);
%     base_Zooscan(i).Townb=str2num(cell2mat(unique(S.sample_tow_nb)));
%     base_Zooscan(i).Towtype=str2num(cell2mat(unique(S.sample_tow_type)));
%     base_Zooscan(i).Nettype=unique(S.sample_net_type);
%     base_Zooscan(i).Netmesh=str2num(cell2mat(unique(S.sample_net_mesh)));
%     base_Zooscan(i).Netsurf=str2num(cell2mat(unique(S.sample_net_surf)));
%     base_Zooscan(i).Zmax=str2num(cell2mat(unique(S.sample_zmax)));
%     base_Zooscan(i).Zmin=str2num(cell2mat(unique(S.sample_zmin)));
%     base_Zooscan(i).Vol=str2num(cell2mat(unique(S.sample_tot_vol)));
%     base_Zooscan(i).Sample_comments=unique(S.sample_comment);
    
    
    %% getting pixel size in micrometer and converting in mm
    %unique(process_particle_pixel_size__m)
    %pixelsize=(str2num((cell2mat(unique(S.process_pixel))))).^(0.5);  % in mm/pixel
    %base_Zooscan(i).pixelsize=pixelsize;
    %pixelsize=0.174^(0.5);
    %%
    I = find(contains(samplebase,temp));
    volumes = base(I(1)).histnb.data.SampledVolume_L_;   %basetot(1,I).hisnb(:,3).*basetot(1,I).volimg0;   
    depthUVP = base(I(1)).histnb.data.Depth_m_; %basetot(1,I).hisnb(:,1);
%    if isempty(basetot(I).zoopuvp5)==0;
%       pixelsize = basetot(I).zoopuvp5.pixel;
%    end
    base_Zooscan(i).pixelsize = pixelsize(I(1));
    base_Zooscan(i).tot = table(NaN(ceil(max(depthUVP)/up_sz_depthbin),1),...
        NaN(ceil(max(depthUVP)/up_sz_depthbin),1), cell(ceil(max(depthUVP)/up_sz_depthbin),1),...
        cell(ceil(max(depthUVP)/up_sz_depthbin),1), cell(ceil(max(depthUVP)/up_sz_depthbin),1),...
        NaN(ceil(max(depthUVP)/up_sz_depthbin),1), NaN(ceil(max(depthUVP)/up_sz_depthbin),1),...
        cell(ceil(max(depthUVP)/up_sz_depthbin),1), cell(ceil(max(depthUVP)/up_sz_depthbin),1),...
        cell(ceil(max(depthUVP)/up_sz_depthbin),1), cell(ceil(max(depthUVP)/up_sz_depthbin),1),...
        cell(ceil(max(depthUVP)/up_sz_depthbin),1), cell(ceil(max(depthUVP)/up_sz_depthbin),1),...
        cell(ceil(max(depthUVP)/up_sz_depthbin),1), 'VariableNames',{'higher_bin_depth','lower_bin_depth','vol','conver',...
        'depthstrata','totvol','totconver','object_annotation_hierarchy',...
        'depth','major','minor','area_exc','area','perimferet'});
    % bin total depth
    bin_higher_depth = 0;
    if up_sz_depthbin <= MLD.mld(i)
        bin_lower_depth = up_sz_depthbin + 2.5;
    else
        bin_lower_depth = down_sz_depthbin + 2.5;
    end
    pp = 1;
    while bin_higher_depth < max(depthUVP)
        I = depthUVP > bin_higher_depth & depthUVP <= bin_lower_depth;
%         return
        base_Zooscan(i).tot.higher_bin_depth(pp) = bin_higher_depth;
        base_Zooscan(i).tot.lower_bin_depth(pp) = max(depthUVP(I));
        base_Zooscan(i).tot.vol{pp} = volumes(I)/1000; %in m3
        base_Zooscan(i).tot.conver{pp} = 1./(volumes(I)/1000);
        base_Zooscan(i).tot.depthstrata{pp} = depthUVP; %in m3

        base_Zooscan(i).tot.totvol(pp) = nansum(volumes(I)/1000); %in m3
        base_Zooscan(i).tot.totconver(pp) = 1/nansum(volumes(I)/1000);
%         return
%         base_Zooscan(i).tot.Scanfilename=
%         base_Zooscan(i).FracIds=unique(S.acq_id);
%         [nfrac, no_use]=size(base_Zooscan(i).FracIds);
        S.object_depth_max = S.object_depth_max+offset_depth;
        I = S.object_depth_max > bin_higher_depth & S.object_depth_max <= bin_lower_depth;
        %base_Zooscan(i).tot.Fracmin{pp} = unique(S.acq_min_mesh(I));
        %base_Zooscan(i).tot.Fracsup{pp} = unique(S.acq_max_mesh(I));
        %base_Zooscan(i).tot.Fracnb{pp} = unique(S.acq_sub_part(I));
        %base_Zooscan(i).tot.Scanned_Objects{pp} =
        %base_Zooscan(i).tot.Resolution{pp} = unique(S.process_img_resolution(I));

%         base_Zooscan(i).tot.object_annotation_hierarchy = unique(S.object_annotation_hierarchy(I));
        Idlist = unique([Idlist; S.object_annotation_hierarchy(I)]);
        base_Zooscan(i).tot.object_annotation_hierarchy{pp} = S.object_annotation_hierarchy(I);
        base_Zooscan(i).tot.depth{pp} = S.object_depth_max(I); %,S.object_depth_max
        %base_Zooscan(i).tot.object_annotation_hierarchy{pp} = S.object_annotation_hierarchy;
        base_Zooscan(i).tot.major{pp} = S.object_major(I)*base_Zooscan(i).pixelsize; %object_perimmajor
        base_Zooscan(i).tot.minor{pp} = S.object_minor(I)*base_Zooscan(i).pixelsize;
        base_Zooscan(i).tot.area_exc{pp} = S.object_area_exc(I)*(base_Zooscan(i).pixelsize^2);
        base_Zooscan(i).tot.area{pp} = S.object_area(I)*(base_Zooscan(i).pixelsize^2);    %object__area
        base_Zooscan(i).tot.perimferet{pp} = S.object_feret(I)*base_Zooscan(i).pixelsize;   % object_perimferet   %object_feretareaexc
        %base_Zooscan(i).tot.conver{pp} =  str2num(cell2mat(base_Zooscan(i).d.Fracnb))./base_Zooscan(i).Vol; %per cubic meter
        if bin_higher_depth <= MLD.mld(i)
            bin_lower_depth = bin_lower_depth + up_sz_depthbin;
        else
            bin_lower_depth = bin_lower_depth + down_sz_depthbin;
        end
        if bin_higher_depth <= MLD.mld(i)
            bin_higher_depth = bin_lower_depth - up_sz_depthbin;
        else
            bin_higher_depth = bin_lower_depth - down_sz_depthbin;
        end
        pp = pp + 1;
    end
    base_Zooscan(i).tot(isnan(base_Zooscan(i).tot.higher_bin_depth),:) = [];
    fprintf('done\n')
    waitbar(i/m)
end
close(h)

%save base_temporary base_Zooscan
 %load base_temporary
%return
%% now working on spectra
%% -------------- setting for calculus of the size spectra ------------------------------
% Area, Minor and Major are in mm in base_zooscan
smin = 0.000000000001;     %set the lower limit of the biovolume spectra that will be calculated
smax = 10000;      %set the upper limit of the biovolume spectra that will be calculated
k = 2^(1/4);      %set logarithmic base used to calculate bins of the size spectra
%according to platt & denman's theory (1978). scaling
%exponent set to 0.25.
uu = 1; % change here the first size class for the regression (and put option = 1 for performing on the spectra from uu to end)
%zoo_groups = sort(unique(zoo_groups));
%% starting to process the base
[~,m] = size(base_Zooscan);
zoo_groups = Idlist;
%return
nsamp = 0;
for i=1:m
    fprintf('Processing base %s ... ', cell2mat(base_Zooscan(i).SampleID))
    for j = 1:size(base_Zooscan(i).tot,1)
        nsamp = nsamp + 1;
        %[nfrac, no_use]=size(base_Zooscan(i).FracIds);
        fracnb = 1;
        conver1 = base_Zooscan(i).tot.conver{j};
        depth = base_Zooscan(i).tot.depth{j};
        depthstrata = base_Zooscan(i).tot.depthstrata{j};

        depthstrata = depthstrata(1:size(conver1,1));
%         conver2=NaN(size(depth,1),1);
        %% caution only works with 5m depth intervals
%         toto = 0.5+depthstrata/5;
        toto2 = ceil(depth/5);
        I = toto2 > base_Zooscan(i).tot.lower_bin_depth(j)/5;% for the few cases where depth observations sligthly overpass max depth (by max 2.5m...)
        toto2(I) = floor(base_Zooscan(i).tot.lower_bin_depth(j)/5);
        toto2 = toto2 - (min(toto2)-1);
%         toto2 = toto2 - sz_depthbin/5*(j-1);
        conver2 = conver1(toto2);
        base_Zooscan(i).tot.conver{j} = conver2;
        base_Zooscan(i).tot.converorig{j} = conver1;
        base_Zooscan(i).tot.conver{j} = conver2;
    end
    [base_Zooscan, SStot] = process_abundances_spectres_multiples_UVP(base_Zooscan(i).tot,...
        base_Zooscan, smin, smax, k, uu, zoo_groups, i, fracnb);
    fprintf('done\n')
end

%save base_temporary base_Zooscan
%load base_temporary

%return

%% now correcting from the fact that the data are in ind.m3 for each strata observed.... and cumulated
% thus abundances/biovolumes needs to be divided by the number of depth
% strata accumulated
for i=1:m
    for j = 1:size(base_Zooscan(i).tot,1)
        ndepth = size(base_Zooscan(i).tot.depthstrata{j},1);
        base_Zooscan(i).tot.Ab{j} = base_Zooscan(i).tot.Ab{j}/ndepth;
        base_Zooscan(i).tot.Yab{j} = base_Zooscan(i).tot.Yab{j}/ndepth;
        base_Zooscan(i).tot.Ybv_Plain_Area_BV_spectra{j} = base_Zooscan(i).tot.Ybv_Plain_Area_BV_spectra{j}/ndepth;
        base_Zooscan(i).tot.Ybv_Riddled_Area_BV_spectra{j} = base_Zooscan(i).tot.Ybv_Riddled_Area_BV_spectra{j}/ndepth;
        base_Zooscan(i).tot.Bv{j} = base_Zooscan(i).tot.Bv{j}/ndepth;
        base_Zooscan(i).tot.Ybv_Ellipsoid_BV_spectra{j} = base_Zooscan(i).tot.Ybv_Ellipsoid_BV_spectra{j}/ndepth;
    end
end

%% plus producing regrouped groups
table_groupage=readtable('zooregroup_zooscan.xlsx','ReadVariableNames',false);  %all copoda as omnivorous
table_groupage=table2cell(table_groupage);
    
    %% if needing to add new functional/trophic finction, mofify and add your desired within the excell file
% currently trophic groups includes
% -1= do not feed
% 1 phototrophs
% 1.5 mixotrophs
% 2 grazers
% 2.5 omnivorous
% 3 predators
% 3.5 unknown trophic group
% (note a 0.5 place is available for bacteria living from dissolved
% matter and potentially a 0 place is possible for viruses... but the
% "placement" is still subject to debate)

%% now working on regrouping taxa per functional/trophic groups
Zoo_groups=table_groupage(:,1);
% [n,p]=size(zoo_groups);

%return
%% finding if temporary groups are used
istemporary=0;
test = contains(zoo_groups,'temporary_');

if sum(test)>0
    answer = questdlg('Your files includes one or several temporary "t00X" categories. Do you have any "functional/trophic" mapping existing for those', ...
        'temporary categories mapping', ...
        'Yes please load them','No please create them','No please ignore them (not recommended)','Yes please load them');
    switch answer
    case 'No please ignore them (not recommended)'
        zoo_groups(test)=[];
%         [n,p]=size(zoo_groups);
        %% updating to remove temporary groups
        for i=1:m
            for j = 1:size(base_Zooscan(i).tot,1)
                base_Zooscan(i).tot.Zoo_groups{j}(test)=[];
                base_Zooscan(i).tot.Ab{j}(test)=[];              % abundance per fraction rapportée au volume (#/m3)
                base_Zooscan(i).tot.Bv{j}(test)=[];               % abundance per fraction rapportée au volume (#/m3)

                base_Zooscan(i).tot.Yab{j}(:,test)=[];
                base_Zooscan(i).tot.Ybv_Plain_Area_BV_spectra{j}(:,test)=[];
                base_Zooscan(i).tot.Ybv_Riddled_Area_BV_spectra{j}(:,test)=[];
                base_Zooscan(i).tot.Ybv_Ellipsoid_BV_spectra{j}(:,test)=[];
            end
        end
    case 'Yes please load them'
        [file,path] = uigetfile('*.xlsx')
        addontemp = readtable([path file],'ReadVariableNames',false);  %all copoda as herbivorous
        addontemp = table2cell(addontemp);
        Zoo_groups = [Zoo_groups; addontemp(:,1)];
        istemporary = 1;
    case 'No please create them'

        groups=table_groupage(:,2:end);
%         [n,p]=size(groups);

        group1=unique(groups(:,1));
        group2=unique(groups(:,2));
        group3=unique(cellstr(num2str(cell2mat(groups(:,3)))));

        %return
        new_taxa = zoo_groups(test);
        p = length(new_taxa);
        newfunctional = [];

        for i=1:p
            settings = settingsdlg('Description', ['A new temporary taxonomic group have been found ' char(new_taxa(i))],...
                'title' , 'New taxa functional mapping',...
                'Alive' , group1 ,...
                'functional group' , group2 , ...
                'trophic group' , group3 ,...
                'WindowWidth' , 800)
            newfunctional = [newfunctional settings];
        end
        newfunctional=struct2table(newfunctional);
        newfunctional(:,4)=[];
        %% updating the xls reference list
        addontemp=[new_taxa table2cell(newfunctional)];
        tosave=array2table(addontemp);
        [file,path] = uiputfile('temporarymapping_instrument_net_location.xlsx');
        filename = fullfile(path,file);
        writetable(tosave,filename,'WriteVariableNames',0);
        Zoo_groups=[Zoo_groups; addontemp(:,1)];           
        istemporary=1;
    end
end
%% checking if no "new" groups are present
[n,~]=size(zoo_groups);
new_taxa={};p=0;

for i=1:n
    %J=strcmp(char(zoo_groups(i)),Zoo_groups(1:I-1,:));
    J=strcmp(char(zoo_groups(i)),Zoo_groups);
    if sum(J)==0
        p=p+1;
        new_taxa(p,1)=zoo_groups(i);
    end
end

clear Zoo_group
%% proposing a mapping for the new groups
groups = table_groupage(:,2:end);
[n,~] = size(groups);

group1=unique(groups(:,1));
group2=unique(groups(:,2));
group3=unique(cellstr(num2str(cell2mat(groups(:,3)))));

%return
p = length(new_taxa);
newfunctional=[];

for i=1:p
    settings = settingsdlg('Description', ['A new taxonomic group have been found ' char(new_taxa(i))],...
        'title' , 'New taxa functional mapping',...
        'Alive' , group1 ,...
        'functional group' , group2 , ...
        'trophic group' , group3 ,...
        'WindowWidth' , 800)
    newfunctional=[newfunctional settings];
end
if p>0
newfunctional=struct2table(newfunctional);
newfunctional(:,4)=[];
%% updating the xls reference list
%
addon=[new_taxa table2cell(newfunctional)];

table_groupage=[table_groupage; addon];

tosave=array2table(table_groupage);
cd(directoryoftoolbox);
writetable(tosave,'zooregroup_zooscan.xlsx','WriteVariableNames',0);
cd(folder);
end
if istemporary==1
    table_groupage=[table_groupage; addontemp];
end

%% producing the regrouped groups
for i=1:m
    base_regroup = f_regroup_all(table_groupage,base_Zooscan(i).tot);
    base_Zooscan(i).regroupped = base_regroup;
end
    
%% preparing resume files on abundance / BV per taxa per sample
Zoo_groups = base_Zooscan(i).regroupped.Zoo_groups{1};
Ab_resume = NaN(nsamp,size(Zoo_groups,1));
Bv_resume = NaN(nsamp,size(Zoo_groups,1));
samplelist = cell(nsamp,1);
depthslice = cell(nsamp,1);
nsa = 0;
for i=1:m
    for j = 1:size(base_Zooscan(i).tot,1)
        nsa = nsa + 1;
        Ab_resume(nsa,:) = base_Zooscan(i).regroupped.Ab{j};
        Bv_resume(nsa,:) = base_Zooscan(i).regroupped.Bv{j};
        samplelist(nsa) = base_Zooscan(i).SampleID;
        depthslice{nsa} = num2str(j);
    end
end

%% saving the final bases and resume files

instrument=char(list2(indx2));

prompt = {'instrument:','project:'};
title = 'save base under the name:';
dims = [1 35];
definput = {instrument,'pointB_Regent_1995_2019'};
answer = inputdlg(prompt,title,dims,definput)

save(['base_instrument_' char(answer(1)) '_' char(answer(2))],'base_Zooscan','-v7.3')
%save base_spectre_zooscan_regent_point_B base_Zooscan
%save base_spectre_flowcam_168b20 base_spectres

colname = strcat(samplelist, cellfun(@(c)['_depthslice_' c], depthslice, 'uni',false));

%Abtable=table(Ab_resume,'VariableNames',Zoo_groups,'RowNames',samplelist);    % do not work because the name of taxa are TOO LONG
Abtable=array2table(Ab_resume','VariableNames', colname, 'rownames', Zoo_groups);
Bvtable=array2table(Bv_resume','VariableNames', colname, 'rownames',Zoo_groups);
writetable(Abtable,'Abundance_resume.csv','WriteRowNames',true)
writetable(Bvtable,'Biovolume_resume.csv','WriteRowNames',true)
    
  
