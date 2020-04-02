% fonction de groupage du zooplankton pour calculer abondance et spectres
% des PFT

%variables d'entrée:
% 1 table_groupage : un tableau xls avec 3 colonnes. 
%       1iere colonne = groupes originaux
%       2ieme colonne = nouveaux noms des groupes living vs non_living, 
%       3ieme colonne = nouveaux noms des groupes PFT
% 2 abundances : tableau d'abondances (ind/m3)
%       # lignes = # de groupes de zoo identifiés à la base 
%       # colonnes = # d'échantillons
% 3 biovolumes : tableau de biovolumes (mm3/m3)
%       # lignes = # de groupes de zoo identifiés à la base 
%       # colonnes = # d'échantillons
% 4 size_spectra : SStot en sortie de la routine de creation de la
% base_spectres

function [base_regroup] = f_regroup_all(table_groupage,base)

tab=table_groupage;
%tab_orig=id;

base_regroup = table(base.higher_bin_depth, base.lower_bin_depth, cell(size(base,1),1),...
    cell(size(base,1),1), cell(size(base,1),1), cell(size(base,1),1),...
    cell(size(base,1),1), cell(size(base,1),1), cell(size(base,1),1),...
    NaN(size(base,1),1), NaN(size(base,1),1), NaN(size(base,1),1),...
    cell(size(base,1),1), cell(size(base,1),1),...
    'VariableNames',{'higher_bin_depth','lower_bin_depth','Ab','Ybv_Plain_Area_BV_spectra',...
    'Ybv_Riddled_Area_BV_spectra','Bv','Ybv_Ellipsoid_BV_spectra','Zoo_groups',...
    'originalplace','allplace','livingplace','notlivingplace','pftplace','trophicplace'});

for dep = 1:size(base,1)
    tab_orig=base.Zoo_groups{dep};

    g_pft2=table_groupage(:,3);

    to_keep = zeros(size(tab,1),1);
    for i=1:size(tab_orig,1)
        a=strcmp(tab_orig(i),tab(:,1));
        if sum(a)>1
           warning(['warning' tab_orig(i) 'is in double, please check reference excel file before going further']);
        end
        to_keep = to_keep + double(a);
    end

    if max(to_keep)>1
       warning('warning there is a double in zoo_groups, please check reference excel file before going further');
    end

    tab=tab(logical(to_keep),:);

    g_orig=tab(:,1);
    g_living=tab(:,2);% idée: recherche des strings "living" et "not-living"
    g_pft=tab(:,3);
    g_trophic=tab(:,4);

    %données de spectres de taille
    % Ab
    % Yab
    % Ybv_Plain_Area_BV_spectra
    % Ybv_Riddled_Area_BV_spectra
    % Bv
    % Ybv_Ellipsoid_BV_spectra

    Ab=base.Ab{dep};
    %Yab=base.Yab;
    Ybv_Plain_Area_BV_spectra=base.Ybv_Plain_Area_BV_spectra{dep};
    Ybv_Riddled_Area_BV_spectra=base.Ybv_Riddled_Area_BV_spectra{dep};
    Bv=base.Bv{dep};
    ss=base.Ybv_Ellipsoid_BV_spectra{dep};

    clear SA
    %% regroupement total = 'all'
    Ab_t=sum(Ab);
    %Yab_t=sum(Yab,2);
    Ybv_Plain_Area_BV_spectra_t=sum(Ybv_Plain_Area_BV_spectra,2);
    Ybv_Riddled_Area_BV_spectra_t=sum(Ybv_Riddled_Area_BV_spectra,2);
    Bv_t=sum(Bv);

    sst=sum(ss,2);

    %% regroupement living vs non living
    new_groups = unique(g_living);
    Ab_n1 = NaN(1, size(new_groups,1));
    Ybv_Plain_Area_BV_spectra_n1 = NaN(size(Ybv_Plain_Area_BV_spectra,1), size(new_groups,1));
    Ybv_Riddled_Area_BV_spectra_n1 = NaN(size(Ybv_Riddled_Area_BV_spectra,1), size(new_groups,1));
    Bv_n1 = NaN(1,size(new_groups,1));
    ssn1 = NaN(size(ss,1),size(new_groups,1));
    for i=1:size(new_groups,1)
        ng = new_groups(i);
        f_ng = strcmp(g_living, ng);

        % spectres bv
        ss_r = ss(:,f_ng);
        Ab_r = Ab(f_ng);
        %Yab_r=Yab(:,f_ng);
        Ybv_Plain_Area_BV_spectra_r = Ybv_Plain_Area_BV_spectra(:,f_ng);
        Ybv_Riddled_Area_BV_spectra_r = Ybv_Riddled_Area_BV_spectra(:,f_ng);
        Bv_r=Bv(f_ng);

        Ab_n1(i) = sum(Ab_r);
        %Yab_n1(:,i)=sum(Yab_r,2);
        Ybv_Plain_Area_BV_spectra_n1(:,i) = sum(Ybv_Plain_Area_BV_spectra_r,2);
        Ybv_Riddled_Area_BV_spectra_n1(:,i) = sum(Ybv_Riddled_Area_BV_spectra_r,2);
        Bv_n1(i) = sum(Bv_r);
        ssn1(:,i) = sum( ss_r,2);

        clear SSJ ssj ss_r f_ng ng Bv_r Ybv_Riddled_Area_BV_spectra_r Ybv_Plain_Area_BV_spectra_r  Ab_r
    end

    %% regroupement pft
    new_groups=unique(g_pft2);
    Ab_n2 = NaN(1,size(new_groups,1));
    Ybv_Plain_Area_BV_spectra_n2 = NaN(size(Ybv_Plain_Area_BV_spectra,1), size(new_groups,1));
    Ybv_Riddled_Area_BV_spectra_n2 = NaN(size(Ybv_Riddled_Area_BV_spectra,1), size(new_groups,1));
    Bv_n2 = NaN(1,size(new_groups,1));
    ssn2 = NaN(size(ss,1),size(new_groups,1));
    for i=1:size(new_groups,1)
        ng = new_groups(i);
        f_ng = strcmp(g_pft, ng);

        %spectres bv
        ss_r = ss(:,f_ng);
        Ab_r =  Ab(f_ng);
        %Yab_r=Yab(:,f_ng);
        Ybv_Plain_Area_BV_spectra_r = Ybv_Plain_Area_BV_spectra(:,f_ng);
        Ybv_Riddled_Area_BV_spectra_r = Ybv_Riddled_Area_BV_spectra(:,f_ng);
        Bv_r = Bv(f_ng);

        Ab_n2(i) = sum(Ab_r);
        %Yab_n2(:,i)=sum(Yab_r,2);
        Ybv_Plain_Area_BV_spectra_n2(:,i) = sum(Ybv_Plain_Area_BV_spectra_r,2);
        Ybv_Riddled_Area_BV_spectra_n2(:,i) = sum(Ybv_Riddled_Area_BV_spectra_r,2);
        Bv_n2(i) = sum(Bv_r);
        ssn2(:,i) = sum( ss_r,2) ;

        clear SSJ ssj ss_r  f_ng ng Bv_r Ybv_Riddled_Area_BV_spectra_r Ybv_Plain_Area_BV_spectra_r Ab_r
    end

    %% regroupement trophic
    new_groups=unique(cell2mat(g_trophic));
    Ab_n3 = NaN(1,size(new_groups,1));
    Ybv_Plain_Area_BV_spectra_n3 = NaN(size(Ybv_Plain_Area_BV_spectra,1), size(new_groups,1));
    Ybv_Riddled_Area_BV_spectra_n3 = NaN(size(Ybv_Riddled_Area_BV_spectra,1), size(new_groups,1));
    Bv_n3 = NaN(1,size(new_groups,1));
    ssn3 = NaN(size(ss,1),size(new_groups,1));
    for i=1:size(new_groups,1)
        ng = new_groups(i);
        f_ng = ng==cell2mat(g_trophic);
        %spectres bv
        ss_r = ss(:,f_ng);
        Ab_r = Ab(f_ng);
        %Yab_r=Yab(:,f_ng);
        Ybv_Plain_Area_BV_spectra_r = Ybv_Plain_Area_BV_spectra(:,f_ng);
        Ybv_Riddled_Area_BV_spectra_r = Ybv_Riddled_Area_BV_spectra(:,f_ng);
        Bv_r=Bv(f_ng);

        Ab_n3(i)=sum(Ab_r);
        %Yab_n3(:,i)=sum(Yab_r,2);
        Ybv_Plain_Area_BV_spectra_n3(:,i)=sum(Ybv_Plain_Area_BV_spectra_r,2);
        Ybv_Riddled_Area_BV_spectra_n3(:,i)=sum(Ybv_Riddled_Area_BV_spectra_r,2);
        Bv_n3(i)=sum(Bv_r);
        ssn3(:,i)=sum( ss_r,2);

        clear SSJ ssj ss_r  f_ng ng Bv_r Ybv_Riddled_Area_BV_spectra_r Ybv_Plain_Area_BV_spectra_r Ab_r
    end
    %%
    %Spec_regroup=[ss sst ssn1 ssn2 ssn3];

    Ab_regroup=[Ab' Ab_t Ab_n1 Ab_n2 Ab_n3];
    %Yab_regroup=[Yab Yab_t Yab_n1 Yab_n2 Yab_n3];
    Ybv_Plain_Area_BV_spectra_regroup=[Ybv_Plain_Area_BV_spectra Ybv_Plain_Area_BV_spectra_t Ybv_Plain_Area_BV_spectra_n1 Ybv_Plain_Area_BV_spectra_n2 Ybv_Plain_Area_BV_spectra_n3];
    Ybv_Riddled_Area_BV_spectra_regroup=[Ybv_Riddled_Area_BV_spectra Ybv_Riddled_Area_BV_spectra_t Ybv_Riddled_Area_BV_spectra_n1 Ybv_Riddled_Area_BV_spectra_n2 Ybv_Riddled_Area_BV_spectra_n3];
    Bv_regroup=[Bv' Bv_t Bv_n1 Bv_n2 Bv_n3];
    Ybv_Ellipsoid_BV_spectra_regroup=[ss sst ssn1 ssn2 ssn3];

    %base_regroup.Abtot=base.Abtot; % ne sert a rien
    base_regroup.Ab{dep} = Ab_regroup;
    %base_regroup.Yab=Yab_regroup;
    base_regroup.Ybv_Plain_Area_BV_spectra{dep} = Ybv_Plain_Area_BV_spectra_regroup;
    base_regroup.Ybv_Riddled_Area_BV_spectra{dep} = Ybv_Riddled_Area_BV_spectra_regroup;
    base_regroup.Bv{dep} = Bv_regroup;
    %base_regroup.Bvtot=base.Bvtot; % ne sert a rien
    base_regroup.Ybv_Ellipsoid_BV_spectra{dep} = Ybv_Ellipsoid_BV_spectra_regroup;

    base_regroup.Zoo_groups{dep} = [g_orig;{'all'};unique(g_living);unique(g_pft2);cellstr(num2str(unique(cell2mat(g_trophic))))];
    base_regroup.originalplace{dep} = [1 length(g_orig)];
    base_regroup.allplace(dep) = length(g_orig)+1;
    base_regroup.livingplace(dep) = length(g_orig)+2;
    base_regroup.notlivingplace(dep) = length(g_orig)+3;
    base_regroup.pftplace{dep} = [length(g_orig)+4 length(g_orig)+3+length(unique(g_pft2))];
    base_regroup.trophicplace{dep} = [length(g_orig)+4+length(unique(g_pft2)) length(g_orig)+3+length(unique(g_pft2))+length(unique(cell2mat(g_trophic)))];
    % groups=[g_orig(:,1);{'all'};unique(g_living);unique(g_pft)];
end
