function [array_C_p,array_rho_p, B_rowin_p, B_colin_p,neigh_index_p]=NeighD_AllDpr(obsloc_sorted,m,obsloce)
%indi=1;
Dim_x=size(obsloc_sorted,2);
Ne = size(obsloce,1);

% obs_pred_dist_all=cell(1,Dim_x);
% dist_neigh_mat_all=cell(1,Dim_x);
% array_C_p =cell(Dim_x,1);
% array_rho_p=cell(Dim_x,1);
% for i=1:Dim_x
%     obs_pred_dist_all{i}=zeros([m Ne]);
%     dist_neigh_mat_all{i}=zeros([m m Ne]);
%     array_C_p{i}=zeros([m m (Ne)]);
%     array_rho_p{i}=zeros([m (Ne)]);
% end


array_C_p = zeros([m, m, Ne, Dim_x]);
array_rho_p = zeros([m, 1, Ne, Dim_x]);
obs_pred_dist_all = zeros([m, Dim_x]);

neigh_index_p=zeros([m Ne]);

% B_rowin_p=[];
% B_colin_p=[];
B_rowin_p = zeros(1, m*Ne);
B_colin_p = zeros(1, m*Ne);

for indi=(1):(Ne)
   % indi

    % compute pred-obs distances
    new_loc=obsloce(indi,:);
    %obs_loc1=p_locs(obsloc_sorted, new_loc);
    obs_loc1=obsloc_sorted;
    obs_pred_dist=pdist2(obs_loc1,new_loc);%D_dif2r(obs_loc1,new_loc);%
    
    [sorted_dis,sorted_ind]=sort(obs_pred_dist);
    % for i=1:Dim_x
    %     obs_pred_dist_X{i}=pdist2(obs_loc1(:,i),new_loc(:,i));%D_dif2r(obs_loc1(:,i),new_loc(:,i));    
    %     sorted_dis_X{i}=obs_pred_dist_X{i}(sorted_ind);
    %     obs_pred_dist_all{i}(1:m,indi)=sorted_dis_X{i}(1:m); % THIS is a key diference
    % end

    for i=1:Dim_x
        dtemp = pdist2(obs_loc1(:,i),new_loc(:,i));
        obs_pred_dist_all(1:m,i) = dtemp(sorted_ind(1:m));
    end

    neigh_index_p(1:m,indi)=sorted_ind(1:m);     
    % compute neighbors distance matrix
    obs_loc2=obs_loc1(neigh_index_p(1:m,indi),:);%[xloc1(neigh_index_p(1:m,indi)) yloc1(neigh_index_p(1:m,indi))];
    % for i=1:Dim_x
    %     dist_neigh_mat_all{i}(1:m, 1:m,indi)=pdist2(obs_loc2(:,i),obs_loc2(:,i));%
    %     array_C_p{i}(:,:,indi)=dist_neigh_mat_all{i}(1:(m),1:(m),indi);
    %     array_rho_p{i}(:,indi)=(obs_pred_dist_all{i}(1:(m),indi))';  
    % end

    for i=1:Dim_x
        array_C_p(:,:,indi,i) = pdist2(obs_loc2(:,i), obs_loc2(:,i));
        array_rho_p(:,1,indi,i) = obs_pred_dist_all(1:m,i);
    end

    %dist_neigh_mat_all
    
    % B_rowin_p=[B_rowin_p indi*ones([1 m])];
    % B_colin_p=[B_colin_p (neigh_index_p(1:m,indi))'];
    ind_start = (indi-1)*m+1;
    ind_end = indi*m;
    B_rowin_p(ind_start:ind_end) = indi*ones([1 m]);
    B_colin_p(ind_start:ind_end) = (neigh_index_p(1:m,indi))';

end

end