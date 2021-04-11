function [array_C,array_rho, B_rowin, B_colin,neigh_index]=NeighD_AllD(loc_sorted,m)
indi=1;
Dim_x=size(loc_sorted,2);
nsub=size(loc_sorted, 1);
neigh_index=zeros([m nsub]);

%obs_pred_dist_all=cell(Dim_x,1);
%dist_neigh_mat_all=cell(Dim_x,1);
% array_C = cell(Dim_x,1);
% array_rho = cell(Dim_x,1);
% for i=1:Dim_x
%     obs_pred_dist_all{i} = zeros([m, nsub]);
%     array_C{i} = zeros([m, m, nsub]);
%     array_rho{i} = zeros([m, m, nsub]);
% end

obs_pred_dist_all = zeros([m, Dim_x]);
array_C = zeros([m, m, nsub, Dim_x]);
array_rho = zeros([m, 1, nsub, Dim_x]);


nnzment = m*(m-1)/2 + (nsub-m)*m;
B_rowin = zeros(1, nnzment);
B_colin = zeros(1, nnzment);

for indi=2:m
    % the first m locations, neighbors are previous locations
    %indi

    obs_loc=loc_sorted(1:indi-1,:);
   
    % for i=1:Dim_x
    %     array_C{i}(1:(indi-1),1:(indi-1),indi)=pdist2(obs_loc(:,i), obs_loc(:,i));   
    %     array_rho{i}(1:(indi-1),indi)=pdist2(obs_loc(:,i), loc_sorted(indi,i));
    % end

    for i=1:Dim_x
        array_C(1:(indi-1),1:(indi-1),indi,i) = pdist2(obs_loc(:,i), obs_loc(:,i));
        array_rho(1:(indi-1),1,indi,i) = pdist2(obs_loc(:,i), loc_sorted(indi,i));
    end

    neigh_index(1:(indi-1),indi)=1:(indi-1);

    ind_start = (indi-2)*(indi-1)/2 + 1;
    ind_end = (indi-1)*indi/2;
    B_rowin(ind_start:ind_end) = indi*ones(1, indi-1);
    B_colin(ind_start:ind_end) = (neigh_index(1:(indi-1),indi))';
end

for indi=(m+1): (nsub)

    obs_loc=loc_sorted(1:indi-1,:);
    
    obs_pred_dist1=pdist2(obs_loc, loc_sorted(indi,:));
    
    [sorted_dis1,sorted_ind]=sort(obs_pred_dist1);
    neigh_index(1:m,indi)=sorted_ind(1:m);    

    % for i=1:Dim_x   
    %     obs_pred_dist{i}=pdist2(obs_loc(:,i), loc_sorted(indi,i));%((obs_loc(:,i)-loc_sorted(indi,i)).^2).^(0.5); %CHANGE THIS yloc_sorted
    %     sorted_dist{i}=obs_pred_dist{i}(sorted_ind);
    %     obs_pred_dist_all{i}(1:m,indi)=sorted_dist{i}(1:m); % THIS is a key diference
    % end
    
    for i=1:Dim_x
        dtemp = pdist2(obs_loc(:,i), loc_sorted(indi,i));
        obs_pred_dist_all(1:m,i) = dtemp(sorted_ind(1:m));
    end

    % compute neighbors distance matrix
    obs_loc=loc_sorted(neigh_index(1:m,indi),:);%[xloc_sorted(neigh_index(1:m,indi)) yloc_sorted(neigh_index(1:m,indi))];
    
    % for i=1:Dim_x
    %     %dist_neigh_mat_all{i}(1:m, 1:m,indi)=pdist2(obs_loc(:,i), obs_loc(:,i));
    %     array_C{i}(1:m,1:m,indi)=pdist2(obs_loc(:,i), obs_loc(:,i));
    %     array_rho{i}(1:m,indi)=(obs_pred_dist_all{i}(1:(m),indi))';
    % end
    
    for i=1:Dim_x
        array_C(:,:,indi,i) = pdist2(obs_loc(:,i), obs_loc(:,i));
        array_rho(:,1,indi,i) = obs_pred_dist_all(1:m,i);
    end

    ind_start = (m-1)*m/2 + m*(indi-m-1) + 1;
    ind_end = ind_start + m -1;
    B_rowin(ind_start:ind_end) = indi*ones(1, m);
    B_colin(ind_start:ind_end) = neigh_index(1:m,indi)';
end



end


