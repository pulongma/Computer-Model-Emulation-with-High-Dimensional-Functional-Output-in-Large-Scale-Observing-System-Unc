function [RInv, dvec, Omega]=build_RInv(B_rowin,B_colin,theta,array_rho,array_C,funname)

n = size(array_C, 3); % number of inputs
m = size(array_C, 1);  % number of nearest neighbors

tau2 = theta.tau2;
sig2 = theta.sig2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
dvec=zeros([n 1]);
dvec(1)=sig2;    

BSall_vec = zeros(length(B_rowin), 1);

for indd=2:m    

    rhos=NNGP_Corr_AllD(theta, array_rho(1:(indd-1),1,indd,:), funname);
    Ctemp=NNGP_Corr_AllD(theta, array_C(1:(indd-1),1:(indd-1),indd,:), funname);
    % Ctemp = Ctemp + tau2*eye(size(Ctemp,1));

    L = chol(Ctemp, 'lower');
    Lrho = L\rhos;
    dvec(indd) = sig2 - sig2*(Lrho.'*Lrho);
    ind_start = indd*(indd-1)/2 - indd + 2;
    ind_end = indd*(indd-1)/2;
    BSall_vec(ind_start:ind_end) = Lrho.'/L;
end
     
for indd=(m+1):n

    rhos=NNGP_Corr_AllD(theta, array_rho(:,1,indd,:), funname);
    Ctemp=NNGP_Corr_AllD(theta, array_C(:,:,indd,:), funname);
    % Ctemp = Ctemp + tau2*eye(size(Ctemp,1));
    
    L = chol(Ctemp, 'lower');
    Lrho = L\rhos;
    dvec(indd) = sig2 - sig2*(Lrho.'*Lrho);
    ind_start = m*(m-1)/2 + m*(indd-m-1) + 1;
    ind_end = ind_start + m - 1;
    BSall_vec(ind_start:ind_end) = Lrho.'/L;
end

%% use C mex code
%[dvec, BSall_vec] = buildA(array_C, array_rho, theta); % memory leak, not sure why

% construct I - A in the NNGP paper
A = sparse(B_rowin, B_colin, BSall_vec, n, n);
B = speye(n) - A;
D = spdiags(dvec, 0, n, n);
RInv = B'*spdiags(dvec.^(-1), 0, n, n)*B;

%% add nugget if necessary
Omega = 1/tau2*speye(n,n) + RInv;


end
    