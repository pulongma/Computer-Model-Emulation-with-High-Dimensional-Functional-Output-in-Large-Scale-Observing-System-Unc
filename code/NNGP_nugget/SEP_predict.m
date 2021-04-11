function Zhat = SEP_predict(output, H, param, B_rowin, B_colin, array_rho, array_C, ...
				funname, H_new, neigh_index_p, array_rho_p, array_C_p)



n = size(output,1);
q = size(output,2);
p = size(H,2);
Dim_x = size(array_C,4);
np = size(H_new,1);

m = size(array_rho, 1); % number of nearest neighbors

df = n-p; % degrees of freedom

phi = param.phi;
tau2 = param.tau2;

nsample = length(tau2); % number of MCMC samples


theta.sig2 = 1;
sig2 = theta.sig2;

Zhat = zeros(np, q, nsample);
w = zeros(np, q);
rho_hat = zeros(np,1);

for iter=1:nsample

	theta.phi = phi(iter, :);
	theta.tau2 = tau2(iter);

    [RInv, dvec, Omega]=build_RInv(B_rowin,B_colin,theta,array_rho,array_C,funname);
	
	RInvH = RInv*H;
	[LOmega,~,po] = chol(Omega, 'lower', 'vector');
	OInvH(po,:) = LOmega'\(LOmega\H(po,:));
	HRH = OInvH'*RInvH * (1/theta.tau2);
	L = chol(HRH, 'lower');

	RInvZ = RInv*output*(1/theta.tau2);
	OmegaInvZ(po, :) = LOmega'\(LOmega\output(po,:));
	ZRZ = OmegaInvZ'*RInvZ;
	HRZ = RInvH'*OmegaInvZ * (1/theta.tau2);
	LHRZ = L\HRZ;

	Beta_gls = L'\LHRZ;
	ZGZ = ZRZ - LHRZ'*LHRZ;

	Sigma_gls = ZGZ/df;
	Sigma_chol = chol(Sigma_gls, 'lower');

	output_detrend = output - H*Beta_gls;

	%% perform prediction 
	for k=1:np
		rhos = NNGP_Corr_AllD(theta, array_rho_p(:,1,k,:), funname);
		Ctemp = NNGP_Corr_AllD(theta, array_C_p(:,:,k,:), funname);
		Ctemp = Ctemp + theta.tau2*eye(size(Ctemp,1));
		Crho = Ctemp\rhos;
		rho_hat(k) = sig2 - sig2*rhos'*Crho + theta.tau2;
		%w(k,:) = Crho'*output_detrend(neigh_index_p(:,k), :);
		zhat_mean = H_new(k,:)*Beta_gls + Crho'*output_detrend(neigh_index_p(:,k), :);
		Zhat(k, :, iter) = rmvt(zhat_mean(:), sqrt(rho_hat(k))*Sigma_chol, df);
	end

end

end


