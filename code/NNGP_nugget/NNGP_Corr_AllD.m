function [Cmat]=NNGP_Corr_AllD(theta, d, funname)


Dim=size(d); 
Dim_x = Dim(end);


switch funname
	case 'exp'  % exponential correlation function
		phi = theta.phi;
		Q = 0;
		for k=1:Dim_x
			Q = Q + d(:,:,k).^2/phi(k)^2;
		end

		Cmat = exp(-sqrt(Q)); 

	case 'sqexp'  % squared exponential or Gaussian correlation function 
		phi = theta.phi;
		Q = 0;
		for k=1:Dim_x
			Q = Q + d(:,:,k).^2/phi(k)^2;
		end

		Cmat = exp(-Q); 
  
	case 'matern_5_2'
		phi = theta.phi;

		Cmat = ones(size(d(:,:,1)));
		for k=1:Dim_x
			dtemp = d(:,:,k)/phi(k);
			Cmat = Cmat .* (1+sqrt(5)*dtemp+(5/3)*dtemp.^2).*exp(-sqrt(5)*dtemp);
		end

	case 'aniso_matern_5_2'
		phi = theta.phi;

		Cmat = ones(size(d(:,:,1)));
		dtemp = 0;
		for k=1:Dim_x
			dtemp = dtemp + d(:,:,k).^2/phi(k).^2;
		end
		dtemp = dtemp.^0.5;
		Cmat = (1+sqrt(5)*dtemp+(5/3)*dtemp.^2).*exp(-sqrt(5)*dtemp);

	case 'aniso_CH'
		phi = theta.phi;
		alpha = theta.alpha;
		nu = theta.nu;
		dtemp = 0;
		for k=1:Dim_x
			dtemp = dtemp + d(:,:,k).^2/phi(k).^2;
		end
		con = gammaln(nu+alpha) - gammaln(nu);
		Cmat = exp(con)*kummerU(alpha, 1-nu, dtemp);


end



end