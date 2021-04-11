function [theta_curr, loglik_curr] = SEP_MwH_Gibbs_O21(output, H, theta_curr, B_rowin, B_colin, prior, Delta, array_rho, array_C, funname, nugget_est)

[loglik_curr]=SEP_loglik(output, H, theta_curr, B_rowin, B_colin, array_rho, array_C, funname);


Dim=size(array_C);
Dim_x = Dim(end);

% extract parameters in prior distributions
% parameters in the prior density of tau2
tau2_a = prior.tau2_a;
tau2_b = prior.tau2_b;


% bounds on the parameters
tau2_lb = prior.tau2_bd(1);
tau2_ub = prior.tau2_bd(2);

phi_lb = prior.phi_lb;
phi_ub = prior.phi_ub;


%%%%%%%% update tau2
theta_prop = theta_curr;
% if(nugget_est)

	tau2_curr = theta_curr.tau2;
	tau2_prop=lognrnd(log(tau2_curr), Delta.tau2,1);
	theta_prop.tau2 = tau2_prop;

	% loglikelihood evaluation
	[loglik_prop]=SEP_loglik(output, H, theta_prop, B_rowin, B_colin, array_rho, array_C, funname);

	% log of prior density of cauchy dist
	log_prior_curr = -log(1+tau2_curr^2);
	log_prior_prop = -log(1+tau2_prop^2);


	% from tau2_curr --> tau2_prop.
	log_tau2_tran = -log(tau2_prop)-(log(tau2_prop)-log(tau2_curr))^2/(2*Delta.tau2^2);
	% from tau2_prop --> tau2_curr.
	log_tau2_rtran = -log(tau2_curr)-(log(tau2_curr)-log(tau2_prop))^2/(2*Delta.tau2^2);

	log_diff=(loglik_prop-loglik_curr)+ (log_prior_prop-log_prior_curr) + (log_tau2_rtran-log_tau2_tran);

	a_unif=rand(1); 
	if a_unif<=exp(log_diff)
	    tau2_curr=tau2_prop;
	    loglik_curr=loglik_prop;
	end

	theta_curr.tau2 = tau2_curr;
% end 

%%%%%%%% update phi
phi_prop = theta_curr.phi;

for k=1:(Dim_x)
    theta_prop = theta_curr;

    phi_curr = theta_curr.phi;
    % random walk for phi_x 
    phi_prop(k)=lognrnd(log(phi_curr(k)), Delta.phi(k), 1);
    if phi_prop(k)<=phi_lb(k) || phi_prop(k)>phi_ub(k)
        phi_prop(k) = phi_curr(k);
    else

        theta_prop.phi(k) = phi_prop(k);

        % loglikelihood evaluation
        [loglik_prop]=SEP_loglik(output, H, theta_prop, B_rowin, B_colin, array_rho, array_C, funname);

        % % log prior density of uniform dist
        log_prior_curr = 0;
        log_prior_prop = 0;

        % from phi_curr_X --> phi_prop_X.
        log_pr_tran = -log(phi_prop(k))-(log(phi_prop(k))-log(phi_curr(k)))^2/(2*Delta.phi(k)^2);
        % from phi_prop_X --> phi_curr_X.
        log_pr_rtran = -log(phi_curr(k))-(log(phi_curr(k))-log(phi_prop(k)))^2/(2*Delta.phi(k)^2);


        log_diff=(loglik_prop-loglik_curr)+ (log_prior_prop-log_prior_curr) + (log_pr_rtran-log_pr_tran);

        a_unif=rand(1); 

        if a_unif<=exp(log_diff)
            phi_curr(k)=phi_prop(k);
            loglik_curr=loglik_prop;
        end
    end
    theta_curr.phi(k) = phi_curr(k);

end


end	