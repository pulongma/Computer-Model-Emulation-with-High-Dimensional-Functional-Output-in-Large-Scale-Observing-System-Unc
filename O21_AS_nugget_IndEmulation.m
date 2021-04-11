clear all;
close all;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The following routine implements the independent emulator for O2 band with 
% the first PC
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



cd ~/Documents/OCO2/OCO2code


% load NNGP operations
addpath ./code/NNGP_nugget; 
% load Data
addpath ./data;
%%

aaa=csvread('all_together.csv',1); % input data (view geometry+active variables) and output data (FPCA scores)
Id=aaa(:,1); 
Inst_AziA=aaa(:,2);
Inst_ZenA=aaa(:,3);
Sol_AziA=aaa(:,4);
Sol_ZenA=aaa(:,5);
test.id=aaa(:,6); % selecting training index (randomly selected)
s1=aaa(:,7);
s2=aaa(:,8);
s3=aaa(:,9);
s4=aaa(:,10);

O2_1=aaa(:,11);
O2_2=aaa(:,12);
O2_3=aaa(:,13);
SCO2=aaa(:,14);
WCO2=aaa(:,15);

OCO2_scores=[O2_1];
q = size(OCO2_scores, 2);

%all_loc=[Inst_AziA Inst_ZenA Sol_AziA Sol_ZenA s1 s2 s3 s4];
angles = load("standardized_angles.txt");
all_loc = [angles, s1, s2, s3, s4];
Dim_x=size(all_loc,2);

index_p = (find(aaa(:,6) == 1)); % index of prediction


% For MSPE:
OCO2_scores_testing=OCO2_scores(index_p,:);
input_testing=all_loc(index_p,:);

OCO2_training = OCO2_scores;
input_training = all_loc;
OCO2_training(index_p,:)=[];
input_training(index_p,:)=[];

N=size(OCO2_training,1);

%% NNGP specification

coord_order = 5;
[output_sorted, input_sorted, ind] = sort_data(OCO2_training, input_training, coord_order);


%%

Ho=ones(N,1);
Hp=ones(size(input_testing,1),1);
p = size(Ho,2);

m=20; %
%Ne=size(obs_loce,1);

[array_C, array_rho, B_rowin, B_colin, neigh_index]=NeighD_AllD(input_sorted,m);
neigh_index = int32(neigh_index);


[array_C_p, array_rho_p, B_rowin_p, B_colin_p, neigh_index_p]=NeighD_AllDpr(input_sorted,m,input_testing);
neigh_index_p = int32(neigh_index_p);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
Zdata=output_sorted;
Beta = (Ho'*Ho)\(Ho'*Zdata);

Sigma_gls = (Zdata-(Ho*Beta))'*(Zdata-(Ho*Beta))/(N-p);
theta.tau2 = mean(diag(Sigma_gls))/20;
theta.phi = (max(input_sorted) - min(input_sorted))/3;
funname='aniso_matern_5_2';  % used in MATLAB code, see NNGP_Corr_AllD.m
theta.sig2 = 1; % ensure correlation function for input space

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% tau2\sim IG(a1,b1)
a1=2;
b1=0.1;

% phi \sim Unif (a3,b3)
a3=.00001*ones(1,Dim_x);
b3=(max(input_sorted) - min(input_sorted))*10;



prior.tau2_a = a1;
prior.tau2_b = b1;
prior.phi_a = a3;
prior.phi_b = b3;

prior.tau2_bd = [1e-6, 1];
prior.phi_lb = a3;
prior.phi_ub = b3;



%%%%% pre-allocated quantities to hold MCMC samples

niter=5000;

phi_samp=zeros(niter, Dim_x);
tau2_samp=zeros(niter, 1);
loglik_samp=zeros(niter, 1);

tau2_samp(1) = theta.tau2;
phi_samp(1,:) = theta.phi;

%%%% proposal variance
Delta.tau2=0.1; 
Delta.phi = [0.1, 0.15, 0.09, 0.1, 0.05, 0.1, 0.08, 0.08]; 


% initial values for theta
theta_curr = theta;


%% estimate nugget
nugget_est = true;
theta_curr.tau2=theta.tau2;


rng(12345);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%    
tStart = tic;

for iter=1:niter

    if mod(iter,200)==0
        display(['Iteration: ' int2str(iter)]);
        disp([theta_curr.tau2, theta_curr.phi])
    end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
    %%%%% update w
    [theta_curr, loglik_curr] = SEP_MwH_Gibbs_O21(Zdata, Ho, theta_curr, B_rowin, B_colin, prior, Delta, array_rho, ...
                                array_C, funname, nugget_est);

    tau2_samp(iter) = theta_curr.tau2;
    phi_samp(iter, :) = theta_curr.phi;
    loglik_samp(iter) = loglik_curr;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   

end

tEnd = toc(tStart);


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%   
%%%% update predictions
burnin=2000;
phi_post = phi_samp((burnin+1):niter, :);
tau2_post = tau2_samp((burnin+1):niter);

param.phi = phi_post;
param.tau2 = tau2_post;

Zhat = SEP_predict(Zdata, Ho, param, B_rowin, B_colin, array_rho, array_C, ...
        funname, Hp, neigh_index_p, array_rho_p, array_C_p);

mu_Zhat = mean(Zhat, 3);
SD_Zhat = std(Zhat, 0, 3);
quant_Zhat = quantile(Zhat, [0.025, 0.5, 0.975], 3);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

RMSPE = zeros(q,1);
PCI95 = zeros(q,1);
crps = zeros(q,1);

for i =1:q

RMSPE(i) = sqrt(mean((OCO2_scores_testing(:,i)-mu_Zhat(:,i)).^2));
CI025 = quant_Zhat(:,i,1);
CI975 = quant_Zhat(:,i,3);
PCI95(i) = mean(CI025<=OCO2_scores_testing(:,i)& CI975>=OCO2_scores_testing(:,i));
crps(i) = mean(CRPS(OCO2_scores_testing(:,i), mu_Zhat(:,i), SD_Zhat(:,i)));

end 


dlmwrite("./Results/OCO2_O21_AS_IndEmulation.txt", Zhat)

