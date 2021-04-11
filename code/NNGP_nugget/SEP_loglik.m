function prob = SEP_loglik(output, H, theta, B_rowin, B_colin, array_rho, array_C, funname)

tau2 = theta.tau2;
[RInv, dvec, Omega]=build_RInv(B_rowin,B_colin,theta,array_rho,array_C,funname);

n = size(output, 1);
q = size(output, 2);
p = size(H, 2);


RInvH = RInv*H;

[LOmega,~,po] = chol(Omega, 'lower', 'vector');
log_detOmega = 2*sum(log(diag(LOmega)));
OInvH(po,:) = LOmega'\(LOmega\H(po,:));
HRH = OInvH'*RInvH * (1/tau2);

log_detR = n*log(tau2) + sum(log(dvec)) + log_detOmega;

%HRH = H'*RInv*H;

L = chol(HRH, 'lower');
log_detHRH = 2*sum(log(diag(L)));


RInvZ = RInv*output;
OmegaInvZ(po,:) = LOmega'\(LOmega\output(po,:));
ZRZ = OmegaInvZ'*RInvZ*(1/tau2);
HRZ = RInvH'*OmegaInvZ*(1/tau2);
LHRZ = L\HRZ;
ZGZ = ZRZ - LHRZ'*LHRZ;

log_detZGZ = 2*sum(log(diag(chol(ZGZ, 'lower'))));

prob = -0.5*q*log_detR - 0.5*q*log_detHRH - 0.5*(n-p)*log_detZGZ;

end