function X = rmvt(mu, A, df)
% sampling from a multivariate t distribution with X ~ MVT(mu, A*A', df)
d = size(A,1);
Z = randn(d,1);
W = df/chi2rnd(df);

X = mu + sqrt(W)*A*Z;

end
