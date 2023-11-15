
L = [1 5 10 50 100];
B_zero = zeros(length(L),1);
B_zero_tera = zeros (length(L),1);

beta_2 = -20;
beta_2_ps = beta_2*1e-24;

C = 1/sqrt(3)

X_1 = -abs(C)+sqrt(1+C^2);

for c = 1:length(L)
    B_zero(c) = 1/sqrt(abs(beta_2_ps)*L(c)*X_1);
end