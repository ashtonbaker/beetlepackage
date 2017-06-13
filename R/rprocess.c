double e1 = rnorm(0, sigma_1);
double e2 = rnorm(0, sigma_2);
double e3 = rnorm(0, sigma_3);

L = (sqrt(b * A * exp(-cel * L - cea * A)) + e1) *
    (sqrt(b * A * exp(-cel * L - cea * A)) + e1);

P = (sqrt(L * (1 - ul)) + e2) * (sqrt(L * (1 - ul)) + e2);

A = (sqrt(P * exp(-cpa * A) + A * (1 - ua)) + e3) *
    (sqrt(P * exp(-cpa * A) + A * (1 - ua)) + e3);'
