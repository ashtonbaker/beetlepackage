double eps = 0.000001;
if((abs(L_obs - L) > eps) ||
   (abs(P_obs - P) > eps) ||
   (abs(A_obs - A) > eps)) {
  lik = 0;
} else {
  lik = 1;
}
