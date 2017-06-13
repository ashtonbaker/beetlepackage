step.model <- function(data, params, statenames, paramnames, globals,
                       initializer, rprocess, dmeasure, rmeasure,
                       fromEstimationScale, toEstimationScale,
                       times = 'weeks',
                       t0 = 0,
                       delta.t = 1/7,
                       nstages = c(E=7, L=7, P=7)) {

  if (missing(statenames)) {
    statenames <- c(
      sprintf('E%d', 1:nstages['E']),
      sprintf('L%d', 1:nstages['L']),
      sprintf('P%d', 1:nstages['P']),
      'A', 'A_prev', 'P_prev'
      )
  }

  if (missing(paramnames)) {
    paramnames <- c(
      'b', 'cea', 'cel', 'cpa', 'cpa_force',
      'mu_A', 'mu_L', 'mu_A_force',
      'tau_E', 'tau_L', 'tau_P',
      'od'
      )
  }

  if (missing(globals)) {
    globals <- Csnippet(
      sprintf(
        "
        #include <math.h>
        #define ESTAGES %d
        #define LSTAGES %d
        #define PSTAGES %d

        #define L_0 250
        #define P_0 5
        #define A_0 100
        ",
        nstages['E'],
        nstages['L'],
        nstages['P'])
      )
  }

  if (missing(initializer)) {
    initializer <- Csnippet(
      "
      double *E = &E1;
      double *L = &L1;
      double *P = &P1;

      double gamma_L = (LSTAGES / tau_L) * (1 - mu_L);
      double gamma_P = (PSTAGES / tau_P) * exp((-cpa * A_0) / ESTAGES);

      double mu_l = (LSTAGES / tau_L) * mu_L;
      double mu_p = (PSTAGES / tau_P) * (1 - exp((-cpa * A_0) / ESTAGES));

      double L_rate[LSTAGES] = {0};
      double P_rate[PSTAGES] = {0};

      int k;
      double sum;
      for (k = 0, sum = 0; k < LSTAGES; k++){
      L_rate[k] = pow(gamma_L/(gamma_L + mu_l), k);
      sum += L_rate[k];
      }
      for (k = 0; k < LSTAGES; k++) L_rate[k] /= sum;
      for (k = LSTAGES - 1, sum = 0; k >=0; k--){
      sum += L_rate[k];
      L_rate[k] /= sum;
      }

      for (k = 0, sum = 0; k < PSTAGES; k++){
      P_rate[k] = pow(gamma_P/(gamma_P + mu_p), k);
      sum += P_rate[k];
      }

      for (k = 0; k < PSTAGES; k++) P_rate[k] /= sum;
      for (k = PSTAGES - 1, sum = 0; k >=0; k--){
      sum += P_rate[k];
      P_rate[k] /= sum;
      }


      for (k = 0; k < ESTAGES; k++) E[k] = 0;

      int L_count = L_0;
      for (k = 0; k < LSTAGES - 1; k++){
      L[k] = rbinom(L_count, L_rate[k]);
      L_count -= L[k];
      }
      L[LSTAGES - 1] = L_count;

      int P_count = P_0;
      for (k = 0; k < PSTAGES - 1; k++){
      P[k] = rbinom(P_count, P_rate[k]);
      P_count -= P[k];
      }

      P[PSTAGES - 1] = P_count;

      A = 100;
      "
      )
  }

  if (missing(rprocess)) {
    rprocess <- discrete.time.sim(
      step.fun = Csnippet(
        "
        double *E = &E1;
        double *L = &L1;
        double *P = &P1;

        int time = round(t * 7);

        int k;
        double L_tot = 0;
        for (k = 0; k < LSTAGES; k++) L_tot += L[k];

        double gamma_E = (ESTAGES / tau_E) *
        exp((-cel * L_tot - cea * A) / ESTAGES);
        double gamma_L = (LSTAGES / tau_L) * (1 - mu_L);
        double gamma_P = (PSTAGES / tau_P) * exp((-cpa * A) / PSTAGES);

        double mu_e = (ESTAGES / tau_E) - gamma_E;
        double mu_l = (LSTAGES / tau_L) - gamma_L;
        double mu_p = (PSTAGES / tau_P) - gamma_P;

        double etrans[2*ESTAGES], ltrans[2*LSTAGES], ptrans[2*PSTAGES], adeath;

        // Calculate who goes where
        for (k = 0; k < ESTAGES; k++) {
        // Eggs growing to next stage
        etrans[2*k]   = rbinom(E[k], gamma_E);

        // Eggs dying
        etrans[2*k+1] = rbinom(E[k]-etrans[2*k], mu_e/(1 - gamma_E) );
        }

        for (k = 0; k < LSTAGES; k++) {
        // Larvae growing to next stage
        ltrans[2*k]   = rbinom(L[k], gamma_L);

        // Larvae dying
        ltrans[2*k+1] = rbinom(L[k]-ltrans[2*k], mu_l/(1 - gamma_L));
        }

        for (k = 0; k < PSTAGES; k++) {
        // Pupae growing to next stage
        ptrans[2*k]   = rbinom(P[k], gamma_P);

        // Pupae dying
        ptrans[2*k+1] = rbinom(P[k]-ptrans[2*k], mu_p/(1 - gamma_P) );
        }

        adeath = rbinom(A, mu_A);

        // Bookkeeping
        E[0] += rpois(b*A); // oviposition

        for (k = 0; k < ESTAGES; k++) {
        // Subtract eggs that die or progress
        E[k] -= (etrans[2*k]+etrans[2*k+1]);

        // Add eggs that arrive from previous E stage.
        E[k+1] += etrans[2*k]; // E[ESTAGES] == L[0]!!
        }

        for (k = 0; k < LSTAGES; k++) {
        // Subtract larvae that die or progress
        L[k] -= (ltrans[2*k]+ltrans[2*k+1]);

        // Add larvae that arrive from previous E stage.
        L[k+1] += ltrans[2*k]; // L[LSTAGES] == P[0]!!
        }

        for (k = 0; k < PSTAGES; k++) {
        // Subtract pupae that die or progress
        P[k] -= (ptrans[2*k]+ptrans[2*k+1]);

        // Add pupae that arrive from previous E stage.
        P[k+1] += ptrans[2*k]; // P[PSTAGES] == A[0]!!
        }

        A -= adeath;

        if ((time % 14 == 0) && (time != 0) && (mu_A_force > 0.00001)) {
        double P_tot = 0;

        for (k = 0; k < PSTAGES; k++) P_tot += P[k];

        double A_pred = round((1 - mu_A_force) * A_prev) +
        round(P_prev * exp(-cpa_force * A));

        if (A_pred < A) {
        double A_sub = fmin(A - A_pred, A_prev);
        A = fmax(A - A_sub, 0);
        }

        P_prev = P_tot;
        A_prev = A;
        }

        "),
      delta.t = delta.t
      )
  }

  if (missing(dmeasure)) {
    dmeasure <- Csnippet(
      "
      const double *L = &L1;
      const double *P = &P1;
      double fudge = 1e-9;

      int k;
      double L_tot = 0;
      double P_tot = 0;
      for (k = 0; k < LSTAGES; k++) L_tot += L[k];
      for (k = 0; k < PSTAGES; k++) P_tot += P[k];

      lik = dnbinom_mu(L_obs, 1/od, L_tot + fudge, 1) +
      dnbinom_mu(P_obs, 1/od, P_tot + fudge, 1) +
      dnbinom_mu(A_obs, 1/od, A + fudge,     1);

      lik = (give_log) ? lik : exp(lik);
      "
      )
  }

  if (missing(rmeasure)) {
    rmeasure <- Csnippet(
        "
        const double *L = &L1;
        const double *P = &P1;
        double fudge = 1e-9;

        int k;
        double L_tot = 0;
        double P_tot = 0;
        for (k = 0; k < LSTAGES; k++) L_tot += L[k];
        for (k = 0; k < PSTAGES; k++) P_tot += P[k];

        L_obs = rnbinom_mu(1/od, L_tot + fudge);
        P_obs = rnbinom_mu(1/od, P_tot + fudge);
        A_obs = rnbinom_mu(1/od, A + fudge);
        "
      )
  }

  if (missing(fromEstimationScale)) {
    fromEstimationScale <- Csnippet(
      "
      Tb = exp(b);
      Tcea = expit(cea);
      Tcel = expit(cel);
      Tcpa = expit(cpa);
      Tmu_A = expit(mu_A);
      Tmu_L = expit(mu_L);
      Ttau_E = ESTAGES+exp(tau_E);
      Ttau_L = LSTAGES+exp(tau_L);
      Ttau_P = PSTAGES+exp(tau_P);
      Tod = exp(od);
      "
      )
  }

  if (missing(toEstimationScale)) {
    toEstimationScale <- Csnippet(
      "
      Tb = log(b);
      Tcea = logit(cea);
      Tcel = logit(cel);
      Tcpa = logit(cpa);
      Tmu_A = logit(mu_A);
      Tmu_L = logit(mu_L);
      Ttau_E = log(tau_E-ESTAGES);
      Ttau_L = log(tau_L-LSTAGES);
      Ttau_P = log(tau_P-PSTAGES);
      Tod = log(od);
      "
      )
  }

  if (missing(data)) {
    data <- beetledata
  }

  pomp(
    data  = data,
    times = times,
    t0 = t0,
    statenames = statenames,
    paramnames = paramnames,
    globals = globals,
    initializer = initializer,
    rprocess = rprocess,
    dmeasure = dmeasure,
    rmeasure = rmeasure,
    toEstimationScale   = toEstimationScale,
    fromEstimationScale = fromEstimationScale,
    params = params
    )
}
