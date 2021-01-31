#ifndef SET_UP_G_H
#define SET_UP_G_H

using namespace NEdyson;

// Fills 2x2 GF and SE with functions which obey 
// the relations:
// C^<(t,t') = - C^<(t',t)^\dagger
// C^R(t,t') = - C^R(t',t)^\dagger
// C^R(t,t') = C^A(t',t)^\dagger
// C^\rceil(t, \tau) = -\xi C^\lceil(\beta-\tau, t)^\dagger
// These relations are used explicitly in the energy_conv function
void set_up_G_herm(DYSON &D, GREEN &G, GREEN &S, double dt, double beta) {
  int Nt = G.nt();
  int Ntau = G.ntau();
  int sig = G.sig();
  cplx cplxi = cplx(0., 1.);
  
  // S^R(t,tp) = ((i t^2 tp^2, i t tp), (i t tp, i t^3 tp^3))
  for(int t = 0; t <= Nt; t++) {
    for(int tp = 0; tp <= t; tp++) {
      double dt2 = dt * dt;
      int t2 = t * t;
      int tp2 = tp * tp;
      S.retptr(t,tp)[0] = cplxi * (dt2 * dt2 * t2 * tp2);
      S.retptr(t,tp)[1] = cplxi * (dt2 * t * tp);
      S.retptr(t,tp)[2] = cplxi * (dt2 * t * tp);
      double dt3 = dt2 * dt;
      int t3 = t2 * t;
      int tp3 = tp2 * tp;
      S.retptr(t,tp)[3] = cplxi * (dt3 * dt3 * t3 * tp3);
    }
  }

  // G^<(t,tp) = ((i exp(i(t-tp)), t-tp), (t-tp, i exp(3i(t-tp))))
  for(int tp = 0; tp <= Nt; tp++) {
    for(int t = 0; t <= tp; t++) {
      G.lesptr(t,tp)[0] = cplxi * std::exp(cplxi * (double)(t-tp) * dt);
      G.lesptr(t,tp)[1] = (t-tp) * dt;
      G.lesptr(t,tp)[2] = (t-tp) * dt;
      G.lesptr(t,tp)[3] = cplxi * std::exp(3. * cplxi * (double)(t-tp) * dt);
    }
  }

  // S^<(t,tp) = ((i t tp, i exp(2 i (t-tp))), (i exp(2 i exp(2 i (t-tp))), t-tp))
  for(int tp = 0; tp <= Nt; tp++) {
    for(int t = 0; t <= tp; t++) {
      double dt2 = dt * dt;
      S.lesptr(t,tp)[0] = cplxi * (dt2 * t * tp);
      S.lesptr(t,tp)[1] = cplxi * std::exp(2. * cplxi * (dt * (t-tp)));
      S.lesptr(t,tp)[2] = cplxi * std::exp(2. * cplxi * (dt * (t-tp)));
      S.lesptr(t,tp)[3] = dt * (t-tp);
    }
  }

  // G^R(t,tp) = ((-i cos(3(t-tp)), -sin(2(t-tp))), (-sin(2(t-tp)), -icos(t-tp)))
  for(int t = 0; t <= Nt; t++) {
    for(int tp = 0; tp <= t; tp++) {
      G.retptr(t,tp)[0] = -cplxi * std::cos(3. * dt * (t-tp));
      G.retptr(t,tp)[1] = -std::sin(2. * dt * (t-tp));
      G.retptr(t,tp)[2] = -std::sin(2. * dt * (t-tp));
      G.retptr(t,tp)[3] = -cplxi * std::cos(dt * (t-tp));
    }
  }

  // G^\rceil(t,tau) = ((\tau^2, i t), (exp(i tau), t^2 - itau))
  for(int t = 0; t <= Nt; t++) {
    for(int tau = 0; tau <= Ntau; tau++) {
      double x = D.Convolution().collocation().x_i()(tau);
      double taupt = beta * (x + 1) / 2.;
      G.tvptr(t,tau)[0] = taupt * taupt;
      G.tvptr(t,tau)[1] = cplxi * (dt * t);
      G.tvptr(t,tau)[2] = std::exp(cplxi * taupt);
      G.tvptr(t,tau)[3] = dt * dt * t * t - cplxi * taupt;
    }
  }

  // S^\rceil(t,tau) = ((t \tau, t^2 \tau), (sin(t), i(t - tau)))
  for(int t = 0; t <= Nt; t++) {
    for(int tau = 0; tau <= Ntau; tau++) {
      double x = D.Convolution().collocation().x_i()(tau);
      double taupt = beta * (x + 1) / 2.;
      S.tvptr(t,tau)[0] = dt * t * taupt;
      S.tvptr(t,tau)[1] = dt * dt * t * t * taupt;
      S.tvptr(t,tau)[2] = std::sin(dt * t);
      S.tvptr(t,tau)[3] = cplxi * (dt * t - taupt);
    }
  }

  // G^M(tau) = ((tau, i tau), (tau^2, exp(itau)))
  for(int tau = 0; tau <= Ntau; tau++) {
    double x = D.Convolution().collocation().x_i()(tau);
    double taupt = beta * (x + 1) / 2.;
    G.matptr(tau)[0] = taupt;
    G.matptr(tau)[1] = cplxi * taupt;
    G.matptr(tau)[2] = taupt * taupt;
    G.matptr(tau)[3] = std::exp(cplxi * taupt);
  }
  
  // S^M(tau) = ((i tau, tau), (sin(tau^2), exp(itau)))
  for(int tau = 0; tau <= Ntau; tau++) {
    double x = D.Convolution().collocation().x_i()(tau);
    double taupt = beta * (x + 1) / 2.;
    S.matptr(tau)[0] = taupt * cplxi;
    S.matptr(tau)[1] = taupt;
    S.matptr(tau)[2] = std::sin(taupt * taupt);
    S.matptr(tau)[3] = std::exp(cplxi * taupt);
  }
}


// Not every two time quantity on the contour obeys the \ddag symmetry defined in the nessi
// paper.  Therefore we must sometimes use two GREEN objects to store the extra information
// needed in these cases
void set_up_G(DYSON &D, GREEN &G, GREEN &Gcc, GREEN &S, GREEN &Scc, double dt, double beta) {
  int Nt = G.nt();
  int Ntau = G.ntau();
  int Gsig = G.sig();
  int Ssig = S.sig();
  cplx cplxi = cplx(0., 1.);
  
  // S^R(t,tp) = ((t tp^2, it-tp), (sin(t), cos(tp))
  for(int t = 0; t <= Nt; t++) {
    for(int tp = 0; tp <= t; tp++) {
      S.retptr(t,tp)[0] = dt * dt * dt * t * tp * tp;
      S.retptr(t,tp)[1] = cplxi * (dt * t) - dt * tp;
      S.retptr(t,tp)[2] = std::sin(dt * t);
      S.retptr(t,tp)[3] = std::cos(dt * tp);

      Scc.retptr(t,tp)[0] = -std::conj(dt * dt * dt * tp * t * t);
      Scc.retptr(t,tp)[2] = -std::conj(cplxi * (dt * tp) - dt * t);
      Scc.retptr(t,tp)[1] = -std::conj(std::sin(dt * tp));
      Scc.retptr(t,tp)[3] = -std::conj(std::cos(dt * t));
    }
  }

  // S^<(t,tp) = (t^2 tp, t+2tp), (exp(it), exp(i(t-tp))))
  for(int tp = 0; tp <= Nt; tp++) {
    for(int t = 0; t <= tp; t++) {
      S.lesptr(t,tp)[0] = dt * dt * dt * t * t * tp;
      S.lesptr(t,tp)[1] = dt * t + 2. * dt * tp;
      S.lesptr(t,tp)[2] = std::exp(cplxi * (dt * t));
      S.lesptr(t,tp)[3] = std::exp(cplxi * (double)(t-tp) * dt);

      Scc.lesptr(t,tp)[0] = -std::conj(dt * dt * dt * tp * tp * t);
      Scc.lesptr(t,tp)[2] = -std::conj(dt * tp + 2. * dt * t);
      Scc.lesptr(t,tp)[1] = -std::conj(std::exp(cplxi * (dt * tp)));
      Scc.lesptr(t,tp)[3] = -std::conj(std::exp(cplxi * (double)(tp-t) * dt));
    }
  }

  // S^\rceil(t,tau) = ((t tau, cos(t tau)), (exp(i tau), t-tau))
  // S^\lceil(tau,t) = ((i - t tau, sin(2 tau)), (3tau+t, it))
  for(int t = 0; t <= Nt; t++) {
    for(int tau = 0; tau <= Ntau; tau++) {
      double x = D.Convolution().collocation().x_i()(tau);
      double taupt = beta * (x + 1) / 2.;

      S.tvptr(t,tau)[0] = dt * t * taupt;
      S.tvptr(t,tau)[1] = std::cos(dt * t * taupt);
      S.tvptr(t,tau)[2] = std::exp(cplxi * taupt);
      S.tvptr(t,tau)[3] = dt * t - taupt;

      Scc.tvptr(t,tau)[0] = - (double)Ssig * std::conj(cplxi - dt * t * (beta-taupt));
      Scc.tvptr(t,tau)[2] = - (double)Ssig * std::conj(std::sin(2 * (beta-taupt)));
      Scc.tvptr(t,tau)[1] = - (double)Ssig * std::conj(3*(beta-taupt) + dt * t);
      Scc.tvptr(t,tau)[3] = - (double)Ssig * std::conj(dt * t * cplxi);
    }
  }

  // S^M(tau) = ((1, tau), (2, i tau^2))
  for(int tau = 0; tau <= Ntau; tau++) {
    double x = D.Convolution().collocation().x_i()(tau);
    double taupt = beta * (x + 1) / 2.;

    S.matptr(tau)[0] = 1;
    S.matptr(tau)[1] = taupt;
    S.matptr(tau)[2] = 2;
    S.matptr(tau)[3] = cplxi * taupt * taupt;
  }


  // G^R(t,tp) = ((t+tp, t+2tp), (t+3tp, icos(tp)))
  for(int t = 0; t <= Nt; t++) {
    for(int tp = 0; tp <= t; tp++) {
      G.retptr(t,tp)[0] = dt * (t + tp);
      G.retptr(t,tp)[1] = dt * (t + 2 * tp);
      G.retptr(t,tp)[2] = dt * (t + 3 * tp);
      G.retptr(t,tp)[3] = cplxi * std::cos(dt * tp);

      Gcc.retptr(t,tp)[0] = -std::conj(dt * (tp + t));
      Gcc.retptr(t,tp)[2] = -std::conj(dt * (tp + 2 * t));
      Gcc.retptr(t,tp)[1] = -std::conj(dt * (tp + 3 * t));
      Gcc.retptr(t,tp)[3] = -std::conj(cplxi * std::cos(dt * t));
    }
  }

  // G^<(t,tp) = ((i t, i tp), (t-tp, 4))
  for(int tp = 0; tp <= Nt; tp++) {
    for(int t = 0; t <= tp; t++) {
      G.lesptr(t,tp)[0] = t * dt * cplxi;
      G.lesptr(t,tp)[1] = tp * dt * cplxi;
      G.lesptr(t,tp)[2] = dt * (t-tp);
      G.lesptr(t,tp)[3] = 4;

      Gcc.lesptr(t,tp)[0] = -std::conj(tp * dt * cplxi);
      Gcc.lesptr(t,tp)[2] = -std::conj(t * dt * cplxi);
      Gcc.lesptr(t,tp)[1] = -std::conj(dt * (tp-t));
      Gcc.lesptr(t,tp)[3] = -std::conj(4);
    }
  }

  // G^\rceil(t,tau) = ((t, i2), (tau, tau^2))
  // G^\lceil(tau,t) = ((t tau, itau), (i-tau, t+tau))
  for(int t = 0; t <= Nt; t++) {
    for(int tau = 0; tau <= Ntau; tau++) {
      double x = D.Convolution().collocation().x_i()(tau);
      double taupt = beta * (x + 1) / 2.;

      G.tvptr(t,tau)[0] = dt * t;
      G.tvptr(t,tau)[1] = cplxi * 2.;
      G.tvptr(t,tau)[2] = taupt;
      G.tvptr(t,tau)[3] = taupt * taupt;

      Gcc.tvptr(t,tau)[0] = - (double)Gsig * std::conj(dt * t * (beta - taupt));
      Gcc.tvptr(t,tau)[2] = - (double)Gsig * std::conj(cplxi * (beta - taupt));
      Gcc.tvptr(t,tau)[1] = - (double)Gsig * std::conj(cplxi - (beta - taupt));
      Gcc.tvptr(t,tau)[3] = - (double)Gsig * std::conj(dt * t + (beta - taupt));
    }
  }

  // G^M(tau) = ((exp(i tau), exp(itau)), (cos(tau), i))
  for(int tau = 0; tau <= Ntau; tau++) {
    double x = D.Convolution().collocation().x_i()(tau);
    double taupt = beta * (x + 1) / 2.;

    G.matptr(tau)[0] = std::exp(cplxi * taupt);
    G.matptr(tau)[1] = std::exp(cplxi * taupt);
    G.matptr(tau)[2] = std::cos(taupt);
    G.matptr(tau)[3] = cplxi;
  }
}


void set_up_G_herm(DYSON &D, TTI_GREEN &G, TTI_GREEN &S, double dt, double beta) {
  int Nt = G.nt();
  int Ntau = G.ntau();
  int Gsig = G.sig();
  int Ssig = S.sig();
  cplx cplxi = cplx(0., 1.);
  
  // S^R(t) = ((sint, t), (t, ie^(it)))
  for(int t = 0; t <= Nt; t++) {
    S.retptr(t)[0] = std::sin(dt * t);
    S.retptr(t)[1] = dt * t;
    S.retptr(t)[2] = dt * t;
    S.retptr(t)[3] = cplxi * std::exp(dt * t * cplxi);
  }

  // G^R(t) = ((icost, it^2), (it^2, iexp(i3tp)))
  for(int t = 0; t <= Nt; t++) {
    G.retptr(t)[0] = cplxi * std::cos(dt * t);
    G.retptr(t)[1] = dt * dt * t * t * cplxi;
    G.retptr(t)[2] = dt * dt * t * t * cplxi;
    G.retptr(t)[3] = cplxi * std::exp(3. * dt * t * cplxi);
  }

  // S^L(t) = ((t exp(it),  i), (i, sint))
  for(int t = 0; t >= -Nt; t--) {
    S.lesptr(t)[0] = dt * t * std::exp(dt * t * cplxi);
    S.lesptr(t)[1] = cplxi;
    S.lesptr(t)[2] = cplxi;
    S.lesptr(t)[3] = std::sin(dt * t);
  }

  // G^L(t) = ((i expit, 2), (2, t))
  for(int t = 0; t >= -Nt; t--) {
    G.lesptr(t)[0] = cplxi * std::exp(dt * t * cplxi);
    G.lesptr(t)[1] = 2;
    G.lesptr(t)[2] = 2;
    G.lesptr(t)[3] = dt * t;
  }

  // S^\rceil(t,tau) = ((t tau, cos(t tau)), (exp(i tau), t-tau))
  for(int t = 0; t <= Nt; t++) {
    for(int tau = 0; tau <= Ntau; tau++) {
      double x = D.Convolution().collocation().x_i()(tau);
      double taupt = beta * (x + 1) / 2.;

      S.tvptr(t,tau)[0] = dt * t * taupt;
      S.tvptr(t,tau)[1] = std::cos(dt * t * taupt);
      S.tvptr(t,tau)[2] = std::exp(cplxi * taupt);
      S.tvptr(t,tau)[3] = dt * t - taupt;
    }
  }

  // S^M(tau) = ((1, tau), (2, i tau^2))
  for(int tau = 0; tau <= Ntau; tau++) {
    double x = D.Convolution().collocation().x_i()(tau);
    double taupt = beta * (x + 1) / 2.;

    S.matptr(tau)[0] = 1;
    S.matptr(tau)[1] = taupt;
    S.matptr(tau)[2] = 2;
    S.matptr(tau)[3] = cplxi * taupt * taupt;
  }


  // G^\rceil(t,tau) = ((t, i2), (tau, tau^2))
  for(int t = 0; t <= Nt; t++) {
    for(int tau = 0; tau <= Ntau; tau++) {
      double x = D.Convolution().collocation().x_i()(tau);
      double taupt = beta * (x + 1) / 2.;

      G.tvptr(t,tau)[0] = dt * t;
      G.tvptr(t,tau)[1] = cplxi * 2.;
      G.tvptr(t,tau)[2] = taupt;
      G.tvptr(t,tau)[3] = taupt * taupt;
    }
  }

  // G^M(tau) = ((exp(i tau), exp(itau)), (cos(tau), i))
  for(int tau = 0; tau <= Ntau; tau++) {
    double x = D.Convolution().collocation().x_i()(tau);
    double taupt = beta * (x + 1) / 2.;

    G.matptr(tau)[0] = std::exp(cplxi * taupt);
    G.matptr(tau)[1] = std::exp(cplxi * taupt);
    G.matptr(tau)[2] = std::cos(taupt);
    G.matptr(tau)[3] = cplxi;
  }
}

void set_up_G(DYSON &D, TTI_GREEN &G, TTI_GREEN &Gcc, TTI_GREEN &S, TTI_GREEN &Scc, double dt, double beta) {
  int Nt = G.nt();
  int Ntau = G.ntau();
  int Gsig = G.sig();
  int Ssig = S.sig();
  cplx cplxi = cplx(0., 1.);
  
  // S^R(t) = ((sint, t), (t^2, e^(it)))
  for(int t = 0; t <= Nt; t++) {
    S.retptr(t)[0] = std::sin(dt * t);
    S.retptr(t)[1] = dt * t;
    S.retptr(t)[2] = dt * t * dt * t;
    S.retptr(t)[3] = std::exp(dt * t * cplxi);

    Scc.retptr(t)[0] = -std::conj(std::sin(dt * -t));
    Scc.retptr(t)[2] = -std::conj(-dt * t);
    Scc.retptr(t)[1] = -std::conj(dt * -t * dt * -t);
    Scc.retptr(t)[3] = -std::conj(std::exp(dt * -t * cplxi));
  }

  // G^R(t) = ((icost, i), (sin(t)cps(t), iexp(i3tp)))
  for(int t = 0; t <= Nt; t++) {
    G.retptr(t)[0] = cplxi * std::cos(dt * t);
    G.retptr(t)[1] = cplxi;
    G.retptr(t)[2] = std::sin(dt * t) * std::cos(dt * t);
    G.retptr(t)[3] = cplxi * std::exp(3. * dt * t * cplxi);

    Gcc.retptr(t)[0] = -std::conj( cplxi * std::cos(dt * -t));
    Gcc.retptr(t)[2] = -std::conj( cplxi);
    Gcc.retptr(t)[1] = -std::conj( std::sin(dt * -t) * std::cos(dt * -t));
    Gcc.retptr(t)[3] = -std::conj( cplxi * std::exp(3. * dt * -t * cplxi));
  }

  // S^L(t) = ((exp(it),  i-t), (t, t^2 cost))
  for(int t = 0; t >= -Nt; t--) {
    S.lesptr(t)[0] = std::exp(dt * t * cplxi);
    S.lesptr(t)[1] = cplxi - dt * t;
    S.lesptr(t)[2] = dt * t;
    S.lesptr(t)[3] = dt * dt * t * t * std::cos(dt * t);

    Scc.lesptr(t)[0] = -std::conj( std::exp(dt * -t * cplxi));
    Scc.lesptr(t)[2] = -std::conj( cplxi - dt * -t);
    Scc.lesptr(t)[1] = -std::conj( dt * -t);
    Scc.lesptr(t)[3] = -std::conj( dt * dt * -t * -t * std::cos(dt * -t));
  }

  // G^L(t) = ((i expit, 2+it), (2, (t-i)^2))
  for(int t = 0; t >= -Nt; t--) {
    G.lesptr(t)[0] = cplxi * std::exp(dt * t * cplxi);
    G.lesptr(t)[1] = 2. + dt * t * cplxi;
    G.lesptr(t)[2] = 2.;
    G.lesptr(t)[3] = (dt * t - cplxi) * (dt * t - cplxi);

    Gcc.lesptr(t)[0] = -std::conj( cplxi * std::exp(dt * -t * cplxi));
    Gcc.lesptr(t)[2] = -std::conj( 2. - dt * t * cplxi);
    Gcc.lesptr(t)[1] = -std::conj( 2.);
    Gcc.lesptr(t)[3] = -std::conj( (-dt * t - cplxi) * (-dt * t - cplxi));
  }

  // S^\rceil(t,tau) = ((t tau, cos(t tau)), (exp(i tau), t-tau))
  // S^\lceil(tau,t) = ((i - t tau, sin(2 tau)), (3tau+t, it))
  for(int t = 0; t <= Nt; t++) {
    for(int tau = 0; tau <= Ntau; tau++) {
      double x = D.Convolution().collocation().x_i()(tau);
      double taupt = beta * (x + 1) / 2.;

      S.tvptr(t,tau)[0] = dt * t * taupt;
      S.tvptr(t,tau)[1] = std::cos(dt * t * taupt);
      S.tvptr(t,tau)[2] = std::exp(cplxi * taupt);
      S.tvptr(t,tau)[3] = dt * t - taupt;

      Scc.tvptr(t,tau)[0] = - (double)Ssig * std::conj(cplxi - dt * t * (beta-taupt));
      Scc.tvptr(t,tau)[2] = - (double)Ssig * std::conj(std::sin(2. * (beta-taupt)));
      Scc.tvptr(t,tau)[1] = - (double)Ssig * std::conj(3.*(beta-taupt) + dt * t);
      Scc.tvptr(t,tau)[3] = - (double)Ssig * std::conj(dt * t * cplxi);
    }
  }

  // S^M(tau) = ((1, tau), (2, i tau^2))
  for(int tau = 0; tau <= Ntau; tau++) {
    double x = D.Convolution().collocation().x_i()(tau);
    double taupt = beta * (x + 1) / 2.;

    S.matptr(tau)[0] = 1;
    S.matptr(tau)[1] = taupt;
    S.matptr(tau)[2] = 2;
    S.matptr(tau)[3] = cplxi * taupt * taupt;
  }

  // G^\rceil(t,tau) = ((t, i2), (tau, tau^2))
  // G^\lceil(tau,t) = ((t tau, itau), (i-tau, t+tau))
  for(int t = 0; t <= Nt; t++) {
    for(int tau = 0; tau <= Ntau; tau++) {
      double x = D.Convolution().collocation().x_i()(tau);
      double taupt = beta * (x + 1) / 2.;

      G.tvptr(t,tau)[0] = dt * t;
      G.tvptr(t,tau)[1] = cplxi * 2.;
      G.tvptr(t,tau)[2] = taupt;
      G.tvptr(t,tau)[3] = taupt * taupt;

      Gcc.tvptr(t,tau)[0] = - (double)Gsig * std::conj(dt * t * (beta - taupt));
      Gcc.tvptr(t,tau)[2] = - (double)Gsig * std::conj(cplxi * (beta - taupt));
      Gcc.tvptr(t,tau)[1] = - (double)Gsig * std::conj(cplxi - (beta - taupt));
      Gcc.tvptr(t,tau)[3] = - (double)Gsig * std::conj(dt * t + (beta - taupt));
    }
  }

  // G^M(tau) = ((exp(i tau), exp(itau)), (cos(tau), i))
  for(int tau = 0; tau <= Ntau; tau++) {
    double x = D.Convolution().collocation().x_i()(tau);
    double taupt = beta * (x + 1) / 2.;

    G.matptr(tau)[0] = std::exp(cplxi * taupt);
    G.matptr(tau)[1] = std::exp(cplxi * taupt);
    G.matptr(tau)[2] = std::cos(taupt);
    G.matptr(tau)[3] = cplxi;
  }
}
#endif
