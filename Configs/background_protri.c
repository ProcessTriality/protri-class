
/* =========================================================================================
   background_protri.c
   Erweiterung für CLASS - Integration der Prozess-Trialität-Felder
   Autor: Maximilian Rupp (2025)
   ========================================================================================= */

#include "common.h"
#include <math.h>

/* =========================================================================================
   Parameterstruktur für ProTri-Felder
   ========================================================================================= */
struct protri {
  double epsilon;     /* Kopplungsparameter */
  double Phi_T, Phi_S, Phi_P; /* Feldwerte */
  double dPhi_T, dPhi_S, dPhi_P; /* Zeitableitungen */
  double Hesse[3][3]; /* Hesse-Matrix */
  double signH[3][3]; /* signierte Hesse (effektiv) */
};

/* =========================================================================================
   Hilfsfunktionen
   ========================================================================================= */

/* Potenzial */
static double V(double Phi_T, double Phi_S, double Phi_P, double eps) {
  return 0.5 * (Phi_T*Phi_T + Phi_S*Phi_S - Phi_P*Phi_P)
       + eps * Phi_T * Phi_S * Phi_P;
}

/* Gradienten des Potentials */
static void dV(double Phi_T, double Phi_S, double Phi_P, double eps,
               double *dVdT, double *dVdS, double *dVdP) {
  *dVdT = Phi_T + eps * Phi_S * Phi_P;
  *dVdS = Phi_S + eps * Phi_T * Phi_P;
  *dVdP = -Phi_P + eps * Phi_T * Phi_S;
}

/* Cantor-Projektion (Selbstähnliche Kopplung) */
static void cantor_projection(double *Phi_T, double *Phi_S, double *Phi_P) {
  double PT = *Phi_T, PS = *Phi_S, PP = *Phi_P;
  *Phi_T = exp(-0.1 * PT) * cosh(PS);
  *Phi_S = tanh(PT + PP);
  *Phi_P = sin(PS - PT);
}

/* Signiertes Quadrat des Logarithmus */
static double signed_square(double x) {
  if (x == 0) return 0.0;
  double s = (x > 0) ? 1.0 : -1.0;
  return s * pow(log(fabs(x)), 2);
}

/* =========================================================================================
   Hesse-Matrix (effektiv)
   ========================================================================================= */
static void compute_hesse(struct protri *p) {
  double PT = p->Phi_T, PS = p->Phi_S, PP = p->Phi_P, eps = p->epsilon;
  double H[3][3] = {
    {1.0, eps*PP, eps*PS},
    {eps*PP, 1.0, eps*PT},
    {eps*PS, eps*PT, -1.0}
  };
  for (int i=0;i<3;i++){
    for (int j=0;j<3;j++){
      p->Hesse[i][j] = H[i][j];
      p->signH[i][j] = signed_square(H[i][j]);
    }
  }
}

/* =========================================================================================
   Integration eines Zeitschritts (Runge-Kutta 4)
   ========================================================================================= */
static void protri_step(struct protri *p, double dt) {
  double k1_T, k1_S, k1_P;
  double k2_T, k2_S, k2_P;
  double k3_T, k3_S, k3_P;
  double k4_T, k4_S, k4_P;
  double dVdT, dVdS, dVdP;

  dV(p->Phi_T, p->Phi_S, p->Phi_P, p->epsilon, &dVdT, &dVdS, &dVdP);
  k1_T = p->dPhi_T;
  k1_S = p->dPhi_S;
  k1_P = p->dPhi_P;

  k2_T = p->dPhi_T - 0.5*dt*dVdT;
  k2_S = p->dPhi_S - 0.5*dt*dVdS;
  k2_P = p->dPhi_P - 0.5*dt*dVdP;

  k3_T = p->dPhi_T - 0.5*dt*dVdT;
  k3_S = p->dPhi_S - 0.5*dt*dVdS;
  k3_P = p->dPhi_P - 0.5*dt*dVdP;

  k4_T = p->dPhi_T - dt*dVdT;
  k4_S = p->dPhi_S - dt*dVdS;
  k4_P = p->dPhi_P - dt*dVdP;

  p->Phi_T += (dt/6.0) * (k1_T + 2*k2_T + 2*k3_T + k4_T);
  p->Phi_S += (dt/6.0) * (k1_S + 2*k2_S + 2*k3_S + k4_S);
  p->Phi_P += (dt/6.0) * (k1_P + 2*k2_P + 2*k3_P + k4_P);

  p->dPhi_T -= dt * dVdT;
  p->dPhi_S -= dt * dVdS;
  p->dPhi_P += dt * dVdP;

  cantor_projection(&p->Phi_T, &p->Phi_S, &p->Phi_P);
  compute_hesse(p);
}

/* =========================================================================================
   Hauptfunktion
   ========================================================================================= */
int evolve_protri(double a_start, double a_end, int Nsteps, double eps,
                  double *Phi_T_out, double *Phi_S_out, double *Phi_P_out) {
  struct protri p;
  p.epsilon = eps;
  p.Phi_T = 1.0; p.Phi_S = 1.0; p.Phi_P = 1.0;
  p.dPhi_T = 0.0; p.dPhi_S = 0.0; p.dPhi_P = 0.0;
  double dt = (a_end - a_start)/Nsteps;

  for (int i=0;i<Nsteps;i++) {
    protri_step(&p, dt);
  }

  *Phi_T_out = p.Phi_T;
  *Phi_S_out = p.Phi_S;
  *Phi_P_out = p.Phi_P;
  return _SUCCESS_;
}
