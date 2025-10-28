/****************************************************************************
* Copyright (c) 2025, CEA
* All rights reserved.
*
* Redistribution and use in source and binary forms, with or without modification, are permitted provided that the following conditions are met:
* 1. Redistributions of source code must retain the above copyright notice, this list of conditions and the following disclaimer.
* 2. Redistributions in binary form must reproduce the above copyright notice, this list of conditions and the following disclaimer in the documentation and/or other materials provided with the distribution.
* 3. Neither the name of the copyright holder nor the names of its contributors may be used to endorse or promote products derived from this software without specific prior written permission.
*
* THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.
* IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS;
* OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
*
*****************************************************************************/

#include <Solveur_Newmark.h>
#include <Equation_Navier_Cauchy.h>
#include <Solveur_Masse_base.h>
#include <Schema_Temps_base.h>
#include <SolveurSys.h>
#include <TRUSTTrav.h>
#include <Matrice_Morse.h>
#include <Probleme_base.h>
#include <Discretisation_base.h>
#include <Domaine_dis_base.h>
#include <Process.h>
#include <algorithm>

Implemente_instanciable(Solveur_Newmark, "Newmark", Simple);

Sortie& Solveur_Newmark::printOn(Sortie& os) const { return Simple::printOn(os); }
Entree& Solveur_Newmark::readOn(Entree& is) { return Simple::readOn(is); }

Entree& Solveur_Newmark::lire(const Motcle& motlu,Entree& is)
{
  Motcles les_mots(3);
  {
    les_mots[0] = "beta";
    les_mots[1] = "gamma";
    les_mots[2] = "alpha"; // C = alpha * M damping
  }

  int rang = les_mots.search(motlu);
  switch(rang)
    {
    case 0:
      {
        is >> beta_;
        break;
      }
    case 1:
      {
        is >> gamma_;
        break;
      }
    case 2:
      {
        is >> alpha_;
        break;
      }
    default :
      {
        Cerr << "Keyword : " << motlu << " is not understood in " << que_suis_je() << finl;
        exit();
      }
    }
  return is;
}

bool Solveur_Newmark::iterer_eqn(Equation_base& eqn, const DoubleTab& inut, DoubleTab& current, double dt, int nb_iter, int& ok)
{
  // Newmark-β average-acceleration (β=1/4, γ=1/2) single implicit step
  // u: unknown (displacement), v: velocity, a: acceleration
  // K_eff = K + a0*M + a1*C, with a0 = 1/(β*dt^2) and a1 = γ/(β*dt)
  // Here C = alpha*M (mass-proportional Rayleigh damping).
  // RHS  = f_{n+1} + M * (a0 * u_pred) + C * (a1 * u_pred - v_pred)

  if (!sub_type(Equation_Navier_Cauchy, eqn))
    {
      Cerr << que_suis_je() << " cannot be used with equation type " << eqn.que_suis_je() << ". Expected Equation_Navier_Cauchy." << finl;
      Process::exit();
    }
//   eqn.solv_masse().set_name_of_coefficient_temporel("masse_volumique");

  Equation_Navier_Cauchy& eq = ref_cast(Equation_Navier_Cauchy, eqn);
  SolveurSys& solveur = get_and_set_parametre_implicite(eqn).solveur();

  // Storage for v and a across calls (shape-matched to 'current')
  // (Kept local-static to avoid touching headers; reset if size changes)
  if ((v_old_.size_totale() != current.size_totale()) || (a_old_.size_totale() != current.size_totale()))
    {
      v_old_ = current;
      v_old_ = 0.;
      a_old_ = current;
      a_old_ = 0.;
    }

  // Previous displacement (u_n) from the equation unknown
  const DoubleTab& u_old = eq.inconnue().passe();

  // Newmark coefficients
  const double a0 = 1.0 / (beta_ * dt * dt);         // scales M on LHS
  const double a1 = gamma_ / (beta_ * dt);           // scales C on LHS (Rayleigh with C = alpha_*M)

  // Predictors
  DoubleTrav u_pred(current);
  u_pred = u_old; // start from previous displacement
  // u_pred = u_n + dt*v_n + dt^2*(0.5 - beta_)*a_n
  const int N = current.size_totale();
  for (int i = 0; i < N; i++)
    u_pred.addr()[i] += dt * v_old_.addr()[i] + dt * dt * (0.5 - beta_) * a_old_.addr()[i];

  // v_pred = v_n + dt*(1 - gamma_)*a_n
  DoubleTrav v_pred(current);
  v_pred = v_old_; // start from previous velocity
  for (int i = 0; i < N; i++)
    v_pred.addr()[i] += dt * (1.0 - gamma_) * a_old_.addr()[i];

// Effective mass scaling handled by Solveur_Masse: ajouter_masse(dt_eff, ...)
// In TRUST, ajouter_masse(dt_eff, M, ...) adds (M / dt_eff) to the matrix.
// We need a0*M on the LHS, hence choose dt_eff = 1 / a0 = beta_ * dt^2.
  const double dt_eff = 1.0 / a0; // = beta_ * dt^2

// Assemble stiffness and external forces at the Newmark-predicted state
// Use 'current' as linearization point set to u_pred before assembly
  current = u_pred;
  Matrice_Morse matrice;
  eq.dimensionner_matrice(matrice);
  matrice.get_set_coeff() = 0;

  DoubleTrav rhs(current);
  rhs = 0.;
  statistics().begin_count(STD_COUNTERS::matrix_assembly,statistics().get_last_opened_counter_level()+1);
  eq.assembler(matrice, current, rhs); // K and external terms at u_pred
  statistics().end_count(STD_COUNTERS::matrix_assembly);

// Add effective mass to matrix and RHS (M*a0*u_pred)
  eq.solv_masse().ajouter_masse(dt_eff, matrice, 0 /*implicit*/);
  const bool use_old_volumes = false; //eq.domaine_dis().domaine().deformable();
  eq.solv_masse().ajouter_masse(dt_eff, rhs, u_pred, 0 /*implicit*/, use_old_volumes);

// Add damping contributions if alpha_ != 0: K_eff += a1 * C = a1 * alpha_ * M
// and RHS += C * (a1 * u_pred - v_pred) = alpha_ * M * (a1 * u_pred - v_pred)
  if (alpha_ != 0.)
    {
      // Matrix: add (a1 * alpha_) * M  -> achieved with dt_eff_C_mat = 1 / (a1 * alpha_)
      if (a1 != 0.)
        {
          const double dt_eff_C_mat = 1.0 / (a1 * alpha_);
          eq.solv_masse().ajouter_masse(dt_eff_C_mat, matrice, 0 /*implicit*/);
        }
      // RHS: add alpha_ * M * (a1 * u_pred - v_pred) -> use dt_eff_C_rhs = 1 / alpha_
      DoubleTrav w(current);
      // w = a1 * u_pred - v_pred
      for (int i = 0; i < N; i++)
        w.addr()[i] = a1 * u_pred.addr()[i] - v_pred.addr()[i];
      const double dt_eff_C_rhs = 1.0 / alpha_;
      eq.solv_masse().ajouter_masse(dt_eff_C_rhs, rhs, w, 0 /*implicit*/, use_old_volumes);
    }

// Apply boundary conditions
  eq.modifier_pour_Cl(matrice, rhs);

// Solve for u_{n+1}
  solveur->reinit();
  solveur.resoudre_systeme(matrice, rhs, current); // current := u_{n+1}

// Recover a_{n+1} and v_{n+1}
// a_{n+1} = a0*(u_{n+1} - u_pred)
  DoubleTrav a_new(current);
  for (int i = 0; i < N; i++)
    a_new.addr()[i] = a0 * (current.addr()[i] - u_pred.addr()[i]);
// v_{n+1} = v_pred + gamma_*dt*a_{n+1}
  DoubleTrav v_new(current);
  for (int i = 0; i < N; i++)
    v_new.addr()[i] = v_pred.addr()[i] + gamma_ * dt * a_new.addr()[i];

// Store for next step
  v_old_ = v_new;
  a_old_ = a_new;

// Project/adjust solution if mass solver needs it (kept from original code)
// eq.solv_masse().corriger_solution(current, current);

  eq.valider_iteration();
  ok = 1;
  solveur->reinit();
  return ok;
}
