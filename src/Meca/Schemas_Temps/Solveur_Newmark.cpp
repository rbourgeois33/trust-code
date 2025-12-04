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

Implemente_instanciable(Solveur_Newmark, "Newmark", Simpler);

Sortie& Solveur_Newmark::printOn(Sortie& os) const { return Simpler::printOn(os); }
Entree& Solveur_Newmark::readOn(Entree& is) { return Simpler::readOn(is); }

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
  if (!sub_type(Equation_Navier_Cauchy, eqn))
    {
      Cerr << que_suis_je() << " cannot be used with equation type " << eqn.que_suis_je() << ". Expected Equation_Navier_Cauchy." << finl;
      Process::exit();
    }

  DoubleTab u_old(current); // previous implicit iterate (u^k-1)

  Equation_Navier_Cauchy& eq = ref_cast(Equation_Navier_Cauchy, eqn);
  SolveurSys& solveur = get_and_set_parametre_implicite(eqn).solveur();

  // Storage for v and a across calls (shape-matched to 'current')
  // (Kept local-static to avoid touching headers; reset if size changes)
  if ((v_n_.size_totale() != current.size_totale()) || (a_n_.size_totale() != current.size_totale()))
    {
      v_n_ = current;
      v_n_ = 0.;
      a_n_ = current;
      a_n_ = 0.;
      v_kp1_ = current;
      v_kp1_ = 0.;
      a_kp1_ = current;
      a_kp1_ = 0.;
      u_pred_ = current;
      v_pred_ = current;
    }

  // Newmark coefficients
  const double a0 = 1.0 / (beta_ * dt * dt);         // scales M on LHS
  const double a1 = gamma_ / (beta_ * dt);           // scales C on LHS (Rayleigh with C = alpha_*M)

  const int N = current.size_totale();
  if (nb_iter == 1)
    {
      v_n_ = v_kp1_;
      a_n_ = a_kp1_;
      // u_pred = u_n + dt*v_n + dt^2*(0.5 - beta_)*a_n
      for (int i = 0; i < N; i++)
        u_pred_.addr()[i] = eqn.inconnue().passe().addr()[i] + dt * v_n_.addr()[i] + dt * dt * (0.5 - beta_) * a_n_.addr()[i];

      for (int i = 0; i < N; i++)
        v_pred_.addr()[i] = v_n_.addr()[i] + dt * (1.0 - gamma_) * a_n_.addr()[i];
    }

// Effective mass scaling handled by Solveur_Masse: ajouter_masse(dt_eff, ...)
// In TRUST, ajouter_masse(dt_eff, M, ...) adds (M / dt_eff) to the matrix.
// We need a0*M on the LHS, hence choose dt_eff = 1 / a0 = beta_ * dt^2.
  const double dt_eff = 1.0 / a0; // = beta_ * dt^2

// Assemble stiffness and external forces at the Newmark-predicted state
// Use 'current' as linearization point set to u_pred before assembly
  if (current.dimension_tot(0) > 0) current = u_pred_;
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
  eq.solv_masse().ajouter_masse(dt_eff, rhs, u_pred_, 0 /*implicit*/, use_old_volumes);

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
        w.addr()[i] = a1 * u_pred_.addr()[i] - v_pred_.addr()[i];
      const double dt_eff_C_rhs = 1.0 / alpha_;
      eq.solv_masse().ajouter_masse(dt_eff_C_rhs, rhs, w, 0 /*implicit*/, use_old_volumes);
    }

// Apply boundary conditions
  eq.modifier_pour_Cl(matrice, rhs);

// Solve for u_{n+1}
  solveur->reinit();
  solveur.resoudre_systeme(matrice, rhs, current); // current := u_{n+1}

// Under-relaxation on displacement for fixed-point coupled problems
  const std::string eqn_str = Motcle(eqn.que_suis_je()).getString();
  const bool needs_relaxation = (eqn.probleme().is_coupled() && sub_type(Probleme_Couple_Point_Fixe, eqn.probleme().get_pb_couple()) && relax_factors_.count(eqn_str));
  if (needs_relaxation)
    {
      const double omega = relax_factors_.at(eqn_str);
      for (int i = 0; i < N; i++)
        current.addr()[i] = (1.0 - omega) * u_old.addr()[i] + omega * current.addr()[i];
    }

// Recover a_{n+1} and v_{n+1}
// a_{n+1} = a0*(u_{n+1} - u_pred)
  for (int i = 0; i < N; i++)
    a_kp1_.addr()[i] = a0 * (current.addr()[i] - u_pred_.addr()[i]);
// v_{n+1} = v_pred + gamma_*dt*a_{n+1}
  for (int i = 0; i < N; i++)
    v_kp1_.addr()[i] = v_pred_.addr()[i] + gamma_ * dt * a_kp1_.addr()[i];


  DoubleTab dudt(current);
  dudt -= u_old; // dudt = u^k - u^{k-1}
  double dudt_norme = mp_norme_vect(dudt);
  const double seuil_convg = get_and_set_parametre_implicite(eqn).seuil_convergence_implicite();

  const bool converge = (dudt_norme < seuil_convg);
  Cout << eqn.que_suis_je() << (converge ? " is " : " is not ") << "converged at the implicit iteration " << nb_iter << " ( ||uk-uk-1|| = " << dudt_norme << (converge ? "<" : ">") << " implicit threshold " << seuil_convg << " )" << finl;

  eq.valider_iteration();
  ok = 1;
  solveur->reinit();
  return converge ? 1 : 0;
}
