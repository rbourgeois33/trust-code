/****************************************************************************
* Copyright (c) 2026, CEA
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

Implemente_instanciable(Solveur_Newmark, "Newmark", Solveur_non_lineaire);
// XD newmark solveur_implicite_base solveur_newmark 1 Newmark implicit solver for the resolution of the linear elastodynamic equation.
// XD attr seuil_convergence_implicite floattant seuil_convergence_implicite 0 Keyword to set the value of the convergence criteria for the resolution of the implicit system build to solve either the Navier_Stokes equation (only for Simple and Simpler algorithms) or a scalar equation. It is adviced to use the default value (1e6) to solve the implicit system only once by time step. This value must be decreased when a coupling between problems is considered.
// XD attr seuil_convergence_solveur floattant seuil_convergence_solveur 1 value of the convergence criteria for the resolution of the implicit system build by solving several times per time step the Navier_Stokes equation and the scalar equations if any. This value MUST be used when a coupling between problems is considered (should be set to a value typically of 0.1 or 0.01).
// XD attr seuil_generation_solveur floattant seuil_generation_solveur 1 Option to create a GMRES solver and use vrel as the convergence threshold (implicit linear system Ax=B will be solved if residual error ||Ax-B|| is lesser than vrel).
// XD attr seuil_verification_solveur floattant seuil_verification_solveur 1 Option to check if residual error ||Ax-B|| is lesser than vrel after the implicit linear system Ax=B has been solved.
// XD attr seuil_test_preliminaire_solveur floattant seuil_test_preliminaire_solveur 1 Option to decide if the implicit linear system Ax=B should be solved by checking if the residual error ||Ax-B|| is bigger than vrel.
// XD attr solveur solveur_sys_base solveur 1 Method (different from the default one, Gmres with diagonal preconditioning) to solve the linear system.
// XD attr nb_it_max entier nb_it_max 1 Keyword to set the maximum iterations number for the Gmres.
// XD attr controle_residu rien controle_residu 1 Keyword of Boolean type (by default 0). If set to 1, the convergence occurs if the residu suddenly increases.
// XD attr alpha floattant alpha 1 Damping coefficient for Rayleigh damping: C = alpha * M
// XD attr beta floattant beta 1 Newmark parameter beta (default value 0.25 for average acceleration method).
// XD attr gamma floattant gamma 1 Newmark parameter gamma (default value 0.5 for average acceleration method).

Sortie& Solveur_Newmark::printOn(Sortie& os) const { return Solveur_non_lineaire::printOn(os); }
Entree& Solveur_Newmark::readOn(Entree& is) { return Solveur_non_lineaire::readOn(is); }

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
  eqn.dimensionner_matrice(matrice);
  matrice.get_set_coeff() = 0;

  DoubleTrav rhs(current);
  rhs = 0.;
  statistics().begin_count(STD_COUNTERS::matrix_assembly,statistics().get_last_opened_counter_level()+1);
  eqn.assembler(matrice, current, rhs); // K and external terms at u_pred
  statistics().end_count(STD_COUNTERS::matrix_assembly);

// Add effective mass to matrix and RHS (M*a0*u_pred)
  eqn.solv_masse().ajouter_masse(dt_eff, matrice, 0 /*implicit*/);
  eqn.solv_masse().ajouter_masse(dt_eff, rhs, u_pred_, 0 /*implicit*/);

// Add damping contributions if alpha_ != 0: K_eff += a1 * C = a1 * alpha_ * M
// and RHS += C * (a1 * u_pred - v_pred) = alpha_ * M * (a1 * u_pred - v_pred)
  if (alpha_ != 0.)
    {
      // Matrix: add (a1 * alpha_) * M  -> achieved with dt_eff_C_mat = 1 / (a1 * alpha_)
      if (a1 != 0.)
        {
          const double dt_eff_C_mat = 1.0 / (a1 * alpha_);
          eqn.solv_masse().ajouter_masse(dt_eff_C_mat, matrice, 0 /*implicit*/);
        }
      // RHS: add alpha_ * M * (a1 * u_pred - v_pred) -> use dt_eff_C_rhs = 1 / alpha_
      DoubleTrav w(current);
      // w = a1 * u_pred - v_pred
      for (int i = 0; i < N; i++)
        w.addr()[i] = a1 * u_pred_.addr()[i] - v_pred_.addr()[i];
      const double dt_eff_C_rhs = 1.0 / alpha_;
      eqn.solv_masse().ajouter_masse(dt_eff_C_rhs, rhs, w, 0 /*implicit*/);
    }

// Apply boundary conditions
  eqn.modifier_pour_Cl(matrice, rhs);

// Solve for u_{n+1}
  solveur->reinit();
  solveur.resoudre_systeme(matrice, rhs, current); // current := u_{n+1}

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

  eqn.valider_iteration();
  ok = 1;
  solveur->reinit();
  return converge ? 1 : 0;
}
