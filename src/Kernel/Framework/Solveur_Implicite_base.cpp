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

#include <Solveur_Implicite_base.h>
#include <Probleme_base.h>
#include <SChaine.h>
#include <EChaine.h>

Implemente_base(Solveur_Implicite_base, "Solveur_Implicite_base", Objet_U);
// XD solveur_implicite_base objet_u solveur_implicite_base -1 Class for solver in the situation where the time scheme is the implicit scheme. Solver allows equation diffusion and convection operators to be set as implicit terms.

Sortie& Solveur_Implicite_base::printOn(Sortie& os) const { return os; }

Entree& Solveur_Implicite_base::readOn(Entree& is)
{
  int is_seuil_resol_lu = 0;
  double seuil_convergence_solveur_prov = -1.;
  Motcle accferme = "}";
  Motcle accouverte = "{";
  int le_solveur_est_lu = 0;
  Motcle motlu;
  is >> motlu;
  if (motlu != accouverte) Process::exit("We are waiting a { !");
  Nom le_solveur_lu_;
  is >> motlu;
  while (motlu != accferme)
    {
      if (motlu == "max_iter_implicite")
        {
          Cerr << motlu << finl;
          Cerr << motlu << " is no more longer understood by Simpler but by the time scheme" << finl;
          Process::exit();
        }
      else if (motlu == "seuil_convergence_implicite")
        {
          Cerr << motlu << finl;
          is >> param_defaut_.seuil_convergence_implicite();
        }
      else if (motlu == "seuil_generation_solveur")
        {
          Cerr << motlu << finl;
          is_seuil_resol_lu = 1;
          is >> param_defaut_.seuil_generation_solveur();
        }
      else if (motlu == "seuil_verification_solveur")
        {
          Cerr << motlu << finl;
          is >> param_defaut_.seuil_verification_solveur();
        }
      else if (motlu == "seuil_test_preliminaire_solveur")
        {
          Cerr << motlu << finl;
          is >> param_defaut_.seuil_test_preliminaire_solveur();
        }
      else if (motlu == "seuil_convergence_solveur")
        {
          Cerr << motlu << finl;
          is_seuil_resol_lu = 1;
          is >> seuil_convergence_solveur_prov;
        }
      else if (motlu == "nb_it_max")
        {
          is >> param_defaut_.nb_it_max();
        }
      else if (motlu == "seuil_convergence_variable")
        {
          Cerr << "Option seuil_convergence_variable is not yet understood." << finl;
          Process::exit();
          Cerr << motlu << finl;
          facteur_convg_ = 100;
          is_seuil_convg_variable = 1;
        }
      else if (motlu == "controle_residu")
        {
          controle_residu_ = 1;
        }
      else if (motlu == "solveur")
        {
          le_solveur_est_lu = 1;
          SChaine toto;
          // redefinition de motlu pour garder les minuscules/majuscules
          Nom motlubis;
          int nb_acc = 0;
          int ok = 0;
          while (nb_acc != 0 || !ok)
            {
              is >> motlubis;
              toto << " " << motlubis;
              if (motlubis == "}") nb_acc--;
              else if (motlubis == "{")
                {
                  ok = 1;
                  nb_acc++;
                }
            }
          le_solveur_lu_ = Nom(toto.get_str());
        }
      else
        {
          lire(motlu, is);
        }
      is >> motlu;
    }

  if (seuil_convergence_solveur_prov > 0)
    {
      param_defaut_.set_seuil_solveur_avec_seuil_convergence_solveur(seuil_convergence_solveur_prov);
    }
  if ((is_seuil_resol_lu == 0) && (le_solveur_est_lu == 0))
    {
      Cerr << "Neither the solving object nor the threshold seuil_convergence_solveur has been indicated." << finl;
      Cerr << "At least one of them (or both) must be specified." << finl;
      Process::exit();
    }

  if (param_defaut_.seuil_convergence_implicite() < 0)
    param_defaut_.seuil_convergence_implicite() = DMAXFLOAT;

  if (param_defaut_.seuil_verification_solveur() < 0)
    param_defaut_.seuil_verification_solveur() = DMAXFLOAT;
  // on laisse seuil_test_preliminaire <0

  if (le_solveur_est_lu == 0)
    {
      SChaine toto;
      toto << "Gmres { diag seuil " << param_defaut_.seuil_generation_solveur()
           << " nb_it_max " << param_defaut_.nb_it_max()
           << " controle_residu " << get_controle_residu() << " } " << finl;
      le_solveur_lu_ = Nom(toto.get_str());
    }

  {
    EChaine titi(le_solveur_lu_);
    titi >> param_defaut_.solveur();
  }
  param_defaut_.solveur().nommer("solveur_implicite");
  return is;
}

bool Solveur_Implicite_base::iterer_eqs(LIST(OBS_PTR(Equation_base)) eqs, int n, int& ok)
{
  Process::exit("Iterer_eqs non code");
  return false;
}

Entree& Solveur_Implicite_base::lire(const Motcle& motlu, Entree& is)
{
  Cerr << "Keyword : " << motlu << " is not undertood in " << que_suis_je() << finl;
  Process::exit();
  return is;
}

/*! @brief retourne le parametre_implicte de l'equation si il existe si il n'existe pas le cree.
 *
 * .. si les params sont vides on copie ceux du simpler
 *
 */
OWN_PTR(Parametre_equation_base)& Solveur_Implicite_base::get_and_set_parametre_equation(Equation_base& eqn)
{
  OWN_PTR(Parametre_equation_base)& param = eqn.parametre_equation();
  if (param.est_nul())
    {
      param.typer("Parametre_implicite");
    }
  if (!sub_type(Parametre_implicite,param.valeur()))
    {
      Cerr<<eqn.que_suis_je()<<" has parameters of type "<<param.que_suis_je()<<" not coherent with "<<que_suis_je()<<finl;
      exit();
    }

  Parametre_implicite& param_impl = ref_cast(Parametre_implicite,param.valeur());
  // on regarde si il y a des valeurs par defaut a recopier
  if (param_impl.seuil_convergence_implicite()<0)
    param_impl.seuil_convergence_implicite() = param_defaut_.seuil_convergence_implicite();
  if (param_impl.seuil_verification_solveur()<0)
    param_impl.seuil_verification_solveur() = param_defaut_.seuil_verification_solveur();
  if (param_impl.seuil_test_preliminaire_lu()==0)
    param_impl.seuil_test_preliminaire_solveur() = param_defaut_.seuil_test_preliminaire_solveur();
  if (param_impl.solveur().est_nul())
    param_impl.solveur() = param_defaut_.solveur();

  // Some checks:
  if (eqn.probleme().is_coupled())
    {
      if (param_impl.seuil_convergence_implicite()>0.1*DMAXFLOAT)
        {
          Cerr << finl << "Error!" << finl;
          Cerr << "seuil_convergence_implicite option should be defined in your time scheme." << finl;
          Cerr << "It is a mandatory option for a calculation with coupled problems." << finl;
          exit();
        }
    }
  /*
  else if (param_impl.seuil_convergence_implicite()<0.1*DMAXFLOAT)
  {
        Cerr << finl << "Error!" << finl;
  Cerr << "seuil_convergence_implicite option should be NOT be defined in your time scheme" << finl;
  Cerr << "during calculation of a single problem." << finl;
  Cerr << "Remove the option." << finl;
  exit();
  } */
  return param;
}
