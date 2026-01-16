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

#include <Milieu_Elasticite.h>
#include <Discretisation_base.h>
#include <Schema_Temps_base.h>
#include <Probleme_base.h>
#include <Domaine_dis_base.h>
#include <Champ_Uniforme.h>
#include <Process.h>
#include <Param.h>

Implemente_instanciable(Milieu_Elasticite,"Milieu_Elasticite",Milieu_base);
// XD Milieu_Elasticite milieu_base Milieu_Elasticite -1 Milieu_Elasticite.
// XD attr rho field_base rho 1 Density (kg.m-3).
// XD attr E field_base E 1 Field E
// XD attr nu field_base nu 1 Field nu.
// XD attr alpha field_base alpha 1 Coeff dilatation.

Sortie& Milieu_Elasticite::printOn(Sortie& os) const { return Milieu_base::printOn(os); }

Entree& Milieu_Elasticite::readOn(Entree& is)
{
  Milieu_base::readOn(is);
  ensure_rho_field();
  creer_champs_non_lus();
  return is;
}

void Milieu_Elasticite::set_param(Param& param)
{
  Milieu_base::set_param(param);
  param.ajouter("E", &ch_E_, Param::REQUIRED);
  param.ajouter("nu", &ch_nu_, Param::REQUIRED);
  param.ajouter("alpha", &ch_coeff_dilatation_th_);
}

void Milieu_Elasticite::creer_champs_non_lus()
{
  Milieu_base::creer_champs_non_lus();

  if (ch_lambda_lame_.est_nul())
    {
      ch_lambda_lame_ = ch_E_;
      ch_lambda_lame_->nommer("lambda_lame");
    }
  if (ch_mu_.est_nul())
    {
      ch_mu_ = ch_E_;
      ch_mu_->nommer("mu_lame");
    }
  if (ch_K_.est_nul())
    {
      ch_K_ = ch_E_;
      ch_K_->nommer("bulk_modulus");
    }
}

void Milieu_Elasticite::discretiser(const Probleme_base& pb, const Discretisation_base& dis)
{
  Milieu_base::discretiser(pb, dis);
  if (zdb_.est_nul()) zdb_ = pb.domaine_dis();

  const Domaine_dis_base& domaine_dis = pb.domaine_dis();
  const int nc = pb.equation(0).inconnue().nb_valeurs_temporelles();

  dis.nommer_completer_champ_physique(domaine_dis, "module_de_Young", "Pa", ch_E_.valeur(), pb);
  dis.nommer_completer_champ_physique(domaine_dis, "coefficient_de_Poisson", "", ch_nu_.valeur(), pb);
  dis.nommer_completer_champ_physique(domaine_dis, "lambda_lame", "Pa", ch_lambda_lame_.valeur(), pb);
  dis.nommer_completer_champ_physique(domaine_dis, "mu_lame", "Pa", ch_mu_.valeur(), pb);
  dis.nommer_completer_champ_physique(domaine_dis, "module_volumique", "Pa", ch_K_.valeur(), pb);
  dis.discretiser_champ("champ_elem", domaine_dis, "masse_volumique_lagrangienne", "kg/m^3", 1, nc, pb.schema_temps().temps_courant(), ch_rho_lag_);
  if (ch_coeff_dilatation_th_.non_nul())
    dis.nommer_completer_champ_physique(domaine_dis, "coeff_dilatation_thermique", "K-1", ch_coeff_dilatation_th_.valeur(), pb);

  ch_rho_lag_->associer_eqn(pb.equation(0));

  champs_compris_.ajoute_champ(ch_E_.valeur());
  champs_compris_.ajoute_champ(ch_nu_.valeur());
  champs_compris_.ajoute_champ(ch_lambda_lame_.valeur());
  champs_compris_.ajoute_champ(ch_mu_.valeur());
  champs_compris_.ajoute_champ(ch_K_.valeur());
  champs_compris_.ajoute_champ(ch_rho_lag_.valeur());
  if (ch_coeff_dilatation_th_.non_nul())
    champs_compris_.ajoute_champ(ch_coeff_dilatation_th_.valeur());

  eq_ = pb.equation(0);
}

int Milieu_Elasticite::initialiser(const double temps)
{
  const int ok = Milieu_base::initialiser(temps);

  ch_E_->initialiser(temps);
  ch_nu_->initialiser(temps);
  ch_lambda_lame_->initialiser(temps);
  ch_mu_->initialiser(temps);
  ch_K_->initialiser(temps);
  if (ch_coeff_dilatation_th_.non_nul()) ch_coeff_dilatation_th_->initialiser(temps);
  ch_rho_lag_->passe() = ch_rho_->valeurs()(0, 0);
  ch_rho_lag_->valeurs() = ch_rho_->valeurs()(0, 0);
  ch_rho_lag_->mettre_a_jour(temps);
  update_fields(temps, false);

  return ok;
}

bool Milieu_Elasticite::initTimeStep(double dt)
{
  if (eq_.est_nul()) throw;
  const Schema_Temps_base& sch = eq_->schema_temps(); //on recupere le schema en temps par la 1ere equation

  for (int i = 1; i <= sch.nb_valeurs_futures(); i++)
    {
      ch_rho_lag_->changer_temps_futur(sch.temps_futur(i), i);
      ch_rho_lag_->futur(i) = ch_rho_lag_->valeurs();
    }
  return true;
}

void Milieu_Elasticite::mettre_a_jour(double temps)
{
  if (ch_E_.est_nul() || ch_nu_.est_nul() || ch_lambda_lame_.est_nul() || ch_mu_.est_nul() || ch_K_.est_nul())
    {
      Cerr << que_suis_je() << " cannot update without E, nu, lambda, mu or K fields." << finl;
      Process::exit();
    }

  Milieu_base::mettre_a_jour(temps);

  ch_E_->mettre_a_jour(temps);
  ch_nu_->mettre_a_jour(temps);
  if (ch_coeff_dilatation_th_.non_nul()) ch_coeff_dilatation_th_->mettre_a_jour(temps);

  ch_lambda_lame_->changer_temps(temps);
  ch_mu_->changer_temps(temps);
  ch_K_->changer_temps(temps);
  ch_rho_lag_->mettre_a_jour(temps);

  update_fields(temps, false);
}

void Milieu_Elasticite::verifier_coherence_champs(int& err, Nom& message)
{
  Milieu_base::verifier_coherence_champs(err, message);

  if (ch_E_->valeurs()(0,0) <= 0.)
    {
      message += "The Young modulus must be strictly positive.\n";
      err = 1;
    }

  const double nu_val = ch_nu_->valeurs()(0,0);
  if (nu_val <= -1.0 || nu_val >= 0.5)
    {
      message += "The Poisson ratio must belong to (-1, 0.5).\n";
      err = 1;
    }

  if (ch_rho_->valeurs()(0,0) <= 0.)
    {
      Cerr << que_suis_je() << " expects rho > 0 but received " << ch_rho_->valeurs()(0,0) << finl;
      Process::exit();
    }
}

void Milieu_Elasticite::ensure_rho_field()
{
  if (ch_rho_.est_nul())
    {
      Cerr << que_suis_je() << " requires a rho field to be specified." << finl;
      Process::exit();
    }
  if (!sub_type(Champ_Uniforme, ch_rho_.valeur()))
    {
      Cerr << que_suis_je() << " currently expects rho to be provided as Champ_Uniforme." << finl;
      Process::exit();
    }

  if (ch_rho_->valeurs().size_totale() == 0)
    {
      Cerr << que_suis_je() << " received an empty rho field." << finl;
      Process::exit();
    }

  if (ch_rho_->nb_comp() != 1)
    {
      Cerr << que_suis_je() << " expects rho to be a scalar field." << finl;
      Process::exit();
    }
}

void Milieu_Elasticite::update_fields(double temps, bool update_rho)
{
  const double E_val = ch_E_->valeurs()(0, 0);
  const double nu_val = ch_nu_->valeurs()(0, 0);
  const double denom_mu = 2. * (1. + nu_val);
  const double denom_lambda = (1. + nu_val) * (1. - 2. * nu_val);
  const double denom_K = 3. * (1. - 2. * nu_val);

  if (denom_mu == 0. || denom_lambda == 0. || denom_K == 0.)
    {
      Cerr << que_suis_je() << " cannot compute elastic constants with nu=" << nu_val << finl;
      Process::exit();
    }

  ch_mu_->valeurs() = E_val / denom_mu;
  ch_lambda_lame_->valeurs() = E_val * nu_val / denom_lambda;
  ch_K_->valeurs() = E_val / denom_K;

  if (update_rho)
    {
      Cerr << "Updating rho_lagrangien field based on current volume scaling at time " << temps << finl;
      ch_rho_lag_->valeurs() = ch_rho_lag_->passe();
      zdb_->domaine().apply_old_to_new_volume_scaling(ch_rho_lag_->valeurs(), zdb_.valeur());
      const Schema_Temps_base& sch = eq_->schema_temps(); //on recupere le schema en temps par la 1ere equation
      for (int i = 1; i <= sch.nb_valeurs_futures(); i++)
        ch_rho_lag_->futur(i) = ch_rho_lag_->valeurs();
    }
}
