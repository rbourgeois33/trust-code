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

#include <Equation_Navier_Cauchy.h>
#include <Discretisation_base.h>
#include <Schema_Temps_base.h>
#include <Probleme_base.h>
#include <Domaine_dis_base.h>
#include <Process.h>

Implemente_instanciable(Equation_Navier_Cauchy,"Equation_Navier_Cauchy",Equation_base);

Sortie& Equation_Navier_Cauchy::printOn(Sortie& os) const { return Equation_base::printOn(os); }

Entree& Equation_Navier_Cauchy::readOn(Entree& is)
{
  Equation_base::readOn(is);

  // terme_diffusif.set_fichier("Contrainte_visqueuse");
  // terme_diffusif.set_description("Friction drag exerted by the fluid=Integral(-mu*(grad(u) +grad(u)^T)*ndS) [N] if SI units used");

  return is;
}

void Equation_Navier_Cauchy::set_param(Param& param)
{
  Equation_base::set_param(param);
  param.ajouter_non_std("diffusion",(this));
}

int Equation_Navier_Cauchy::lire_motcle_non_standard(const Motcle& mot, Entree& is)
{
  if (mot=="diffusion")
    {
      Cerr << "Reading and typing of the diffusion operator : " << finl;
      terme_diffusif.associer_diffusivite(diffusivite_pour_transport());
      is >> terme_diffusif;
      terme_diffusif.associer_diffusivite_pour_pas_de_temps(diffusivite_pour_pas_de_temps());
      terme_diffusif.associer_diffusivite(milieu_->mu_lame());
      terme_diffusif.associer_diffusivite_volumique(milieu_->lambda_lame());
      return 1;
    }
  else
    return Equation_base::lire_motcle_non_standard(mot,is);
}

const Operateur& Equation_Navier_Cauchy::operateur(int i) const
{
  if (i != 0) Process::exit();
  return terme_diffusif;
}

Operateur& Equation_Navier_Cauchy::operateur(int i)
{
  if (i != 0) Process::exit();
  return terme_diffusif;
}

void Equation_Navier_Cauchy::associer_milieu_base(const Milieu_base& mil)
{
  if (!sub_type(Milieu_Elastic, mil))
    {
      Cerr << que_suis_je() << " expects a Milieu_Elastic but received a " << mil.que_suis_je() << finl;
      Process::exit();
    }
  milieu_ = ref_cast(Milieu_Elastic, mil);
}

void Equation_Navier_Cauchy::discretiser()
{
  const Discretisation_base& dis = discretisation();
  const Domaine_dis_base& dom = domaine_dis();
  const Schema_Temps_base& sch = schema_temps();

  const int nb_val_temp = sch.nb_valeurs_temporelles();
  const double temps = sch.temps_courant();

  dis.discretiser_champ("vitesse", dom, "deplacement", "m", dimension, nb_val_temp, temps, deplacement_);
  dis.discretiser_champ("champ_elem", dom, "von_mises", "Pa", 1, temps, von_mises_);
  dis.discretiser_champ("champ_elem", dom, "contraintes", "Pa", 3, temps, contraintes_);
  dis.discretiser_champ("champ_elem", dom, "deformations", "", 3, temps, deformations_);
  dis.discretiser_champ("vitesse", dom, "vitesse_noeuds", "m/s", dimension, temps, vitesse_noeuds_);
  deformations_->fixer_nom_compo(0, bidim_axi ? "eps_r" : "eps_xx");
  deformations_->fixer_nom_compo(1, bidim_axi ? "eps_z" : "eps_yy");
  deformations_->fixer_nom_compo(2, bidim_axi ? "eps_theta" : "eps_zz");
  contraintes_->fixer_nom_compo(0, bidim_axi ? "sigma_r" : "sigma_xx");
  contraintes_->fixer_nom_compo(1, bidim_axi ? "sigma_z" : "sigma_yy");
  contraintes_->fixer_nom_compo(2, bidim_axi ? "sigma_theta" : "sigma_zz");
  champs_compris_.ajoute_champ(deplacement_);
  champs_compris_.ajoute_champ(von_mises_);
  champs_compris_.ajoute_champ(contraintes_);
  champs_compris_.ajoute_champ(deformations_);
  champs_compris_.ajoute_champ(vitesse_noeuds_);
  terme_diffusif.associer_eqn(*this);

  Equation_base::discretiser();
}

const Motcle& Equation_Navier_Cauchy::domaine_application() const
{
  static Motcle domaine = "Mecanique";
  return domaine;
}

void Equation_Navier_Cauchy::valider_iteration()
{
  Equation_base::valider_iteration();
  update_velocity();
}

void Equation_Navier_Cauchy::update_velocity()
{
  const double dt = schema_temps().pas_de_temps();
  const DoubleTab& disp_n = deplacement_->valeurs();
  const DoubleTab& disp_nm1 = deplacement_->passe();
  DoubleTab& vit_n = vitesse_noeuds_->valeurs();

  for (int i = 0; i < vit_n.dimension_tot(0); i++)
    for (int j = 0; j < dimension; j++)
      vit_n(i, j) = (disp_n(i, j) - disp_nm1(i, j)) / dt;
}

void Equation_Navier_Cauchy::mettre_a_jour(double temps)
{
  Equation_base::mettre_a_jour(temps);

  von_mises_->changer_temps(temps);
  contraintes_->changer_temps(temps);
  deformations_->changer_temps(temps);
  vitesse_noeuds_->changer_temps(temps);
  update_velocity();
  ref_cast(Operateur_Diff_base, terme_diffusif.l_op_base()).calculer_von_mises(deplacement_->valeurs(), deformations_->valeurs(), contraintes_->valeurs(), von_mises_->valeurs());
}
