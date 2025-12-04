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

#ifndef Equation_Navier_Cauchy_included
#define Equation_Navier_Cauchy_included

#include <Operateur_Diff.h>
#include <Equation_base.h>
#include <Milieu_Elastic.h>
#include <TRUST_Ref.h>

class Milieu_Elastic;

/*! @brief Equation de Navier-Cauchy pour l'elasticite lineaire.
 *
 *  Implementation minimale : toutes les operations numeriques sont a definir.
 */
class Equation_Navier_Cauchy : public Equation_base
{
  Declare_instanciable(Equation_Navier_Cauchy);

public:

  int nombre_d_operateurs() const override { return 1; }
  const Operateur& operateur(int) const override;
  Operateur& operateur(int) override;

  const Champ_Inc_base& inconnue() const override { return deplacement_; }
  Champ_Inc_base& inconnue() override  { return deplacement_; }
  void mettre_a_jour(double temps) override;

  void associer_milieu_base(const Milieu_base&) override;
  const Milieu_base& milieu() const override { return milieu_.valeur(); }
  Milieu_base& milieu() override { return milieu_.valeur(); }
  void discretiser() override;
  void set_param(Param& param) override;
  int lire_motcle_non_standard(const Motcle& mot, Entree& is) override;
  virtual const Champ_Don_base& diffusivite_pour_transport() const { return milieu_->mu_lame(); }
  virtual const Champ_base& diffusivite_pour_pas_de_temps() const { return milieu_->mu_lame(); }
  const Motcle& domaine_application() const override;
  void valider_iteration() override;

private:
  void update_velocity();

  Operateur_Diff terme_diffusif;
  OWN_PTR(Champ_Inc_base) deplacement_;
  OBS_PTR(Milieu_Elastic) milieu_;
  OWN_PTR(Champ_Fonc_base) von_mises_, contraintes_, deformations_, vitesse_noeuds_;
};

#endif
