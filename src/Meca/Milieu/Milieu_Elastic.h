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

#ifndef Milieu_Elastic_included
#define Milieu_Elastic_included

#include <Milieu_base.h>

/*! @brief Milieu elastique isotrope lineaire.
 *
 *  Stocke les proprietes elastiques au sein de champs pour rester coherent avec la hierarchie TRUST.
 */
class Milieu_Elastic : public Milieu_base
{
  Declare_instanciable(Milieu_Elastic);

public:
  void set_param(Param& param) override;
  void creer_champs_non_lus() override;
  void discretiser(const Probleme_base& pb, const Discretisation_base& dis) override;
  int initialiser(const double temps) override;
  void mettre_a_jour(double temps) override;

  const Champ_Don_base& young_modulus() const { return ch_E_.valeur(); }
  Champ_Don_base& young_modulus() { return ch_E_.valeur(); }

  const Champ_Don_base& poisson_ratio() const { return ch_nu_.valeur(); }
  Champ_Don_base& poisson_ratio() { return ch_nu_.valeur(); }

  const Champ_Don_base& lambda_lame() const { return ch_lambda_lame_.valeur(); }
  Champ_Don_base& lambda_lame() { return ch_lambda_lame_.valeur(); }

  const Champ_Don_base& mu_lame() const { return ch_mu_.valeur(); }
  Champ_Don_base& mu_lame() { return ch_mu_.valeur(); }

  const Champ_Don_base& bulk_modulus() const { return ch_K_.valeur(); }
  Champ_Don_base& bulk_modulus() { return ch_K_.valeur(); }

  const Champ_Don_base& thermal_expansion() const { return ch_coeff_dilatation_th_.valeur(); }
  Champ_Don_base& thermal_expansion() { return ch_coeff_dilatation_th_.valeur(); }

  const Champ_Don_base& rho_lagrangien() const { return ch_rho_lag_.valeur(); }
  Champ_Don_base& rho_lagrangien() { return ch_rho_lag_.valeur(); }

protected:
  void verifier_coherence_champs(int& err, Nom& message) override;

private:
  void ensure_rho_field();
  void update_fields(double temps);
  double last_update_ = -1.0;

  OWN_PTR(Champ_Don_base) ch_E_;
  OWN_PTR(Champ_Don_base) ch_rho_lag_;
  OWN_PTR(Champ_Don_base) ch_nu_;
  OWN_PTR(Champ_Don_base) ch_lambda_lame_;
  OWN_PTR(Champ_Don_base) ch_mu_;
  OWN_PTR(Champ_Don_base) ch_K_;
  OWN_PTR(Champ_Don_base) ch_coeff_dilatation_th_;
};

#endif
