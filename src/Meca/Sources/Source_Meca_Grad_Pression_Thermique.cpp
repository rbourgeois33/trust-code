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

#include <Source_Meca_Grad_Pression_Thermique.h>
#include <Milieu_Elastic.h>
#include <Domaine_EF.h>
#include <Equation_base.h>
#include <Champ_Uniforme.h>
#include <Param.h>

Implemente_instanciable(Source_Meca_Grad_Pression_Thermique, "Source_Meca_Grad_Pression_Thermique_EF", Source_base);

Sortie& Source_Meca_Grad_Pression_Thermique::printOn(Sortie& os) const { return Source_base::printOn(os); }

Entree& Source_Meca_Grad_Pression_Thermique::readOn(Entree& is)
{
  Param param(que_suis_je());
  param.ajouter("reference_temperature_field", &T_ref_, Param::REQUIRED);
  param.ajouter("temperature_field", &T_, Param::REQUIRED);
  param.lire_avec_accolades_depuis(is);
  return is;
}

int Source_Meca_Grad_Pression_Thermique::initialiser(double temps)
{
  T_ref_->initialiser(temps);
  equation().discretisation().nommer_completer_champ_physique(equation().domaine_dis(), "temperature_field", "", T_.valeur(), equation().probleme());
  T_->initialiser(temps);
  return Source_base::initialiser(temps);
}

void Source_Meca_Grad_Pression_Thermique::mettre_a_jour(double temps)
{
  T_ref_->mettre_a_jour(temps);
  T_->mettre_a_jour(temps);
}

DoubleTab& Source_Meca_Grad_Pression_Thermique::ajouter(DoubleTab& resu) const
{
  const Domaine_EF& domaine_ef = ref_cast(Domaine_EF, equation().domaine_dis());
  const IntTab& elems = domaine_ef.domaine().les_elems();
  const DoubleTab& Bij_thilde = domaine_ef.Bij_thilde();
  const DoubleTab& IPhi_thilde = domaine_ef.IPhi_thilde();

  const int nb_elem_tot = domaine_ef.nb_elem_tot();
  const int nb_som_elem = domaine_ef.domaine().nb_som_elem();
  const int D = dimension;

  const DoubleTab& val_alpha = ref_cast(Milieu_Elastic, equation().milieu()).thermal_expansion().valeurs();
  const DoubleTab& val_K     = ref_cast(Milieu_Elastic, equation().milieu()).bulk_modulus().valeurs();

  const DoubleTab& val_T     = T_->valeurs();
  const DoubleTab& val_Tref  = T_ref_->valeurs();

  const int cT = sub_type(Champ_Uniforme, T_.valeur()) ? 1 : 0;
  const int cTref = sub_type(Champ_Uniforme, T_ref_.valeur()) ? 1 : 0;

  // Assemble: for each element e, scalar s_e = 3 K α (T - Tref),
  // add to each node i and component d: resu(node,d) += s_e * ∫_Ω ∂N_i/∂x_d dV = s_e * Bij_thilde(e,i,d)
  for (int e = 0; e < nb_elem_tot; e++)
    {
      const double alpha_e = val_alpha(0);
      const double K_e     = val_K(0);
      const double dT_e    = val_T(e * !cT) - val_Tref(e * !cTref);
      const double s_e     = 3.0 * K_e * alpha_e * dT_e; // [Pa]

      if (s_e == 0.) continue;

      for (int i = 0; i < nb_som_elem; i++)
        {
          const int s = elems(e, i);
          // Base strong-form (all components): + s_e * ∫ ∂N_i/∂x_d dV
          for (int d = 0; d < D; d++)
            resu(s, d) += s_e * Bij_thilde(e, i, d);

          // Axisymmetric fix: add + s_e * ∫ (N_i / r) dV on the radial component
          if (bidim_axi)
            resu(s, 0) += s_e * IPhi_thilde(e, i) / domaine_ef.xp(e, 0);
        }
    }
  return resu;
}
