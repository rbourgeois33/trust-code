/****************************************************************************
* Copyright (c) 2024, CEA
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

#ifndef TRUSTSchema_RK_TPP_included
#define TRUSTSchema_RK_TPP_included

template <Ordre_RK _ORDRE_ > template<Ordre_RK _O_>
std::enable_if_t<_O_ == Ordre_RK::DEUX_WILLIAMSON || _O_ == Ordre_RK::TROIS_WILLIAMSON || _O_ == Ordre_RK::QUATRE_WILLIAMSON, int>
TRUSTSchema_RK<_ORDRE_>::faire_un_pas_de_temps_eqn_base_generique(Equation_base& eqn)
{
  static constexpr int NB_PTS = (_ORDRE_ == Ordre_RK::DEUX_WILLIAMSON) ? 2 : 3;

  // Warning sur les 100 premiers pas de temps si facsec est egal a 1 pour faire reflechir l'utilisateur
  if (nb_pas_dt() >= 0 && nb_pas_dt() <= NW && facsec_ == 1) print_warning(NW);

  DoubleTab& xi = eqn.inconnue().valeurs(), &xip1 = eqn.inconnue().futur();
  DoubleTrav present(xi), qi(xi);
  present = xi;
  qi = xi;
  const auto& coeff_a = get_a<_ORDRE_,NB_PTS>();
  const auto& coeff_b = get_b<_ORDRE_,NB_PTS>();
  init_rk_flux_accumulators(eqn);

  for (int i = 0; i < NB_PTS; i++)
    {
      // on fait ca : q_i = a_{i-1} * q_{i-1} + dt * f(x_{i-1})
      qi *= coeff_a[i];
      qi.ajoute(dt_, eqn.derivee_en_temps_inco(xip1));
      accumulate_flux_low_storage(eqn, coeff_a[i], coeff_b[i]);

      // on fait ca : x_i = x_{i-1} + b_i * q_i
      xi.ajoute(coeff_b[i], qi);
    }

  xip1 = xi;
  xi -= present;
  xi /= dt_;
  update_critere_statio(xi, eqn);

  // Update boundary condition on futur:
  eqn.domaine_Cl_dis().imposer_cond_lim(eqn.inconnue(), temps_courant() + pas_de_temps());
  xip1.echange_espace_virtuel();

  // xi = x0;
  xi = present;
  finalize_rk_flux_accumulators(eqn);

  return 1;
}

template <Ordre_RK _ORDRE_ > template<Ordre_RK _O_>
std::enable_if_t<_O_ == Ordre_RK::DEUX_CLASSIQUE || _O_ == Ordre_RK::TROIS_CLASSIQUE || _O_ == Ordre_RK::QUATRE_CLASSIQUE || _O_ == Ordre_RK::QUATRE_CLASSIQUE_3_8, int>
TRUSTSchema_RK<_ORDRE_>::faire_un_pas_de_temps_eqn_base_generique(Equation_base& eqn)
{
  // Warning sur les 100 premiers pas de temps si facsec est egal a 1 pour faire reflechir l'utilisateur
  if (nb_pas_dt() >= 0 && nb_pas_dt() <= NW && facsec_ == 1) print_warning(NW);

  static constexpr bool IS_DEUX = (_O_ == Ordre_RK::DEUX_CLASSIQUE) , IS_TROIS = (_O_ == Ordre_RK::TROIS_CLASSIQUE), IS_QUATRE = (_O_ == Ordre_RK::QUATRE_CLASSIQUE);
  static constexpr int NB_PTS = IS_DEUX ? 2 : ( IS_TROIS ? 3 : 4);
  static constexpr int NB_BUTCHER = IS_DEUX ? 0 : ( IS_TROIS ? 1 : ( IS_QUATRE ? 2 : 3 ));

  DoubleTab& present = eqn.inconnue().valeurs(), &futur = eqn.inconnue().futur();
  if (ki_.size()!=NB_PTS)
    {
      ki_.resize(NB_PTS);
      // just for initializing the array structure ...
      for (int i = 0; i < NB_PTS; i++) ki_[i] = present;
    }

  DoubleTrav sauv(present);
  sauv = present; // sauv = y0
  const auto& butcher = BUTCHER_TAB.at(NB_BUTCHER);
  init_rk_flux_accumulators(eqn);

  // Step 1
  eqn.derivee_en_temps_inco(ki_[0]); // ki[0] = f(y0)
  ki_[0] *= dt_; // ki[0] = h * f(y0)
  accumulate_flux_classical(eqn, butcher.at(0));

  for (int step = 1; step < NB_PTS; step++ ) // ATTENTION : ne touche pas !
    {
      present = sauv; // back to y0
      for (int i = 0; i < step; i++) present.ajoute(get_butcher_coeff<_ORDRE_,NB_PTS-1>().at(step-1).at(i), ki_[i]); // Et ouiiiiiiii :-)

      eqn.derivee_en_temps_inco(ki_[step]);
      ki_[step] *= dt_;
      accumulate_flux_classical(eqn, butcher.at(step));
    }

  futur = sauv; // futur = y1 = y0
  for (int i = 0; i < NB_PTS; i++) futur.ajoute(BUTCHER_TAB.at(NB_BUTCHER).at(i), ki_[i]);

  update_critere_statio(futur, eqn);

  // Update boundary condition on futur:
  eqn.domaine_Cl_dis().imposer_cond_lim(eqn.inconnue(), temps_courant() + pas_de_temps());
  futur.echange_espace_virtuel();

  present = sauv; // back to y0
  finalize_rk_flux_accumulators(eqn);

  return 1;
}

template <Ordre_RK _ORDRE_>
void TRUSTSchema_RK<_ORDRE_>::init_rk_flux_accumulators(Equation_base& eq) const
{
  for (int op = 0; op < eq.nombre_d_operateurs(); op++)
    eq.operateur(op).l_op_base().rk_reset_flux_accumulators();
}

template <Ordre_RK _ORDRE_>
void TRUSTSchema_RK<_ORDRE_>::finalize_rk_flux_accumulators(Equation_base& eq) const
{
  for (int op = 0; op < eq.nombre_d_operateurs(); op++)
    eq.operateur(op).l_op_base().rk_finalize_flux_accumulators();
}

template <Ordre_RK _ORDRE_>
void TRUSTSchema_RK<_ORDRE_>::accumulate_flux_classical(Equation_base& eq, double weight) const
{
  for (int op = 0; op < eq.nombre_d_operateurs(); op++)
    eq.operateur(op).l_op_base().rk_accumulate_flux(weight);
}

template <Ordre_RK _ORDRE_>
void TRUSTSchema_RK<_ORDRE_>::accumulate_flux_low_storage(Equation_base& eq, double ai, double bi) const
{
  for (int op = 0; op < eq.nombre_d_operateurs(); op++)
    eq.operateur(op).l_op_base().rk_low_storage_update_flux(ai, bi);
}

#endif /* TRUSTSchema_RK_TPP_included */
