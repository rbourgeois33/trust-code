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

#ifndef Source_Meca_Grad_Pression_Thermique_included
#define Source_Meca_Grad_Pression_Thermique_included

#include <Source_base.h>

/*! @brief Source couplant gradient de pression et dilatation thermique pour la formulation thermomécanique.
 *
 *  La classe ne contient pour l'instant que le squelette des méthodes nécessaires. Les contributions
 *  physiques seront ajoutées dans une étape ultérieure.
 */
class Source_Meca_Grad_Pression_Thermique : public Source_base
{
  Declare_instanciable(Source_Meca_Grad_Pression_Thermique);

public:
  int initialiser(double temps) override;
  void mettre_a_jour(double temps) override;
  DoubleTab& ajouter(DoubleTab&) const override;
  void associer_domaines(const Domaine_dis_base&, const Domaine_Cl_dis_base&) override {}
  void associer_pb(const Probleme_base&) override {}

private:
  OWN_PTR(Champ_Don_base) T_ref_;
  OWN_PTR(Champ_Don_base) T_;
};

#endif /* Source_Meca_Grad_Pression_Thermique_included */

