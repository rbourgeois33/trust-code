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

#include <CL_Contrainte_Imposee.h>
#include <Domaine_Cl_dis_base.h>
#include <Equation_base.h>
#include <Motcle.h>
#include <Process.h>

Implemente_instanciable(CL_Contrainte_Imposee, "paroi_pression_imposee", Neumann);

// XD paroi_pression_imposee condlim_base CL_Contrainte_Imposee -1 CL_Contrainte_Imposee/paroi_pression_imposee
// XD attr ch front_field_base ch 0 Boundary field type.

Sortie& CL_Contrainte_Imposee::printOn(Sortie& os) const { return Neumann::printOn(os); }
Entree& CL_Contrainte_Imposee::readOn(Entree& is)
{
  if (app_domains.size() == 0) app_domains = { Motcle("Mecanique") };
  return Neumann::readOn(is);
}

void CL_Contrainte_Imposee::verifie_ch_init_nb_comp() const
{
  if (le_champ_front.non_nul())
    {
      const int nb_comp = le_champ_front->nb_comp();
      if (nb_comp != 1)
        {
          Cerr << que_suis_je() << " expects a pressure field with 1 component but received " << nb_comp << finl;
          Process::exit();
        }
    }
}
