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

#include <Navier_Stokes_std.h>
#include <Assembleur_base.h>
#include <Probleme_base.h>
#include <Simpler_Base.h>
#include <Domaine_VF.h>
#include <TRUSTTrav.h>

Implemente_base_sans_constructeur(Simpler_Base,"Simpler_Base",Solveur_non_lineaire);

Sortie& Simpler_Base::printOn(Sortie& os ) const
{
  return os;
}

Entree& Simpler_Base::readOn(Entree& is )
{
  is_seuil_convg_variable = 0;
  controle_residu_ = 0;
  no_qdm_ = 0;
  return Solveur_Implicite_base::readOn(is);
}

Entree& Simpler_Base::lire(const Motcle& motlu, Entree& is)
{
  if (motlu == "no_qdm")
    {
      Cerr << motlu << finl;
      no_qdm_ = 1;
    }
  else if (motlu == "facsec_diffusion_for_sets")
    {
      Cerr << motlu << finl;
      is >> facsec_diffusion_for_sets_;
    }
  else
    {
      Cerr << "Keyword : " << motlu << " is not undertood in " << que_suis_je() << finl;
      Process::exit();
    }
  return is;
}

void Simpler_Base::assembler_matrice_pression_implicite(Equation_base& eqn_NS,const Matrice_Morse& matrice,Matrice& matrice_en_pression_2)
{
  Navier_Stokes_std& eqnNS = ref_cast(Navier_Stokes_std,eqn_NS);
  const IntVect& tab1 = matrice.get_tab1();
  const DoubleVect& coeff = matrice.get_coeff();
  const DoubleTab& present = eqn_NS.inconnue().valeurs();
  int nb_comp = 1;
  int deux_entrees = 0;

  if (present.nb_dim()==2)
    {
      deux_entrees = 1;
      nb_comp = present.dimension(1);
    }

  const Domaine_VF& le_dom = ref_cast(Domaine_VF,eqnNS.domaine_dis());
  if (deux_entrees==0)
    {
      DoubleVect vol2 = le_dom.volumes_entrelaces();
      int ns = vol2.size();
      for (int i=0; i<ns; i++)
        {
          int idiag = tab1[i*nb_comp]-1;
          double ref = coeff[idiag];
          for (int c=1; c<nb_comp; c++)
            if (!est_egal(ref,coeff[tab1[i*nb_comp+c]-1]))
              {
                Cerr<<"Pb dans Piso sur la diagonale case"<<i<<" comp "<< c<<" ref "<<ref<<" valeurs "<<coeff[tab1[i*nb_comp+c]-1]<<finl;
                exit();
              }
          vol2[i] = coeff[idiag];
        }
      vol2.echange_espace_virtuel();
      eqnNS.assembleur_pression()->assembler_mat(matrice_en_pression_2,vol2,1,1);
    }
  else
    {
      DoubleTab vol2(present);
      int ns = vol2.dimension_tot(0);
      for (int i=0; i<ns; i++)
        for (int c=0; c<nb_comp; c++)
          vol2(i,c) = matrice(i*nb_comp+c,i*nb_comp+c);

      vol2.echange_espace_virtuel();
      eqnNS.assembleur_pression()->assembler_mat(matrice_en_pression_2,vol2,1,1);
    }

  SolveurSys& solveur_pression_ = eqnNS.solveur_pression();
  solveur_pression_->reinit();
}
