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
#include <Frontiere_deplacement_impose.h>
#include <Equation_base.h>
#include <Motcle.h>

Implemente_instanciable(Frontiere_deplacement_impose, "Frontiere_deplacement_impose", Scalaire_impose_paroi);
// XD Frontiere_deplacement_impose dirichlet Frontiere_deplacement_impose 0 Neutron flux prescribed at inlet condition
// XD attr ch front_field_base ch 0 Boundary field type.

// XXX TODO FIXME Yo Yannick faut mettre ca dans TRUST dans TRAD_2.org mais j'ai la flemme :/
// faut retirer attr criteres_convergence .... bla bla bla

// XD sets piso sets -1 Stability-Enhancing Two-Step solver which is useful for a multiphase problem.

Sortie& Frontiere_deplacement_impose::printOn(Sortie& s) const { return Scalaire_impose_paroi::printOn(s); }
Entree& Frontiere_deplacement_impose::readOn(Entree& s)
{
  if (app_domains.size() == 0) app_domains = { Motcle("Mecanique") };
  return Scalaire_impose_paroi::readOn(s);
}
