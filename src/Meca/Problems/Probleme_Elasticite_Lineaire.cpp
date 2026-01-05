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

#include <Probleme_Elasticite_Lineaire.h>
#include <Milieu_Elasticite.h>
#include <Process.h>

Implemente_instanciable(Probleme_Elasticite_Lineaire,"Probleme_Elasticite_Lineaire",Probleme_base);


Sortie& Probleme_Elasticite_Lineaire::printOn(Sortie& os) const { return os; }
Entree& Probleme_Elasticite_Lineaire::readOn(Entree& is) { return Probleme_base::readOn(is); }

const Equation_base& Probleme_Elasticite_Lineaire::equation(int i) const
{
  if (i != 0) Process::exit();
  return equation_mecanique_;
}

Equation_base& Probleme_Elasticite_Lineaire::equation(int i)
{
  if (i != 0) Process::exit();
  return equation_mecanique_;
}

void Probleme_Elasticite_Lineaire::abort_if_unimplemented(const char* method) const
{
  Cerr << que_suis_je() << "::" << method << " is not implemented yet." << finl;
  Process::exit();
}

