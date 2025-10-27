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

#ifndef Probleme_Elasticite_Lineaire_included
#define Probleme_Elasticite_Lineaire_included

#include <Probleme_base.h>
#include <Equation_Navier_Cauchy.h>

class Milieu_Elastic;

/*! @brief Probleme minimal d'elasticite lineaire.
 *
 *  L'implementation numerique est a completer.
 */
class Probleme_Elasticite_Lineaire : public Probleme_base
{
  Declare_instanciable(Probleme_Elasticite_Lineaire);

public:

  int nombre_d_equations() const override { return 1; }
  const Equation_base& equation(int) const override;
  Equation_base& equation(int) override;

  void associer_milieu_base(const Milieu_base& mil) override { equation_mecanique_.associer_milieu_base(mil); }
  const Milieu_base& milieu() const override { return equation_mecanique_.milieu(); }
  Milieu_base& milieu() override { return equation_mecanique_.milieu(); }

private:
  void abort_if_unimplemented(const char* method) const;

  Equation_Navier_Cauchy equation_mecanique_;
};

#endif
