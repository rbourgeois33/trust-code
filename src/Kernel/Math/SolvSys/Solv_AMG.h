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

#ifndef Solv_AMG_included
#define Solv_AMG_included

#include <SolveurSys.h>
#include <Solv_Petsc.h>
#include <Perf_counters.h>
#include <EChaine.h>

/*! @brief AMD solver wrapper to switch to the more robust/performant AMG preconditioner on CPU/GPU Nvidia/GPU AMD
 *  and have datafile unicity, as changing solver/preconditioner for different architecture is painful
 */
class Solv_AMG : public SolveurSys_base
{
  Declare_instanciable(Solv_AMG);
public :
  virtual int solveur_direct() const override { return 0; };
  virtual int resoudre_systeme(const Matrice_Base& mat, const DoubleVect& b, DoubleVect& x) override;
  virtual int resoudre_systeme(const Matrice_Base& mat, const DoubleVect& b, DoubleVect& x, int niter_max) override
  {
    statistics().end_count(STD_COUNTERS::system_solver,-1,0);
    int res = solveur_.resoudre_systeme(mat, b, x);
    statistics().begin_count(STD_COUNTERS::system_solver,statistics().get_last_opened_counter_level()+1);
    return res;
  };
private :
  void create_amg();
  void create_block_amg(int,Nom);
  SolveurSys solveur_;
  Nom library_="", solver_="", options_="";
  double rtol_=-1, atol_=-1, st_=-1;
  bool impr_ = false;
  bool petsc_cg_issue_ = false;
};

#endif


