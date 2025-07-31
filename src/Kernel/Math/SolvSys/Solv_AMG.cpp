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

#include <Solv_AMG.h>
#include <EChaine.h>
#include <Motcle.h>
#include <Solv_AMGX.h>
#include <Solv_Petsc_GPU.h>
#ifdef TRUST_USE_ROCM
#include <rocm-core/rocm_version.h>
#endif
#include <comm_incl.h> // Mandatory to have MPIX_CUDA_AWARE_SUPPORT defined or not
#include <MD_Vector_composite.h>

Implemente_instanciable(Solv_AMG,"Solv_AMG",SolveurSys_base);
// XD amg solveur_sys_base amg 0 Wrapper for AMG preconditioner-based solver which switch for the best one on CPU/GPU Nvidia/GPU AMD
// XD attr solveur chaine solveur 0 not_set
// XD attr option_solveur bloc_lecture option_solveur 0 not_set

// printOn
Sortie& Solv_AMG::printOn(Sortie& s ) const
{
  s << chaine_lue_;
  return s;
}

/**
 * @brief Reads the configuration for the AMG solver from the input stream.
 *
 * This function parses the input stream to configure the AMG solver parameters,
 * including the relative tolerance (RTOL), absolute tolerance (ATOL), and whether
 * to print additional information (IMPR). It supports different solver libraries
 * based on the available hardware (CPU or GPU).
 *
 * The expected input format is:
 *   amg solver { rtol value [impr] }
 *
 * @param is The input stream from which to read the solver configuration.
 * @return The input stream after reading the configuration.
 *
 * @throws Process::exit if the input syntax is incorrect or if an unsupported
 * library is specified.
 */
Entree& Solv_AMG::readOn(Entree& is)
{
  // amg GCP|BISGTSTAB|GMRES { atol|rtol doublee [st double] [impr]  }
  is >> solver_;
  if ((Motcle)solver_!="GCP")
    {
      Cerr << solver_ << " not supported yet for AMG !" << finl;
      Process::exit();
    }
  Motcle motcle;
  is >> motcle;
  while (motcle != "}")
    {
      if (motcle=="{") {}
      else if (motcle=="RTOL") is >> rtol_;
      else if (motcle=="ATOL") is >> atol_;
      else if (motcle=="ST") is >> st_;
      else if (motcle=="IMPR") impr_ = true;
      else if (motcle=="READ_MATRIX") set_read_matrix(true);
      else if (motcle=="SAVE_MATRIX_PETSC_FORMAT") set_save_matrix(2);
      else if (motcle=="SEUIL") Process::exit("Use atol 'absolute tolerance' instead of seuil.");
      else
        {
          options_+=" ";
          options_+=motcle;
          Cerr << motcle << " is not an option of AMG solver." << finl;
          Process::exit();
        }
      is >> motcle;
    }
  if (atol_<0 && rtol_<0) Process::exit("atol or rtol should be defined in AMG solver.");
  return is;
}

void Solv_AMG::create_block_amg(int n, Nom precond)
{
  // ToDo: not efficient on P0P1Pa (n==3)
  chaine_lue_="cli { -ksp_type ";
  chaine_lue_+=petsc_cg_issue_ ? "bcgs" : "cg";
  chaine_lue_+=rtol_>0 ? Nom(rtol_, " -ksp_rtol %e") : "";
  chaine_lue_+=atol_>0 ? Nom(atol_, " -ksp_atol %e") : "";
  chaine_lue_+=" -ksp_norm_type UNPRECONDITIONED \
-pc_type fieldsplit \
-pc_fieldsplit_type additive";
  if (precond=="gamg")
    {
      // ToDo: fix crash on multi-GPU (issue sent to PETSc support)
      chaine_lue_+=" -fieldsplit_P0_ksp_type preonly \
-fieldsplit_P0_pc_type gamg \
-fieldsplit_P0_pc_gamg_threshold 0.01 \
-fieldsplit_P0_pc_gamg_square_graph 1 \
-fieldsplit_P1_ksp_type preonly \
-fieldsplit_P1_pc_type gamg \
-fieldsplit_P1_pc_gamg_threshold 0.01 \
-fieldsplit_P1_pc_gamg_square_graph 1";
      if (n==3)
        {
          chaine_lue_+=" -fieldsplit_Pa_ksp_type preonly \
-fieldsplit_Pa_ksp_type preonly \
-fieldsplit_Pa_pc_type gamg \
-fieldsplit_Pa_pc_gamg_threshold 0.01 \
-fieldsplit_Pa_pc_gamg_square_graph 1";
        }
    }
  else if (precond=="boomeramg")
    {
      chaine_lue_+=" -fieldsplit_P0_ksp_type preonly \
-fieldsplit_P0_pc_type hypre \
-fieldsplit_P0_pc_hypre_type boomeramg \
-fieldsplit_P0_pc_hypre_boomeramg_strong_threshold 0.1 \
-fieldsplit_P0_pc_hypre_boomeramg_print_statistics 1 \
-fieldsplit_P1_ksp_type preonly \
-fieldsplit_P1_pc_type hypre \
-fieldsplit_P1_pc_hypre_type boomeramg \
-fieldsplit_P1_pc_hypre_boomeramg_strong_threshold 0.1 \
-fieldsplit_P1_pc_hypre_boomeramg_print_statistics 1";
      if (n==3)
        {
          chaine_lue_+=" -fieldsplit_Pa_ksp_type preonly \
-fieldsplit_Pa_pc_type hypre \
-fieldsplit_Pa_pc_hypre_type boomeramg \
-fieldsplit_Pa_pc_hypre_boomeramg_strong_threshold 0.1 \
-fieldsplit_Pa_pc_hypre_boomeramg_print_statistics 1";
        }
      // To avoid this issue on Nvidia:  CUSPARSE ERROR (code = 11, insufficient resources) at csr_spgemm_device_cusparse.c:152
#ifdef TRUST_USE_CUDA
      if (n==2) chaine_lue_+=" -fieldsplit_P0_pc_mg_galerkin_mat_product_algorithm hypre";
      if (n==2) chaine_lue_+=" -fieldsplit_P1_pc_mg_galerkin_mat_product_algorithm hypre";
      if (n==3) chaine_lue_+=" -fieldsplit_Pa_pc_mg_galerkin_mat_product_algorithm hypre";
#endif
    }
  else if (precond=="amgx")
    {
      Cerr << "Warning! PETSc with AmgX preconditioner was not tested yet for nnz>2^31 !" << finl;
      chaine_lue_+=" -fieldsplit_P0_ksp_type preonly \
-fieldsplit_P0_pc_type amgx \
-fieldsplit_P0_pc_amgx_verbose 1 \
-fieldsplit_P0_pc_amgx_print_grid_stats 1 \
-fieldsplit_P1_ksp_type preonly \
-fieldsplit_P1_pc_type amgx \
-fieldsplit_P1_pc_amgx_verbose 1 \
-fieldsplit_P1_pc_amgx_print_grid_stats 1";
      if (n==3)
        {
          chaine_lue_+=" -fieldsplit_Pa_ksp_type preonly \
-fieldsplit_Pa_pc_type amgx \
-fieldsplit_Pa_pc_amgx_verbose 1 \
-fieldsplit_Pa_pc_amgx_print_grid_stats 1";
        }
    }
  else
    Process::exit("Error in Solv_AMG::create_block_amg");
  chaine_lue_ +=" }";
}

Nom boomeramg(double st)
{
  Nom chaine(" { precond boomeramg { }");
  if (st>=0)
    {
      chaine += " cli { -pc_hypre_boomeramg_strong_threshold";
      chaine += Nom(st, "%e");
      chaine += " }";
    }
  return chaine;
}
void Solv_AMG::create_amg()
{
  // We select the more efficient/robust one:
  chaine_lue_ = solver_;
#if defined(TRUST_USE_CUDA)
  library_ = "petsc_gpu";
  chaine_lue_ += boomeramg(st_); // Best GPU solver
#if defined(MPIX_CUDA_AWARE_SUPPORT)
  // KSP divergence with cg+boomeramg/amgx on multi-node with MPI Cuda Aware so we switch to bcgs !
  // Strangely KSPSolve is 2x-3x slower on A100X vs MI250X... rocsparse better than cusparse ? And Kokkos-Kernels ?
  // With Cuda backend, BCGS+AmgX seems better than BCGS+Boomeramg (A100)
  if (Process::nproc()>4)
    {
      petsc_cg_issue_ = true;
      library_ = "amgx";
      chaine_lue_ = solver_;
      chaine_lue_ += " { precond c-amg {";
      if (st_>=0) chaine_lue_ += Nom(st_, " p:strength_threshold %e");
      chaine_lue_ += " }";
    }
#endif
#elif defined(TRUST_USE_ROCM)
  library_ = "petsc_gpu";
  const char* value = std::getenv("ROCM_ARCH");
  if (value != nullptr && std::string(value) == "gfx1100")
    {
      if (st_>=0) Process::exit("st option not supported yet in Solv_AMG");
      if (Process::is_parallel())
        chaine_lue_ += " { precond ua-amg { }";  // Converge mais plus lent que sa-amg
      else
        chaine_lue_ += " { precond sa-amg { }";  // Crash en parallele
    }
  else
    chaine_lue_ += boomeramg(st_); // Best GPU solver (// sa-amg is slow...)
#else
  library_ = "petsc";
  chaine_lue_ += boomeramg(st_); // Best CPU solver
#endif
  chaine_lue_ += rtol_>0 ? Nom(rtol_, " rtol %e") : Nom(atol_, " atol %e");
  if (impr_) chaine_lue_ += " impr";
  if (options_!="") chaine_lue_ += options_;
  chaine_lue_ += " }";
}

int Solv_AMG::resoudre_systeme(const Matrice_Base& mat, const DoubleVect& b, DoubleVect& x)
{
  // We don't create solver during readOn as usual but just before solve to get more infos about matrix/vectors to fine tune
  if (solveur_.est_nul())
    {
      create_amg();
      int nb_blocks = sub_type(MD_Vector_composite, b.get_md_vector().valeur()) ? ref_cast(MD_Vector_composite, b.get_md_vector().valeur()).nb_parts() : 1;
      if (nb_blocks>1)
        {
          // Block matrix : we use PCFieldsplit (eg: VEF) for preconditioner
          // Much better convergence for P0P1 for instance
          Cerr << "Detecting " << nb_blocks << "x" << nb_blocks << " blocks into the matrix. Creating a specific block preconditioning:" << finl;
          if (chaine_lue_.contient("gamg"))
            create_block_amg(nb_blocks, "gamg");
          else if (chaine_lue_.contient("boomeramg"))
            create_block_amg(nb_blocks, "boomeramg");
          else if (library_=="amgx")
            {
              library_ = "petsc_gpu";
              create_block_amg(nb_blocks, "amgx");
            }
        }
      Cerr << "====================================================================" << finl;
      Cerr << "Creating AMG solver: " << library_ << " " << chaine_lue_ << finl;
      Cerr << "====================================================================" << finl;
      EChaine entree(chaine_lue_);
      Nom nom_solveur("Solv_");
      nom_solveur+=library_;
      solveur_.typer(nom_solveur);
      solveur_.nommer("solveur_pression");
      if (library_=="amgx")
        ref_cast(Solv_AMGX, solveur_.valeur()).create_solver(entree);
      else if (library_=="petsc")
        ref_cast(Solv_Petsc, solveur_.valeur()).create_solver(entree);
      else if (library_=="petsc_gpu")
        ref_cast(Solv_Petsc_GPU, solveur_.valeur()).create_solver(entree);
      else
        Process::exit("Unsupported case in Solv_AMG::readOn");
      solveur_->set_save_matrix(save_matrix());
      solveur_->set_read_matrix(read_matrix());
    }
  statistics().end_count(STD_COUNTERS::system_solver,-1,0);
  int res = solveur_.resoudre_systeme(mat, b, x);
  statistics().begin_count(STD_COUNTERS::system_solver,statistics().get_last_opened_counter_level()+1);
  return res;
}
