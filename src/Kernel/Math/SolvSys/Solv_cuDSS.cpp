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

#include <Solv_cuDSS.h>
#include <TRUSTVect.h>
#include <EChaine.h>
#include <Motcle.h>
#include <SFichier.h>
#include <Comm_Group_MPI.h>
#include <MD_Vector_std.h>
#include <MD_Vector_composite.h>
#include <Device.h>
#include <Perf_counters.h>
#include <Array_tools.h>
#include <TRUSTTrav.h>

#define CUDSS_CALL_AND_CHECK(call, status, msg) \
    do { \
        status = call; \
        if (status != CUDSS_STATUS_SUCCESS) { \
            printf("Example FAILED: CUDSS call ended unsuccessfully with status = %d, details: " #msg "\n", status); \
            Process::exit(); \
        } \
    } while(0);

Implemente_instanciable_sans_constructeur_ni_destructeur(Solv_cuDSS, "Solv_cuDSS", Solv_Externe);
// XD cuDSS petsc cuDSS 0 Solver via cuDSS API
// XD attr solveur chaine solveur 0 not_set
// XD attr option_solveur bloc_lecture option_solveur 0 not_set

// printOn
Sortie& Solv_cuDSS::printOn(Sortie& s ) const
{
  s << chaine_lue_;
  return s;
}

// Read parameters
Entree& Solv_cuDSS::readOn(Entree& entree)
{
#ifdef cuDSS_
  Cerr <<"[cuDSS] reading jdd  "<< finl;
  lecture(entree);
  EChaine is(get_chaine_lue());
  Motcle accolade_ouverte("{"), accolade_fermee("}");
  Motcle solver, motlu;
  is >> solver;   // On lit le solveur en premier puis les options du solveur
  is >> motlu; // On lit l'accolade
  if (motlu != accolade_ouverte)
    {
      Cerr << "[cuDSS] Error while reading the parameters of the solver " << solver << " { ... }" << finl;
      Cerr << "We expected " << accolade_ouverte << " instead of " << motlu << finl;
      exit();
    }

  if (solver==(Motcle)"LU")
    {
      std::cout<<"[cuDSS] LU is read, matrix is assumed non symmetric"<<std::endl;
      mtype=CUDSS_MTYPE_GENERAL;
      mview=CUDSS_MVIEW_FULL;
    }
  else if (solver==(Motcle)"Cholesky")
    {
      std::cout<<"[cuDSS] Cholesky is read, matrix is assumed symmetric"<<std::endl;
      mtype=CUDSS_MTYPE_SYMMETRIC;
      mview=CUDSS_MVIEW_UPPER;
    }
  // Lecture des parametres du solver (LU non symetric, cholesky symmetric)
  // LU|Cholesky { algo name [impr] }
  is >> motlu;
  while (motlu!=accolade_fermee)
    {
      if (motlu==(Motcle)"impr")
        {
          fixer_limpr(1);
        }
      else if (motlu==(Motcle)"algo")
        {
          is >> motlu; // Lecture de l'algo

          if (motlu==(Motcle)"0") {  reorder_alg = CUDSS_ALG_DEFAULT;}
          else if (motlu==(Motcle)"1") {  reorder_alg = CUDSS_ALG_1;}
          else if (motlu==(Motcle)"2") {  reorder_alg = CUDSS_ALG_2;}
          else if (motlu==(Motcle)"3") {  reorder_alg = CUDSS_ALG_3;}
          else if (motlu==(Motcle)"4") {  reorder_alg = CUDSS_ALG_4;}
          else if (motlu==(Motcle)"5") {  reorder_alg = CUDSS_ALG_5;}
          else
            {
              Cerr << "[cuDSS] algo number "<<motlu<<" unrecognized (0-5) only"<<finl;
              Process::exit();
            }
          Cerr << "[cuDSS] algo number "<<motlu<<" read. "<<finl;
          Cerr <<"[cuDSS] Note from cuDSS doc: Different values represent different algorithms (for reordering, factorization, etc.) and can lead to significant differences in accuracy and performance. It is currently recommended to use CUDSS_ALG_DEFAULT (0) and only in case accuracy or performance are not sufficient, one can experiment with other values."<<finl;

        }
      else
        {
          Cerr << motlu << "[cuDSS] keyword not recognized for cuDSS solver " << solver << finl;
          Process::exit();
        }
      is >> motlu;
    }
#else
  Process::exit("Sorry, cuDSS solvers not available with this build (cuDSS readon).");
#endif
  return entree;
}


Solv_cuDSS::Solv_cuDSS(const Solv_cuDSS& org)
{
  // on relance la lecture .... Necessary for some reasons
  EChaine recup(org.get_chaine_lue());
  readOn(recup);
}

Solv_cuDSS::~Solv_cuDSS()
{
#ifdef cuDSS_

  if (Axb_are_built)
    {
      cudssMatrixDestroy(A);
      cudssMatrixDestroy(b);
      cudssMatrixDestroy(x);
    }

  if (solver_is_built)
    {
      cudssDataDestroy(handle, solverData);
      cudssConfigDestroy(solverConfig);
      cudssDestroy(handle);
    }
#endif
}

#include <Debog.h>
int Solv_cuDSS::resoudre_systeme(const Matrice_Base& a, const DoubleVect& bvect, DoubleVect& xvect)
{
#ifdef cuDSS_
  if ((nouvelle_matrice()||first_solve))
    {
      /* conversion matric base to csr */
      /*build the csr matrix on host*/
      construit_matrice_morse_intermediaire(a, csr_);
      const Matrice_Morse& csr = csr_.nb_lignes() ? csr_ : ref_cast(Matrice_Morse, a);

      /* create cudss matrix, vector and solver and size them */
      /* Does very few things if the matrix has not changed*/
      Create_objects(csr);

      /* give the device data / csr pointers to the cudss matrix */
      set_pointers_A(csr);

      /* sizes should never change */
      assert(a.nb_lignes() == n);
      assert(bvect.size_totale() == n);
      assert(xvect.size_totale() == n);
    }

  /* analysis and facto can be only done when the matrix changes */
  statistics().begin_count(STD_COUNTERS::gpu_library);
  if ((nouvelle_matrice()||first_solve))
    {
      /* Symbolic factorization x, and b are unused*/
      CUDSS_CALL_AND_CHECK(cudssExecute(handle, CUDSS_PHASE_ANALYSIS, solverConfig, solverData,
                                        A, x, b), status, "cudssExecute for analysis");

      /* Factorization x, and b are unused*/
      CUDSS_CALL_AND_CHECK(cudssExecute(handle, CUDSS_PHASE_FACTORIZATION, solverConfig,
                                        solverData, A, x, b), status, "cudssExecute for facto");
    }


  //The pointers to x and b should never change between calls to resoudre
// assert(x_values_d==const_cast<double*>(xvect.view_rw<1>().data()));
// assert(b_values_d==const_cast<double*>(bvect.view_ro<1>().data()));

  /* give the device data / csr pointers to the cudss vectors x and b */
  /* This has to be done at each solve to ensure sync mechanism*/
  set_pointers_xb(bvect, xvect);

  /*some checks*/
  /* sizes should never change */
  assert(a.nb_lignes()==n);
  assert(bvect.size_totale()==n);
  assert(xvect.size_totale()==n);


  /* Solving */
  CUDSS_CALL_AND_CHECK(cudssExecute(handle, CUDSS_PHASE_SOLVE, solverConfig, solverData,
                                    A, x, b), status, "cudssExecute for solve");
  cudaStreamSynchronize(nullptr); // Important factorization and solve phases are asynchronous
  statistics().end_count(STD_COUNTERS::gpu_library);

  /*compute error in debug mode */
#ifndef NDEBUG
  DoubleVect test(bvect);
  test*=-1;
  a.ajouter_multvect(xvect,test);
  double vrai_residu = mp_norme_vect(test);
  assert(vrai_residu<1e-5);
  Cout << "||Ax-b||=" << vrai_residu << finl;
#endif
  first_solve = false;
  fixer_nouvelle_matrice(0);
  return 1;
#else
  Process::exit("Sorry, cuDSS solvers not available with this build (resoudre_systeme).");
  return -1;
#endif
}

void Solv_cuDSS::Create_objects(const Matrice_Morse& csr)
{
#ifdef cuDSS_
  /* Create the solver and matrix / vector objects */
  /* Solver and vectors are created only in the first call */
  /*Matrix can be destroyed and re-created if nnz changes*/
  if (Process::is_parallel())
    {
      Process::exit("Sorry, cuDSS is a sequential solver.");
    }

  /* get dimensions and nnz */
  /* check that n does not change between solves */
  int new_n = csr.get_tab1().size()-1;
  int new_nnz = csr.get_coeff().size();
#ifndef NDEBUG
  assert(((new_n==n)||(first_solve))); // n should never change between solves
#endif


  //Re-build A only if new nnz is different from old one, possible if matrix changes
  if((first_solve)||(new_nnz != nnz))
    {
      nnz=new_nnz;
      n = new_n;

      /* destroy if A is built */
      if (Axb_are_built)
        cudssMatrixDestroy(A);

      /* create matrix, we give nullptrs*/
      CUDSS_CALL_AND_CHECK(cudssMatrixCreateCsr(&A, n, n, nnz, nullptr, nullptr,
                                                nullptr, nullptr, CUDA_R_32I, CUDA_R_64F, mtype, mview,
                                                base), status, "cudssMatrixCreateCsr");
    }

  if (first_solve)
    {
      /* create rhs and vector we give nullptrs*/
      // n should never change between solves, no need to resize

      CUDSS_CALL_AND_CHECK(cudssMatrixCreateDn(&b, n, nrhs, n, nullptr, CUDA_R_64F,
                                               CUDSS_LAYOUT_COL_MAJOR), status, "cudssMatrixCreateDn for b");
      CUDSS_CALL_AND_CHECK(cudssMatrixCreateDn(&x, n, nrhs, n, nullptr, CUDA_R_64F,
                                               CUDSS_LAYOUT_COL_MAJOR), status, "cudssMatrixCreateDn for x");

      Axb_are_built=true;

      /* Create the solver */

      /* Create handle */
      CUDSS_CALL_AND_CHECK(cudssCreate(&handle), status, "cudssCreate");
      /* Create config */
      CUDSS_CALL_AND_CHECK(cudssConfigCreate(&solverConfig), status, "cudssConfigCreate");
      /* set options to config */
      CUDSS_CALL_AND_CHECK(cudssConfigSet(solverConfig, CUDSS_CONFIG_REORDERING_ALG,
                                          &reorder_alg, sizeof(cudssAlgType_t)), status, "cudssConfigSet");
      /* create solver data */
      CUDSS_CALL_AND_CHECK(cudssDataCreate(handle, &solverData), status, "cudssDataCreate");

      solver_is_built = true;
    }
#endif
}

void Solv_cuDSS::set_pointers_A(const Matrice_Morse& csr)
{
#ifdef cuDSS_
  /* get pointers */
  csr_offsets_d = const_cast<int*>(csr.get_tab1().view_ro<1>().data());
  csr_columns_d = const_cast<int*>(csr.get_tab2().view_ro<1>().data());
  csr_values_d = const_cast<double*>(csr.get_coeff().view_ro<1>().data());

  /* give pointers to A*/
  CUDSS_CALL_AND_CHECK(cudssMatrixSetCsrPointers(A,
                                                 csr_offsets_d, /*rowStart*/
                                                 nullptr,/*rowEnd*/
                                                 csr_columns_d,
                                                 csr_values_d), status, "cudssMatrixSetCsrPointers")
#endif
}

void Solv_cuDSS::set_pointers_xb(const DoubleVect& bvect, DoubleVect& xvect)
{
#ifdef cuDSS_

  /* get pointers */
  x_values_d = const_cast<double*>(xvect.view_wo<1>().data());
  b_values_d = const_cast<double*>(bvect.view_ro<1>().data());

  /* give pointers to vectors*/
  CUDSS_CALL_AND_CHECK(cudssMatrixSetValues(x, x_values_d), status, "cudssMatrixSetValues for x");
  CUDSS_CALL_AND_CHECK(cudssMatrixSetValues(b, b_values_d), status, "cudssMatrixSetValues for b");
#endif
}

