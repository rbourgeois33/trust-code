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

#include <Comm_Group_MPI.h>
#include <My_Comm_Group.h>
#include <Comm_Group.h>
#include <PE_Groups.h>
#include <Motcle.h>
#include <Param.h>
//#include <Probleme_base.h>

Implemente_instanciable(My_Comm_Group, "My_Comm_Group", Interprete);
// XD My_Comm_Group interprete My_Comm_Group 1 This keyword allows to create a user MPI Comm Group of size N using the processors allocated to TRUST. The set of processors is split in N subsets.

static OWN_PTR(Comm_Group) my_comm_group;

Entree& My_Comm_Group::readOn(Entree& is) { return Interprete::readOn(is); }

Sortie& My_Comm_Group::printOn(Sortie& os) const { return Interprete::printOn(os); }

//void test()
//{
//  Cerr << "TESTING CLASS My_Comm_Group" << finl;
//  Cerr << endl;
//  Process::barrier();
//
//  // get pb
//  if (Interprete::objet_existant("pb"))
//    {
//      const Probleme_base& pb=ref_cast(Probleme_base, Interprete::objet("pb"));
//      const int nb_elem = pb.domaine().nb_elem();
//      const int nb_procs = Process::nproc();
//
//      if (PE_Groups::enter_group(PE_Groups::get_user_defined_group()))
//        {
//          Cerr << "Process::mp_sum(nb_elem) on group = " << Process::mp_sum(nb_elem) << finl;
//          Process::barrier();
//          std::cerr << "Proc local : " << Process::me() << " , proc global : " <<  PE_Groups::groupe_TRUST().rank() << " , nb_elem : " << nb_elem << std::endl;
//
//          PE_Groups::exit_group();
//        }
//
//      Cerr << "@@@@@@@@@@ GLOBAL : Process::mp_sum(nb_elem) = " << Process::mp_sum(nb_elem) << finl;
//      Cerr << endl;
//      Process::barrier();
//
//      std::vector<int> global_nb_elem;
//      global_nb_elem.assign(nb_procs, -123 /* default */);
//      const Comm_Group_MPI& comm = ref_cast(Comm_Group_MPI, PE_Groups::groupe_TRUST());
//
//      MPI_Allgather(&nb_elem, 1, MPI_ENTIER, global_nb_elem.data(), 1, MPI_ENTIER, comm.get_mpi_comm());
//
//      Cerr << "@@@@@@@@@@ GLOBAL : gather all nb_elem in global_nb_elem : "<< finl;
//      for (auto& itr : global_nb_elem)
//        Cerr << "   - " << itr << endl;
//
//      Cerr << endl;
//      Process::barrier();
//
//      std::vector<int> global_nb_elem2;
//
//      Cerr << "Testing local gather all nb_elem on group " << finl;
//      const auto& grp = PE_Groups::get_user_defined_group();
//      if (PE_Groups::enter_group(grp))
//        {
//          global_nb_elem2.assign(Process::nproc(), -123 /* default */);
//          const Comm_Group_MPI& comm_loc = ref_cast(Comm_Group_MPI, grp);
//
//          MPI_Allgather(&nb_elem, 1, MPI_ENTIER, global_nb_elem2.data(), 1, MPI_ENTIER, comm_loc.get_mpi_comm());
//
//          for (auto& itr : global_nb_elem2)
//            Cerr << "   - " << itr << endl;
//
//          PE_Groups::exit_group();
//        }
//      Cerr << endl;
//
//      Process::barrier();
//      Process::exit();
//    }
//}

Entree& My_Comm_Group::interpreter(Entree& is)
{
  int nb_groups = -123;
  Param param(que_suis_je());
  param.ajouter("Group_nb", &nb_groups, Param::REQUIRED); // XD_ADD_P entier Number of groups to define in your Comm Group.
  param.lire_avec_accolades_depuis(is);

  if (Process::is_sequential())
    return is; /* rien */

#ifndef MPI_
  Process::exit("What !!! You need an MPI TRUST version to use My_Comm_Group !!!");
#endif

  assert(my_comm_group.est_nul());

  const int nb_procs = Process::nproc();
  const int rank = Process::me();

  if (nb_groups > nb_procs)
    nb_groups = nb_procs; // sinon

  const int base_size = nb_procs / nb_groups;
  const int extra = nb_procs % nb_groups; // nombre des procs avec 1 de plus !

  //  ex : 0,1,2,3   4,5,6,7   8,9,10
  //  Total = 11 procs

  int count = 0;
  for (int i = 0; i < nb_groups; i++)
    {
      const int group_size = base_size + (i < extra ? 1 : 0); // le 1 de plus envoyer aux premiers !!!

      if (rank >= count && rank < count + group_size)
        {
          ArrOfInt tab(group_size);
          for (int j = 0; j < group_size; j++)
            tab[j] = count + j;

          PE_Groups::create_group(tab, my_comm_group);
          PE_Groups::initialize_user_defined_group(my_comm_group.valeur());
          break; // le groupe a été trouvé, pas besoin de continuer
        }

      count += group_size;
    }

//  test();

  return is;
}

const Comm_Group& My_Comm_Group::get_my_comm_group()
{
  assert(my_comm_group.non_nul());
  return my_comm_group.valeur();
}
