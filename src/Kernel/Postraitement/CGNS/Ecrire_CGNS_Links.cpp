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
#include <communications.h>
#include <unordered_set>
#include <Ecrire_CGNS.h>
#include <Domaine.h>
#include <unistd.h>

#ifdef HAS_CGNS

/*
 * ***************** *
 * METHODS POUR LINK *
 * ***************** *
 */

void Ecrire_CGNS::cgns_open_grid_base_link_file()
{
  assert(Option_CGNS::USE_LINKS && !postraiter_domaine_);
  std::string fn;

  if (Process::is_parallel() && Option_CGNS::FILE_PER_COMM_GROUP && PE_Groups::has_user_defined_group())
    {
      const auto& grp = PE_Groups::get_user_defined_group();
      if (PE_Groups::enter_group(grp))
        {
          proc_maitre_local_comm_ = PE_Groups::groupe_TRUST().rank();
          envoyer_broadcast(proc_maitre_local_comm_, 0); // XXX should do this !
          fn = (Nom(baseFile_name_)).nom_me(proc_maitre_local_comm_).getString() + ".grid.cgns"; // file name

          unlink(fn.c_str());

          cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::PAR>(fn, fileId_, false);
          PE_Groups::exit_group();
        }
      Cerr << "**** Multiple parallel CGNS files " << baseFile_name_ << "_XXXX.grid.cgns opened !" << finl;
    }
  else
    {
      fn = baseFile_name_ + ".grid.cgns"; // file name

      unlink(fn.c_str());

      if (Process::is_parallel())
        cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::PAR>(fn, fileId_);
      else
        cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_);
    }
}

void Ecrire_CGNS::cgns_close_grid_or_solution_link_file(const double t, const TYPE_LINK_CGNS type, bool is_cerr)
{
  assert(Option_CGNS::USE_LINKS && !postraiter_domaine_);
  std::string fn; // file name

  if (type == TYPE_LINK_CGNS::GRID)
    {
      if (Process::is_parallel() && Option_CGNS::FILE_PER_COMM_GROUP && PE_Groups::has_user_defined_group())
        fn =  baseFile_name_ + "_XXXX.grid.cgns";
      else
        fn = baseFile_name_ + ".grid.cgns";
    }
  else if (type == TYPE_LINK_CGNS::SOLUTION)
    {
      if (Process::is_parallel() && Option_CGNS::FILE_PER_COMM_GROUP && PE_Groups::has_user_defined_group())
        fn = baseFile_name_ + "_XXXX_" + ".solution." + cgns_helper_.convert_double_to_string(t) + ".cgns";
      else
        fn = baseFile_name_ + ".solution." + cgns_helper_.convert_double_to_string(t) + ".cgns";
    }
  else if (type == TYPE_LINK_CGNS::FINAL_LINK)
    fn = baseFile_name_ + ".cgns";
  else
    Process::exit("Error in Ecrire_CGNS::cgns_close_grid_or_solution_link_file !!! \n");

  if (Process::is_parallel())
    {
      if ( Option_CGNS::FILE_PER_COMM_GROUP && PE_Groups::has_user_defined_group())
        {
          cgns_helper_.cgns_close_file<TYPE_RUN_CGNS::PAR>(fn /* inutile */, fileId_, false);
          Cerr << "**** Multiple parallel CGNS files " << fn << " closed !" << finl;
        }
      else
        cgns_helper_.cgns_close_file<TYPE_RUN_CGNS::PAR>(fn, fileId_, is_cerr);
    }
  else
    cgns_helper_.cgns_close_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_, is_cerr);
}

void Ecrire_CGNS::cgns_fill_info_grid_link_file(const char* basename, const CGNS_TYPE& cgns_type_elem, const int icelldim, const int nb_som, const int nb_elem, const bool is_polyedre)
{
  cellDim_.push_back(icelldim);
  baseZone_name_.push_back(std::string(basename));
  sizeId_.push_back( { nb_som, nb_elem } );

  if (cgns_type_elem == CGNS_ENUMV(NGON_n)) // cas polyedre
    {
      if (is_polyedre)
        connectname_.push_back( { "NGON_n", "NFACE_n" });
      else
        connectname_.push_back({ "NGON_n" });
    }
  else
    connectname_.push_back({ "Elem" });
}

// Attention : ind 0 => ELEM et SOM, ind 1 => FACES (si besoin pour faces) !
void Ecrire_CGNS::gather_local_sizeId_multi_loc(std::vector<std::vector<cgsize_t>>& sizeId_som_local_comm_tmp, std::vector<std::vector<cgsize_t>>& sizeId_elem_local_comm_tmp) const
{
  const bool has_elem_field = fld_loc_map_.count("ELEM");
  const bool has_som_field = fld_loc_map_.count("SOM");
  const bool has_faces_field = fld_loc_map_.count("FACES");

  sizeId_som_local_comm_tmp.push_back(std::vector<cgsize_t>());
  sizeId_elem_local_comm_tmp.push_back(std::vector<cgsize_t>());

  if (has_elem_field || has_som_field)
    {
      sizeId_som_local_comm_tmp.back().assign(Process::nproc(), -123 /* default */);
      sizeId_elem_local_comm_tmp.back().assign(Process::nproc(), -123 /* default */);

      int ind_base = -1;
      if (has_elem_field)
        ind_base = TRUST_2_CGNS::get_index_nom_vector(doms_written_, fld_loc_map_.at("ELEM"));
      else
        ind_base = TRUST_2_CGNS::get_index_nom_vector(doms_written_, fld_loc_map_.at("SOM"));

      MPI_Allgather(&sizeId_[ind_base][0], 1, MPI_LONG, sizeId_som_local_comm_tmp[0].data(), 1, MPI_LONG, Comm_Group_MPI::get_trio_u_world());
      MPI_Allgather(&sizeId_[ind_base][1], 1, MPI_LONG, sizeId_elem_local_comm_tmp[0].data(), 1, MPI_LONG, Comm_Group_MPI::get_trio_u_world());
    }

  if (has_faces_field)
    {
      sizeId_som_local_comm_tmp.push_back(std::vector<cgsize_t>(Process::nproc(), -123 /* default */));
      sizeId_elem_local_comm_tmp.push_back(std::vector<cgsize_t>(Process::nproc(), -123 /* default */));

      const int ind_base = TRUST_2_CGNS::get_index_nom_vector(doms_written_, fld_loc_map_.at("FACES"));

      MPI_Allgather(&sizeId_[ind_base][0], 1, MPI_LONG, sizeId_som_local_comm_tmp[1].data(), 1, MPI_LONG, Comm_Group_MPI::get_trio_u_world());
      MPI_Allgather(&sizeId_[ind_base][1], 1, MPI_LONG, sizeId_elem_local_comm_tmp[1].data(), 1, MPI_LONG, Comm_Group_MPI::get_trio_u_world());
    }
}

void Ecrire_CGNS::cgns_write_final_link_file_comm_group()
{
  if (vec_proc_maitre_local_comm_.empty())
    {
      vec_proc_maitre_local_comm_.assign(Process::nproc(), -123 /* default */);
      MPI_Allgather(&proc_maitre_local_comm_, 1, MPI_ENTIER, vec_proc_maitre_local_comm_.data(), 1, MPI_ENTIER, Comm_Group_MPI::get_trio_u_world());

      std::unordered_set<int> seen;

      for (int val : vec_proc_maitre_local_comm_)
        if (seen.insert(val).second)
          unique_vec_proc_maitre_local_comm_.push_back(val); // si val pas dedans

      const bool has_elem_field = fld_loc_map_.count("ELEM");
      const bool has_som_field = fld_loc_map_.count("SOM");
      const bool has_faces_field = fld_loc_map_.count("FACES");

      // Attention : ind 0 => ELEM et SOM, ind 1 => FACES (si besoin pour faces) !
      std::vector<std::vector<cgsize_t>> sizeId_som_local_comm_tmp, sizeId_elem_local_comm_tmp;
      gather_local_sizeId_multi_loc(sizeId_som_local_comm_tmp, sizeId_elem_local_comm_tmp);

      const int nb_grps = static_cast<int>(unique_vec_proc_maitre_local_comm_.size());

      // pour elem/som => ind 0
      sizeId_som_local_comm_.push_back(std::vector<cgsize_t>());
      sizeId_elem_local_comm_.push_back(std::vector<cgsize_t>());

      if (has_elem_field || has_som_field)
        {
          sizeId_som_local_comm_.back().assign(nb_grps, -123 /* default */);
          sizeId_elem_local_comm_.back().assign(nb_grps, -123 /* default */);

          for (int i = 0; i < nb_grps; i++)
            {
              int proc_grp = unique_vec_proc_maitre_local_comm_[i];
              sizeId_som_local_comm_[0][i] = sizeId_som_local_comm_tmp[0][proc_grp];
              sizeId_elem_local_comm_[0][i] = sizeId_elem_local_comm_tmp[0][proc_grp];
            }
        }

      // pour faces et si ca existe => ind 1
      if (has_faces_field)
        {
          sizeId_som_local_comm_.push_back(std::vector<cgsize_t>(Process::nproc(), -123 /* default */));
          sizeId_elem_local_comm_.push_back(std::vector<cgsize_t>(Process::nproc(), -123 /* default */));

          for (int i = 0; i < nb_grps; i++)
            {
              int proc_grp = unique_vec_proc_maitre_local_comm_[i];
              sizeId_som_local_comm_[1][i] = sizeId_som_local_comm_tmp[1][proc_grp];
              sizeId_elem_local_comm_[1][i] = sizeId_elem_local_comm_tmp[1][proc_grp];
            }
        }
    }

  if (!Process::me())
    {
      std::string fn = baseFile_name_ + ".cgns";

      unlink(fn.c_str());
      cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_, true);


      const int nb_grps = static_cast<int>(unique_vec_proc_maitre_local_comm_.size());

      for (auto& itr : fld_loc_map_)
        {
          const std::string& LOC = itr.first;
          /*
           * Quel ind a linker avec ??
           */
          int ind_base = 0; // par defaut c'est la base 0
          if (fld_loc_map_.count(LOC))
            {
              const auto& nom_dom = fld_loc_map_.at(LOC);
              ind_base = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);
            }


          if (cg_base_write(fileId_, baseZone_name_[ind_base].c_str(), cellDim_[ind_base], Objet_U::dimension, &baseId_[ind_base]) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_comm_group : cg_base_write !" << finl, TRUST_CGNS_ERROR();

          if (cg_biter_write(fileId_, baseId_[ind_base], "TimeIterValues", static_cast<int>(time_post_.size())) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_comm_group : cg_biter_write !" << finl, TRUST_CGNS_ERROR();

          if (cg_goto(fileId_, baseId_[ind_base], "BaseIterativeData_t", 1, "end") != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_comm_group : cg_goto BaseIterativeData_t !" << finl, TRUST_CGNS_ERROR();

          cgsize_t nuse = static_cast<cgsize_t>(time_post_.size());
          if (cg_array_write("TimeValues", CGNS_DOUBLE_TYPE, 1, &nuse, time_post_.data()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_comm_group : cg_array_write TimeValues !" << finl, TRUST_CGNS_ERROR();

          if (cg_simulation_type_write(fileId_, baseId_[ind_base], CGNS_ENUMV(TimeAccurate)) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_comm_group : cg_simulation_type_write !" << finl, TRUST_CGNS_ERROR();

          for (int gid = 0; gid < nb_grps; gid++)
            {
              int proc_grp = unique_vec_proc_maitre_local_comm_[gid];
              std::string file_group_id = Nom(baseFile_name_).nom_me(proc_grp).getString();
              file_group_id = TRUST_2_CGNS::remove_slash_linkfile(file_group_id);

              std::string zone_name = Nom("Zone").nom_me(proc_grp).getString();
              std::string linkfile = file_group_id + ".grid.cgns";

              cgsize_t isize[3][1];
              const int ind_som_elem_local_comm = (LOC == "FACES") ? 1 : 0;
              isize[0][0] = sizeId_som_local_comm_[ind_som_elem_local_comm][gid];
              isize[1][0] = sizeId_elem_local_comm_[ind_som_elem_local_comm][gid];
              isize[2][0] = 0;

              int zoneId_tmp = -1;
              if (cg_zone_write(fileId_, baseId_[ind_base], zone_name.c_str(), isize[0], CGNS_ENUMV(Unstructured), &zoneId_tmp) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_comm_group : cg_zone_write !" << finl, TRUST_CGNS_ERROR();

              std::string linkpath = "/" + baseZone_name_[ind_base] + "/" + baseZone_name_[ind_base] + "/GridCoordinates/";

              if (cg_goto(fileId_, baseId_[ind_base], "Zone_t", gid + 1, "end") != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_comm_group : cg_goto Zone_t !" << finl, TRUST_CGNS_ERROR();

              if (cg_link_write("GridCoordinates", linkfile.c_str(), linkpath.c_str()) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_comm_group : cg_link_write GridCoordinates !" << finl, TRUST_CGNS_ERROR();

              for (auto& con : connectname_[ind_base])
                {
                  linkpath = "/" + baseZone_name_[ind_base] + "/" + baseZone_name_[ind_base] + "/" + con + "/";
                  if (cg_link_write(con.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
                    Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_comm_group : cg_link_write connectivity " << con << finl, TRUST_CGNS_ERROR();
                }

              for (auto& itr_t : time_post_)
                {
                  std::string solname = "FlowSolution" + cgns_helper_.convert_double_to_string(itr_t) + "_" + LOC;
                  linkfile = file_group_id + ".solution." + cgns_helper_.convert_double_to_string(itr_t) + ".cgns";
                  linkpath = "/" + baseZone_name_[ind_base] + "/" + baseZone_name_[ind_base] + "/" + solname + "/";

                  if (cg_link_write(solname.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
                    Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_comm_group : cg_link_write FlowSolution " << solname << finl, TRUST_CGNS_ERROR();
                }

              cgsize_t idata[2] = {CGNS_STR_SIZE, static_cast<cgsize_t>(time_post_.size())};
              if (cg_ziter_write(fileId_, baseId_[ind_base], zoneId_tmp, "ZoneIterativeData") != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_comm_group : cg_ziter_write !" << finl, TRUST_CGNS_ERROR();

              if (cg_goto(fileId_, baseId_[ind_base], "Zone_t", gid + 1, "ZoneIterativeData_t",  1, "end") != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_comm_group : cg_goto ZoneIterativeData_t !" << finl, TRUST_CGNS_ERROR();

              const char* solname = (LOC == "SOM") ? solname_som_.c_str() : ( (LOC == "FACES") ? solname_faces_.c_str() : solname_elem_.c_str());
              if (cg_array_write("FlowSolutionPointers", CGNS_ENUMV(Character), 2, idata, solname) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_comm_group : cg_array_write FlowSolutionPointers !" << finl, TRUST_CGNS_ERROR();
            }

        }
      cgns_helper_.cgns_close_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_, true);
    }
}

void Ecrire_CGNS::cgns_open_solution_link_file(const double t, bool is_link)
{
  assert(Option_CGNS::USE_LINKS && !postraiter_domaine_);

  const bool enter_group_comm = Process::is_parallel() && Option_CGNS::FILE_PER_COMM_GROUP && PE_Groups::has_user_defined_group();

  std::string fn;

  if (is_link)
    fn = baseFile_name_ + ".cgns"; // file name
  else
    {
      if (enter_group_comm)
        fn = (Nom(baseFile_name_)).nom_me(proc_maitre_local_comm_).getString() + ".solution." + cgns_helper_.convert_double_to_string(t) + ".cgns"; // file name
      else
        fn = baseFile_name_ + ".solution." + cgns_helper_.convert_double_to_string(t) + ".cgns"; // file name
    }

  unlink(fn.c_str());

  if (Process::is_parallel())
    {
      cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::PAR>(fn, fileId_, enter_group_comm ? false : true);

      if (enter_group_comm)
        Cerr << "**** Multiple parallel CGNS files " << baseFile_name_ << "_XXXX.solution." + cgns_helper_.convert_double_to_string(t) + ".cgns opened !" << finl;
    }
  else
    cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_, true);


  const bool multi_loc = static_cast<int>(fld_loc_map_.size() > 1);

  for (auto &itr : fld_loc_map_)
    {
      const std::string& LOC = itr.first;
      const int ind_base = TRUST_2_CGNS::get_index_nom_vector(doms_written_, fld_loc_map_.at(LOC));

      const std::string BZname = multi_loc ? baseZone_name_[ind_base] + "_" + LOC : baseZone_name_[ind_base];

      if (cg_base_write(fileId_, BZname.c_str(), cellDim_[ind_base], Objet_U::dimension, &baseId_[ind_base]) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_open_solution_file : cg_base_write !" << finl, TRUST_CGNS_ERROR();

      cgsize_t isize[3][1];
      isize[0][0] = sizeId_[ind_base][0];
      isize[1][0] = sizeId_[ind_base][1];
      isize[2][0] = 0;

      if (cg_zone_write(fileId_, baseId_[ind_base], baseZone_name_[ind_base].c_str() /* Dom name */, isize[0], CGNS_ENUMV(Unstructured), &zoneId_[ind_base]) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_domaine_seq : cgns_open_solution_file !" << finl, TRUST_CGNS_ERROR();

      std::string linkfile = baseFile_name_ + ".grid.cgns"; // file name

      if (enter_group_comm)
        linkfile = (Nom(baseFile_name_)).nom_me(proc_maitre_local_comm_).getString() + ".grid.cgns"; // file name

      linkfile = TRUST_2_CGNS::remove_slash_linkfile(linkfile);

      std::string linkpath = "/" + baseZone_name_[ind_base] + "/" + baseZone_name_[ind_base] + "/GridCoordinates/";

      if (cg_goto(fileId_, baseId_[ind_base], "Zone_t", 1, "end") != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_open_solution_file : cg_goto !" << finl, TRUST_CGNS_ERROR();

      if (cg_link_write("GridCoordinates", linkfile.c_str(), linkpath.c_str()) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_open_solution_file : cg_link_write !" << finl, TRUST_CGNS_ERROR();

      for (auto &itr_conn : connectname_[ind_base])
        {
          linkpath = "/" + baseZone_name_[ind_base] + "/" + baseZone_name_[ind_base] + "/" + itr_conn + "/";

          if (cg_link_write(itr_conn.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_open_solution_file : cg_link_write !" << finl, TRUST_CGNS_ERROR();
        }
    }
}



void Ecrire_CGNS::cgns_write_final_link_file()
{
  if (Process::is_parallel() && Option_CGNS::FILE_PER_COMM_GROUP && PE_Groups::has_user_defined_group())
    {
      cgns_write_final_link_file_comm_group();
      return;
    }

  cgns_init_MPI(true); // set self mpi

  if (!Process::me())
    {
      // Fichier link maintenant
      cgns_open_solution_link_file( -123., true /* dernier fichier => link */);

      for (auto& itr : fld_loc_map_)
        {
          const std::string& LOC = itr.first;

          /*
           * Quel ind a linker avec ??
           */
          int ind_base = 0; // par defaut c'est la base 0
          if (fld_loc_map_.count(LOC))
            {
              const auto& nom_dom = fld_loc_map_.at(LOC);
              ind_base = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);
            }

          // link solutions
          for (auto& itr_t : time_post_)
            {
              std::string solname = "FlowSolution" + cgns_helper_.convert_double_to_string(itr_t) + "_" + LOC;

              std::string linkfile = baseFile_name_ + ".solution." + cgns_helper_.convert_double_to_string(itr_t) + ".cgns"; // file name
              linkfile = TRUST_2_CGNS::remove_slash_linkfile(linkfile);

              std::string linkpath = "/" + baseZone_name_[ind_base] + "/" + baseZone_name_[ind_base] + "/" + solname + "/";

              if (cg_link_write(solname.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file : cg_link_write !" << finl, TRUST_CGNS_ERROR();
            }

          cgns_helper_.cgns_write_iters<TYPE_ECRITURE_CGNS::SEQ>(true /* has_field */, 1, fileId_, baseId_[ind_base], ind_base /* 1st Zone */,
                                                                 zoneId_, LOC, solname_som_, solname_elem_, solname_faces_, time_post_);

        }
      cgns_close_grid_or_solution_link_file(-123. /* inutile*/, TYPE_LINK_CGNS::FINAL_LINK, true); // on ferme
    }

  cgns_init_MPI(); // back to COMM_WORLD
}

void Ecrire_CGNS::cgns_write_link_file_for_multiple_files()
{
  if (postraiter_domaine_) return; /* Do nothing */

  // Attention : ind 0 => ELEM et SOM, ind 1 => FACES (si besoin pour faces) !
  std::vector<std::vector<cgsize_t>> sizeId_som_local, sizeId_elem_local;
  gather_local_sizeId_multi_loc(sizeId_som_local, sizeId_elem_local);

  if (!Process::me()) // Only master proc writes !
    {
      std::string fn = baseFile_name_ + ".cgns"; // file name
      Cerr << "Option_CGNS::MULTIPLE_FILES is used ... so we write a unique link file " << fn << " ..." << finl;

      unlink(fn.c_str());

      cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_, true);

      for (auto& itr : fld_loc_map_)
        {
          const std::string& LOC = itr.first;
          const int ind_base = TRUST_2_CGNS::get_index_nom_vector(doms_written_, fld_loc_map_.at(LOC));

          if (cg_base_write(fileId_, baseZone_name_[ind_base].c_str(), cellDim_[ind_base], Objet_U::dimension, &baseId_[ind_base]) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_base_write !" << finl, TRUST_CGNS_ERROR();

          if (cg_biter_write(fileId_, baseId_[ind_base], "TimeIterValues", static_cast<int>(time_post_.size())) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_biter_write !" << finl, TRUST_CGNS_ERROR();

          if (cg_goto(fileId_, baseId_[ind_base], "BaseIterativeData_t", 1, "end") != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_goto BaseIterativeData_t !" << finl, TRUST_CGNS_ERROR();

          cgsize_t nuse = static_cast<cgsize_t>(time_post_.size());
          if (cg_array_write("TimeValues", CGNS_DOUBLE_TYPE, 1, &nuse, time_post_.data()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_array_write TimeValues !" << finl, TRUST_CGNS_ERROR();

          if (cg_simulation_type_write(fileId_, baseId_[ind_base], CGNS_ENUMV(TimeAccurate)) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_simulation_type_write !" << finl, TRUST_CGNS_ERROR();

          for (int proc = 0; proc < Process::nproc(); proc++)
            {
              std::string file_group_id = Nom(baseFile_name_).nom_me(proc).getString();
              std::string linkfile = file_group_id + ".cgns";
              std::string zone_name = Nom("Zone").nom_me(proc).getString();

              const int ind_som_elem_local_comm = (LOC == "FACES") ? 1 : 0;
              cgsize_t isize[3][1];
              isize[0][0] = sizeId_som_local[ind_som_elem_local_comm][proc];
              isize[1][0] = sizeId_elem_local[ind_som_elem_local_comm][proc];
              isize[2][0] = 0;

              int zoneId_tmp = -1;
              if (cg_zone_write(fileId_, baseId_[ind_base], zone_name.c_str(), isize[0], CGNS_ENUMV(Unstructured), &zoneId_tmp) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_zone_write !" << finl, TRUST_CGNS_ERROR();

              std::string linkpath = "/" + baseZone_name_[ind_base] + "/" + baseZone_name_[ind_base] + "/GridCoordinates/";

              if (cg_goto(fileId_, baseId_[ind_base], "Zone_t", proc + 1, "end") != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_goto Zone_t !" << finl, TRUST_CGNS_ERROR();

              if (cg_link_write("GridCoordinates", linkfile.c_str(), linkpath.c_str()) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_link_write GridCoordinates !" << finl, TRUST_CGNS_ERROR();

              for (auto& con : connectname_[ind_base])
                {
                  linkpath = "/" + baseZone_name_[ind_base] + "/" + baseZone_name_[ind_base] + "/" + con + "/";
                  if (cg_link_write(con.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
                    Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_link_write connectivity " << con << finl, TRUST_CGNS_ERROR();
                }

              for (auto& itr_t : time_post_)
                {
                  std::string solname = "FlowSolution" + cgns_helper_.convert_double_to_string(itr_t) + "_" + LOC;
//                  linkfile = file_group_id + "_" + LOC + ".solution." + cgns_helper_.convert_double_to_string(itr_t) + ".cgns";
                  linkpath = "/" + baseZone_name_[ind_base] + "/" + baseZone_name_[ind_base] + "/" + solname + "/";

                  if (cg_link_write(solname.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
                    Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_link_write FlowSolution " << solname << finl, TRUST_CGNS_ERROR();
                }

              cgsize_t idata[2] = {CGNS_STR_SIZE, static_cast<cgsize_t>(time_post_.size())};
              if (cg_ziter_write(fileId_, baseId_[ind_base], zoneId_tmp, "ZoneIterativeData") != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_ziter_write !" << finl, TRUST_CGNS_ERROR();

              if (cg_goto(fileId_, baseId_[ind_base], "Zone_t", proc + 1, "ZoneIterativeData_t",  1, "end") != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_goto ZoneIterativeData_t !" << finl, TRUST_CGNS_ERROR();

              const char* solname = (LOC == "SOM") ? solname_som_.c_str() : ( (LOC == "FACES") ? solname_faces_.c_str() : solname_elem_.c_str());
              if (cg_array_write("FlowSolutionPointers", CGNS_ENUMV(Character), 2, idata, solname) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_array_write FlowSolutionPointers !" << finl, TRUST_CGNS_ERROR();

            }

        }
      cgns_helper_.cgns_close_file<TYPE_RUN_CGNS::SEQ>(baseFile_name_ + ".cgns", fileId_, true);
    }
}

#endif /* HAS_CGNS */
