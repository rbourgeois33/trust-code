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

#include <Ecrire_CGNS.h>
#include <unordered_set>
#include <Domaine.h>
#include <unistd.h>

#ifdef HAS_CGNS

void Ecrire_CGNS::link_multi_loc_support_pb_deformable()
{
  // TODO FIXME : a factoriser avec 3 methodes ...
  // loop and write linked supports !
  for (auto &itr : fld_loc_map_)
    {
      const std::string& LOC = itr.first;
      assert (LOC != "FACES" && has_elem_som_loc_);

      const Nom nom_dom = fld_loc_map_.at(LOC);
      const int index_glob = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);

      const Nom nom_dom_mod = TRUST_2_CGNS::modify_domaine_name_for_link(nom_dom, LOC);
      const int ind_base = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom_mod);

      if (cg_base_write(fileId_, nom_dom.getChar(), cellDim_[ind_base], Objet_U::dimension, &baseId_[index_glob]) != CG_OK)
        Cerr << "Error Ecrire_CGNS::link_multi_loc_support_pb_deformable : cg_base_write !" << finl, TRUST_CGNS_ERROR();

      cgsize_t isize[3][1];
      isize[0][0] = sizeId_[ind_base][0];
      isize[1][0] = sizeId_[ind_base][1];
      isize[2][0] = 0;

      if (cg_zone_write(fileId_, baseId_[index_glob], nom_dom.getChar() /* Dom name */, isize[0], CGNS_ENUMV(Unstructured), &zoneId_[index_glob]) != CG_OK)
        Cerr << "Error Ecrire_CGNS::link_multi_loc_support_pb_deformable : cgns_open_solution_file !" << finl, TRUST_CGNS_ERROR();

      std::string linkfile = ""; // XXX this file

      std::string linkpath = "/" + baseZone_name_[ind_base] + "/" + baseZone_name_[ind_base] + "/GridCoordinates/";

      if (cg_goto(fileId_, baseId_[index_glob], "Zone_t", 1, "end") != CG_OK)
        Cerr << "Error Ecrire_CGNS::link_multi_loc_support_pb_deformable : cg_goto !" << finl, TRUST_CGNS_ERROR();

      if (cg_link_write("GridCoordinates", linkfile.c_str(), linkpath.c_str()) != CG_OK)
        Cerr << "Error Ecrire_CGNS::link_multi_loc_support_pb_deformable : cg_link_write !" << finl, TRUST_CGNS_ERROR();

      for (auto &itr_conn : connectname_[ind_base])
        {
          linkpath = "/" + baseZone_name_[ind_base] + "/" + baseZone_name_[ind_base] + "/" + itr_conn + "/";

          if (cg_link_write(itr_conn.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::link_multi_loc_support_pb_deformable : cg_link_write !" << finl, TRUST_CGNS_ERROR();
        }
    }

  multi_loc_deformable_support_linked_ = true; // of course !
}

void Ecrire_CGNS::cgns_write_final_link_file_comm_group_pb_deformable()
{
  if (!Process::me())
    {
      std::string fn = baseFile_name_ + ".cgns";

      unlink(fn.c_str());
      cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_, true);

      const int nb_grps = static_cast<int>(unique_vec_proc_maitre_local_comm_.size());
      for (auto &itr : fld_loc_map_)
        {
          const std::string& LOC = itr.first;
          const Nom nom_dom = itr.second;
          const int index_glob = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);

          int ind_base = index_glob;
          if (has_elem_som_loc_ && LOC != "FACES")
            {
              const Nom nom_dom_mod = TRUST_2_CGNS::modify_domaine_name_for_link(nom_dom, LOC);
              ind_base = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom_mod);
            }

          if (cg_base_write(fileId_, nom_dom.getChar(), cellDim_[ind_base], Objet_U::dimension, &baseId_[index_glob]) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_base_write !" << finl, TRUST_CGNS_ERROR();

          if (cg_biter_write(fileId_, baseId_[index_glob], "TimeIterValues", static_cast<int>(time_post_.size())) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_biter_write !" << finl, TRUST_CGNS_ERROR();

          if (cg_goto(fileId_, baseId_[index_glob], "BaseIterativeData_t", 1, "end") != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_goto BaseIterativeData_t !" << finl, TRUST_CGNS_ERROR();

          cgsize_t nuse = static_cast<cgsize_t>(time_post_.size());
          if (cg_array_write("TimeValues", CGNS_DOUBLE_TYPE, 1, &nuse, time_post_.data()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_array_write TimeValues !" << finl, TRUST_CGNS_ERROR();

          if (cg_simulation_type_write(fileId_, baseId_[index_glob], CGNS_ENUMV(TimeAccurate)) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_simulation_type_write !" << finl, TRUST_CGNS_ERROR();

          for (int gid = 0; gid < nb_grps; gid++)
            {
              int proc_grp = unique_vec_proc_maitre_local_comm_[gid];
              std::string zone_name = Nom("Zone").nom_me(proc_grp).getString();
              //          cgsize_t isize[3] = { sizeId_[ind_base][0], sizeId_[ind_base][1], 0 };
              cgsize_t isize[3][1];
              isize[0][0] = sizeId_som_local_comm_[0][gid];
              isize[1][0] = sizeId_elem_local_comm_[0][gid];
              isize[2][0] = 0;
              int zoneId_tmp = -1;
              if (cg_zone_write(fileId_, baseId_[index_glob], zone_name.c_str(), isize[0], CGNS_ENUMV(Unstructured), &zoneId_tmp) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_zone_write !" << finl, TRUST_CGNS_ERROR();

              std::string grid_name, grid_name_loc, linkfile, linkpath;
              bool conn_written = false;
              std::string file_group_id = Nom(baseFile_name_).nom_me(proc_grp).getString();

              for (auto &itr_t : time_post_)
                {
                  linkfile = file_group_id + ".solution." + cgns_helper_.convert_double_to_string(itr_t) + ".cgns";
                  TRUST_2_CGNS::remove_slash_linkfile(linkfile);

                  grid_name_loc = "GridCoordinates";

                  if (conn_written) // Pas la premiere fois
                    grid_name_loc += cgns_helper_.convert_double_to_string(itr_t) + "_" + LOC;

                  grid_name_loc.resize(CGNS_STR_SIZE, ' ');
                  grid_name += grid_name_loc;

                  linkpath = "/" + baseZone_name_[ind_base] + "/" + baseZone_name_[ind_base] + "/GridCoordinates/";

                  if (cg_goto(fileId_, baseId_[index_glob], "Zone_t", gid + 1, "end") != CG_OK)
                    Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_goto Zone_t !" << finl, TRUST_CGNS_ERROR();

                  if (cg_link_write(grid_name_loc.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
                    Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_link_write GridCoordinates !" << finl, TRUST_CGNS_ERROR();

                  // Lier toutes les sections / connectivites (une seule fois par zone)
                  if (!conn_written)
                    for (auto &itr_conn : connectname_[ind_base])
                      {
                        linkpath = "/" + baseZone_name_[ind_base] + "/" + baseZone_name_[ind_base] + "/" + itr_conn + "/";
                        if (cg_link_write(itr_conn.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
                          Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_link_write connectivity !" << finl, TRUST_CGNS_ERROR();

                        conn_written = true;
                      }

                  std::string solname = "FlowSolution" + cgns_helper_.convert_double_to_string(itr_t) + "_" + LOC;
                  linkpath = "/" + nom_dom.getString() + "/" + nom_dom.getString() + "/" + solname + "/";
                  if (cg_link_write(solname.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
                    Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_link_write FlowSolution !" << finl, TRUST_CGNS_ERROR();
                }

              cgsize_t idata[2] = { CGNS_STR_SIZE, static_cast<cgsize_t>(time_post_.size()) };

              if (cg_ziter_write(fileId_, baseId_[index_glob], zoneId_tmp, "ZoneIterativeData") != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_ziter_write !" << finl, TRUST_CGNS_ERROR();

              if (cg_goto(fileId_, baseId_[index_glob], "Zone_t", gid + 1, "ZoneIterativeData_t", 1, "end") != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_goto ZoneIterativeData_t !" << finl, TRUST_CGNS_ERROR();

              if (cg_array_write("GridCoordinatesPointers", CGNS_ENUMV(Character), 2, idata, grid_name.c_str()) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_array_write GridCoordinatesPointers !" << finl, TRUST_CGNS_ERROR();

              const char *solname = (LOC == "SOM") ? solname_som_.c_str() : ((LOC == "FACES") ? solname_faces_.c_str() : solname_elem_.c_str());
              if (cg_array_write("FlowSolutionPointers", CGNS_ENUMV(Character), 2, idata, solname) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_array_write FlowSolutionPointers !" << finl, TRUST_CGNS_ERROR();
            }
        }

      cgns_helper_.cgns_close_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_, true);
    }
}

void Ecrire_CGNS::cgns_write_final_link_file_pb_deformable()
{
  if (Process::is_parallel() && Option_CGNS::FILE_PER_COMM_GROUP && PE_Groups::has_user_defined_group())
    {
      if (vec_proc_maitre_local_comm_.empty())
        cgns_write_final_link_file_comm_group();

      cgns_write_final_link_file_comm_group_pb_deformable();

      return;
    }

  cgns_init_MPI(true); // set self mpi

  if (!Process::me()) // seul le proc 0 ecrit le fichier link
    {
      cgns_open_solution_link_file(-123., true);

      for (auto& itr : fld_loc_map_)
        {
          const std::string& LOC = itr.first;
          const Nom nom_dom = itr.second;
          const int index_glob = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);

          int ind_base = index_glob;
          if (has_elem_som_loc_ && LOC != "FACES")
            {
              const Nom nom_dom_mod = TRUST_2_CGNS::modify_domaine_name_for_link(nom_dom, LOC);
              ind_base = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom_mod);
            }

          if (cg_base_write(fileId_, nom_dom.getChar(), cellDim_[ind_base], Objet_U::dimension, &baseId_[index_glob]) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_base_write !" << finl, TRUST_CGNS_ERROR();

          if (cg_biter_write(fileId_, baseId_[index_glob], "TimeIterValues", static_cast<int>(time_post_.size())) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_biter_write !" << finl, TRUST_CGNS_ERROR();

          if (cg_goto(fileId_, baseId_[index_glob], "BaseIterativeData_t", 1, "end") != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_goto BaseIterativeData_t !" << finl, TRUST_CGNS_ERROR();

          cgsize_t nuse = static_cast<cgsize_t>(time_post_.size());
          if (cg_array_write("TimeValues", CGNS_DOUBLE_TYPE, 1, &nuse, time_post_.data()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_array_write TimeValues !" << finl, TRUST_CGNS_ERROR();

          if (cg_simulation_type_write(fileId_, baseId_[index_glob], CGNS_ENUMV(TimeAccurate)) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_simulation_type_write !" << finl, TRUST_CGNS_ERROR();

//          cgsize_t isize[3] = { sizeId_[ind_base][0], sizeId_[ind_base][1], 0 };
          cgsize_t isize[3][1];
          isize[0][0] = sizeId_[ind_base][0];
          isize[1][0] = sizeId_[ind_base][1];
          isize[2][0] = 0;

          if (cg_zone_write(fileId_, baseId_[index_glob], nom_dom.getChar(), isize[0], CGNS_ENUMV(Unstructured), &zoneId_[index_glob]) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_zone_write !" << finl, TRUST_CGNS_ERROR();

          std::string grid_name, grid_name_loc, linkfile, linkpath;
          bool conn_written = false;

          for (auto& itr_t : time_post_)
            {
              linkfile = baseFile_name_ + ".solution." + cgns_helper_.convert_double_to_string(itr_t) + ".cgns";
              TRUST_2_CGNS::remove_slash_linkfile(linkfile);

              grid_name_loc = "GridCoordinates";

              if (conn_written) // Pas la premiere fois
                grid_name_loc += cgns_helper_.convert_double_to_string(itr_t) + "_" + LOC;

              grid_name_loc.resize(CGNS_STR_SIZE, ' ');
              grid_name += grid_name_loc;

              linkpath = "/" + baseZone_name_[ind_base] + "/" + baseZone_name_[ind_base] + "/GridCoordinates/";

              if (cg_goto(fileId_, baseId_[index_glob], "Zone_t", zoneId_[index_glob], "end") != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_goto Zone_t !" << finl, TRUST_CGNS_ERROR();

              if (cg_link_write(grid_name_loc.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_link_write GridCoordinates !" << finl, TRUST_CGNS_ERROR();

              // Lier toutes les sections / connectivites (une seule fois par zone)
              if (!conn_written)
                for (auto& itr_conn : connectname_[ind_base])
                  {
                    linkpath = "/" + baseZone_name_[ind_base] + "/" + baseZone_name_[ind_base] + "/" + itr_conn + "/";
                    if (cg_link_write(itr_conn.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
                      Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_link_write connectivity !" << finl, TRUST_CGNS_ERROR();

                    conn_written = true;
                  }

              std::string solname = "FlowSolution" + cgns_helper_.convert_double_to_string(itr_t) + "_" + LOC;
              linkpath = "/" + nom_dom.getString() + "/" + nom_dom.getString() + "/" + solname + "/";
              if (cg_link_write(solname.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_link_write FlowSolution !" << finl, TRUST_CGNS_ERROR();
            }

          cgsize_t idata[2] = {CGNS_STR_SIZE, static_cast<cgsize_t>(time_post_.size())};

          if (cg_ziter_write(fileId_, baseId_[index_glob], zoneId_[index_glob], "ZoneIterativeData") != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_ziter_write !" << finl, TRUST_CGNS_ERROR();

          if (cg_goto(fileId_, baseId_[index_glob], "Zone_t", zoneId_[index_glob], "ZoneIterativeData_t", 1, "end") != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_goto ZoneIterativeData_t !" << finl, TRUST_CGNS_ERROR();

          if (cg_array_write("GridCoordinatesPointers", CGNS_ENUMV(Character), 2, idata, grid_name.c_str()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_array_write GridCoordinatesPointers !" << finl, TRUST_CGNS_ERROR();

          const char* solname = (LOC == "SOM") ? solname_som_.c_str() : ( (LOC == "FACES") ? solname_faces_.c_str() : solname_elem_.c_str());
          if (cg_array_write("FlowSolutionPointers", CGNS_ENUMV(Character), 2, idata, solname) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_array_write FlowSolutionPointers !" << finl, TRUST_CGNS_ERROR();
        }

      cgns_close_grid_or_solution_link_file(-123., TYPE_LINK_CGNS::FINAL_LINK, true);
    }

  cgns_init_MPI(); // back to COMM_WORLD
}

void Ecrire_CGNS::cgns_write_domaine_deformable_seq(const Domaine * domaine,const Nom& nom_dom, const DoubleTab& les_som, const IntTab& les_elem, const Motcle& type_elem)
{
  const int ind = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);
  TRUST_2_CGNS& TRUST2CGNS = T2CGNS_[ind];
  TRUST2CGNS.associer_domaine_TRUST(domaine, domaine_dis_.non_nul() ? &(domaine_dis_.valeur()) : nullptr, les_som, les_elem, postraiter_domaine_);

  CGNS_TYPE cgns_type_elem = TRUST2CGNS.convert_elem_type(type_elem);
  const bool is_polyedre = (type_elem == "POLYEDRE" || type_elem == "PRISME" || type_elem == "PRISME_HEXAG");
  const int icelldim = TRUST2CGNS.topo_dim_from_elem(cgns_type_elem, is_polyedre); // avant ca : icelldim = les_som.dimension(1)
  const int iphysdim = Objet_U::dimension, nb_som = les_som.dimension(0), nb_elem = les_elem.dimension(0);

  std::vector<double> xCoords, yCoords, zCoords;
  TRUST2CGNS.fill_coords(xCoords, yCoords, zCoords);

  int coordsId;

  char basename[CGNS_STR_SIZE];
  strcpy(basename, nom_dom.getChar()); // dom name

  if (cg_base_write(fileId_, basename, icelldim, iphysdim, &baseId_[ind]) != CG_OK)
    Cerr << "Error Ecrire_CGNS::cgns_write_domaine_seq : cg_base_write !" << finl, TRUST_CGNS_ERROR();

  cgsize_t isize[3][1];
  isize[0][0] = nb_som;
  isize[1][0] = nb_elem;
  isize[2][0] = 0; /* boundary vertex size (zero if elements not sorted) */

  if (nb_elem)
    {
      /* Create zone & grid coords */
      cgns_helper_.cgns_write_zone_grid_coord<TYPE_ECRITURE_CGNS::SEQ>(icelldim, fileId_, baseId_[ind], basename /* Dom name */, isize[0],
                                                                       zoneId_[ind], xCoords, yCoords, zCoords, coordsId, coordsId, coordsId);

      /* Set element connectivity */
      std::string linkfile = baseFile_name_ + ".solution." + cgns_helper_.convert_double_to_string(time_post_[0]) + ".cgns";
      TRUST_2_CGNS::remove_slash_linkfile(linkfile);

      if (cg_goto(fileId_, baseId_[ind], "Zone_t", zoneId_[ind], "end") != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_seq : cg_goto Zone_t !" << finl, TRUST_CGNS_ERROR();

      for (auto &itr_conn : connectname_[ind])
        {
          std::string linkpath = "/" + baseZone_name_[ind] + "/" + baseZone_name_[ind] + "/" + itr_conn + "/";

          if (cg_link_write(itr_conn.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_seq : cg_link_write connectivity !" << finl, TRUST_CGNS_ERROR();
        }
    }
}

void Ecrire_CGNS::cgns_write_domaine_deformable_par_in_zone(const Domaine * domaine,const Nom& nom_dom, const DoubleTab& les_som, const IntTab& les_elem, const Motcle& type_elem)
{
#ifdef MPI_
  const int ind = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);
  TRUST_2_CGNS& TRUST2CGNS = T2CGNS_[ind];
  TRUST2CGNS.associer_domaine_TRUST(domaine, domaine_dis_.non_nul() ? &(domaine_dis_.valeur()) : nullptr, les_som, les_elem, postraiter_domaine_);

  CGNS_TYPE cgns_type_elem = TRUST2CGNS.convert_elem_type(type_elem);
  const bool is_polyedre = (type_elem == "POLYEDRE" || type_elem == "PRISME" || type_elem == "PRISME_HEXAG");
  const int icelldim = TRUST2CGNS.topo_dim_from_elem(cgns_type_elem, is_polyedre); // avant ca : icelldim = les_som.dimension(1)
  const int nb_elem = les_elem.dimension(0), iphysdim = Objet_U::dimension;

  std::vector<double> xCoords, yCoords, zCoords;
  TRUST2CGNS.fill_coords(xCoords, yCoords, zCoords);

  char basename[CGNS_STR_SIZE];
  strcpy(basename, nom_dom.getChar()); // dom name

  if (cg_base_write(fileId_, basename, icelldim, iphysdim, &baseId_[ind]) != CG_OK)
    Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_par_in_zone : cg_base_write !" << finl, TRUST_CGNS_ERROR();

  const int ns_tot = TRUST2CGNS.get_ns_tot(), ne_tot = TRUST2CGNS.get_ne_tot();
  const bool enter_group_comm = Option_CGNS::FILE_PER_COMM_GROUP && PE_Groups::has_user_defined_group() && !postraiter_domaine_;
  const int proc_me = enter_group_comm ? TRUST2CGNS.get_proc_me_local_comm() : Process::me();

  cgsize_t isize[3][1];
  isize[0][0] = (ns_tot == 0 && enter_group_comm) ? 1 : ns_tot; // si ns_tot = 0, on va juste creer une zone vide
  isize[1][0] = (ne_tot == 0 && enter_group_comm) ? 1 : ne_tot; // si ne_tot = 0, on va juste creer une zone vide
  isize[2][0] = 0; /* boundary vertex size (zero if elements not sorted) */

  int coordsIdx = -123, coordsIdy = -123, coordsIdz = -123;

  cgns_helper_.cgns_write_zone_grid_coord<TYPE_ECRITURE_CGNS::PAR_IN>(icelldim, fileId_, baseId_[ind], basename /* Dom name */, isize[0],
                                                                      zoneId_[ind], xCoords, yCoords, zCoords, coordsIdx, coordsIdy, coordsIdz);

  if (ne_tot == 0 && ns_tot == 0) return; // XXX Elie Saikali : zone vide creer, rien a faire de plus ... (cas FILE_PER_COMM_GROUP !!!)

  if (nb_elem > 0) // seulement si le proc a qlq chose a ecrire
    {
      const std::vector<int>& incr_max_som = TRUST2CGNS.get_global_incr_max_som(),
                              &incr_min_som = TRUST2CGNS.get_global_incr_min_som();

      cgsize_t min = incr_min_som[proc_me], max = incr_max_som[proc_me];
      assert (min < max);

      cgns_helper_.cgns_write_grid_coord_data<TYPE_ECRITURE_CGNS::PAR_IN>(icelldim, fileId_, baseId_[ind], zoneId_[ind],
                                                                          coordsIdx, coordsIdy, coordsIdz, min, max, xCoords, yCoords, zCoords);
    }

  /* Set element connectivity */
  std::string linkfile = baseFile_name_ + ".solution." + cgns_helper_.convert_double_to_string(time_post_[0]) + ".cgns";
  TRUST_2_CGNS::remove_slash_linkfile(linkfile);

  if (cg_goto(fileId_, baseId_[ind], "Zone_t", zoneId_[ind], "end") != CG_OK)
    Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_par_in_zone : cg_goto Zone_t !" << finl, TRUST_CGNS_ERROR();

  for (auto &itr_conn : connectname_[ind])
    {
      std::string linkpath = "/" + baseZone_name_[ind] + "/" + baseZone_name_[ind] + "/" + itr_conn + "/";

      if (cg_link_write(itr_conn.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_domaine_deformable_par_in_zone : cg_link_write connectivity !" << finl, TRUST_CGNS_ERROR();
    }

#endif /*MPI_*/
}

#endif /* HAS_CGNS */
