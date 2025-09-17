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

#ifdef HAS_CGNS

void Ecrire_CGNS::cgns_write_final_link_file_pb_deformable()
{
  if (Process::is_parallel() && Option_CGNS::FILE_PER_COMM_GROUP && PE_Groups::has_user_defined_group())
    {
      throw; // FIXME
      cgns_write_final_link_file_comm_group();
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
              linkfile = TRUST_2_CGNS::remove_slash_linkfile(linkfile);

              grid_name_loc = "GridCoordinates";

              if (conn_written) // Pas la premiere fois
                grid_name_loc += cgns_helper_.convert_double_to_string(itr_t) + "_" + LOC;

              grid_name_loc.resize(CGNS_STR_SIZE, ' ');
              grid_name += grid_name_loc;

              linkpath = "/" + nom_dom.getString() + "/" + nom_dom.getString() + "/GridCoordinates/";

              if (cg_goto(fileId_, baseId_[index_glob], "Zone_t", zoneId_[index_glob], "end") != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_goto Zone_t !" << finl, TRUST_CGNS_ERROR();

              if (cg_link_write(grid_name_loc.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file_pb_deformable : cg_link_write GridCoordinates !" << finl, TRUST_CGNS_ERROR();

              // Lier toutes les sections / connectivites (une seule fois par zone)
              if (!conn_written)
                for (auto& itr_conn : connectname_[ind_base])
                  {
                    linkpath = "/" + nom_dom.getString() + "/" + nom_dom.getString()+ "/" + itr_conn + "/";
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

  /* Fill coords */
  std::vector<double> xCoords, yCoords, zCoords;
  TRUST2CGNS.fill_coords(xCoords, yCoords, zCoords);

  const int icelldim = les_som.dimension(1), iphysdim = Objet_U::dimension,
            nb_som = les_som.dimension(0), nb_elem = les_elem.dimension(0);

  int coordsId;

  /* Base write */
  char basename[CGNS_STR_SIZE];
  strcpy(basename, nom_dom.getChar()); // dom name

  if (cg_base_write(fileId_, basename, icelldim, iphysdim, &baseId_.back()) != CG_OK)
    Cerr << "Error Ecrire_CGNS::cgns_write_domaine_seq : cg_base_write !" << finl, TRUST_CGNS_ERROR();

  /* Vertex, cell & boundary vertex sizes */
  cgsize_t isize[3][1];
  isize[0][0] = nb_som;
  isize[1][0] = nb_elem;
  isize[2][0] = 0; /* boundary vertex size (zero if elements not sorted) */

  const bool is_polyedre = (type_elem == "POLYEDRE" || type_elem == "PRISME" || type_elem == "PRISME_HEXAG");

  /* Write all */
  if (nb_elem) // XXX cas // mais MULTIPLE_FILES
    {
      /* Create zone & grid coords */
      cgns_helper_.cgns_write_zone_grid_coord<TYPE_ECRITURE_CGNS::SEQ>(icelldim, fileId_, baseId_, basename /* Dom name */, isize[0],
                                                                       zoneId_, xCoords, yCoords, zCoords, coordsId, coordsId, coordsId);

      /* Set element connectivity */
      int sectionId;
      cgsize_t start = 1, end;

      if (cgns_type_elem == CGNS_ENUMV(NGON_n)) // cas polyedre
        {
          throw;
          assert (domaine != nullptr);
          std::vector<cgsize_t> sf, sf_offset;

          end = start + static_cast<cgsize_t>(TRUST2CGNS.convert_connectivity_ngon(sf, sf_offset, is_polyedre)) -1;

          if (cg_poly_section_write(fileId_, baseId_.back(), zoneId_.back(), "NGON_n", CGNS_ENUMV(NGON_n), start, end, 0, sf.data(), sf_offset.data(), &sectionId))
            TRUST_CGNS_ERROR();

          if (is_polyedre) // Pas pour polygone
            {
              std::vector<cgsize_t> ef, ef_offset;

              start = end + 1;
              end = start + static_cast<cgsize_t>(TRUST2CGNS.convert_connectivity_nface(ef, ef_offset)) -1;

              if (cg_poly_section_write(fileId_, baseId_.back(), zoneId_.back(), "NFACE_n", CGNS_ENUMV(NFACE_n), start, end, 0, ef.data(), ef_offset.data(), &sectionId))
                TRUST_CGNS_ERROR();
            }
        }
      else
        {
          int nodes_per_elem = -1;

          switch(cgns_type_elem)
            {
            case CGNS_ENUMV(HEXA_8):
              nodes_per_elem = 8;
              break;
            case CGNS_ENUMV(QUAD_4):
              nodes_per_elem = 4;
              break;
            case CGNS_ENUMV(TETRA_4):
              nodes_per_elem = 4;
              break;
            case CGNS_ENUMV(TRI_3):
              nodes_per_elem = 3;
              break;
            case CGNS_ENUMV(BAR_2):
              nodes_per_elem = 2;
              break;
            default:
              Cerr << "Type not yet coded in TRUST_2_CGNS::convert_connectivity ! Call the 911 !" << finl;
              Process::exit();
            }

          const std::vector<cgsize_t>& elems = TRUST2CGNS.get_connectivity_elem();

          end = start + static_cast<cgsize_t>(elems.size()) / nodes_per_elem - 1;

          if (cg_section_write(fileId_, baseId_.back(), zoneId_.back(), "Elem", cgns_type_elem, start, end, 0, elems.data(), &sectionId) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_domaine_seq : cg_section_write !" << finl, TRUST_CGNS_ERROR();
        }
    }
}

#endif /* HAS_CGNS */
