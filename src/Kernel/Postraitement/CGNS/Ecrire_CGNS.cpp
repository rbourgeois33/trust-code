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

#include <Format_Post_CGNS.h>
#include <Comm_Group_MPI.h>
#include <communications.h>
#include <TRUST_2_MED.h>
#include <Ecrire_CGNS.h>
#include <Domaine_VF.h>
#include <Polyedre.h>
#include <Polygone.h>
#include <Domaine.h>
#include <unistd.h>

#ifdef HAS_CGNS

#include <cgns_io.h>

void Ecrire_CGNS::cgns_set_base_name(const Nom& fn)
{
  baseFile_name_ = fn.getString();
}

void Ecrire_CGNS::cgns_associer_domaine_dis(const Domaine_dis_base& domaine_dis_base)
{
  domaine_dis_ = domaine_dis_base;
}

void Ecrire_CGNS::cgns_init_MPI(bool is_self) const
{
#ifdef MPI_

  if (Option_CGNS::FILE_PER_COMM_GROUP && PE_Groups::has_user_defined_group() && !postraiter_domaine_)
    {
      const Comm_Group_MPI& comm_loc = ref_cast(Comm_Group_MPI, PE_Groups::get_user_defined_group());
      if (cgp_mpi_comm(comm_loc.get_mpi_comm()) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_init_MPI : cgp_mpi_comm -- Comm_Group_MPI !" << finl, TRUST_CGNS_ERROR();
    }
  else
    {
      if (cgp_mpi_comm(/* XXX MPI_COMM_WORLD */ Comm_Group_MPI::get_trio_u_world()) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_init_MPI : cgp_mpi_comm -- trio_u_world !" << finl, TRUST_CGNS_ERROR();
    }

  /*
   * Elie Saikali XXX NOTA BENE XXX : We need sometimes to write planes, boundaries where some processors have no elements/vertices ...
   * This wont work because the default PIO mode of CGNS is COLLECTIVE. We want it to be INDEPENDENT !
   */
  if (cgp_pio_mode((CGNS_ENUMT(PIOmode_t)) CGP_INDEPENDENT) != CG_OK)
    Cerr << "Error Ecrire_CGNS::cgns_init_MPI : cgp_pio_mode !" << finl, TRUST_CGNS_ERROR();
#endif
}

void Ecrire_CGNS::cgns_open_file()
{
  if (Process::is_parallel()) cgns_init_MPI(); // 1er truc a faire

  fill_infos_loc();

  if (is_deformable_) /* Si deformable => force to use links ! */
    {
      if (Process::is_parallel() && Option_CGNS::FILE_PER_COMM_GROUP && PE_Groups::has_user_defined_group())
        init_proc_maitre_local_comm();
      else if (!Option_CGNS::USE_LINKS) /* Si deformable et pas FILE_PER_COMM_GROUP/USE_LINKS => force to use links ! */
        {
          Option_CGNS::USE_LINKS = true;
          Option_CGNS::SINGLE_SAFE_FILE = false;
        }
    }

  if (Option_CGNS::USE_LINKS && !postraiter_domaine_)
    return; /* rien a faire si USE_LINKS ou FILE_PER_COMM_GROUP */

  const std::string fn = baseFile_name_ + ".cgns"; // file name

  if (Process::is_parallel())
    cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::PAR>(fn, fileId_);
  else
    cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_);
}

void Ecrire_CGNS::fill_infos_loc()
{
  if (postraiter_domaine_)
    {
      has_elem_som_loc_ = true;
      has_som_field_ = true;
      has_elem_field_ = true;
      return;
    }

  assert (loc_vect_.non_nul());
  assert (static_cast<int>(loc_vect_->size()) <= 3 && static_cast<int>(loc_vect_->size()) > 0);

  for (auto& itr : loc_vect_.valeur())
    {
      if (itr == "FACES") has_faces_field_ = true;
      else if (itr == "SOM") has_som_field_ = true;
      else if (itr == "ELEM") has_elem_field_ = true;
      else
        throw std::runtime_error("Ecrire_CGNS::fill_infos_loc => Unsupported LOC : " + itr);
    }

  if (has_som_field_ && has_elem_field_)
    has_elem_som_loc_ = true;
}

void Ecrire_CGNS::finir_ecriture(double temps)
{
  if (postraiter_domaine_) return; /* rien à faire */

  if (Option_CGNS::USE_LINKS)
    {
      cgns_close_grid_or_solution_link_file(temps, TYPE_LINK_CGNS::SOLUTION);

      /* rewrite the link file so you can visualize during simulation !!! */
      if (is_deformable_)
        cgns_write_final_link_file_pb_deformable();
      else
        cgns_write_final_link_file();
    }
  else
    {
      if (Option_CGNS::SINGLE_SAFE_FILE && singlefile_open_)
        {
          const bool will_flush = (Option_CGNS::FLUSH_EVERY_N > 0) &&
                                  (step_single_file_counter_ % Option_CGNS::FLUSH_EVERY_N == 0);

          const bool will_close = (Option_CGNS::CLOSE_EVERY_N > 0) &&
                                  (step_single_file_counter_ % Option_CGNS::CLOSE_EVERY_N == 0);

          /* single but SAFE file => update iterateurs + close to force flush on disc so you can visualize during simulation !!! */
          if (!first_time_post_ && (will_flush || will_close))
            {
              cgns_update_iterative_singlefile();

              if (!will_close)
                cgns_flush_to_disk();
            }

          if (will_close)
            {
              if (Process::is_parallel())
                cgns_helper_.cgns_close_file<TYPE_RUN_CGNS::PAR>(baseFile_name_ /* inutile */, fileId_, false /* print */);
              else
                cgns_helper_.cgns_close_file<TYPE_RUN_CGNS::SEQ>(baseFile_name_ /* inutile */, fileId_, false /* print */);

              singlefile_open_ = false;
            }
        }
    }
}

void Ecrire_CGNS::cgns_finir()
{
  if (Option_CGNS::USE_LINKS && !postraiter_domaine_)
    return; /* All done */

  if ( Option_CGNS::SINGLE_SAFE_FILE && !singlefile_open_)
    return; /* All done */

  const std::string fn = baseFile_name_ + ".cgns"; // file name

  if (Process::is_parallel())
    {
      if (!postraiter_domaine_)
        {
          if (!Option_CGNS::SINGLE_SAFE_FILE)
            {
              if (Option_CGNS::PARALLEL_OVER_ZONE)
                cgns_write_iters_par_over_zone();
              else
                cgns_write_iters_par_in_zone();
            }
          else
            cgns_update_iterative_singlefile();
        }

      cgns_helper_.cgns_close_file<TYPE_RUN_CGNS::PAR>(fn, fileId_);
    }
  else
    {
      if (!postraiter_domaine_)
        {
          if (!Option_CGNS::SINGLE_SAFE_FILE)
            cgns_write_iters_seq();
          else
            cgns_update_iterative_singlefile();
        }

      cgns_helper_.cgns_close_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_);
    }
}

void Ecrire_CGNS::cgns_add_time(const double t)
{
  if (first_time_post_ && !time_post_.empty())
    first_time_post_ = false;

  if (Option_CGNS::USE_LINKS && !postraiter_domaine_)
    if (!first_time_post_ || is_deformable_) /* Si pas deformable, la 1er fois dans cgns_write_field (fill field_loc_map) */
      cgns_open_solution_link_file(t);

  if (Option_CGNS::SINGLE_SAFE_FILE && !postraiter_domaine_)
    {
      step_single_file_counter_++;

      if (!singlefile_open_)
        {
          if (!first_time_post_)
            {
              std::string fn = baseFile_name_ + ".cgns";

              if (Process::is_parallel())
                cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::PAR, TYPE_MODE_CGNS::MODIFY>(fn, fileId_, /*print*/false);
              else
                cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::SEQ, TYPE_MODE_CGNS::MODIFY>(fn, fileId_, /*print*/false);
            }

          singlefile_open_ = true;
        }
    }

  time_post_.push_back(t); // add time_post
  fieldName_dumped_.clear();
  flowId_elem_++, flowId_som_++, flowId_faces_++; // increment
  fieldId_elem_ = 0, fieldId_som_ = 0, fieldId_faces_ = 0; // reset
  solname_elem_written_ = false, solname_som_written_ = false, solname_faces_written_ = false; // reset
  multi_loc_deformable_support_linked_ = false; // reset
}

void Ecrire_CGNS::cgns_flush_to_disk() const
{
  int cgio_num = 0;

  if (cg_get_cgio(fileId_, &cgio_num) != CG_OK)
    Cerr << "Error Ecrire_CGNS::cgns_flush_to_disk : cg_get_cgio !" << finl, TRUST_CGNS_ERROR();

  if (cgio_flush_to_disk(cgio_num) != 0)
    Cerr << "Error Ecrire_CGNS::cgns_flush_to_disk : cgio_flush_to_disk !" << finl, TRUST_CGNS_ERROR();
}

void Ecrire_CGNS::ensure_modify_open_singlefile()
{
  if (ensure_modify_done_ || Option_CGNS::USE_LINKS || postraiter_domaine_) return;

  const std::string fn = baseFile_name_ + ".cgns";

  /* Close file for first time */
  if (Process::is_parallel())
    cgns_helper_.cgns_close_file<TYPE_RUN_CGNS::PAR>(fn, fileId_, /*print*/ false);
  else
    cgns_helper_.cgns_close_file<TYPE_RUN_CGNS::SEQ>(fn, fileId_, /*print*/ false);

  /* Reopen file for first time with MODIFY mode */
  if (Process::is_parallel())
    cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::PAR, TYPE_MODE_CGNS::MODIFY>(fn, fileId_, /*print*/ false);
  else
    cgns_helper_.cgns_open_file<TYPE_RUN_CGNS::SEQ, TYPE_MODE_CGNS::MODIFY>(fn, fileId_, /*print*/ false);

  ensure_modify_done_ = true;
}

/*  MAJ des iterateurs : pour le mode SINGLE_SAFE_FILE */
void Ecrire_CGNS::cgns_update_iterative_singlefile()
{
  if (Option_CGNS::USE_LINKS || postraiter_domaine_)
    return;

  if (!ensure_modify_done_)
    ensure_modify_open_singlefile();

  const cgsize_t nb_dt_post = static_cast<cgsize_t>(time_post_.size());
  if (!nb_dt_post)
    return;

  std::string grid_ptrs;

  if (is_deformable_)
    {
      std::string one = "GridCoordinates";
      one.resize(CGNS_STR_SIZE, ' ');
      grid_ptrs.reserve(static_cast<size_t>(CGNS_STR_SIZE) * nb_dt_post);
      for (cgsize_t i = 0; i < nb_dt_post; ++i)
        grid_ptrs += one;
    }

  cgsize_t idims[2] = { CGNS_STR_SIZE, nb_dt_post };
  const bool over_zones = Option_CGNS::PARALLEL_OVER_ZONE;

  // helper local : supprimer TOUTES les occurrences d'un DataArray_t <name> sous le noeud courant
  auto delete_all_arrays_named = [](const char *arrname)
  {
    // On est deja positionne sur le parent (BaseIterativeData_t / ZoneIterativeData_t) => Supprime toutes les occurrences
    for (;;)
      {
        if (cg_delete_node(const_cast<char*>(arrname)) == CG_OK) continue;
        break; // plus de noeud de ce nom
      }
  };

  // Pour chaque localisation ecrite (ELEM/SOM/FACES)
  std::vector<int> ind_doms_dumped;
  for (const auto &itr : fld_loc_map_)
    {
      const std::string& LOC = itr.first;
      const Nom& dom = itr.second;
      int ind_base = -123;
      int index_glob = TRUST_2_CGNS::get_index_nom_vector(doms_written_, dom);
      ind_doms_dumped.push_back(index_glob);

      if (has_elem_som_loc_ && LOC != "FACES")
        {
          const Nom dom_mod = TRUST_2_CGNS::modify_domaine_name_for_link(dom, LOC);
          ind_base = TRUST_2_CGNS::get_index_nom_vector(doms_written_, dom_mod);
        }
      else
        ind_base = index_glob;


      // index glob
      // Re-ecrire TimeIterValues avec la TAILLE COURANTE nb_dt_post
      if (cg_biter_write(fileId_, baseId_[index_glob], "TimeIterValues", static_cast<int>(nb_dt_post)) != CG_OK)
        Cerr << "Error cgns_update_iterative_singlefile: cg_biter_write !" << finl, TRUST_CGNS_ERROR();

      // Se positionner sur BaseIterativeData_t
      if (cg_goto(fileId_, baseId_[index_glob], "BaseIterativeData_t", 1, "end") != CG_OK)
        Cerr << "Error cgns_update_iterative_singlefile: cg_goto BaseIterativeData_t !" << finl, TRUST_CGNS_ERROR();

      // Supprimer + Re-ecrire TimeValues
      delete_all_arrays_named("TimeValues");
      if (cg_array_write("TimeValues", CGNS_DOUBLE_TYPE, 1, &nb_dt_post, time_post_.data()) != CG_OK)
        Cerr << "Error cgns_update_iterative_singlefile: cg_array_write TimeValues !" << finl, TRUST_CGNS_ERROR();

      if (cg_simulation_type_write(fileId_, baseId_[index_glob], CGNS_ENUMV(TimeAccurate)) != CG_OK)
        Cerr << "Error cgns_update_iterative_singlefile: cg_simulation_type_write !" << finl, TRUST_CGNS_ERROR();

      const char* solname = (LOC == "SOM") ? solname_som_.c_str() : (LOC == "FACES") ? solname_faces_.c_str() : solname_elem_.c_str();

      if (over_zones)
        {
          const int nb_zones = static_cast<int>(zoneId_par_[ind_base].size());
          for (int iz = 0; iz < nb_zones; ++iz)
            {
              const int zone_ord = iz + 1; // ordinal (1-based)

              if (cg_goto(fileId_, baseId_[ind_base], "Zone_t", zone_ord, "ZoneIterativeData_t", 1, "end") != CG_OK)
                {
                  if (cg_ziter_write(fileId_, baseId_[ind_base], zoneId_par_[ind_base][iz], "ZoneIterativeData") != CG_OK)
                    Cerr << "Error cgns_update_iterative_singlefile: cg_ziter_write (PAR_OVER) !" << finl, TRUST_CGNS_ERROR();
                  if (cg_goto(fileId_, baseId_[ind_base], "Zone_t", zone_ord, "ZoneIterativeData_t", 1, "end") != CG_OK)
                    Cerr << "Error cgns_update_iterative_singlefile: cg_goto ZoneIterativeData_t (PAR_OVER) !" << finl, TRUST_CGNS_ERROR();
                }

              if (is_deformable_)
                {
                  delete_all_arrays_named("GridCoordinatesPointers");
                  if (cg_array_write("GridCoordinatesPointers", CGNS_ENUMV(Character), 2, idims, grid_ptrs.c_str()) != CG_OK)
                    Cerr << "Error cgns_update_iterative_singlefile: cg_array_write GridCoordinatesPointers !" << finl, TRUST_CGNS_ERROR();
                }

              delete_all_arrays_named("FlowSolutionPointers");
              if (cg_array_write("FlowSolutionPointers", CGNS_ENUMV(Character), 2, idims, solname) != CG_OK)
                Cerr << "Error cgns_update_iterative_singlefile: cg_array_write FlowSolutionPointers !" << finl, TRUST_CGNS_ERROR();
            }
        }
      else
        {
          if (cg_goto(fileId_, baseId_[index_glob], "Zone_t", 1, "ZoneIterativeData_t", 1, "end") != CG_OK)
            {
              if (cg_ziter_write(fileId_, baseId_[index_glob], zoneId_[ind_base], "ZoneIterativeData") != CG_OK)
                Cerr << "Error cgns_update_iterative_singlefile: cg_ziter_write !" << finl, TRUST_CGNS_ERROR();
              if (cg_goto(fileId_, baseId_[index_glob], "Zone_t", 1, "ZoneIterativeData_t", 1, "end") != CG_OK)
                Cerr << "Error cgns_update_iterative_singlefile: cg_goto ZoneIterativeData_t !" << finl, TRUST_CGNS_ERROR();
            }

          if (is_deformable_)
            {
              delete_all_arrays_named("GridCoordinatesPointers");
              if (cg_array_write("GridCoordinatesPointers", CGNS_ENUMV(Character), 2, idims, grid_ptrs.c_str()) != CG_OK)
                Cerr << "Error cgns_update_iterative_singlefile: cg_array_write GridCoordinatesPointers !" << finl, TRUST_CGNS_ERROR();
            }

          delete_all_arrays_named("FlowSolutionPointers");
          if (cg_array_write("FlowSolutionPointers", CGNS_ENUMV(Character), 2, idims, solname) != CG_OK)
            Cerr << "Error cgns_update_iterative_singlefile: cg_array_write FlowSolutionPointers !" << finl, TRUST_CGNS_ERROR();
        }
    }

  /* 2 : on iter sur les autres domaines; ie: domaine dis */
  for (int i = 0; i < static_cast<int>(doms_written_.size()); i++)
    {
      if (std::find(ind_doms_dumped.begin(), ind_doms_dumped.end(), i) == ind_doms_dumped.end()) // indice pas dans ind_doms_dumped
        {
          const Nom& nom_dom = doms_written_[i];
          const int ind = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);
          assert(ind > -1);

          cgns_helper_.cgns_write_iters<TYPE_ECRITURE_CGNS::SEQ>(false /* has_field */, 1 /* nb_zones_to_write */, fileId_, baseId_[ind], ind, zoneId_, "rien",
                                                                 solname_som_, solname_elem_,solname_faces_, time_post_);
        }
      else { /* Do Nothing */ }
    }
}

void Ecrire_CGNS::cgns_write_domaine(const Domaine * dom,const Nom& nom_dom, const DoubleTab& som, const IntTab& elem, const Motcle& type_e)
{
  std::string nom_dom_modifie = TRUST_2_CGNS::modify_domaine_name_for_post(nom_dom);

  if (Option_CGNS::USE_LINKS && !postraiter_domaine_ && !is_deformable_)
    if (!grid_file_opened_)
      {
        cgns_open_grid_base_link_file();
        grid_file_opened_ = true; // On ouvre pour .grid.cgns
      }

  if (Process::is_parallel())
    {
      if (Option_CGNS::PARALLEL_OVER_ZONE || postraiter_domaine_)
        cgns_write_domaine_par_over_zone(dom, Nom(nom_dom_modifie), som, elem, type_e);
      else
        cgns_write_domaine_par_in_zone(dom, Nom(nom_dom_modifie), som, elem, type_e);
    }
  else
    cgns_write_domaine_seq(dom, Nom(nom_dom_modifie), som, elem, type_e);
}

void Ecrire_CGNS::cgns_write_field(const Domaine& domaine, const Noms& noms_compo, double temps,
                                   const Nom& id_du_champ, const Nom& id_du_domaine, const Nom& localisation,
                                   const DoubleTab& valeurs)
{
  const std::string LOC = Motcle(localisation).getString();

  if (is_deformable_ && LOC == "FACES")
    {
      has_faces_field_ = false;
      Cerr << "Field " << id_du_champ << " located at " << LOC << " is skipped !! " << finl;
      return; // TODO FIXME
    }

  /* 1 : if first time called ... build different links to support mixed locations */
  if (first_time_post_)
    cgns_fill_field_loc_map(domaine, LOC);

  if (is_deformable_ && has_elem_som_loc_ && !multi_loc_deformable_support_linked_)
    link_multi_loc_support_pb_deformable();

  /* 2 : on ecrit */
  const int nb_cmp = valeurs.dimension(1);

  if (Process::is_parallel())
    {
      for (int i = 0; i < nb_cmp; i++)
        {
          const Motcle field_name = nb_cmp > 1 ? Motcle(noms_compo[i]) : id_du_champ;

          if (std::find(fieldName_dumped_.begin(), fieldName_dumped_.end(), field_name) == fieldName_dumped_.end()) // pas dedans => faut ecrire !
            {
              if (Option_CGNS::PARALLEL_OVER_ZONE || postraiter_domaine_)
                cgns_write_field_par_over_zone(i /* compo */, temps, field_name, id_du_domaine, localisation, fld_loc_map_.at(LOC), valeurs);
              else
                cgns_write_field_par_in_zone(i /* compo */, temps, field_name, id_du_domaine, localisation, fld_loc_map_.at(LOC), valeurs);

              fieldName_dumped_.push_back(field_name);
            }
          else
            Cerr << "Field " << field_name << " is already written => we skip it ..." << finl;
        }
    }
  else
    for (int i = 0; i < nb_cmp; i++)
      {
        const Motcle field_name = nb_cmp > 1 ? Motcle(noms_compo[i]) : id_du_champ;
        if (std::find(fieldName_dumped_.begin(), fieldName_dumped_.end(), field_name) == fieldName_dumped_.end()) // pas dedans => faut ecrire !
          {
            cgns_write_field_seq(i /* compo */, temps, field_name, id_du_domaine, localisation, fld_loc_map_.at(LOC), valeurs);
            fieldName_dumped_.push_back(field_name);
          }
        else
          Cerr << "Field " << field_name << " is already written => we skip it ..." << finl;
      }
}

/*
 * *********************************** *
 * METHODES PRIVEES CLASSE Ecrire_CGNS *
 * *********************************** *
 */

void Ecrire_CGNS::cgns_fill_field_loc_map(const Domaine& domaine, const std::string& LOC)
{
  assert (static_cast<int>(time_post_.size()) == 1 && first_time_post_);

  /* pour les champs aux faces, il faut un support ! */
  if (has_faces_field_ && !is_dual_)
    {
      Nom nom_dom = domaine.le_nom();
      nom_dom += "_FACES";
      Cerr << "###  Building a new CGNS zone to host the fields located at FACES !" << finl;
      cgns_write_domaine_dual(domaine, 0 /* pas premier post ... mais inutile */, nom_dom);
    }

  if (!Option_CGNS::USE_LINKS || postraiter_domaine_)
    {
      Nom nom_dom = domaine.le_nom();
      if (LOC == "FACES" || has_elem_som_loc_)
        {
          nom_dom += "_";
          nom_dom += LOC;
        }

      if (!fld_loc_map_.count(LOC))
        {
          fld_loc_map_.insert( { LOC, nom_dom } );

          if (LOC != "FACES" && has_elem_som_loc_)
            add_new_linked_base(LOC, nom_dom);
        }
    }
  else // Option_CGNS::USE_LINKS
    {
      assert (Option_CGNS::USE_LINKS);
      if (grid_file_opened_ && !is_deformable_)
        {
          cgns_close_grid_or_solution_link_file(0. /* inutile */, TYPE_LINK_CGNS::GRID, false);
          grid_file_opened_ = false;
        }

      if (!solution_file_opened_)
        {
          Nom nom_dom;
          std::string loc_link;

          if (has_elem_field_)
            {
              loc_link = "ELEM";
              assert (!fld_loc_map_.count(loc_link));
              nom_dom = domaine.le_nom();
              if (has_elem_som_loc_)
                nom_dom += "_ELEM";
              fld_loc_map_.insert( { loc_link, nom_dom } );

              if (has_elem_som_loc_)
                cgns_init_solution_link_file(loc_link, nom_dom);
            }

          if (has_som_field_)
            {
              loc_link = "SOM";
              assert (!fld_loc_map_.count(loc_link));
              nom_dom = domaine.le_nom();
              if (has_elem_som_loc_)
                nom_dom += "_SOM";
              fld_loc_map_.insert( { loc_link, nom_dom } );

              if (has_elem_som_loc_)
                cgns_init_solution_link_file(loc_link, nom_dom);
            }

          if (has_faces_field_)
            {
              loc_link = "FACES";
              assert (!fld_loc_map_.count(loc_link));
              nom_dom = domaine.le_nom();
              nom_dom += "_FACES";
              fld_loc_map_.insert( { loc_link, nom_dom } );
            }

          if (!is_deformable_)
            cgns_open_solution_link_file(time_post_.back()); // 1ere ouverture sol file ici ! puis dans cgns_add_time

          solution_file_opened_ = true;
        }
    }
}

/*
 * ******************** *
 * VERSION SEQUENTIELLE *
 * ******************** *
 */
void Ecrire_CGNS::cgns_write_domaine_seq(const Domaine * domaine,const Nom& nom_dom, const DoubleTab& les_som, const IntTab& les_elem, const Motcle& type_elem)
{
  if (is_deformable_ && !first_time_post_)
    {
      cgns_write_domaine_deformable_seq(domaine, nom_dom, les_som, les_elem, type_elem);
      return;
    }

  /* 1 : Instance of TRUST_2_CGNS */
  T2CGNS_.push_back(TRUST_2_CGNS());
  TRUST_2_CGNS& TRUST2CGNS = T2CGNS_.back();
  TRUST2CGNS.associer_domaine_TRUST(domaine, domaine_dis_.non_nul() ? &(domaine_dis_.valeur()) : nullptr, les_som, les_elem, postraiter_domaine_);
  if (is_dual_ && Objet_U::dimension == 3)
    {
      assert(fs_dual_.size() > 0 && ef_dual_.size() > 0);
      TRUST2CGNS.associer_connec_pour_dual(fs_dual_, ef_dual_);
    }

  doms_written_.push_back(nom_dom);

  CGNS_TYPE cgns_type_elem = TRUST2CGNS.convert_elem_type(type_elem);
  const bool is_polyedre = (type_elem == "POLYEDRE" || type_elem == "PRISME" || type_elem == "PRISME_HEXAG");
  const int icelldim = TRUST2CGNS.topo_dim_from_elem(cgns_type_elem, is_polyedre); // avant ca : icelldim = les_som.dimension(1)
  const int iphysdim = Objet_U::dimension, nb_som = les_som.dimension(0), nb_elem = les_elem.dimension(0);

  /* 2 : Fill coords */
  std::vector<double> xCoords, yCoords, zCoords;
  TRUST2CGNS.fill_coords(xCoords, yCoords, zCoords);

  int coordsId;

  /* 3 : Base write */
  baseId_.push_back(-123); // pour chaque dom, on a une baseId
  char basename[CGNS_STR_SIZE];
  strcpy(basename, nom_dom.getChar()); // dom name

  if (cg_base_write(fileId_, basename, icelldim, iphysdim, &baseId_.back()) != CG_OK)
    Cerr << "Error Ecrire_CGNS::cgns_write_domaine_seq : cg_base_write !" << finl, TRUST_CGNS_ERROR();

  /* 4 : Vertex, cell & boundary vertex sizes */
  cgsize_t isize[3][1];
  isize[0][0] = nb_som;
  isize[1][0] = nb_elem;
  isize[2][0] = 0; /* boundary vertex size (zero if elements not sorted) */


  cgns_fill_info_grid_link_file(basename, cgns_type_elem, icelldim, nb_som, nb_elem, is_polyedre);

  zoneId_.push_back(-123);

  /* 5 : Write all */
  if (nb_elem)
    {
      /* 5.1 : Create zone & grid coords */
      cgns_helper_.cgns_write_zone_grid_coord<TYPE_ECRITURE_CGNS::SEQ>(icelldim, fileId_, baseId_.back(), basename /* Dom name */, isize[0],
                                                                       zoneId_.back(), xCoords, yCoords, zCoords, coordsId, coordsId, coordsId);

      /* 5.2 : Set element connectivity */
      int sectionId;
      cgsize_t start = 1, end;

      if (cgns_type_elem == CGNS_ENUMV(NGON_n)) // cas polyedre
        {
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
          std::vector<cgsize_t> elems;
          const int nsom = TRUST2CGNS.convert_connectivity(cgns_type_elem, elems);

          end = start + static_cast<cgsize_t>(elems.size()) / nsom - 1;

          if (cg_section_write(fileId_, baseId_.back(), zoneId_.back(), "Elem", cgns_type_elem, start, end, 0, elems.data(), &sectionId) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_domaine_seq : cg_section_write !" << finl, TRUST_CGNS_ERROR();
        }
    }
  TRUST2CGNS.clear_vectors();
}

void Ecrire_CGNS::cgns_write_field_seq(const int comp, const double temps, const Nom& id_du_champ, const Nom& id_du_domaine, const Nom& localisation, const Nom& nom_dom, const DoubleTab& valeurs)
{
  std::string LOC = Motcle(localisation).getString();
  Motcle id_du_champ_modifie = TRUST_2_CGNS::modify_field_name_for_post(id_du_champ, id_du_domaine, LOC, fieldId_som_, fieldId_elem_, fieldId_faces_);
  Nom& id_champ = id_du_champ_modifie;

  /* 2 : Get corresponding domain index */
  const int ind = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);
  assert(ind > -1);

  const int nb_vals = valeurs.dimension(0);

  if (nb_vals)
    {
      /* 3 : Write solution names for iterative data later */
      cgns_helper_.cgns_sol_write<TYPE_ECRITURE_CGNS::SEQ>(1 /* nb_zones_to_write */, fileId_, baseId_[ind], ind, temps, zoneId_, LOC,
                                                           solname_som_, solname_elem_, solname_faces_,
                                                           solname_som_written_, solname_elem_written_, solname_faces_written_,
                                                           flowId_som_, flowId_elem_, flowId_faces_);

      /* 4 : Fill field values & dump to cgns file */
      if (LOC == "FACES")
        {
          const Domaine_VF& dom_vf = ref_cast(Domaine_VF, domaine_dis_.valeur());
          DoubleTrav new_vals;
          TRUST_2_CGNS::map_face_values(dom_vf, valeurs, new_vals);

          cgns_helper_.cgns_field_write_data<TYPE_ECRITURE_CGNS::SEQ>(fileId_, baseId_[ind], ind, zoneId_, LOC,
                                                                      flowId_som_, flowId_elem_, flowId_faces_, comp,
                                                                      id_champ, new_vals, fieldId_som_, fieldId_elem_, fieldId_faces_);

        }
      else
        cgns_helper_.cgns_field_write_data<TYPE_ECRITURE_CGNS::SEQ>(fileId_, baseId_[ind], ind, zoneId_, LOC,
                                                                    flowId_som_, flowId_elem_, flowId_faces_, comp,
                                                                    id_champ, valeurs, fieldId_som_, fieldId_elem_, fieldId_faces_);
    }
}

void Ecrire_CGNS::cgns_write_iters_seq()
{
  assert(static_cast<int>(baseId_.size()) == static_cast<int>(zoneId_.size()));
  std::vector<int> ind_doms_dumped;

  /* 1 : on iter juste sur le map fld_loc_map_; ie: pas domaine dis ... */
  for (auto &itr : fld_loc_map_)
    {
      const std::string& LOC = itr.first;
      const Nom& nom_dom = itr.second;
      const int ind = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);
      ind_doms_dumped.push_back(ind);
      assert(ind > -1);

      cgns_helper_.cgns_write_iters<TYPE_ECRITURE_CGNS::SEQ>(true /* has_field */, 1 /* nb_zones_to_write */, fileId_, baseId_[ind], ind, zoneId_, LOC,
                                                             solname_som_, solname_elem_, solname_faces_, time_post_);
    }

  /* 2 : on iter sur les autres domaines; ie: domaine dis */
  for (int i = 0; i < static_cast<int>(doms_written_.size()); i++)
    {
      if (std::find(ind_doms_dumped.begin(), ind_doms_dumped.end(), i) == ind_doms_dumped.end()) // indice pas dans ind_doms_dumped
        {
          const Nom& nom_dom = doms_written_[i];
          const int ind = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);
          assert(ind > -1);

          cgns_helper_.cgns_write_iters<TYPE_ECRITURE_CGNS::SEQ>(false /* has_field */, 1 /* nb_zones_to_write */, fileId_, baseId_[ind], ind, zoneId_, "rien",
                                                                 solname_som_, solname_elem_,solname_faces_, time_post_);
        }
      else { /* Do Nothing */ }
    }
}

/*
 * *************************** *
 * VERSION PARALLELE OVER ZONE *
 * *************************** *
 */
void Ecrire_CGNS::cgns_write_domaine_par_over_zone(const Domaine * domaine,const Nom& nom_dom, const DoubleTab& les_som, const IntTab& les_elem, const Motcle& type_elem)
{
#ifdef MPI_
  assert (!Option_CGNS::USE_LINKS || postraiter_domaine_);
  assert (!is_deformable_);
  doms_written_.push_back(nom_dom);

  /* 1 : Instance of TRUST_2_CGNS */
  T2CGNS_.push_back(TRUST_2_CGNS());
  TRUST_2_CGNS& TRUST2CGNS = T2CGNS_.back();
  TRUST2CGNS.associer_domaine_TRUST(domaine, domaine_dis_.non_nul() ? &(domaine_dis_.valeur()) : nullptr,les_som, les_elem, postraiter_domaine_);
  if (is_dual_ && Objet_U::dimension == 3)
    {
      assert(fs_dual_.size() > 0 && ef_dual_.size() > 0);
      TRUST2CGNS.associer_connec_pour_dual(fs_dual_, ef_dual_);
    }
  CGNS_TYPE cgns_type_elem = TRUST2CGNS.convert_elem_type(type_elem);
  const bool is_polyedre = (type_elem == "POLYEDRE" || type_elem == "PRISME" || type_elem == "PRISME_HEXAG");
  const int icelldim = TRUST2CGNS.topo_dim_from_elem(cgns_type_elem, is_polyedre); // avant ca : icelldim = les_som.dimension(1)
  const int iphysdim = Objet_U::dimension, proc_me = Process::me(),
            nb_som = les_som.dimension(0), nb_elem = les_elem.dimension(0);

  /* 2 : Fill coords */
  std::vector<double> xCoords, yCoords, zCoords;
  TRUST2CGNS.fill_coords(xCoords, yCoords, zCoords);


  /* 3 : Base write */
  baseId_.push_back(-123); // pour chaque dom, on a une baseId
  char basename[CGNS_STR_SIZE];
  strcpy(basename, nom_dom.getChar()); // dom name

  if (cg_base_write(fileId_, basename, icelldim, iphysdim, &baseId_.back()) != CG_OK)
    Cerr << "Error Ecrire_CGNS::cgns_write_domaine_par_over_zone : cg_base_write !" << finl, TRUST_CGNS_ERROR();

  /* 4 : We need global nb_elems/nb_soms => MPI_Allgather. Thats the only information required ! */

  cgns_fill_info_grid_link_file(basename, cgns_type_elem, icelldim, nb_som, nb_elem, is_polyedre);

  TRUST2CGNS.fill_global_infos(); // XXX

  if (cgns_type_elem == CGNS_ENUMV(NGON_n)) /*cas polygone/polyedre */
    TRUST2CGNS.fill_global_infos_poly(is_polyedre);

  /* 5 : CREATION OF FILE STRUCTURE : zones, coords & sections
   *
   *  - All processors THAT HAVE nb_elem > 0 write the same information.
   *  - Only zone meta-data is written to the library at this stage ... So no worries ^^
   */
  std::vector<int> coordsIdx, coordsIdy, coordsIdz, sectionId, sectionId2;
  std::string zonename;

  const int nb_zones_to_write = TRUST2CGNS.nb_procs_writing();
  const bool all_write = TRUST2CGNS.all_procs_write(); // all procs will write !

  // on boucle seulement sur les procs qui n'ont pas des nb_elem 0
  zoneId_.clear(); // XXX commencons par ca
  const std::vector<int>& global_nb_elem = TRUST2CGNS.get_global_nb_elem(),
                          &global_nb_som = TRUST2CGNS.get_global_nb_som(),
                           &proc_non_zero_elem = TRUST2CGNS.get_proc_non_zero_elem();

  for (int i = 0; i != nb_zones_to_write; i++)
    {
      const int indZ = all_write ? i : proc_non_zero_elem[i]; // procID
      const int ne_loc = global_nb_elem[indZ], ns_loc = global_nb_som[indZ]; /* nb_elem & nb_som local */
      assert (ne_loc > 0);

      cgsize_t start = 1, end = ne_loc;
      cgsize_t isize[3][1];
      isize[0][0] = ns_loc;
      isize[1][0] = end;
      isize[2][0] = 0; /* boundary vertex size (zero if elements not sorted) */

      zoneId_.push_back(-123);
      zonename = nom_dom.nom_me(indZ).getString();
      zonename.resize(CGNS_STR_SIZE, ' ');

      coordsIdx.push_back(-123), coordsIdy.push_back(-123);
      if (icelldim > 2)
        coordsIdz.push_back(-123);

      /* 5.1 : Create zone & Construct the grid coordinates nodes */
      cgns_helper_.cgns_write_zone_grid_coord<TYPE_ECRITURE_CGNS::PAR_OVER>(icelldim, fileId_, baseId_.back(), zonename.c_str(), isize[0],
                                                                            zoneId_.back(), xCoords, yCoords, zCoords,
                                                                            coordsIdx.back(), coordsIdy.back(), coordsIdz.empty() ? coordsIdy.back() /* inutile */ : coordsIdz.back());

      /* 5.2 : Construct the sections to host connectivity later */
      sectionId.push_back(-123);

      if (cgns_type_elem == CGNS_ENUMV(NGON_n)) // cas polyedre
        {
          if (is_polyedre) // Pas pour polygone
            {
              end = start + static_cast<cgsize_t>(TRUST2CGNS.get_global_nb_face_som()[indZ]) -1;
              cgsize_t maxoffset = static_cast<cgsize_t>(TRUST2CGNS.get_global_nb_face_som_offset()[indZ]);

              if (cgp_poly_section_write(fileId_, baseId_.back(), zoneId_.back(), "NGON_n", CGNS_ENUMV(NGON_n), start, end, maxoffset, 0, &sectionId.back()) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_domaine_par_over_zone : cgp_poly_section_write !" << finl, TRUST_CGNS_ERROR();

              sectionId2.push_back(-123);
              start = end + 1;

              end = start + static_cast<cgsize_t>(TRUST2CGNS.get_global_nb_elem_face()[indZ]) -1;
              maxoffset = static_cast<cgsize_t>(TRUST2CGNS.get_global_nb_elem_face_offset()[indZ]);

              if (cgp_poly_section_write(fileId_, baseId_.back(), zoneId_.back(), "NFACE_n", CGNS_ENUMV(NFACE_n), start, end, maxoffset, 0, &sectionId2.back()) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_domaine_par_over_zone : cgp_poly_section_write !" << finl, TRUST_CGNS_ERROR();
            }
          else // polygon
            {
              end = start + static_cast<cgsize_t>(TRUST2CGNS.get_global_nb_elem()[indZ]) -1; /* ici pareil comme get_global_nb_elem_som ... fais moi confiance ... */
              cgsize_t maxoffset = static_cast<cgsize_t>(TRUST2CGNS.get_global_nb_elem_som_offset()[indZ]);

              if (cgp_poly_section_write(fileId_, baseId_.back(), zoneId_.back(), "NGON_n", CGNS_ENUMV(NGON_n), start, end, maxoffset, 0, &sectionId.back()) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_domaine_par_over_zone : cgp_poly_section_write !" << finl, TRUST_CGNS_ERROR();
            }
        }
      else
        {
          if (cgp_section_write(fileId_, baseId_.back(), zoneId_.back(), "Elem", cgns_type_elem, start, end, 0, &sectionId.back()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_domaine_par_over_zone : cgp_section_write !" << finl, TRUST_CGNS_ERROR();
        }
    }

  zoneId_par_.push_back(zoneId_); // XXX : Dont touch

//  Process::barrier();

  /* 6 : Write grid coordinates & set connectivity */
  if (nb_elem > 0) // this proc will write !
    {
      cgsize_t min = 1, max = nb_som;
      int indx = -123;
      if (all_write) indx = proc_me;
      else
        for (int i = 0; i < nb_zones_to_write; i++)
          if (proc_non_zero_elem[i] == proc_me)
            {
              indx = i;
              break;
            }

      /* 6.1 : Write grid coordinates */
      cgns_helper_.cgns_write_grid_coord_data<TYPE_ECRITURE_CGNS::PAR_OVER>(icelldim, fileId_, baseId_.back(), zoneId_par_.back()[indx],
                                                                            coordsIdx[indx], coordsIdy[indx], coordsIdz.empty() ? -123 : coordsIdz[indx],
                                                                            min, max, xCoords, yCoords, zCoords);

      /* 6.2 : Set element connectivity */
      if (cgns_type_elem == CGNS_ENUMV(NGON_n)) // cas polyedre
        {
          if (is_polyedre) // Pas pour polygone
            {
              const std::vector<cgsize_t>& fs = TRUST2CGNS.get_local_fs(),
                                           &fs_offset = TRUST2CGNS.get_local_fs_offset();

              max = min + TRUST2CGNS.get_nb_fs() - 1;

              if (cgp_poly_elements_write_data(fileId_, baseId_.back(), zoneId_par_.back()[indx], sectionId[indx], min, max, fs.data(), fs_offset.data()) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_domaine_par_over_zone : cgp_poly_elements_write_data !" << finl, TRUST_CGNS_ERROR();

              const std::vector<cgsize_t>& ef = TRUST2CGNS.get_local_ef(),
                                           &ef_offset = TRUST2CGNS.get_local_ef_offset();

              min = max + 1, max = min + TRUST2CGNS.get_nb_ef() - 1;

              if (cgp_poly_elements_write_data(fileId_, baseId_.back(), zoneId_par_.back()[indx], sectionId2[indx], min, max, ef.data(), ef_offset.data()) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_domaine_par_over_zone : cgp_poly_elements_write_data !" << finl, TRUST_CGNS_ERROR();
            }
          else
            {
              const std::vector<cgsize_t>& es = TRUST2CGNS.get_local_es(),
                                           &es_offset = TRUST2CGNS.get_local_es_offset();

              max = min + TRUST2CGNS.get_nb_es() -1;

              if (cgp_poly_elements_write_data(fileId_, baseId_.back(), zoneId_par_.back()[indx], sectionId[indx], min, max, es.data(), es_offset.data()) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_domaine_par_over_zone : cgp_poly_elements_write_data !" << finl, TRUST_CGNS_ERROR();
            }
        }
      else
        {
          std::vector<cgsize_t> elems;
          TRUST2CGNS.convert_connectivity(cgns_type_elem, elems);

          max = nb_elem; /* now we need local elem */
          if (cgp_elements_write_data(fileId_, baseId_.back(), zoneId_par_.back()[indx], sectionId[indx], min, max, elems.data()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_domaine_par_over_zone : cgp_elements_write_data !" << finl, TRUST_CGNS_ERROR();
        }
    }
  TRUST2CGNS.clear_vectors();
#endif
}

void Ecrire_CGNS::cgns_write_field_par_over_zone(const int comp, const double temps, const Nom& id_du_champ, const Nom& id_du_domaine, const Nom& localisation, const Nom& nom_dom, const DoubleTab& valeurs)
{
#ifdef MPI_
  assert (!Option_CGNS::USE_LINKS || postraiter_domaine_);
  std::string LOC = Motcle(localisation).getString();
  Motcle id_du_champ_modifie = TRUST_2_CGNS::modify_field_name_for_post(id_du_champ, id_du_domaine, LOC, fieldId_som_, fieldId_elem_, fieldId_faces_);
  Nom& id_champ = id_du_champ_modifie;

  /* 2 : Get corresponding domain index */
  const int proc_me = Process::me(), nb_vals = valeurs.dimension(0);
  const int ind = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);
  assert(ind > -1);

  int ind_new = ind;
  if (ind > (static_cast<int>(T2CGNS_.size()) -1) )
    {
      const Nom nom_dom_mod = TRUST_2_CGNS::modify_domaine_name_for_link(nom_dom, LOC);
      ind_new = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom_mod);
    }

  /* 3 : CREATION OF FILE STRUCTURE
   *
   *  - All processors THAT HAVE nb_vals > 0 write the same information.
   *  - Only field meta-data is written to the library at this stage ... So no worries ^^
   *  - And just once per dt !
   */
  const TRUST_2_CGNS& TRUST2CGNS = T2CGNS_[ind_new];

  const int nb_zones_to_write = TRUST2CGNS.nb_procs_writing();
  const bool all_write = TRUST2CGNS.all_procs_write(); // all procs will write !

  cgns_helper_.cgns_sol_write<TYPE_ECRITURE_CGNS::PAR_OVER>(nb_zones_to_write, fileId_, baseId_[ind], ind, temps, zoneId_par_[ind], LOC,
                                                            solname_som_, solname_elem_, solname_faces_,
                                                            solname_som_written_, solname_elem_written_, solname_faces_written_,
                                                            flowId_som_, flowId_elem_, flowId_faces_);

  cgns_helper_.cgns_field_write<TYPE_ECRITURE_CGNS::PAR_OVER>(nb_zones_to_write, fileId_, baseId_[ind], ind, zoneId_par_[ind], LOC,
                                                              flowId_som_, flowId_elem_, flowId_faces_, id_champ.getChar(),
                                                              fieldId_som_, fieldId_elem_, fieldId_faces_);

  /* 4 : Fill field values & dump to cgns file */
  if (nb_vals > 0) // this proc will write !
    {
      cgsize_t min = 1, max = nb_vals;
      int indx = -123;
      const std::vector<int>& proc_non_zero_write= TRUST2CGNS.get_proc_non_zero_elem();
      if (all_write) indx = proc_me;
      else
        for (int i = 0; i < nb_zones_to_write; i++)
          if (proc_non_zero_write[i] == proc_me)
            {
              indx = i;
              break;
            }

      if (LOC == "FACES")
        {
          const Domaine_VF& dom_vf = ref_cast(Domaine_VF, domaine_dis_.valeur());
          DoubleTrav new_vals;
          TRUST_2_CGNS::map_face_values(dom_vf, valeurs, new_vals);

          max = new_vals.dimension(0); // XXX

          cgns_helper_.cgns_field_write_data<TYPE_ECRITURE_CGNS::PAR_OVER>(fileId_, baseId_[ind], indx /* XXX */, zoneId_par_[ind], LOC,
                                                                           flowId_som_, flowId_elem_, flowId_faces_,
                                                                           fieldId_som_, fieldId_elem_, fieldId_faces_,
                                                                           comp, min, max, new_vals);
        }
      else
        cgns_helper_.cgns_field_write_data<TYPE_ECRITURE_CGNS::PAR_OVER>(fileId_, baseId_[ind], indx /* XXX */, zoneId_par_[ind], LOC,
                                                                         flowId_som_, flowId_elem_, flowId_faces_,
                                                                         fieldId_som_, fieldId_elem_, fieldId_faces_,
                                                                         comp, min, max, valeurs);
    }
#endif
}

void Ecrire_CGNS::cgns_write_iters_par_over_zone()
{
#ifdef MPI_
//  assert(static_cast<int>(baseId_.size()) == static_cast<int>(zoneId_.size())); // XXX No not for // !!!
  std::vector<int> ind_doms_dumped;

  /* 1 : on iter juste sur le map fld_loc_map_; ie: pas domaine dis ... */
  for (auto &itr : fld_loc_map_)
    {
      const std::string& LOC = itr.first;
      const Nom& nom_dom = itr.second;
      const int ind = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);
      ind_doms_dumped.push_back(ind);
      assert(ind > -1);

      int ind_new = ind;
      if (ind > (static_cast<int>(T2CGNS_.size()) -1) )
        {
          const Nom nom_dom_mod = TRUST_2_CGNS::modify_domaine_name_for_link(nom_dom, LOC);
          ind_new = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom_mod);
        }

      const TRUST_2_CGNS& TRUST2CGNS = T2CGNS_[ind_new];
      const int nb_zones_to_write = TRUST2CGNS.nb_procs_writing();

      cgns_helper_.cgns_write_iters<TYPE_ECRITURE_CGNS::PAR_OVER>(true /* has_field */, nb_zones_to_write, fileId_, baseId_[ind], ind, zoneId_par_[ind], LOC,
                                                                  solname_som_, solname_elem_, solname_faces_, time_post_);
    }

  /* 2 : on iter sur les autres domaines; ie: domaine dis */
  for (int i = 0; i < static_cast<int>(doms_written_.size()); i++)
    {
      if (std::find(ind_doms_dumped.begin(), ind_doms_dumped.end(), i) == ind_doms_dumped.end()) // indice pas dans ind_doms_dumped
        {
          const Nom& nom_dom = doms_written_[i];
          const int ind = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);
          assert(ind > -1);

          int ind_new = ind;
          if (ind > (static_cast<int>(T2CGNS_.size()) -1) )
            {
              Nom nom_dom_mod = nom_dom;
              if (nom_dom.finit_par("_ELEM"))
                nom_dom_mod = TRUST_2_CGNS::modify_domaine_name_for_link(nom_dom, "ELEM");
              else if (nom_dom.finit_par("_SOM"))
                nom_dom_mod = TRUST_2_CGNS::modify_domaine_name_for_link(nom_dom, "SOM");

              ind_new = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom_mod);
            }

          const TRUST_2_CGNS& TRUST2CGNS = T2CGNS_[ind_new];
          const int nb_zones_to_write = TRUST2CGNS.nb_procs_writing();

          cgns_helper_.cgns_write_iters<TYPE_ECRITURE_CGNS::PAR_OVER>(false /* has_field */, nb_zones_to_write, fileId_, baseId_[ind], ind, zoneId_par_[ind], "rien",
                                                                      solname_som_, solname_elem_, solname_faces_, time_post_);
        }
      else { /* Do Nothing */ }
    }
#endif
}

/*
 * ************************* *
 * VERSION PARALLELE IN ZONE *
 * ************************* *
 */
void Ecrire_CGNS::cgns_write_domaine_par_in_zone(const Domaine * domaine,const Nom& nom_dom, const DoubleTab& les_som, const IntTab& les_elem, const Motcle& type_elem)
{
#ifdef MPI_

  if (is_deformable_ && !first_time_post_)
    {
      cgns_write_domaine_deformable_par_in_zone(domaine, nom_dom, les_som, les_elem, type_elem);
      return;
    }

  doms_written_.push_back(nom_dom);

  /* 1 : Instance of TRUST_2_CGNS */
  T2CGNS_.push_back(TRUST_2_CGNS());
  TRUST_2_CGNS& TRUST2CGNS = T2CGNS_.back();
  TRUST2CGNS.associer_domaine_TRUST(domaine, domaine_dis_.non_nul() ? &(domaine_dis_.valeur()) : nullptr, les_som, les_elem, postraiter_domaine_);
  if (is_dual_ && Objet_U::dimension == 3)
    {
      assert(fs_dual_.size() > 0 && ef_dual_.size() > 0);
      TRUST2CGNS.associer_connec_pour_dual(fs_dual_, ef_dual_);
    }
  CGNS_TYPE cgns_type_elem = TRUST2CGNS.convert_elem_type(type_elem);
  const bool is_polyedre = (type_elem == "POLYEDRE" || type_elem == "PRISME" || type_elem == "PRISME_HEXAG");
  const int icelldim = TRUST2CGNS.topo_dim_from_elem(cgns_type_elem, is_polyedre); // avant ca : icelldim = les_som.dimension(1)
  const int nb_elem = les_elem.dimension(0), iphysdim = Objet_U::dimension;

  /* 2 : Fill coords */
  std::vector<double> xCoords, yCoords, zCoords;
  TRUST2CGNS.fill_coords(xCoords, yCoords, zCoords);

  /* 3 : Base write */
  baseId_.push_back(-123); // pour chaque dom, on a une baseId
  char basename[CGNS_STR_SIZE];
  strcpy(basename, nom_dom.getChar()); // dom name

  if (cg_base_write(fileId_, basename, icelldim, iphysdim, &baseId_.back()) != CG_OK)
    Cerr << "Error Ecrire_CGNS::cgns_write_domaine_par_in_zone : cg_base_write !" << finl, TRUST_CGNS_ERROR();

  /* 4 : CREATION OF FILE STRUCTURE : zones, coords & sections
   *
   *  - All processors write the same information.
   *  - XXX XXX XXX Only ONE zone meta-data is written to the library at this stage ...
   */

  TRUST2CGNS.fill_global_infos(); // XXX

  const bool enter_group_comm = Option_CGNS::FILE_PER_COMM_GROUP && PE_Groups::has_user_defined_group() && !postraiter_domaine_;
  const int proc_me = enter_group_comm ? TRUST2CGNS.get_proc_me_local_comm() : Process::me();

  if (cgns_type_elem == CGNS_ENUMV(NGON_n)) /*cas polygone/polyedre */
    TRUST2CGNS.fill_global_infos_poly(is_polyedre);

  const int ns_tot = TRUST2CGNS.get_ns_tot(), ne_tot = TRUST2CGNS.get_ne_tot();

  assert (enter_group_comm || (!enter_group_comm && ns_tot > 0 && ne_tot > 0));

  /* 4.1 : Create zone & grid */
  cgsize_t isize[3][1];
  isize[0][0] = (ns_tot == 0 && enter_group_comm) ? 1 : ns_tot; // si ns_tot = 0, on va juste creer une zone vide
  isize[1][0] = (ne_tot == 0 && enter_group_comm) ? 1 : ne_tot; // si ne_tot = 0, on va juste creer une zone vide
  isize[2][0] = 0; /* boundary vertex size (zero if elements not sorted) */

  cgns_fill_info_grid_link_file(basename, cgns_type_elem, icelldim,
                                (ns_tot == 0 && enter_group_comm) ? 1 : ns_tot,
                                (ne_tot == 0 && enter_group_comm) ? 1 : ne_tot,
                                is_polyedre);

  int coordsIdx = -123, coordsIdy = -123, coordsIdz = -123, sectionId = -123, sectionId2 = -123;
  zoneId_.push_back(-123);

  cgns_helper_.cgns_write_zone_grid_coord<TYPE_ECRITURE_CGNS::PAR_IN>(icelldim, fileId_, baseId_.back(), basename /* Dom name */, isize[0],
                                                                      zoneId_.back(), xCoords, yCoords, zCoords, coordsIdx, coordsIdy, coordsIdz);

  if (ne_tot == 0 && ns_tot == 0) return; // XXX Elie Saikali : zone vide creer, rien a faire de plus ... (cas FILE_PER_COMM_GROUP !!!)

  /* 4.2 : Construct the sections to host connectivity later */
  cgsize_t start = -123, end = -123;
  if (cgns_type_elem == CGNS_ENUMV(NGON_n)) // cas polyedre
    {
      cgsize_t maxoffset = -123;

      if (is_polyedre) // Pas pour polygone
        {
          const int nb_fs = TRUST2CGNS.get_nfs_tot();
          const int nb_fs_offset = TRUST2CGNS.get_nfs_offset_tot();

          start = 1, end = start + nb_fs - 1;
          maxoffset = nb_fs_offset;
          assert(start <= end);

          if (cgp_poly_section_write(fileId_, baseId_.back(), zoneId_.back(), "NGON_n", CGNS_ENUMV(NGON_n), start, end, maxoffset, 0, &sectionId) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_domaine_par_in_zone : cgp_poly_section_write !" << finl, TRUST_CGNS_ERROR();

          const int nb_ef = TRUST2CGNS.get_nef_tot();
          const int nb_ef_offset = TRUST2CGNS.get_nef_offset_tot();

          start = end + 1, end = start + nb_ef - 1;
          maxoffset = nb_ef_offset;
          assert(start <= end);

          if (cgp_poly_section_write(fileId_, baseId_.back(), zoneId_.back(), "NFACE_n", CGNS_ENUMV(NFACE_n), start, end, maxoffset, 0, &sectionId2) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_domaine_par_in_zone : cgp_poly_section_write !" << finl, TRUST_CGNS_ERROR();
        }
      else // polygon
        {
          const int nb_es = ne_tot;
          const int nb_es_offset = TRUST2CGNS.get_nes_offset_tot();

          start = 1, end = start + nb_es - 1;
          maxoffset = nb_es_offset;

          if (cgp_poly_section_write(fileId_, baseId_.back(), zoneId_.back(), "NGON_n", CGNS_ENUMV(NGON_n), start, end, maxoffset, 0, &sectionId) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_domaine_par_in_zone : cgp_poly_section_write !" << finl, TRUST_CGNS_ERROR();
        }
    }
  else
    {
      start = 1, end = ne_tot;
      assert(start <= end);

      if (cgp_section_write(fileId_, baseId_.back(), zoneId_.back(), "Elem", cgns_type_elem, start, end, 0, &sectionId) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_domaine_par_in_zone : cgp_section_write !" << finl, TRUST_CGNS_ERROR();
    }

  /* 5 : Write grid coordinates & set connectivity */
  if (nb_elem > 0) // seulement si le proc a qlq chose a ecrire
    {
      const std::vector<int>& incr_max_som = TRUST2CGNS.get_global_incr_max_som(),
                              &incr_min_som = TRUST2CGNS.get_global_incr_min_som();

      cgsize_t min = incr_min_som[proc_me], max = incr_max_som[proc_me];
      assert (min < max);

      /* 5.1 : Write grid coordinates */
      cgns_helper_.cgns_write_grid_coord_data<TYPE_ECRITURE_CGNS::PAR_IN>(icelldim, fileId_, baseId_.back(), zoneId_.back(),
                                                                          coordsIdx, coordsIdy, coordsIdz, min, max, xCoords, yCoords, zCoords);

      /* 5.2 : Set element connectivity */
      if (cgns_type_elem == CGNS_ENUMV(NGON_n)) // cas polyedre
        {
          if (is_polyedre)
            {
              const std::vector<cgsize_t>& fs = TRUST2CGNS.get_local_fs(),
                                           &fs_offset = TRUST2CGNS.get_local_fs_offset();

              const std::vector<int>& incr_min_face_som = TRUST2CGNS.get_global_incr_min_face_som(),
                                      &incr_max_face_som = TRUST2CGNS.get_global_incr_max_face_som();

              min = incr_min_face_som[proc_me], max = incr_max_face_som[proc_me];
              assert (min < max);

              if (cgp_poly_elements_write_data(fileId_, baseId_.back(), zoneId_.back(), sectionId, min, max, fs.data(), fs_offset.data()) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_domaine_par_in_zone : cgp_poly_elements_write_data !" << finl, TRUST_CGNS_ERROR();

              const std::vector<cgsize_t>& ef = TRUST2CGNS.get_local_ef(),
                                           &ef_offset = TRUST2CGNS.get_local_ef_offset();

              const std::vector<int>& incr_min_elem_face = TRUST2CGNS.get_global_incr_min_elem_face(),
                                      &incr_max_elem_face = TRUST2CGNS.get_global_incr_max_elem_face();

              min = incr_max_face_som.back() + incr_min_elem_face[proc_me]; // BOOM
              max = incr_max_face_som.back() + incr_max_elem_face[proc_me]; // BEEM
              assert (min <= max);

              if (cgp_poly_elements_write_data(fileId_, baseId_.back(), zoneId_.back(), sectionId2, min, max, ef.data(), ef_offset.data()) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_domaine_par_in_zone : cgp_poly_elements_write_data !" << finl, TRUST_CGNS_ERROR();
            }
          else
            {
              const std::vector<cgsize_t>& es = TRUST2CGNS.get_local_es(),
                                           &es_offset = TRUST2CGNS.get_local_es_offset();

              const std::vector<int>& incr_max_elem = TRUST2CGNS.get_global_incr_max_elem(),
                                      &incr_min_elem = TRUST2CGNS.get_global_incr_min_elem();

              min = incr_min_elem[proc_me], max = incr_max_elem[proc_me];
              assert (min <= max);

              if (cgp_poly_elements_write_data(fileId_, baseId_.back(), zoneId_.back(), sectionId, min, max, es.data(), es_offset.data()) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_domaine_par_in_zone : cgp_poly_elements_write_data !" << finl, TRUST_CGNS_ERROR();
            }
        }
      else
        {
          std::vector<cgsize_t> elems;
          TRUST2CGNS.convert_connectivity(cgns_type_elem, elems);

          const std::vector<int>& incr_max_elem = TRUST2CGNS.get_global_incr_max_elem(),
                                  &incr_min_elem = TRUST2CGNS.get_global_incr_min_elem();

          min = incr_min_elem[proc_me], max = incr_max_elem[proc_me];
          assert (min <= max);

          if (cgp_elements_write_data(fileId_, baseId_.back(), zoneId_.back(), sectionId, min, max, elems.data()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_domaine_par_in_zone : cgp_elements_write_data !" << finl, TRUST_CGNS_ERROR();
        }
    }
  TRUST2CGNS.clear_vectors();
#endif
}

void Ecrire_CGNS::cgns_write_field_par_in_zone(const int comp, const double temps, const Nom& id_du_champ, const Nom& id_du_domaine, const Nom& localisation, const Nom& nom_dom, const DoubleTab& valeurs)
{
#ifdef MPI_
  const int nb_vals = valeurs.dimension(0);
  const int ind = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom);
  assert(ind > -1);

  std::string LOC = Motcle(localisation).getString();
  Motcle id_du_champ_modifie = TRUST_2_CGNS::modify_field_name_for_post(id_du_champ, id_du_domaine, LOC, fieldId_som_, fieldId_elem_, fieldId_faces_);
  Nom& id_champ = id_du_champ_modifie;

  /* 1 : CREATION OF FILE STRUCTURE
   *
   *  - All processors write the same information.
   *  - Only field meta-data is written to the library at this stage ... So no worries ^^
   *  - And just once per dt !
   */
  cgns_helper_.cgns_sol_write<TYPE_ECRITURE_CGNS::PAR_IN>(1 /* nb_zones_to_write */, fileId_, baseId_[ind], ind, temps, zoneId_, LOC,
                                                          solname_som_, solname_elem_, solname_faces_,
                                                          solname_som_written_, solname_elem_written_, solname_faces_written_,
                                                          flowId_som_, flowId_elem_, flowId_faces_);

  cgns_helper_.cgns_field_write<TYPE_ECRITURE_CGNS::PAR_IN>(1 /* nb_zones_to_write */, fileId_, baseId_[ind], ind, zoneId_, LOC,
                                                            flowId_som_, flowId_elem_,flowId_faces_,
                                                            id_champ.getChar(), fieldId_som_, fieldId_elem_, fieldId_faces_);

  /* 2 : Fill field values & dump to cgns file */
  if (nb_vals > 0) // this proc will write !
    {
      int ind_new = ind;

      if (ind > (static_cast<int>(T2CGNS_.size()) -1) )
        {
          const Nom nom_dom_mod = TRUST_2_CGNS::modify_domaine_name_for_link(nom_dom, LOC);
          ind_new = TRUST_2_CGNS::get_index_nom_vector(doms_written_, nom_dom_mod);
        }

      const TRUST_2_CGNS& TRUST2CGNS = T2CGNS_[ind_new];
      const bool enter_group_comm = Option_CGNS::FILE_PER_COMM_GROUP && PE_Groups::has_user_defined_group() && !postraiter_domaine_;
      const int proc_me = enter_group_comm ? TRUST2CGNS.get_proc_me_local_comm() : Process::me();

      cgsize_t min = -123, max = -123;

      if (LOC == "SOM")
        {
          const std::vector<int>& incr_max_som = TRUST2CGNS.get_global_incr_max_som(),
                                  &incr_min_som = TRUST2CGNS.get_global_incr_min_som();

          min = incr_min_som[proc_me], max = incr_max_som[proc_me];
        }
      else
        {
          const std::vector<int>& incr_max_elem = TRUST2CGNS.get_global_incr_max_elem(),
                                  &incr_min_elem = TRUST2CGNS.get_global_incr_min_elem();

          min = incr_min_elem[proc_me], max = incr_max_elem[proc_me];
        }

      if (LOC == "FACES")
        {
          const Domaine_VF& dom_vf = ref_cast(Domaine_VF, domaine_dis_.valeur());
          DoubleTrav new_vals;
          TRUST_2_CGNS::map_face_values(dom_vf, valeurs, new_vals);

          cgns_helper_.cgns_field_write_data<TYPE_ECRITURE_CGNS::PAR_IN>(fileId_, baseId_[ind], ind, zoneId_, LOC,
                                                                         flowId_som_, flowId_elem_, flowId_faces_,
                                                                         fieldId_som_, fieldId_elem_, fieldId_faces_,
                                                                         comp, min, max, new_vals);
        }
      else
        cgns_helper_.cgns_field_write_data<TYPE_ECRITURE_CGNS::PAR_IN>(fileId_, baseId_[ind], ind, zoneId_, LOC,
                                                                       flowId_som_, flowId_elem_, flowId_faces_,
                                                                       fieldId_som_, fieldId_elem_, fieldId_faces_,
                                                                       comp, min, max, valeurs);
    }
#endif
}

void Ecrire_CGNS::cgns_write_iters_par_in_zone()
{
  cgns_write_iters_seq();
}

/*
 * *************** *
 * Write Dual Mesh *
 * *************** *
 */
void Ecrire_CGNS::cgns_write_domaine_dual(const Domaine& domaine, const int est_le_premier_post, const Nom& nom_dom_faces)
{
  Cerr << "Writing the Dual mesh of " << domaine.le_nom() << " in a CGNS format ..." << finl;
  assert(domaine_dis_.non_nul());
  if (Objet_U::dimension==0)
    Process::exit("Dimension is not defined. Check your data file.");
  const Domaine_VF& dom_vf = ref_cast(Domaine_VF, domaine_dis_.valeur());
  const auto& dual_m = dom_vf.get_mc_dual_mesh();

  // Check the mesh
#ifndef NDEBUG
  dual_m->checkConsistency();
#endif

  const Nom dom_dual_nom = (nom_dom_faces != "??") ? nom_dom_faces : Nom(dual_m->getName());
  Domaine dom_dual;
  dom_dual.nommer(dom_dual_nom);

  DoubleTab sommets;

  // Get the nodes: size and fill sommets
  int nnodes = static_cast<int>(dual_m->getNumberOfNodes());
  const double *coord = dual_m->getCoords()->begin();
  sommets.resize(nnodes, Objet_U::dimension);
  std::copy(coord, coord+sommets.size_array(), sommets.addr());

  // Get cell connectivity
  int ncells = static_cast<int>(dual_m->getNumberOfCells());

  ArrOfInt conn, connIndex;
  int conn_size = static_cast<int>(dual_m->getNodalConnectivity()->getNbOfElems()),
      conn_indx_size= static_cast<int>(dual_m->getNodalConnectivityIndex()->getNbOfElems());

  const auto *c  = dual_m->getNodalConnectivity()->begin(),
              *cI = dual_m->getNodalConnectivityIndex()->begin();

  conn.resize(conn_size);
  std::copy(c, c + conn_size, conn.addr());
  connIndex.resize(conn_indx_size);
  std::copy(cI, cI + conn_indx_size, connIndex.addr());

  int mesh_type_cell = static_cast<int>(conn[connIndex[0]]);  // type is always an int.
  Motcle type_cell;

  if (mesh_type_cell == INTERP_KERNEL::NORM_TRI3)
    type_cell = "Triangle";
  else if (mesh_type_cell == INTERP_KERNEL::NORM_POLYHED)
    type_cell = "Polyedre";
  else
    {
      Cerr << "Cell type " << mesh_type_cell << " is not supported yet. It should be only triangle (2D) and polyedre (3D). Call the 911 !!" << finl;
      Process::exit();
    }

  Elem_geom type_elem;
  type_elem.typer(type_cell);

  IntTab les_elems;
  // Fill les_elem : Different treatment according type_elem:
  if (sub_type(Polyedre, type_elem.valeur()))
    {
      int marker = 0;
      for (int i = 0; i < conn_size; i++)
        if (conn[i]<0) marker++;
      int num_nodes = conn_size - ncells - marker;
      int nfaces = ncells + marker;
      ArrOfInt nodes(num_nodes), facesIndex(nfaces+1), polyhedronIndex(ncells+1);
      int face=0, node = 0;
      for (int i = 0; i < ncells; i++)
        {
          polyhedronIndex[i] = face; // Index des polyedres

          const int index = connIndex[i] + 1;
          const int nb_som = static_cast<int>(connIndex[i + 1] - index);
          for (int j = 0; j < nb_som; j++)
            {
              if (j==0 || conn[index + j]<0)
                facesIndex[face++] = node; // Index des faces:
              if (conn[index + j]>=0)
                nodes[node++] = conn[index + j]; // Index local des sommets de la face
            }
        }
      facesIndex[nfaces] = node;
      polyhedronIndex[ncells] = face;
      ref_cast(Polyedre,type_elem.valeur()).affecte_connectivite_numero_global(nodes, facesIndex, polyhedronIndex, les_elems);
    }
  else // Tous les autres types
    {
      for (int i = 0; i < ncells; i++)
        {
          const int index = connIndex[i] + 1;
          const int nb_som = static_cast<int>(connIndex[i + 1] - index);
          if (i==0) les_elems.resize(ncells, nb_som); // Size les_elems2
          for (int j = 0; j < nb_som; j++)
            les_elems(i, j) = conn[index + j];
        }
    }

  // Converting from MED to TRUST connectivity
  conn_trust_to_med(les_elems,type_elem->que_suis_je(),false);

  dom_dual.les_sommets() = sommets; // fill sommets
  dom_dual.type_elem() = type_elem;

  dom_dual.type_elem()->associer_domaine(dom_dual);
  dom_dual.les_elems() = les_elems;

  // write dual_mesh
  is_dual_ = true;
  // we fill face/som & elem faces conn ET seulement si poly !!!
  if (Objet_U::dimension == 3)
    fill_connectivity_from_mc_mesh(dual_m, fs_dual_, ef_dual_);

  cgns_write_domaine(&dom_dual, dom_dual_nom, sommets, les_elems, type_cell);
}

#endif /* HAS_CGNS */
