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
#include <Ecrire_CGNS.h>
#include <Domaine.h>
#include <unistd.h>

#ifdef HAS_CGNS

/*
 * ***************** *
 * METHODS POUR LINK *
 * ***************** *
 */
void Ecrire_CGNS::cgns_fill_info_grid_link_file(const char* basename, const CGNS_TYPE& cgns_type_elem, const int icelldim, const int nb_som, const int nb_elem, const bool is_polyedre)
{
  if (connectname_.empty())
    {
      cellDim_ = icelldim;
      baseZone_name_ = std::string(basename);
      sizeId_ = { nb_som, nb_elem };
      if (cgns_type_elem == CGNS_ENUMV(NGON_n)) // cas polyedre
        {
          if (is_polyedre)
            connectname_ = { "NGON_n", "NFACE_n" };
          else
            connectname_.push_back("NGON_n");
        }
      else
        connectname_.push_back("Elem");
    }
}

void Ecrire_CGNS::cgns_open_close_link_files(const double t)
{
  if (grid_file_opened_)
    {
      cgns_close_grid_solution_link_file(0 /* only one index here */, baseFile_name_ + ".grid.cgns", true);
      grid_file_opened_ = false;
    }

  if (!time_post_.empty()) /* 1er fois, on fais dans cgns_write_field => fill field_loc_map */
    {
      assert ((static_cast<int>(fld_loc_map_.size()) <= 2)); // ELEM, SOM au max pour le moment
      for (auto itr = fld_loc_map_.begin(); itr != fld_loc_map_.end(); ++itr)
        {
          const int ind = static_cast<int>(std::distance(fld_loc_map_.begin(), itr));
          cgns_close_grid_solution_link_file(ind, baseFile_name_);
          cgns_open_solution_link_file(ind, itr->first, t);
        }
    }
}

void Ecrire_CGNS::cgns_open_grid_base_link_file()
{
  assert(Option_CGNS::USE_LINKS && !postraiter_domaine_);
  std::string fn = baseFile_name_ + ".grid.cgns"; // file name
  unlink(fn.c_str());
  if (Process::is_parallel())
    {
#ifdef MPI_
      if (cgp_open(fn.c_str(), CG_MODE_WRITE, &fileId_) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_open_grid_file : cgp_open !" << finl, TRUST_CGNS_ERROR();

      Cerr << "**** Parallel CGNS file " << fn << " opened !" << finl;
#endif
    }
  else
    {
      if (cg_open(fn.c_str(), CG_MODE_WRITE, &fileId_) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_open_grid_file : cg_open !" << finl, TRUST_CGNS_ERROR();

      Cerr << "**** CGNS file " << fn << " opened !" << finl;
    }
}

void Ecrire_CGNS::cgns_open_solution_link_file(const int ind, const std::string& LOC, const double t, bool is_link)
{
  assert(Option_CGNS::USE_LINKS && !postraiter_domaine_);
  const bool mult_loc = (static_cast<int>(fld_loc_map_.size()) > 1);

  std::string fn;
  True_int& fileId = (ind == 0 ? fileId_ : fileId2_); // XXX : ref

  if (is_link)
    fn = !mult_loc ? baseFile_name_ + ".cgns" : baseFile_name_ + "_" + LOC + ".cgns"; // file name
  else
    fn = baseFile_name_ + "_" + LOC + ".solution." + cgns_helper_.convert_double_to_string(t) + ".cgns"; // file name

  unlink(fn.c_str());

  if (Process::is_parallel())
    {
#ifdef MPI_
      if (cgp_open(fn.c_str(), CG_MODE_WRITE, &fileId) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_open_solution_file : cgp_open !" << finl, TRUST_CGNS_ERROR();

      if (is_link)
        Cerr << "**** Parallel CGNS file " << fn << " opened !" << finl;
#endif
    }
  else
    {
      if (cg_open(fn.c_str(), CG_MODE_WRITE, &fileId) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_open_solution_file : cg_open !" << finl, TRUST_CGNS_ERROR();

      if (is_link)
        Cerr << "**** CGNS file " << fn << " opened !" << finl;
    }

  if (cg_base_write(fileId, baseZone_name_.c_str(), cellDim_, Objet_U::dimension, &baseId_[0]) != CG_OK)
    Cerr << "Error Ecrire_CGNS::cgns_open_solution_file : cg_base_write !" << finl, TRUST_CGNS_ERROR();

  cgsize_t isize[3][1];
  isize[0][0] = sizeId_[0];
  isize[1][0] = sizeId_[1];
  isize[2][0] = 0;

  if (cg_zone_write(fileId, baseId_[0], baseZone_name_.c_str() /* Dom name */, isize[0], CGNS_ENUMV(Unstructured), &zoneId_[0]) != CG_OK)
    Cerr << "Error Ecrire_CGNS::cgns_write_domaine_seq : cgns_open_solution_file !" << finl, TRUST_CGNS_ERROR();

  std::string linkfile = baseFile_name_ + ".grid.cgns"; // file name
  linkfile = TRUST_2_CGNS::remove_slash_linkfile(linkfile);

  std::string linkpath = "/" + baseZone_name_ + "/" + baseZone_name_ + "/GridCoordinates/";

  if (cg_goto(fileId, baseId_[0], "Zone_t", 1, "end") != CG_OK)
    Cerr << "Error Ecrire_CGNS::cgns_open_solution_file : cg_goto !" << finl, TRUST_CGNS_ERROR();

  if (cg_link_write("GridCoordinates", linkfile.c_str(), linkpath.c_str()) != CG_OK)
    Cerr << "Error Ecrire_CGNS::cgns_open_solution_file : cg_link_write !" << finl, TRUST_CGNS_ERROR();

  for (auto &itr : connectname_)
    {
      linkpath = "/" + baseZone_name_ + "/" + baseZone_name_ + "/" + itr + "/";

      if (cg_link_write(itr.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_open_solution_file : cg_link_write !" << finl, TRUST_CGNS_ERROR();
    }
}

void Ecrire_CGNS::cgns_close_grid_solution_link_file(const int ind, const std::string& fn, bool is_cerr)
{
  assert(Option_CGNS::USE_LINKS && !postraiter_domaine_);
  const True_int fileId = (ind == 0 ? fileId_ : fileId2_);

  if (Process::is_parallel())
    {
#ifdef MPI_
      if (cgp_close(fileId) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_close_solution_file : cgp_close !" << finl, TRUST_CGNS_ERROR();

      if (is_cerr) Cerr << "**** Parallel CGNS file " << fn << " closed !" << finl;
#endif
    }
  else
    {
      if (cg_close(fileId) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_close_solution_file : cg_close !" << finl, TRUST_CGNS_ERROR();

      if (is_cerr) Cerr << "**** CGNS file " << fn << " closed !" << finl;
    }
}

void Ecrire_CGNS::cgns_write_final_link_file()
{
  const bool mult_loc = (static_cast<int>(fld_loc_map_.size()) > 1);

  for (auto itr = fld_loc_map_.begin(); itr != fld_loc_map_.end(); ++itr)
    {
      const int ind = static_cast<int>(std::distance(fld_loc_map_.begin(), itr));

      // XXX a pas oublier, dernier sol fichier ... faut le fermer
      cgns_close_grid_solution_link_file(ind, baseFile_name_);

      // Fichier link maintenant
      const std::string& LOC = itr->first;
      cgns_open_solution_link_file(ind, LOC, -123., true /* dernier fichier => link */);

      const True_int fileId = (ind == 0 ? fileId_ : fileId2_);

      // link solutions
      for (auto& itr_t : time_post_)
        {
          std::string solname = "FlowSolution" + cgns_helper_.convert_double_to_string(itr_t) + "_" + LOC;

          std::string linkfile = baseFile_name_ + "_" + LOC + ".solution." + cgns_helper_.convert_double_to_string(itr_t) + ".cgns"; // file name
          linkfile = TRUST_2_CGNS::remove_slash_linkfile(linkfile);

          std::string linkpath = "/" + baseZone_name_ + "/" + baseZone_name_ + "/" + solname + "/";

          if (cg_link_write(solname.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_final_link_file : cg_link_write !" << finl, TRUST_CGNS_ERROR();
        }

      cgns_helper_.cgns_write_iters<TYPE_ECRITURE::SEQ>(true /* has_field */, 1 /* nb_zones_to_write */, fileId, baseId_[0], 0 /* 1st Zone */, zoneId_, LOC, solname_som_, solname_elem_, time_post_);

      cgns_close_grid_solution_link_file(ind, !mult_loc ? baseFile_name_ + ".cgns" : baseFile_name_ + "_" + LOC + ".cgns", true); // on ferme
    }
}

void Ecrire_CGNS::cgns_write_link_file_for_multiple_files()
{
  if (postraiter_domaine_) return; /* Do nothing */

  Process::barrier();
  if (Process::je_suis_maitre()) // Only master proc writes !
    {
      std::string fn = baseFile_name_ + ".cgns"; // file name
      Cerr << "Option_CGNS::MULTIPLE_FILES is used ... so we write a unique link file " << fn << " ..." << finl;

      True_int fileId_l = -123, baseId_l = -123, zoneId_l = -123;
      True_int fileId = -123, baseId = 1, zoneId = 1, cell_dim = -123, phys_dim = -123;
      True_int nbndry = -123, iparent_flag = -123, nsols = -123, nsections = -123;
      char basename[CGNS_STR_SIZE], zonename[CGNS_STR_SIZE], sectionname[CGNS_STR_SIZE], solname[CGNS_STR_SIZE];

      cgsize_t isize[3][1], istart, iend;
      std::vector<std::string> connectname, sols;
      CGNS_TYPE itype;
      CGNS_LOC loc;

      /* Step 1 : on ouvre baseFile_name_0000.cgns et on lit */
      fn = (Nom(baseFile_name_)).nom_me(0).getString() + ".cgns"; // file name

      if (cg_open(fn.c_str(), CG_MODE_READ, &fileId) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_open !" << finl, TRUST_CGNS_ERROR();

      if (cg_base_read(fileId, baseId, basename, &cell_dim, &phys_dim) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_base_read !" << finl, TRUST_CGNS_ERROR();

      if (cg_zone_read(fileId, baseId, zoneId, zonename, isize[0]) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_zone_read !" << finl, TRUST_CGNS_ERROR();

      if (cg_nsections(fileId, baseId, zoneId, &nsections) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_nsections !" << finl, TRUST_CGNS_ERROR();

      for (int index_sect = 1; index_sect <= nsections; index_sect++)
        {
          if (cg_section_read(fileId, baseId, zoneId, index_sect, sectionname, &itype, &istart, &iend, &nbndry, &iparent_flag) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_section_read !" << finl, TRUST_CGNS_ERROR();

          connectname.push_back(std::string(sectionname));
        }

      if (cg_nsols(fileId, baseId, zoneId, &nsols) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_nsols !" << finl, TRUST_CGNS_ERROR();

      for (int i = 1; i <= nsols; i++)
        {
          if (cg_sol_info(fileId, baseId, zoneId, i, solname, &loc) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_sol_info !" << finl, TRUST_CGNS_ERROR();

          sols.push_back(std::string(solname));
        }

      if (cg_close(fileId) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_close !" << finl, TRUST_CGNS_ERROR();

      /* Step 2 : on ouvre le link file, on laisse ouvert et on ecrit dans la premiere noeud ... */
      fn = baseFile_name_ + ".cgns"; // file name
      unlink(fn.c_str());

      if (cg_open(fn.c_str(), CG_MODE_WRITE, &fileId_l) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_open !" << finl, TRUST_CGNS_ERROR();

      Cerr << "**** CGNS file " << fn << " opened !" << finl;

      if (cg_base_write(fileId_l, basename, cell_dim, phys_dim, &baseId_l) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_base_write !" << finl, TRUST_CGNS_ERROR();

      if (cg_biter_write(fileId_l, baseId_l, "TimeIterValues", static_cast<int>(sols.size())) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_link_file : cg_biter_write !" << finl, TRUST_CGNS_ERROR();

      if (cg_goto(fileId_l, baseId_l, "BaseIterativeData_t", 1, "end") != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_link_file : cg_goto !" << finl, TRUST_CGNS_ERROR();

      std::string linkfile = (Nom(baseFile_name_)).nom_me(0).getString() + ".cgns"; // file name
      linkfile = TRUST_2_CGNS::remove_slash_linkfile(linkfile);

      std::string linkpath = "/" + std::string(basename) + "/TimeIterValues/TimeValues/";

      if (cg_link_write("TimeValues", linkfile.c_str(), linkpath.c_str()) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_link_write !" << finl, TRUST_CGNS_ERROR();

      if (cg_simulation_type_write(fileId_l, baseId_l, CGNS_ENUMV(TimeAccurate)) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_link_file : cg_simulation_type_write !" << finl, TRUST_CGNS_ERROR();

      /* Step 3 : On lit tous les fichiers pour recuperer sizes, et on ecrit les links */
      for (int proc = 0; proc < Process::nproc(); proc++)
        {
          /* on lit les cgns files */
          fn = (Nom(baseFile_name_)).nom_me(proc).getString() + ".cgns"; // file name

          if (cg_open(fn.c_str(), CG_MODE_READ, &fileId) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_open !" << finl, TRUST_CGNS_ERROR();

          if (cg_zone_read(fileId, baseId, zoneId, zonename, isize[0]) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_zone_read !" << finl, TRUST_CGNS_ERROR();

          if (cg_close(fileId) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_close !" << finl, TRUST_CGNS_ERROR();

          /* on ecrit dans le link file : 2eme noeud sous domaine */
          Nom new_zonename = Nom(zonename).nom_me(proc);
          linkfile = fn; // file name

          if (cg_zone_write(fileId_l, baseId_l, new_zonename.getChar(), isize[0], CGNS_ENUMV(Unstructured), &zoneId_l) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_base_write !" << finl, TRUST_CGNS_ERROR();

          linkpath = "/" + std::string(basename) + "/" + std::string(zonename) + "/GridCoordinates/";

          if (cg_goto(fileId_l, baseId_l, "Zone_t", proc + 1, "end") != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_goto !" << finl, TRUST_CGNS_ERROR();

          if (cg_link_write("GridCoordinates", linkfile.c_str(), linkpath.c_str()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_link_write !" << finl, TRUST_CGNS_ERROR();

          for (auto &itr : connectname)
            {
              linkpath = "/" + std::string(basename) + "/" + std::string(zonename) + "/" + itr + "/";
              if (cg_link_write(itr.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_link_write !" << finl, TRUST_CGNS_ERROR();
            }

          // link solutions
          for (auto &itr : sols)
            {
              linkpath = "/" + std::string(basename) + "/" + std::string(zonename) + "/" + itr + "/";
              if (cg_link_write(itr.c_str(), linkfile.c_str(), linkpath.c_str()) != CG_OK)
                Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_link_write !" << finl, TRUST_CGNS_ERROR();
            }

          if (cg_ziter_write(fileId_l, baseId_l, zoneId_l, "ZoneIterativeData") != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_ziter_write !" << finl, TRUST_CGNS_ERROR();

          if (cg_goto(fileId_l, baseId_l, "Zone_t", zoneId_l, "ZoneIterativeData_t", 1, "end") != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_goto !" << finl, TRUST_CGNS_ERROR();

          linkpath = "/" + std::string(basename) + "/" + std::string(zonename) + "/ZoneIterativeData/FlowSolutionPointers/";
          if (cg_link_write("FlowSolutionPointers", linkfile.c_str(), linkpath.c_str()) != CG_OK)
            Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_link_write !" << finl, TRUST_CGNS_ERROR();
        }

      /* Step 4 : on ferme le link file, on laisse ouvert et on ecrit la base */
      if (cg_close(fileId_l) != CG_OK)
        Cerr << "Error Ecrire_CGNS::cgns_write_link_file_for_multiple_files : cg_close !" << finl, TRUST_CGNS_ERROR();

      Cerr << "**** CGNS file " << baseFile_name_ + ".cgns" << " closed !" << finl;
    }
}

#endif /* HAS_CGNS */
