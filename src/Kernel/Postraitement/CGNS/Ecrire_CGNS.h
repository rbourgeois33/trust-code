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

#ifndef Ecrire_CGNS_included
#define Ecrire_CGNS_included

#include <Domaine_forward.h>
#include <TRUST_2_CGNS.h>
#include <map>

class Nom;

class Ecrire_CGNS
{
#ifdef HAS_CGNS
public:
  void cgns_associer_domaine_dis(const Domaine_dis_base& );
  void cgns_init_MPI(bool is_self = false);
  void cgns_set_postraiter_domain() { postraiter_domaine_ = true; }
  void cgns_set_is_dual_domain() { is_dual_ = true; }
  void cgns_set_is_deformable_domain() { is_deformable_ = true; }
  void cgns_set_loc_vector(const std::vector<std::string>& vec) { loc_vect_ = vec; }
  void cgns_set_base_name(const Nom& );
  void cgns_open_file();
  void cgns_finir();
  void cgns_add_time(const double );
  void cgns_write_domaine_dual(const Domaine& , const int , const Nom& nom_dom = "??");
  void cgns_write_domaine(const Domaine * ,const Nom& , const DoubleTab& , const IntTab& , const Motcle& );
  void cgns_write_field(const Domaine&, const Noms&, double, const Nom&, const Nom&, const Nom&, const DoubleTab&);
  void finir_ecriture(double);

private:

  void fill_infos_loc();
  Ecrire_CGNS_helper cgns_helper_;
  OBS_PTR(Domaine_dis_base) domaine_dis_;
  OBS_PTR(std::vector<std::string>) loc_vect_;
  std::vector<TRUST_2_CGNS> T2CGNS_;

  std::map<std::string, Nom> fld_loc_map_; /* { Loc , Nom_dom } */
  std::vector<Nom> doms_written_;
  std::vector<Nom> fieldName_dumped_; /* filled just once to see what fields are already written ! */
  std::string solname_elem_ = "", solname_som_ = "", solname_faces_ = "", baseFile_name_ = "";
  std::vector<double> time_post_;
  std::vector<int> baseId_, zoneId_;

  bool has_elem_field_ = false, has_faces_field_ = false, has_som_field_ = false, has_elem_som_loc_ = false;
  bool solname_elem_written_ = false, solname_som_written_ = false, solname_faces_written_ = false;
  bool postraiter_domaine_ = false;
  bool first_time_post_ = true;

  int fileId_ = -123;
  int flowId_elem_ = 0, fieldId_elem_ = 0;
  int flowId_som_ = 0, fieldId_som_ = 0;
  int flowId_faces_ = 0, fieldId_faces_ = 0;

  // specifique pour link
  bool grid_file_opened_ = false, solution_file_opened_ = false; /* Management of link files */
  std::vector<std::string> baseZone_name_;
  std::vector<std::vector<std::string>> connectname_;
  std::vector<std::vector<cgsize_t>> sizeId_;
  std::vector<int> cellDim_;

  // specifique par in zone
  std::vector<std::vector<int>> zoneId_par_; /* par ordre d'ecriture du domaine */

  // specifique FILE_PER_COMM_GROUP
  int proc_maitre_local_comm_ = -123;
  std::vector<int> vec_proc_maitre_local_comm_, unique_vec_proc_maitre_local_comm_;
  std::vector<std::vector<cgsize_t>> sizeId_som_local_comm_, sizeId_elem_local_comm_; // Attention : ind 0 => ELEM et SOM, ind 1 => FACES !

  // specifique maillage dual pour faces
  IntTab fs_dual_, ef_dual_;
  bool is_dual_ = false;
  inline IntTab& get_fs_dual() { return fs_dual_;}
  inline const IntTab& get_fs_dual() const { return fs_dual_;}
  inline IntTab& get_ef_dual() { return ef_dual_;}
  inline const IntTab& get_ef_dual() const { return ef_dual_;}

  // gestion elem/som/faces
  void cgns_fill_field_loc_map(const Domaine&, const std::string&);

  // Methodes pour Domaine Deformable
  void cgns_write_final_link_file_pb_deformable();
  void link_multi_loc_support_pb_deformable();
  bool is_deformable_ = false, multi_loc_deformable_support_linked_ = false;

  // Methodes pour LINK
  void cgns_fill_info_grid_link_file(const char*, const CGNS_TYPE&, const int, const int, const int, const bool);
  void cgns_open_grid_base_link_file();
  void cgns_init_solution_link_file(const std::string& , const Nom&);
  void cgns_open_solution_link_file(const double, bool is_link = false);
  void cgns_write_final_link_file();
  void cgns_write_final_link_file_comm_group();
  void cgns_write_link_file_for_multiple_files();
  void cgns_close_grid_or_solution_link_file(const double, const TYPE_LINK_CGNS,  bool is_cerr = true);
  void gather_local_sizeId_multi_loc(std::vector<std::vector<cgsize_t>>& , std::vector<std::vector<cgsize_t>>& ) const ;
  void add_new_linked_base(const std::string&, const Nom&);
  void add_new_linked_base_par_over_zone(const std::string&, const Nom&, const Nom&, const int);

  // Version sequentielle
  void cgns_write_domaine_seq(const Domaine * ,const Nom& , const DoubleTab& , const IntTab& , const Motcle& );
  void cgns_write_domaine_deformable_seq(const Domaine * ,const Nom& , const DoubleTab& , const IntTab& , const Motcle& );
  void cgns_write_field_seq(const int, const double, const Nom&, const Nom&, const Nom&, const Nom&, const DoubleTab&);
  void cgns_write_iters_seq();

  // Version parallele over zone
  void cgns_write_domaine_par_over_zone(const Domaine * ,const Nom& , const DoubleTab& , const IntTab& , const Motcle& );
  void cgns_write_field_par_over_zone(const int, const double, const Nom&, const Nom&, const Nom&, const Nom&, const DoubleTab&);
  void cgns_write_iters_par_over_zone();

  // Version parallele in zone
  void cgns_write_domaine_par_in_zone(const Domaine * ,const Nom& , const DoubleTab& , const IntTab& , const Motcle& );
  void cgns_write_domaine_deformable_par_in_zone(const Domaine * ,const Nom& , const DoubleTab& , const IntTab& , const Motcle& );
  void cgns_write_field_par_in_zone(const int, const double, const Nom&, const Nom&, const Nom&, const Nom&, const DoubleTab&);
  void cgns_write_iters_par_in_zone();

#endif /* HAS_CGNS */
};

inline void verify_if_cgns(const char * nom_funct)
{
#ifdef HAS_CGNS
#ifdef MPI_
  return;
#else /* NOT MPI_ */
  if (Process::is_parallel())
    Process::exit("Parallel CGNS files need MPI installed ... ");
#endif /* MPI_ */
#else /* NOT HAS_CGNS */
  Cerr << "Format_Post_CGNS::" <<  nom_funct << " should not be called since TRUST is not compiled with the CGNS library !!! " << finl;
  Process::exit();
#endif /* HAS_CGNS */
}

#endif /* Ecrire_CGNS_included */
