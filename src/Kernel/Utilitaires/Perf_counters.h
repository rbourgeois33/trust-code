/****************************************************************************
 * Copyright (c) 2024, CEA
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

#include <string>
#include <vector>
#include <map>
#include <tuple>
#include <chrono>
#include <array>
// This file contains all of the needed for the description of the counter associated with the tracking of performance in the TRUST code.
// A class Time is introduced in Perf_counters.cpp for extracting the time
class Counter;

enum class STD_COUNTERS : unsigned int
{ temps_total_execution_counter_ ,
  initialisation_calcul_counter_ ,
  timestep_counter_ ,
  solv_sys_counter_,
  solv_sys_petsc_counter_,
  diffusion_implicite_counter_,
  dt_counter_ ,
  nut_counter_ ,
  convection_counter_ ,
  diffusion_counter_ ,
  decay_counter_ ,
  gradient_counter_ ,
  divergence_counter_ ,
  source_counter_ ,
  postraitement_counter_ ,
  sauvegarde_counter_ ,
  temporary_counter_ ,
  assemblage_sys_counter_ ,
  update_vars_counter_ ,
  update_fields_counter_ ,
  mettre_a_jour_counter_  ,
  divers_counter_ ,
  probleme_fluide_ ,
  probleme_combustible_ ,
  echange_vect_counter_ ,
  mpi_sendrecv_counter_  ,
  mpi_send_counter_ ,
  mpi_recv_counter_ ,
  mpi_bcast_counter_ ,
  mpi_alltoall_counter_ ,
  mpi_allgather_counter_ ,
  mpi_gather_counter_ ,
  mpi_partialsum_counter_,
  mpi_sumdouble_counter_ ,
  mpi_mindouble_counter_,
  mpi_maxdouble_counter_ ,
  mpi_sumfloat_counter_ ,
  mpi_minfloat_counter_ ,
  mpi_maxfloat_counter_ ,
  mpi_sumint_counter_ ,
  mpi_minint_counter_ ,
  mpi_maxint_counter_ ,
  mpi_barrier_counter_ ,
  gpu_library_counter_ ,
  gpu_kernel_counter_ ,
  gpu_copytodevice_counter_ ,
  gpu_copyfromdevice_counter_ ,
  IO_EcrireFicPartageMPIIO_counter_ ,
  IO_EcrireFicPartageBin_counter_ ,
  interprete_scatter_counter_
};


/*!  @brief This class store counters in TRUST
 *
 */
class Perf_counters
{
public:


  //	static Perf_counters & GetInstance() {
  //		static Perf_counters * inst = new Perf_counters(); // initialized only once!
  //		return inst;
  //	}

  /*! @brief Create the basic counters of TRUST which are always used
   *
   * Should always be called at the start of any case on every processors simultaneously
   */
  void declare_base_counters();

  /*! @brief Create a new counter and add it to the vector of counters
   *
   * @param counter_level
   * @param counter_description
   * @param counter_family
   * @param is_comm
   * @return create a new counter
   */
  void create_custom_counter(int counter_level, std::string counter_description, std::string counter_family ,bool is_comm);


  /*! @brief Sort the vector of counters based on counter level
   *
   */
  void sort_counters_by_level();

  /*! @brief Start the count of a counter
   *
   */
  //	void begin_count(Counter c);


  /*! Standard counters
   *
   * @param
   */
  void begin_count(const STD_COUNTERS std_cnt);

  /*! Custom counters
   *
   * @param custom_count_name
   */
  void begin_count(const std::string& custom_count_name);


  /*! @brief End the count of a counter and update the counter values
   *
   * @param c is the counter to end the count
   * @param count_increment is the count increment. If not specified, then it is equal to 1
   * @param quantity_increment is the increment of custom variable quantity. If not specified, it is set to 0.
   */
  void end_count(const std::string& custom_count_name, int count_increment=1, double quantity_increment=0);

  /*! @brief End the count of a counter and update the counter values
   *
   * @param c is the counter to end the count
   * @param count_increment is the count increment. If not specified, then it is equal to 1
   * @param quantity_increment is the increment of custom variable quantity. If not specified, it is set to 0.
   */
  void end_count(const STD_COUNTERS std_cnt, int count_increment, double quantity_increment);

  /*! @brief Stop all counters, has to be called on every processor simultaneously
   *
   */
  void stop_counters();

  /*! @brief Restart all counters, has to be called on every processor simultaneously
   *
   */
  void restart_counters();

  /*!
   *
   * @param is_the_three_first_time_steps_elapsed == True, then the thre first time steps are discarded
   */
  inline void set_three_first_steps_elapsed(bool is_the_three_first_time_steps_elapsed)
  {
    three_first_steps_elapsed_ = is_the_three_first_time_steps_elapsed;
  }

  /*! @brief Create the csv.TU file.
   *
   * Some local sub-functions are defined in Perf_counters.cpp for constructing the csv.TU_file
   */
  void print_performance_to_csv(std::string message, bool mode_append);

  /*! @brief Create the global .Tu file with agglomerated stats
   *
   * Some local sub-functions are defined in Perf_counters.cpp for constructing the csv.TU_file
   */
  void print_global_TU(std::string message, bool mode_append);

  /*! Change the maximum counter level to be print in the Tu files
   *
   * @param new_max_counter_lvl_to_print
   */
  void set_max_counter_lvl_to_print(int new_max_counter_lvl_to_print);

  void compute_avg_min_max_var_per_step(int tstep);

  std::string get_os();

  std::string get_cpu();

  std::string get_gpu();

  std::string get_date();

private:
  Perf_counters();
  ~Perf_counters();
  bool three_first_steps_elapsed_;
  int max_counter_lvl_to_print;
  /* std::array <Counter, static_cast<int>(STD_COUNTERS::LAST_COUNT)-1> std_counters_ ;
  std::map <std::string, Counter> custom_counters_ ; */

  std::array <Counter,static_cast<int>(STD_COUNTERS::interprete_scatter_counter_)-1> std_counters_ ; // array of the standard counters of TRUST, always used in practice
  std::map <std::string, Counter> custom_counter_map_str_to_counter_ ; // Link the custom counters descriptions to the counter type

};


