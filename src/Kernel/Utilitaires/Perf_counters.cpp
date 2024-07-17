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
#include <array>
#include <map>
#include <tuple>
#include <chrono>
#include "Perf_counters.h"

/**************************************************************************************************************************
 *
 * 		Introduction of the class counter that described the behavior of a single counter in TRUST
 *
 **************************************************************************************************************************/
class Counter
{
public:
  friend class Perf_counters;
  Counter();

  Counter(int counter_level, std::string counter_name, std::string counter_family = "None", bool is_comm = false);

  inline std::chrono::time_point<std::chrono::high_resolution_clock> get_time_now()
  {
    return std::chrono::high_resolution_clock::now();
  }

  inline bool is_a_comm_counter()
  {
    return is_comm_;
  }

  inline void begin_count_()
  {
    last_open_time_ = Counter::get_time_now();
    is_running_ = true;
  }

  void end_count_(int count_increment, double quantity_increment);
  //check std::chrono

  /*! @brief update variables : avg_time_per_step_ , min_time_per_step_ , max_time_per_step_ , sd_time_per_step_
   *
   */
  void compute_avg_min_max_var_per_step();

protected:
  std::string description_;
  int level_;
  std::string family_ ;
  bool is_comm_;
  bool is_running_;
  int count_;
  double quantity_;
  std::chrono::duration<double> total_time_;
  std::chrono::time_point<std::chrono::high_resolution_clock> last_open_time_;
  double avg_time_per_step_;
  double min_time_per_step_;
  double max_time_per_step_;
  double sd_time_per_step_;
}
;

Counter::Counter(int counter_level, std::string counter_name, std::string counter_family , bool is_comm )
{
  description_ = counter_name;
  level_ = counter_level;
  family_ = counter_family;
  is_comm_ = is_comm;
  is_running_ = false;
  count_ = 0;
  quantity_ = 0.0;
  avg_time_per_step_ = 0.0 ;
  min_time_per_step_ = 0.0 ;
  max_time_per_step_ = 0.0 ;
  sd_time_per_step_ = 0.0 ;
  total_time_ = std::chrono::duration<double>::zero() ;
  last_open_time_ = std::chrono::high_resolution_clock::now();
}

Counter::Counter()
{
  description_ = "";
  level_ = -1;
  family_ = "";
  is_comm_ = false;
  is_running_ = false;
  count_ = 0;
  quantity_ = 0.0;
  avg_time_per_step_ = 0.0 ;
  min_time_per_step_ = 0.0 ;
  max_time_per_step_ = 0.0 ;
  sd_time_per_step_ = 0.0 ;
  total_time_ = std::chrono::duration<double>::zero() ;
  last_open_time_ = std::chrono::high_resolution_clock::now();
}

void Counter::end_count_(int count_increment, double quantity_increment)
{
  assert(is_running_);
  std::chrono::time_point <std::chrono::high_resolution_clock> time_stop_counter = Counter::get_time_now();
  total_time_ += time_stop_counter-last_open_time_;
}

/**************************************************************************************************************************
 *
 * 					Declaration of the class Perf_counters
 *
 **************************************************************************************************************************/
/*! @Brief declare all standard counters of TRUST inside an array
 *
 */

void Perf_counters::declare_base_counters()
{
  // Macro counters
  std_counters_[STD_COUNTERS::temps_total_execution_counter_]= Counter::Counter(1, "Total time");
  std_counters_[STD_COUNTERS::initialisation_calcul_counter_] = Counter::Counter(1, "Prepare computation");
  std_counters_[STD_COUNTERS::timestep_counter_] = Counter::Counter(1, "Number of time steps");
  std_counters_[STD_COUNTERS::solv_sys_counter_] = Counter::Counter(1, "System's solver");
  std_counters_[STD_COUNTERS::solv_sys_petsc_counter_] = Counter::Counter(1, "Petsc solver");
  std_counters_[STD_COUNTERS::diffusion_implicite_counter_] = Counter::Counter(1, "Basic_equation::conjugated gradient implicit diffusion");
  std_counters_[STD_COUNTERS::dt_counter_] = Counter::Counter(1, "Computation of the time step");
  std_counters_[STD_COUNTERS::nut_counter_] = Counter::Counter(1, "Turbulence model::update");
  std_counters_[STD_COUNTERS::convection_counter_] = Counter::Counter(1, "Convection operator::add/compute");
  std_counters_[STD_COUNTERS::diffusion_counter_] = Counter::Counter(1, "Diffusion operator::add/compute");
  std_counters_[STD_COUNTERS::decay_counter_] = Counter::Counter(1, "Decreasing operator::add/compute");
  std_counters_[STD_COUNTERS::gradient_counter_] = Counter::Counter(1, "Gradient operator::add/compute");
  std_counters_[STD_COUNTERS::divergence_counter_] = Counter::Counter(1, "Divergence operator::add/compute");
  std_counters_[STD_COUNTERS::source_counter_] = Counter::Counter(1, "Source_terms::add/compute");
  std_counters_[STD_COUNTERS::postraitement_counter_] = Counter::Counter(1, "Base problem::post-treatment");
  std_counters_[STD_COUNTERS::sauvegarde_counter_] = Counter::Counter(1, "Base problem::save");
  std_counters_[STD_COUNTERS::assemblage_sys_counter_] = Counter::Counter(1, "Assemble matrix");
  // Pb_multiphase counters
  std_counters_[STD_COUNTERS::update_vars_counter_] = Counter::Counter(1, "Implicit scheme 4eqs::update variables");
  std_counters_[STD_COUNTERS::update_fields_counter_] = Counter::Counter(1, "Base two-phase flow problem::updateGivenFields");
  std_counters_[STD_COUNTERS::mettre_a_jour_counter_] = Counter::Counter(1, "Update");
  // MPI communication counters
  std_counters_[STD_COUNTERS::mpi_sendrecv_counter_]  = Counter::Counter(2, "MPI_send_recv", "MPI_sendrecv", true);
  std_counters_[STD_COUNTERS::mpi_send_counter_]     = Counter::Counter(2, "MPI_send",      "MPI_sendrecv", true);
  std_counters_[STD_COUNTERS::mpi_recv_counter_]      = Counter::Counter(2, "MPI_recv",      "MPI_sendrecv", true);
  std_counters_[STD_COUNTERS::mpi_bcast_counter_]    = Counter::Counter(2, "MPI_broadcast", "MPI_sendrecv", true);
  std_counters_[STD_COUNTERS::mpi_alltoall_counter_]  = Counter::Counter(2, "MPI_alltoall",  "MPI_sendrecv", true);
  std_counters_[STD_COUNTERS::mpi_allgather_counter_] = Counter::Counter(2, "MPI_allgather", "MPI_sendrecv", true);
  std_counters_[STD_COUNTERS::mpi_gather_counter_]    = Counter::Counter(2, "MPI_gather", "MPI_sendrecv", true);
  std_counters_[STD_COUNTERS::mpi_partialsum_counter_]= Counter::Counter(2, "MPI_partialsum","MPI_allreduce", true);
  std_counters_[STD_COUNTERS::mpi_sumdouble_counter_] = Counter::Counter(2, "MPI_sumdouble", "MPI_allreduce", true);
  std_counters_[STD_COUNTERS::mpi_mindouble_counter_] = Counter::Counter(2, "MPI_mindouble", "MPI_allreduce", true);
  std_counters_[STD_COUNTERS::mpi_maxdouble_counter_] = Counter::Counter(2, "MPI_maxdouble", "MPI_allreduce", true);
  std_counters_[STD_COUNTERS::mpi_sumfloat_counter_]  = Counter::Counter(2, "MPI_sumfloat", "MPI_allreduce", true);
  std_counters_[STD_COUNTERS::mpi_minfloat_counter_]  = Counter::Counter(2, "MPI_minfloat", "MPI_allreduce", true);
  std_counters_[STD_COUNTERS::mpi_maxfloat_counter_]  = Counter::Counter(2, "MPI_maxfloat", "MPI_allreduce", true);
  std_counters_[STD_COUNTERS::mpi_sumint_counter_]    = Counter::Counter(2, "MPI_sumint",    "MPI_allreduce", true);
  std_counters_[STD_COUNTERS::mpi_minint_counter_]    = Counter::Counter(2, "MPI_minint",    "MPI_allreduce", true);
  std_counters_[STD_COUNTERS::mpi_maxint_counter_]    = Counter::Counter(2, "MPI_maxint",    "MPI_allreduce", true);
  std_counters_[STD_COUNTERS::mpi_barrier_counter_]   = Counter::Counter(2, "MPI_barrier",   "MPI_allreduce", true);
  // GPU
  std_counters_[STD_COUNTERS::gpu_library_counter_]        = Counter::Counter(2, "GPU_library",       "GPU_library");
  std_counters_[STD_COUNTERS::gpu_kernel_counter_]         = Counter::Counter(2, "GPU_kernel",        "GPU_kernel");
  std_counters_[STD_COUNTERS::gpu_copytodevice_counter_  ] = Counter::Counter(2, "GPU_copyToDevice",  "GPU_copy");
  std_counters_[STD_COUNTERS::gpu_copyfromdevice_counter_] = Counter::Counter(2, "GPU_copyFromDevice","GPU_copy");
  // Count the writing time in EcrireFicPartageXXX (big chunk of data in files XYZ or LATA)
  // For those two counters, quantity corresponds to the amount of written bytes
  std_counters_[STD_COUNTERS::IO_EcrireFicPartageMPIIO_counter_] = Counter::Counter(2, "MPI_File_write_all", "IO"); // Call on each processor
  std_counters_[STD_COUNTERS::IO_EcrireFicPartageBin_counter_] = Counter::Counter(2, "write", "IO"); // Call on master
  // Execution of Scatter::interpreter
  std_counters_[STD_COUNTERS::interprete_scatter_counter_] = Counter::Counter(2, "Scatter");
  // Problemes
  std_counters_[STD_COUNTERS::probleme_fluide_] = Counter::Counter(3, "pb_fluide");
  std_counters_[STD_COUNTERS::probleme_combustible_] = Counter::Counter(3, "pb_combustible");

  // Appels a DoubleVect::echange_espace_virtuel()
  std_counters_[STD_COUNTERS::echange_vect_counter_] = Counter::Counter(2, "DoubleVect/IntVect::echange_espace_virtuel", "None", true);

}

void Perf_counters::create_custom_counter(int counter_level, std::string counter_description, std::string counter_family ,bool is_comm)
{
  Counter new_counter = Counter::Counter(counter_level, counter_description, counter_family ,is_comm);
  custom_counter_map_str_to_counter_.insert({counter_description,new_counter});
}


/*! Standard counters
 *
 * @param
 */
void Perf_counters::begin_count(STD_COUNTERS std_cnt)
{
  std_counters_[std_cnt].begin_count();
}

/*! Custom counters
 *
 * @param custom_count_name
 */
void Perf_counters::begin_count(const std::string& custom_count_name)
{
  custom_counter_map_str_to_counter_[custom_count_name].begin_count();
}

/*! @brief End the count of a counter and update the counter values
 *
 * @param c is the counter to end the count
 * @param count_increment is the count increment. If not specified, then it is equal to 1
 * @param quantity_increment is the increment of custom variable quantity. If not specified, it is set to 0.
 */
void Perf_counters::end_count(const STD_COUNTERS std_cnt, int count_increment, double quantity_increment)
{
  assert(std_counters_[std_cnt].is_running_);
  //if (three_first_steps_elapsed_)
  std_counters_[std_cnt].end_count_(count_increment, quantity_increment);
}


/*! @brief End the count of a counter and update the counter values
 *
 * @param c is the counter to end the count
 * @param count_increment is the count increment. If not specified, then it is equal to 1
 * @param quantity_increment is the increment of custom variable quantity. If not specified, it is set to 0.
 */
void Perf_counters::end_count(const std::string& custom_count_name, int count_increment, double quantity_increment)
{
  assert(custom_counter_map_str_to_counter_[custom_count_name].is_running_);
  custom_counter_map_str_to_counter_[custom_count_name].end_count_(count_increment,quantity_increment);
}





