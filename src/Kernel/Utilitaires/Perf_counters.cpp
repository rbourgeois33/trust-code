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
#include <stdio.h>
#include <algorithm>
#include <string.h>
#include <vector>
#include <array>
#include <map>
#include <tuple>
#include <chrono>
#include <assert.h>
#include <sys/utsname.h>
#include <math.h>
#include <cmath>
#include <iostream>
#include <sstream>
#include <fstream>
#include <EntreeSortie.h>
#include <Perf_counters.h>
#include <iomanip>
#include <EcrFicPartage.h>
#include <SChaine.h>
#include <communications.h>
#include <TRUSTArray.h>
#include <Comm_Group_MPI.h>
#include <memory>
#include <iomanip>

/**************************************************************************************************************************
 *
 *      Introduction of the class counter that described the behavior of a single counter in TRUST
 *
 **************************************************************************************************************************/
class Counter
{
public:
  friend class Perf_counters;
  Counter(bool to_print_in_global_TU, int counter_level, std::string counter_name, std::string counter_family = "None", bool is_comm = false);

  void begin_count_(int counter_level, std::chrono::time_point<std::chrono::high_resolution_clock> t);

  void end_count_(int count_increment, int quantity_increment, std::chrono::time_point<std::chrono::high_resolution_clock> t_stop);
  //check std::chrono

  inline void set_parent(Counter * parent_counter) { parent_ = parent_counter;}

  inline void set_to_print(bool to_print) {to_print_in_global_TU_ = to_print;}

  inline double get_time_() {return total_time_.count();}

  /*! @brief update variables : avg_time_per_step_ , min_time_per_step_ , max_time_per_step_ , sd_time_per_step_
   *
   */
  void compute_avg_min_max_var_per_step();

  std::array< std::array<double,4> ,4> compute_min_max_avg_sd_();

  void reset();

protected:
  const std::string description_;
  int level_;
  const std::string family_ ;
  const bool is_comm_;
  int count_;
  int quantity_;
  Counter* parent_;
  std::chrono::duration<double> total_time_;
  std::chrono::duration<double> time_alone_;   // time when the counter is open minus the time where an counter of lower lvl was open
  std::chrono::duration<double> time_ts_;  // total time tracked during the current time_steps
  std::chrono::time_point<std::chrono::high_resolution_clock> open_time_ts_;
  std::chrono::time_point<std::chrono::high_resolution_clock> last_open_time_;
  std::chrono::time_point<std::chrono::high_resolution_clock> last_open_time_alone_;
  double avg_time_per_step_;
  double min_time_per_step_;
  double max_time_per_step_;
  double sd_time_per_step_;
  bool to_print_in_global_TU_;
  bool is_running_;
}
;

Counter::Counter(bool to_print_in_global_TU, int counter_level, std::string counter_name, std::string counter_family , bool is_comm)
  :description_(counter_name), family_(counter_family), is_comm_(is_comm)
{
  level_ = counter_level;
  count_ = 0;
  quantity_ = 0;
  avg_time_per_step_ = 0.0 ;
  min_time_per_step_ = 0.0 ;
  max_time_per_step_ = 0.0 ;
  sd_time_per_step_ = 0.0 ;
  open_time_ts_ = std::chrono::time_point<std::chrono::high_resolution_clock>();
  total_time_ = std::chrono::duration<double>::zero() ;
  time_alone_ = std::chrono::duration<double>::zero() ;
  time_ts_ = std::chrono::duration<double>::zero() ;
  last_open_time_alone_= std::chrono::time_point<std::chrono::high_resolution_clock>();
  last_open_time_ = std::chrono::time_point<std::chrono::high_resolution_clock>();
  parent_ = nullptr;
  to_print_in_global_TU_=to_print_in_global_TU;
  is_running_ = false;
}

void Counter::begin_count_(int counter_level, std::chrono::time_point<std::chrono::high_resolution_clock> t)
{
  if (counter_level != level_)
    {
      level_ = counter_level; ///< You have changed the level of your counter
    }
  is_running_ = true;
  last_open_time_ = t;
  last_open_time_alone_ = t;
  if (parent_!= nullptr)
    {
      parent_->time_alone_ += std::chrono::duration<double> (t - last_open_time_alone_);
      parent_->last_open_time_alone_ =  std::chrono::time_point<std::chrono::high_resolution_clock>();
    }
}

void Counter::end_count_(int count_increment, int quantity_increment, std::chrono::time_point<std::chrono::high_resolution_clock> t_stop)
{
  if (last_open_time_ == std::chrono::time_point<std::chrono::high_resolution_clock>() )
    Process::exit("Last open_time was not properly set"+ description_);
  if (last_open_time_alone_ == std::chrono::time_point<std::chrono::high_resolution_clock>())
    Process::exit("Last open_time alone was not properly set" + description_);
  std::chrono::duration<double> t_tot = t_stop-last_open_time_;
  std::chrono::duration<double> t_alone = t_stop - last_open_time_alone_;
  quantity_ += quantity_increment;
  total_time_ += t_tot;
  time_alone_ += t_alone;
  count_ += count_increment;
  if (parent_!= nullptr)
    {
      parent_-> last_open_time_alone_ = t_stop;
      parent_ = nullptr;
    }
  is_running_ = false;
  last_open_time_ = std::chrono::time_point<std::chrono::high_resolution_clock>();
  last_open_time_alone_ = std::chrono::time_point<std::chrono::high_resolution_clock>();
  open_time_ts_ = std::chrono::time_point<std::chrono::high_resolution_clock>();
}

std::array< std::array<double,4> ,4> Counter::compute_min_max_avg_sd_()
{
  assert(Process::is_parallel());
  double min,max,avg,sd ;

  min = Process::mp_min(total_time_.count());
  max = Process::mp_max(total_time_.count());
  avg = Process::mp_sum(total_time_.count())/Process::nproc();
  sd = sqrt(Process::mp_sum((total_time_.count()-avg)*(total_time_.count()-avg))/Process::nproc());

  std::array<double,4> min_max_avg_sd_time = {min,max,avg,sd};

  min = (double)Process::mp_min(quantity_);
  max = (double)Process::mp_max(quantity_);
  avg = (double)Process::mp_sum(quantity_)/Process::nproc();
  sd = sqrt((double)Process::mp_sum((quantity_-avg)*(quantity_-avg))/Process::nproc());

  std::array<double,4> min_max_avg_sd_quantity = {min,max,avg,sd};

  min = (double)Process::mp_min(count_);
  max = (double)Process::mp_max(count_);
  avg = (double)Process::mp_sum(count_)/Process::nproc();
  sd = sqrt((double)Process::mp_sum((count_-avg)*(count_-avg))/Process::nproc());

  std::array<double,4> min_max_avg_sd_count_ = {min,max,avg,sd};

  min = (double)Process::mp_min(time_alone_.count());
  max = (double)Process::mp_max(time_alone_.count());
  avg = (double)Process::mp_sum(time_alone_.count())/Process::nproc();
  sd = sqrt((double)Process::mp_sum((time_alone_.count()-avg)*(time_alone_.count()-avg))/Process::nproc());

  std::array<double,4> min_max_avg_sd_time_alone_ = {min,max,avg,sd};

  return {min_max_avg_sd_time,min_max_avg_sd_quantity,min_max_avg_sd_count_,min_max_avg_sd_time_alone_ };
}

void Counter::reset()
{
  count_ = 0;
  quantity_ = 0;
  avg_time_per_step_ = 0.0 ;
  min_time_per_step_ = 0.0 ;
  max_time_per_step_ = 0.0 ;
  sd_time_per_step_ = 0.0 ;
  open_time_ts_ = std::chrono::time_point<std::chrono::high_resolution_clock>();
  total_time_ = std::chrono::duration<double>::zero() ;
  last_open_time_alone_= std::chrono::time_point<std::chrono::high_resolution_clock>();
  last_open_time_ = std::chrono::time_point<std::chrono::high_resolution_clock>();
  time_alone_= std::chrono::duration<double>::zero() ;;   // time when the counter is open minus the time where an counter of lower lvl was open
  time_ts_= std::chrono::duration<double>::zero() ;;  // total time tracked during the current time_steps
}

/**************************************************************************************************************************
 *
 * 					Declaration of the class Perf_counters
 *
 **************************************************************************************************************************/
/*! @Brief declare all standard counters of TRUST inside an array
 *
 */

Perf_counters::Perf_counters()
  : std_counters_
{
  // Macro counters
  new Counter(false,-1, "Total time"),
  new Counter(true,0, "Prepare computation"),
  new Counter(false,0, "Time loop"),
  new Counter(true,1, "Number of linear system resolutions Ax=B"),
  new Counter(true,1, "Petsc solver"),
  new Counter(true,1, "Number of linear system resolutions for implicit diffusion:"),
  new Counter(true, 1, "Computation of the time step dt"),
  new Counter(true, 1, "Turbulence model::update"),
  new Counter(true, 1, "Convection operator::add/compute"),
  new Counter(true, 1, "Diffusion operator::add/compute"),
  new Counter(true, 1, "Gradient operator::add/compute"),
  new Counter(true, 1, "Divergence operator::add/compute"),
  new Counter(true, 1, "Source_terms::add/compute"),
  new Counter(true, 1, "Post-treatment"),
  new Counter(true, 0, "Back-up operations"),
  new Counter(true, 1,"Read file for restart"),
  new Counter(true, 1, "Number of matrix assemblies for the implicit scheme:"),
  new Counter(true, 1, "Update"),
  // MPI communication counters
  new Counter(false, 2, "MPI_send_recv", "MPI_sendrecv", true),
  new Counter(false, 2, "MPI_send",      "MPI_sendrecv", true),
  new Counter(false, 2, "MPI_recv",      "MPI_sendrecv", true),
  new Counter(false, 2, "MPI_broadcast", "MPI_sendrecv", true),
  new Counter(false, 2, "MPI_alltoall",  "MPI_sendrecv", true),
  new Counter(false, 2, "MPI_allgather", "MPI_sendrecv", true),
  new Counter(false, 2, "MPI_gather", "MPI_sendrecv", true),
  new Counter(false, 2, "MPI_partialsum","MPI_allreduce", true),
  new Counter(false, 2, "MPI_sumdouble", "MPI_allreduce", true),
  new Counter(false, 2, "MPI_mindouble", "MPI_allreduce", true),
  new Counter(false, 2, "MPI_maxdouble", "MPI_allreduce", true),
  new Counter(false, 2, "MPI_sumfloat", "MPI_allreduce", true),
  new Counter(false, 2, "MPI_minfloat", "MPI_allreduce", true),
  new Counter(false, 2, "MPI_maxfloat", "MPI_allreduce", true),
  new Counter(false, 2, "MPI_sumint",    "MPI_allreduce", true),
  new Counter(false, 2, "MPI_minint",    "MPI_allreduce", true),
  new Counter(false, 2, "MPI_maxint",    "MPI_allreduce", true),
  new Counter(false, 2, "MPI_barrier",   "MPI_allreduce", true),
  // GPU
  new Counter(true, 2, "GPU_library",       "GPU_library"),
  new Counter(true, 2, "GPU_kernel",        "GPU_kernel"),
  new Counter(true, 2, "GPU_copyToDevice",  "GPU_copy"),
  new Counter(true, 2, "GPU_copyFromDevice","GPU_copy"),
  // Count the writing time in EcrireFicPartageXXX (big chunk of data in files XYZ or LATA)
  // For those two counters, quantity corresponds to the amount of written bytes
  new Counter(true, 1, "MPI_File_write_all", "IO"), // Call on each processor
  new Counter(true, 1, "write", "IO"), // Call on master
  // Execution of Scatter::interpreter
  new Counter(true, 2, "Scatter_interprete"),
  new Counter(true, 2, "DoubleVect/IntVect::virtual_swap", "None", true),
  new Counter(true, 2, "Scatter::read_domaine"),
}
{
  /*
  // Macro counters
  std_counters_[static_cast<int>(STD_COUNTERS::total_execution_time)]= new Counter(false,-1, "Total time");
  std_counters_[static_cast<int>(STD_COUNTERS::computation_start_up)] = new Counter(true,0, "Prepare computation");
  std_counters_[static_cast<int>(STD_COUNTERS::timeloop)] = new Counter(false,0, "Time loop");
  std_counters_[static_cast<int>(STD_COUNTERS::system_solver)] =new Counter(true,1, "Number of linear system resolutions Ax=B");
  std_counters_[static_cast<int>(STD_COUNTERS::petsc_solver)] =new Counter(true,1, "Petsc solver");
  std_counters_[static_cast<int>(STD_COUNTERS::implicit_diffusion)] =new Counter(true,1, "Number of linear system resolutions for implicit diffusion:");
  std_counters_[static_cast<int>(STD_COUNTERS::compute_dt)] =new Counter(true, 1, "Computation of the time step dt");
  std_counters_[static_cast<int>(STD_COUNTERS::turbulent_viscosity)] =new Counter(true, 1, "Turbulence model::update");
  std_counters_[static_cast<int>(STD_COUNTERS::convection)] =new Counter(true, 1, "Convection operator::add/compute");
  std_counters_[static_cast<int>(STD_COUNTERS::diffusion)] =new Counter(true, 1, "Diffusion operator::add/compute");
  std_counters_[static_cast<int>(STD_COUNTERS::gradient)] =new Counter(true, 1, "Gradient operator::add/compute");
  std_counters_[static_cast<int>(STD_COUNTERS::divergence)] =new Counter(true, 1, "Divergence operator::add/compute");
  std_counters_[static_cast<int>(STD_COUNTERS::rhs)] =new Counter(true, 1, "Source_terms::add/compute");
  std_counters_[static_cast<int>(STD_COUNTERS::postreatment)] =new Counter(true, 1, "Post-treatment");
  std_counters_[static_cast<int>(STD_COUNTERS::backup_file)] =new Counter(true, 1, "Back-up operations");
  std_counters_[static_cast<int>(STD_COUNTERS::restart)] =new Counter(true, 1,"Read file for restart");
  std_counters_[static_cast<int>(STD_COUNTERS::matrix_assembly)] =new Counter(true, 1, "Number of matrix assemblies for the implicit scheme:");
  std_counters_[static_cast<int>(STD_COUNTERS::update_variables)] =new Counter(true, 1, "Update");
  // MPI communication counters
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_sendrecv)] =new Counter(false, 2, "MPI_send_recv", "MPI_sendrecv", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_send)] =new Counter(false, 2, "MPI_send",      "MPI_sendrecv", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_recv)] =new Counter(false, 2, "MPI_recv",      "MPI_sendrecv", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_bcast)] =new Counter(false, 2, "MPI_broadcast", "MPI_sendrecv", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_alltoall)] =new Counter(false, 2, "MPI_alltoall",  "MPI_sendrecv", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_allgather)] =new Counter(false, 2, "MPI_allgather", "MPI_sendrecv", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_gather)] =new Counter(false, 2, "MPI_gather", "MPI_sendrecv", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_partialsum)]=new Counter(false, 2, "MPI_partialsum","MPI_allreduce", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_sumdouble)] =new Counter(false, 2, "MPI_sumdouble", "MPI_allreduce", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_mindouble)] =new Counter(false, 2, "MPI_mindouble", "MPI_allreduce", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_maxdouble)] =new Counter(false, 2, "MPI_maxdouble", "MPI_allreduce", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_sumfloat)] =new Counter(false, 2, "MPI_sumfloat", "MPI_allreduce", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_minfloat)] =new Counter(false, 2, "MPI_minfloat", "MPI_allreduce", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_maxfloat)] =new Counter(false, 2, "MPI_maxfloat", "MPI_allreduce", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_sumint)] =new Counter(false, 2, "MPI_sumint",    "MPI_allreduce", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_minint)] =new Counter(false, 2, "MPI_minint",    "MPI_allreduce", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_maxint)] =new Counter(false, 2, "MPI_maxint",    "MPI_allreduce", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_barrier)] =new Counter(false, 2, "MPI_barrier",   "MPI_allreduce", true);
  // GPU
  std_counters_[static_cast<int>(STD_COUNTERS::gpu_library)] =new Counter(true, 2, "GPU_library",       "GPU_library");
  std_counters_[static_cast<int>(STD_COUNTERS::gpu_kernel)]   =new Counter(true, 2, "GPU_kernel",        "GPU_kernel");
  std_counters_[static_cast<int>(STD_COUNTERS::gpu_copytodevice)] =new Counter(true, 2, "GPU_copyToDevice",  "GPU_copy");
  std_counters_[static_cast<int>(STD_COUNTERS::gpu_copyfromdevice)] =new Counter(true, 2, "GPU_copyFromDevice","GPU_copy");
  // Count the writing time in EcrireFicPartageXXX (big chunk of data in files XYZ or LATA)
  // For those two counters, quantity corresponds to the amount of written bytes
  std_counters_[static_cast<int>(STD_COUNTERS::IO_EcrireFicPartageMPIIO)] =new Counter(true, 2, "MPI_File_write_all", "IO"); // Call on each processor
  std_counters_[static_cast<int>(STD_COUNTERS::IO_EcrireFicPartageBin)] =new Counter(true, 2, "write", "IO"); // Call on master
  // Execution of Scatter::interpreter
  std_counters_[static_cast<int>(STD_COUNTERS::interprete_scatter)] =new Counter(true, 2, "Scatter_interprete");
  // Problemes
  // std_counters_[STD_COUNTERS::fluid_problem_] =new Counter(true, 3, "Fluid problem"));
  // std_counters_[STD_COUNTERS::fuel_problem_] =new Counter(true, 3, "Fuel problem"));
  // Appels a DoubleVect::echange_espace_virtuel()
  std_counters_[static_cast<int>(STD_COUNTERS::virtual_swap)] =new Counter(true, 2, "DoubleVect/IntVect::virtual_swap", "None", true);
  std_counters_[static_cast<int>(STD_COUNTERS::read_scatter)] =new Counter(true, 2, "Scatter::lire_domaine");
  */
  nb_steps_elapsed_ = 3;
  end_cache_=false;
  time_loop_=false;
  counters_stop_=false;
  time_cache_ = std::chrono::duration<double>::zero();
  computation_time_ = std::chrono::duration<double>::zero();
  last_opened_counter_ = nullptr;
}

Perf_counters::~Perf_counters()
{
  for (Counter* c :std_counters_)
    delete c;
  for (auto it = custom_counter_map_str_to_counter_.begin(); it != custom_counter_map_str_to_counter_.end(); ++it)
    delete it->second;
}
void Perf_counters::create_custom_counter(std::string counter_description , int counter_level,  std::string counter_family , bool is_comm,  bool to_print_in_global_TU)
{
  if (counter_level <=0)
    Process::exit("Custom counters should not be set with a zero or negative level value");
  if (custom_counter_map_str_to_counter_.count(counter_description)==0)
    {
      std::string error_msg = "Another custom counter already has the same name that the one you have given to your new counter: " + counter_description + "\n";
      Process::exit(error_msg);
    }
  auto result =custom_counter_map_str_to_counter_.emplace(counter_description, new Counter(to_print_in_global_TU,counter_level, counter_description, counter_family ,is_comm));
  if (!result.second)
    Process::exit("Failed to insert the new custom counter in the custom counter map");
}

void Perf_counters::print_in_global_TU(const STD_COUNTERS& name, bool to_print_or_not_to_print)
{
  Counter* c = access_std_counter(name);
  c->set_to_print(to_print_or_not_to_print);
}

void Perf_counters::print_in_global_TU(const std::string& name, bool to_print_or_not_to_print)
{
  Counter* c = access_custom_counter(name);
  c->set_to_print(to_print_or_not_to_print);
}

int Perf_counters::get_last_opened_counter_level() const
{
  if (last_opened_counter_ !=nullptr)
    return last_opened_counter_->level_;
  else
    return(-2);
}

void Perf_counters::check_begin(Counter* const c, int counter_lvl, std::chrono::time_point<std::chrono::high_resolution_clock> t)
{
  if (last_opened_counter_ != nullptr)
    {
      if (c->is_running_)
        Process::exit("The counter that you are trying to start is already running:" + c->description_);
      int expected_lvl = last_opened_counter_->level_ +1;
      if (counter_lvl != expected_lvl)
        {
          std::stringstream error_msg ;
          error_msg << "The counter you are trying to start does not have the expected level, counter running: " << last_opened_counter_->description_ << " counter that you try to open: " << c->description_  << " ; expected level: "<<  expected_lvl << std::endl ;
          Process::exit(error_msg.str());
        }
      if (time_loop_)
        {
          c->open_time_ts_ = t;
        }
      c->set_parent(last_opened_counter_);
    }
  last_opened_counter_ =c;
}

void Perf_counters::check_end(Counter* const c, std::chrono::time_point<std::chrono::high_resolution_clock> t)
{
  if (!c->is_running_)
    Process::exit("You are trying to close a counter that is not running: " + c->description_);
  if (last_opened_counter_ != c)
    {
      std::string error_msg = "The counter you are trying to close is not the last opened, counter: " + c->description_;
      Process::exit(error_msg);
    }
  if (time_loop_)
    {
      c->time_ts_ += t - c->open_time_ts_;
      c->open_time_ts_ = std::chrono::time_point<std::chrono::high_resolution_clock>();
    }
  last_opened_counter_ = c->parent_;
}

/*!
 *
 * @param std_cnt name in the enumerate STD_COUNTERS that corresponds to the counter you try to start
 * @param counter_lvl wanted lvl of the counter you try to start. It has to be equal to the lvl of the last called counter +1. counter_lvl become the new level of the counter you try to start.
 */
void Perf_counters::begin_count(const STD_COUNTERS& std_cnt, int counter_lvl)
{
  Counter* c = access_std_counter(std_cnt);
  if (counter_lvl == -100000)
    counter_lvl = c->level_;
  if (!time_loop_ || end_cache_ || c->level_ <= 0)
    {
      std::chrono::time_point<std::chrono::high_resolution_clock> t = std::chrono::high_resolution_clock::now();
      check_begin(c, counter_lvl,t);
      c->begin_count_(counter_lvl,t);
    }
}

/*!
 *
 * @param custom_count_name key of the map (custom_counter_map_str_to_counter_) of the custom counter you try to close
 * @param counter_lvl
 */
void Perf_counters::begin_count(const std::string& custom_count_name, int counter_lvl)
{
  Counter* c = access_custom_counter(custom_count_name);
  if (counter_lvl == -100000)
    counter_lvl = c->level_;
  if (!time_loop_ || end_cache_ )
    {
      std::chrono::time_point<std::chrono::high_resolution_clock> t = std::chrono::high_resolution_clock::now();
      check_begin(c, counter_lvl,t);
      c->begin_count_(counter_lvl,t);
    }
}

/*! @brief End the count of a counter and update the counter values
 *
 * @param c is the counter to end the count
 * @param count_increment is the count increment. If not specified, then it is equal to 1
 * @param quantity_increment is the increment of custom variable quantity. If not specified, it is set to 0.
 */
void Perf_counters::end_count(const STD_COUNTERS& std_cnt, int count_increment, int quantity_increment)
{
  Counter* c = access_std_counter(std_cnt);
  if (!time_loop_ || end_cache_ || c->level_ <= 0)
    {
      std::chrono::time_point<std::chrono::high_resolution_clock> t = std::chrono::high_resolution_clock::now();
      check_end(c, t);
      c->end_count_(count_increment, quantity_increment,t);
      if (c->level_ == -1)
        computation_time_ += c->total_time_;
    }
}


/*! @brief End the count of a counter and update the counter values
 *
 * @param c is the custom counter to end the count
 * @param count_increment is the count increment. If not specified, then it is equal to 1
 * @param quantity_increment is the increment of custom variable quantity. If not specified, it is set to 0.
 */
void Perf_counters::end_count(const std::string& custom_count_name, int count_increment, int quantity_increment)
{
  Counter* c = access_custom_counter(custom_count_name);
  if (!time_loop_ || end_cache_ || c->level_ <= 0)
    {
      std::chrono::time_point<std::chrono::high_resolution_clock> t = std::chrono::high_resolution_clock::now();
      assert(custom_counter_map_str_to_counter_.count(custom_count_name) > 0);
      check_end(c, t);
      c->end_count_(count_increment, quantity_increment,t);
    }
}

void Perf_counters::start_timeloop()
{
  if (last_opened_counter_==nullptr)
    Process::exit("You are trying to start the time loop before the start-up");
  time_loop_=true;
}

void Perf_counters::end_timeloop()
{
  if (!time_loop_)
    Process::exit("The time loop has not started, but you are trying to end it");
  time_loop_=false;
}

void Perf_counters::start_time_step()
{
  assert (last_opened_counter_!=nullptr);
  std::chrono::time_point<std::chrono::high_resolution_clock> t = std::chrono::high_resolution_clock::now();
  if (last_opened_counter_ == nullptr)
    Process::exit("You are trying to start a time step outside the time loop");
  Counter* c = last_opened_counter_;
  while (c->parent_ != nullptr && c->level_>=0)
    {
      c->open_time_ts_ = t;
      c = c->parent_;
    }
  for (Counter* c_std :std_counters_)
    c_std->time_ts_=std::chrono::duration<double>::zero();
  for (auto map_it = custom_counter_map_str_to_counter_.begin(); map_it != custom_counter_map_str_to_counter_.end(); ++map_it)
    map_it->second->time_ts_=std::chrono::duration<double>::zero();
}

double Perf_counters::get_time_since_last_open(const STD_COUNTERS& name)
{
  Counter* c = access_std_counter(name);
  std::chrono::duration<double> t = std::chrono::duration<double>::zero();
  if (!time_loop_ || end_cache_)
    {
      if (!c->is_running_)
        Process::exit("The counter is not running: " + c->description_);
      t = std::chrono::high_resolution_clock::now() - c->last_open_time_;
    }
  else
    Cerr << "Time of this step is set to zero if you are during the cache" << std::endl;
  return (t.count());
}

double Perf_counters::get_time_since_last_open(const std::string& name)
{
  Counter* c = access_custom_counter(name);
  std::chrono::duration<double> t = std::chrono::duration<double>::zero();
  if (!time_loop_ || end_cache_)
    {
      if (!c->is_running_)
        Process::exit("The counter is not running: " + c->description_);
      t = std::chrono::high_resolution_clock::now() - c->last_open_time_;
    }
  else
    Cerr << "Time of this step is set to zero if you are during the cache" << std::endl;
  return (t.count());
}

double Perf_counters::get_total_time(const STD_COUNTERS& name)
{
  Counter* c = access_std_counter(name);
  std::chrono::duration<double> t = c->total_time_;
  if (c->is_running_)
    t += std::chrono::high_resolution_clock::now() - c->last_open_time_;
  return (t.count());
}

double Perf_counters::get_total_time(const std::string& name)
{
  Counter* c = access_custom_counter(name);
  std::chrono::duration<double> t = c->total_time_;
  if (c->is_running_)
    t += std::chrono::high_resolution_clock::now() - c->last_open_time_;
  return (t.count());
}

/*! @brief Compute for each counter open during a time step avg_time_per_step_, min_time_per_step_, max_time_per_step_ and sd_time_per_step_
 *
 * Called at the end of each time step, and only then
 * @param tstep number of time steps elapsed since the start of the computation
 *
 * Nota : if the counter is called before the time step loop, it is not accounted for in the computation.
 */
void Perf_counters::compute_avg_min_max_var_per_step(unsigned int tstep)
{
  Perf_counters::stop_counters(); ///< stop_counters already updated c->tim_ts_
  if (last_opened_counter_ == nullptr)
    Process::exit("You are trying to compute the statistics of a time steps but have not open any counter");
  Counter* c_r=access_std_counter(STD_COUNTERS::total_execution_time);
  if (!time_loop_)
    Process::exit("You are trying to compute time loop statistics outside of the time loop");
  if (!end_cache_)
    {
      end_cache_ = tstep > nb_steps_elapsed_;
      if (end_cache_)
        {
          Process::barrier();
          time_cache_ = std::chrono::high_resolution_clock::now() - c_r->last_open_time_;
        }
    }
  int step = tstep - nb_steps_elapsed_;
  auto compute = [&](Counter* c)
  {
    if (c!=nullptr && time_loop_ && c->level_>0 )
      {
        c->min_time_per_step_ = (c->min_time_per_step_ < (c->time_ts_).count()) ? c->min_time_per_step_ : (c->time_ts_).count();
        c->max_time_per_step_ = (c->min_time_per_step_ > (c->time_ts_).count()) ? c->min_time_per_step_ : (c->time_ts_).count();
        c->avg_time_per_step_ = ((step-1)*c->avg_time_per_step_ + (c->time_ts_).count())/step;
        c->sd_time_per_step_ += ((c->time_ts_).count()* (c->time_ts_).count() - 2*((c->time_ts_).count())* c->avg_time_per_step_ +  c->avg_time_per_step_ *  c->avg_time_per_step_)/step;
        c->sd_time_per_step_ = sqrt(c->sd_time_per_step_);
        if (c->sd_time_per_step_ < 0)
          c->sd_time_per_step_ = 0;
      }
    c->open_time_ts_=std::chrono::time_point<std::chrono::high_resolution_clock>();
  };
  if (end_cache_ || !time_loop_)
    {
      for (Counter *c: std_counters_)
        compute(c);
      if (!custom_counter_map_str_to_counter_.empty())
        {
          for (auto map_it = custom_counter_map_str_to_counter_.begin(); map_it != custom_counter_map_str_to_counter_.end(); ++map_it)
            compute(map_it->second);
        }
    }
  Perf_counters::restart_counters();
}

/*! @brief Return a string that contains os information
 *
 * for example:
 * node name   = is247793
 * system name = Linux
 * machine     = x86_64
 * release     = 6.5.0-44-generic
 * version     = #44~22.04.1-Ubuntu SMP PREEMPT_DYNAMIC + current date and time
 *
 * @return string with is, os info and date
 */
std::string Perf_counters::get_os() const
{
  std::string result;
  struct utsname buffer;
  if (uname(&buffer) == -1)
    return "Error: Unable to retrieve OS info";
  result += std::string(buffer.nodename) + "__";
  result += std::string(buffer.sysname) + "__";
  result += std::string(buffer.machine) + "__";
  result += std::string(buffer.release)+ "__";
  result += std::string(buffer.version);
  return result;
}

/*!
 *
 * @return string that contains cpu model and number of proc
 */

std::string Perf_counters::get_cpu() const
{
  int result;
  result = std::system("lscpu 2>/dev/null | grep 'Model name' > cpu_detail.txt");
  if (result !=0)
    Cerr << "Bash command in Perf_counters::get_cpu failed"<<std::endl;
  result = std::system("lscpu 2>/dev/null | grep 'Thread(s) per core' >> cpu_detail.txt");
  if (result !=0)
    Cerr << "Bash command in Perf_counters::get_cpu failed" <<std::endl;
  std::ifstream file("cpu_detail.txt");
  if (!file.is_open())
    {
      Cerr << "Failed to open file in get_cpu: " << std::endl;
      return "";
    }
  std::string line1, line2;
  std::getline(file, line1); // Read the first line
  std::getline(file, line2); // Read the second line
  file.close();
  // Concatenate the two lines with a space in between
  std::string str = line1 + " ; " + line2;
  result = std::system("rm cpu_detail.txt");
  if (result !=0)
    Cerr << "Bash command in Perf_counters::get_cpu failed"<<std::endl;
  return (str);
}
/*!
 *
 * @return string with gpu model name
 */
std::string Perf_counters::get_gpu() const
{
  std::string gpu_description = "No GPU was used for the computation";
#ifdef TRUST_USE_CUDA
  int result = std::system("nvidia-smi 2>/dev/null | grep NVIDIA > gpu_detail.txt");
  if (result !=0)
    Cerr << "Bash command in Perf_counters::get_gpu failed"<<std::endl;
  std::stringstream gpu_desc;
  gpu_desc << std::ifstream("gpu_detail.txt").rdbuf();
  result = std::system("rm gpu_detail.txt");
  if (result !=0)
    Cerr << "Bash command in Perf_counters::get_gpu failed"<<std::endl;
  gpu_description = gpu_desc.str();
#endif
#ifdef TRUST_USE_HIP
  int result_ = std::system("rocminfo 2>/dev/null | grep Marketing > gpu_detail.txt");
  if (result_ !=0)
    Cerr << "Bash command in Perf_counters::get_gpu failed"<<std::endl;
  std::stringstream gpu_desc;
  gpu_desc << std::ifstream("gpu_detail.txt").rdbuf();
  result_ = std::system("rm gpu_detail.txt");
  if (result_ !=0)
    Cerr << "Bash command in Perf_counters::get_gpu failed"<<std::endl;
  gpu_description = gpu_desc.str();
#endif
  return gpu_description;
}

/*!
 *
 * @return string with complete date and time : DD-MM-YYYY -- hour:minute:second
 */
std::string Perf_counters::get_date() const
{
  time_t now = time(0);
  std::stringstream date;
  struct tm tstruct = *localtime(&now);
  date<< std::put_time(&tstruct, "%d-%m-%Y -- %X");
  return date.str();
}

/*! @brief Stop all counters, has to be called on every processor simultaneously
 *
 */
void Perf_counters::stop_counters()
{
  std::chrono::time_point<std::chrono::high_resolution_clock> t_stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_elapsed_before_stop;
  if (counters_stop_)
    Process::exit("The counter are already stop, you can't stop them two times in a row \n");
  if (last_opened_counter_ != nullptr)
    {
      Counter* c = last_opened_counter_;
      c->time_alone_ += t_stop -c->last_open_time_alone_;
      c->last_open_time_alone_ = std::chrono::time_point<std::chrono::high_resolution_clock>();
      while (c!=nullptr)
        {
          time_elapsed_before_stop= t_stop -c->last_open_time_;
          if (time_loop_ )
            c->time_ts_ += time_elapsed_before_stop;
          c->total_time_ += time_elapsed_before_stop;
          c = c->parent_;
        }
    }
  Counter* c_time = access_std_counter(STD_COUNTERS::total_execution_time);
  computation_time_ += c_time->total_time_;
  counters_stop_=true;
}

/*! @brief Restart all counters, has to be called on every processor simultaneously
 *
 */
void Perf_counters::restart_counters()
{
  std::chrono::time_point<std::chrono::high_resolution_clock> t_restart = std::chrono::high_resolution_clock::now();
  if (counters_stop_==false)
    Process::exit("Try to restart counters but they have never been stopped before");
  if (last_opened_counter_ != nullptr)
    {
      Counter* c = last_opened_counter_;
      c->last_open_time_alone_ = t_restart;
      while (c !=nullptr )
        {
          c->last_open_time_ = t_restart;
          c->open_time_ts_ = t_restart;
          c = c->parent_;
        }
    }
  counters_stop_=false;
}

void Perf_counters::reset_counters()
{
  for (Counter* c :std_counters_)
    {
      if (c!=nullptr)
        c->reset();
    }
  for (auto map_it = custom_counter_map_str_to_counter_.begin(); map_it != custom_counter_map_str_to_counter_.end(); ++map_it)
    {
      Counter* c_custom = map_it->second;
      if (c_custom!=nullptr)
        c_custom->reset();
    }
}

static void build_line_csv(std::stringstream& lines, const std::vector<std::string>& line_items, const std::vector<int>& item_size)
{
  int size_of_str_to_add = 50;
  long unsigned int len_line = line_items.size();
  for (long unsigned int i=0 ; i<len_line ; i++)
    {
      size_of_str_to_add = item_size[i];
      lines << std::setw(size_of_str_to_add) ; ///< Ensure that each item of a column has the same size
      lines << line_items[i];
      if (i == len_line -1)
        lines << std::endl ; ///< if end_line == True, then add a break line but no delimiter
      else
        lines << " \t"; ///< Put the column delimiter if we are not at the end of the line
    }
}

void Perf_counters::print_performance_to_csv(const std::string& message,const bool mode_append)
{
  /*
  stop_counters();
  assert(!message.empty());
  std::stringstream perfs;   ///< Stringstream that contains stats for each processor
  std::stringstream perfs_globales;   ///< Stringstream that contains stats average on the processors : processor number = -1
  std::stringstream File_header;      ///< Stringstream that contains the lines at the start of the file

  unsigned int length_line = 24; ///< number of item of a line of the _csv.Tu file
  std::vector<int> item_size(length_line,20);   ///< Contains the the width of the printed string, 20 for numbers by default
  std::vector<std::string> line_items(length_line,"");   ///< Contains the data of a line that we want to print in the _csv.TU file.

  std::stringstream tmp_item; ///< Create a temporary stringstream for converting wanted line items in string to construct the line_items vector and therefore
  unsigned int nb_procs =  Process::nproc();

  /// We specify the width of large items of lines of the _csv.Tu file for making it readable by human
  item_size[0] = 50;
  item_size[2] = 40;
  item_size[3] = 45;

  if ( (Process::je_suis_maitre()) && (message == "Computation start-up statistics") )
    {
      File_header << "# Detailed performance log file. See the associated validation form for an example of data analysis"<< std::endl;
      File_header << "# Date of the computation:     " << get_date() << std::endl;
      File_header << "# OS used:     " << get_os() << std::endl;
      File_header << "# CPU info:     " << get_cpu() << std::endl;
      File_header << "# GPU info:     " << get_gpu() << std::endl;
      File_header << "# Number of processor used = " << nb_procs << std::endl;
      File_header << "# The time was measured by the following method using std::chrono::high_resolution_clock::now()" << std::endl ;
      File_header << "# By default, only averaged statistics on all processor are printed. For accessing the detail per processor, add 'stat_per_proc_perf_log 1' in the data file"<< std::endl;
      File_header << "# Processor number equal to -1 corresponds to the performance of the calculation averaged on the processors during the simulation step" << std::endl;
      File_header << "# If a counter does not belong in any particular family, then counter family is set to None" << std::endl;
      File_header << "# Count means the number of time the counter is called during the overall calculation step." << std::endl;
      File_header << "# Min, max and SD accounts respectively for the minimum, maximum and Standard Deviation of the quantity of the previous row." << std::endl;
      File_header << "# Quantity is a custom variable that depends on the counter. It is used to compute bandwidth for communication counters for example. See the table at the end of the introduction on statistics in TRUST form for more details." << std::endl;
      File_header << "#" << std::endl << "#" << std::endl;
      /// Then we create a vector line_items that contains each item we want to print
      line_items[0] = "Overall_simulation_step";
      line_items[1] = "Processor_Number";
      line_items[2] = "Counter_family";
      line_items[3] = "Counter_name";
      line_items[4] = "Counter_level";
      line_items[5] = "Is_comm";
      line_items[6] = "%_total_time";
      line_items[7] = "total time_(s)";
      line_items[8] = "t_min";
      line_items[9] = "t_max";
      line_items[10] = "t_SD";
      line_items[11] = "time alone(s)";
      line_items[12] = "t_alone_min";
      line_items[13] = "t_alone_max";
      line_items[14] = "t_alone_SD";
      line_items[15] = "count";
      line_items[16] = "time_per_step";
      line_items[17] = "tps_min";
      line_items[18] = "tps_max";
      line_items[19] = "tps_SD";
      line_items[20] = "Quantity";
      line_items[21] = "q_min";
      line_items[22] = "q_max";
      line_items[23] = "q_SD";
      assert(item_size.size()==length_line);
      assert(line_items.size()==item_size.size());
      /// After filling line_items and item_size, we use the function build_line_csv to build the line at the expected format
      build_line_csv(File_header,line_items,item_size);
    }

  /// Check if all of the processors see the same number of counter, if not print an error message in perfs_globales
  bool skip_globals = Objet_U::disable_TU;
  int total_nb_of_counters = (int)std_counters_.size() + (int)custom_counter_map_str_to_counter_.size();
  int min_total_nb_of_counters = (int) Process::mp_min(total_nb_of_counters);
  int max_total_nb_of_counters = (int) Process::mp_max(total_nb_of_counters);

  if ( (max_total_nb_of_counters - min_total_nb_of_counters)!=0 )
    {
      if (Process::je_suis_maitre())
        {
          perfs_globales << "Unable to collect statistics :" << std::endl
                         << " there is not the same number of counters on all"
                         " processors."<< std::endl;
        }
      skip_globals = true; ///< If min_nb_of_counters != max_nb_of_counters, aggregated stats are not printed
    }
  Counter* c_time = access_std_counter(STD_COUNTERS::total_execution_time);
  std::chrono::duration<double> total_time = c_time->total_time_;
  c_time = access_std_counter(STD_COUNTERS::timeloop);
  int nb_ts = c_time->count_;

  if (nb_ts == 0)
    {
      if (Process::je_suis_maitre())
        {
          perfs_globales << "The computation didn't start" << std::endl;
        }
      skip_globals = true; ///< If min_nb_of_counters != max_nb_of_counters, aggregated stats are not printed
    }


  //*****************************************************************************************************************************
  //                                            Lambda function to build line
  //*****************************************************************************************************************************
  int level; ///< Level of details of the counter
  bool is_comm; ///< Equal to 1 if the counter is a communication counter, 0 otherwise
  int count;  ///< number of time the counter is open and closed
  int quantity, min_quantity=0, max_quantity=0; ///< A custom quantity which depends on the counter. Used for example to compute the bandwidth
  double time,time_alone,min_time_alone=0.,max_time_alone=0.,SD_time_alone=0.0;
  double percent_time, min_time=0.0, max_time=0.0; ///< Percent of the total time used in the method tracked by the counter
  double SD_time=0.0, SD_quantity=0.0; ///< the standard dev of all the prev vars
  double avg_time_per_step, min_time_per_step, max_time_per_step, sd_time_per_step;

  auto fill_items = [&](int proc_number, const std::string& desc, const std::string& familly)
  {
    tmp_item << message; ///< Convert into string the item we want to print in the line, here the overall simulation step
    line_items[0] = tmp_item.str(); ///< Add the item to the vector line_itmes, used to construct the line of the _csv.TU file

    tmp_item.str(""); ///< Empties the temporary stringstream

    tmp_item<< proc_number;
    line_items[1] = tmp_item.str(); ///< Add the processor number to the vector line_items
    tmp_item.str("");

    tmp_item<< familly;
    line_items[2] = tmp_item.str(); ///< Add the counter's family to the vector line_items, null if the counter does not belong in a family
    tmp_item.str("");

    tmp_item << desc;
    line_items[3] = tmp_item.str(); ///< Add the counter's name to the vector line_items
    tmp_item.str("");

    tmp_item<< level;
    line_items[4] = tmp_item.str(); ///< Add the counter's level to the vector line_items
    tmp_item.str("");

    tmp_item<< is_comm;
    line_items[5] = tmp_item.str(); ///< Add 1 if the counter is a communication counter, 0 otherwise
    tmp_item.str("");

    tmp_item<< std::setprecision(4);
    tmp_item<< percent_time;
    line_items[6] = tmp_item.str(); ///< Add the percent of total time used by the operation tracked by counter i to the vector line_items
    tmp_item.str("");

    tmp_item << std::scientific << std::setprecision(7);
    tmp_item<< time;
    line_items[7] = tmp_item.str(); ///< Time elapsed when using the operation tracked by counter i
    tmp_item.str("");

    tmp_item<< min_time;
    line_items[8] = tmp_item.str(); ///< Detail per proc, so the min, max, avg and SD on proc is equal to 0
    tmp_item.str("");

    tmp_item<< max_time;
    line_items[9] = tmp_item.str(); ///< Detail per proc, so the min, max, avg and SD on proc is equal to 0
    tmp_item.str("");

    tmp_item<< SD_time;
    line_items[10] = tmp_item.str(); ///< Detail per proc, so the min, max, avg and SD on proc is equal to 0
    tmp_item.str("");

    tmp_item<< time_alone;
    line_items[11] = tmp_item.str(); ///< Detail per proc, so the min, max, avg and SD on proc is equal to 0
    tmp_item.str("");

    tmp_item<< min_time_alone;
    line_items[11] = tmp_item.str(); ///< Detail per proc, so the min, max, avg and SD on proc is equal to 0
    tmp_item.str("");

    tmp_item<< max_time_alone;
    line_items[11] = tmp_item.str(); ///< Detail per proc, so the min, max, avg and SD on proc is equal to 0
    tmp_item.str("");

    tmp_item<< SD_time_alone;
    line_items[11] = tmp_item.str(); ///< Detail per proc, so the min, max, avg and SD on proc is equal to 0
    tmp_item.str("");

    tmp_item<< count;
    line_items[12] = tmp_item.str();  ///< Number of time the counter was called on the overall simulation step
    tmp_item.str("");

    tmp_item<< avg_time_per_step;
    line_items[13] = tmp_item.str(); ///< Averaged time elapsed by time step for the operation tracked by the counter on the overall simulation step
    tmp_item.str("");

    tmp_item<< min_time_per_step;
    line_items[14] = tmp_item.str(); ///< Minimum time elapsed by time step for the operation tracked by the counter on the overall simulation step
    tmp_item.str("");

    tmp_item<< max_time_per_step;
    line_items[15] = tmp_item.str(); ///< Maximum time elapsed by time step for the operation tracked by the counter on the overall simulation step
    tmp_item.str("");

    tmp_item<< sqrt(sd_time_per_step);
    line_items[16] = tmp_item.str(); ///< Standard Deviation of time elapsed by time step for the operation tracked by the counter on the overall simulation step
    tmp_item.str("");

    tmp_item<< quantity;
    line_items[17] = tmp_item.str(); ///< Custom variable that depends on the counter the overall simulation step
    tmp_item.str("");

    tmp_item<< min_quantity;
    line_items[18] = tmp_item.str(); ///< Detail per proc, so the min, max and SD on proc is equal to 0
    tmp_item.str("");

    tmp_item<< max_quantity;
    line_items[19] = tmp_item.str(); ///< Detail per proc, so the min, max and SD on proc is equal to 0
    tmp_item.str("");

    tmp_item<< SD_quantity;
    line_items[20] = tmp_item.str(); ///< Detail per proc, so the min, max and SD on proc is equal to 0
    tmp_item.str("");
  };

  auto extract_stats = [&](Counter * c_lambda)
  {
    if (c_lambda->count_ > 0)
      {
        level = c_lambda->level_; ///< Level of details of the counter
        is_comm = c_lambda->is_comm_; ///< Equal to 1 if the counter is a communication counter, 0 otherwise
        time = c_lambda->total_time_.count();
        time_alone = c_lambda->time_alone_.count();
        count = c_lambda->count_;
        quantity = c_lambda->quantity_;
        avg_time_per_step = c_lambda->avg_time_per_step_;
        min_time_per_step = c_lambda->min_time_per_step_;
        max_time_per_step = c_lambda->max_time_per_step_;
        sd_time_per_step = c_lambda->sd_time_per_step_;
        min_time = 0.;
        max_time = 0.;
        SD_time = 0.;
        min_quantity = 0.;
        max_quantity = 0.;
        SD_quantity = 0.;
        min_time_alone = 0.;
        max_time_alone = 0.;
        SD_time_alone = 0.;
        if (Objet_U::stat_per_proc_perf_log && Process::is_parallel())
          {
            fill_items(Process::me(),c_lambda->description_, c_lambda->family_);
            build_line_csv(perfs,line_items,item_size);  ///< Build the line of the stats associated on the counter i for a single proc
          }
        if (Process::je_suis_maitre()&& Process::is_parallel())
          {
            std::array< std::array<double,4> ,4> table = c_lambda->compute_min_max_avg_sd_();
            time = table[0][2];
            min_time = table[0][0];
            max_time = table[0][1];
            SD_time = table[0][3];
            quantity = static_cast <int>(std::floor(table[1][2]));
            min_quantity = static_cast <int>(std::floor(table[1][0]));
            max_quantity = static_cast <int>(std::floor(table[1][1]));
            SD_quantity = table[1][3];
            time_alone = table[3][2];
            min_time_alone = table[3][0];
            max_time_alone = table[3][1];
            SD_time_alone = table[3][3];
            if (! skip_globals )
              {
                fill_items(-1,c_lambda->description_,c_lambda->family_);
                build_line_csv(perfs_globales,line_items,item_size);
              }
          }
      }
  };
  for (Counter* c_std : std_counters_)
    {
      if (c_std!=nullptr)
        extract_stats(c_std);
    }

  for (auto map_it = custom_counter_map_str_to_counter_.begin(); map_it != custom_counter_map_str_to_counter_.end(); ++map_it)
    {
      Counter* c_custom = map_it->second;
      if (c_custom!=nullptr)
        extract_stats(c_custom);
    }
  Nom CSV(Objet_U::nom_du_cas());
  CSV+="_csv_new.TU";
  std::string root=Sortie_Fichier_base::root;
  Sortie_Fichier_base::root = "";
  EcrFicPartage file(CSV, mode_append ? (ios::out | ios::app) : (ios::out));
  Sortie_Fichier_base::root = root;
  file << File_header.str();
  file << perfs_globales.str();
  file << perfs.str();
  file.syncfile();
  restart_counters();
  */
}

/*!@brief Function used for computing communication statistics in the global.TU
 *
 * @param time
 * @param quantity
 * @param count
 *
 * The three parameters are updated by their mean value over the processors
 * @return an array of array that contains the min, max, average and standard deviation over the processors of the three parameters of the function
 */
inline std::array< std::array<double,4> ,3> compute_min_max_avg_sd(double& time, int& quantity, int& count)
{
  assert(Process::is_parallel());
  double qty,cnt,min,max,avg,sd ;
  qty=static_cast<double>(quantity);
  cnt = static_cast<double>(count);
  min = Process::mp_min(time);
  max = Process::mp_max(time);
  avg = Process::mp_sum(time)/Process::nproc();
  sd = sqrt(Process::mp_sum((time-avg)*(time-avg))/Process::nproc());
  time =avg;

  std::array<double,4> min_max_avg_sd_time = {min,max,avg,sd};

  min = Process::mp_min(qty);
  max = Process::mp_max(qty);
  avg = Process::mp_sum(qty)/Process::nproc();
  sd = sqrt(Process::mp_sum((qty-avg)*(qty-avg))/Process::nproc());
  quantity = static_cast<int>(std::floor (avg));

  std::array<double,4> min_max_avg_sd_quantity = {min,max,avg,sd};

  min = 1.0*Process::mp_min(cnt);
  max = 1.0*Process::mp_max(cnt);
  avg = 1.0*Process::mp_sum(cnt)/Process::nproc();
  sd = sqrt(Process::mp_sum((cnt-avg)*(cnt-avg))/Process::nproc());
  count = static_cast<int>(std::floor (avg));

  std::array<double,4> min_max_avg_sd_count = {min,max,avg,sd};

  return {min_max_avg_sd_time,min_max_avg_sd_quantity,min_max_avg_sd_count};
}

/*!
 *
 * @param message
 * @param mode_append
 */
void Perf_counters::print_global_TU(const std::string& message, const bool mode_append)
{
  assert(!message.empty());
  stop_counters();
  std::stringstream perfs_TU;   ///< Stringstream that contains algomerated stats that will be printed in the .TU
  std::stringstream perfs_GPU;   ///< Stringstream that contains algomerated stats that will be printed in the .TU
  std::stringstream perfs_IO;   ///< Stringstream that contains algomerated stats that will be printed in the .TU
  std::stringstream captions;
  std::stringstream File_header;      ///< Stringstream that contains the File header

  int nb_procs =  Process::nproc();
  double comm_allreduce_t = 0.0, comm_sendrecv_t = 0.0;
  int comm_allreduce_q = 0.0,comm_sendrecv_q = 0.0,comm_allreduce_c = 0,comm_sendrecv_c = 0;
  std::array< std::array<double,4> ,3> min_max_avg_sd_t_q_c_sendrecv_comm ;
  for (std::array<double,4>& arr: min_max_avg_sd_t_q_c_sendrecv_comm)
    for (double & d : arr)
      d=-1;
  std::array< std::array<double,4> ,3> min_max_avg_sd_t_q_c_allreduce_comm = min_max_avg_sd_t_q_c_sendrecv_comm ;
  std::array< std::array<double,4> ,4> min_max_avg_sd_t_q_c_echange_espace_virtuel ;
  for (std::array<double,4>& arr: min_max_avg_sd_t_q_c_echange_espace_virtuel)
    for (double & d : arr)
      d=-1;
  std::array< std::array<double,4> ,4> min_max_avg_sd_t_q_c_sendrecv = min_max_avg_sd_t_q_c_echange_espace_virtuel;
  Counter* c = access_std_counter(STD_COUNTERS::timeloop);
  int nb_ts = (Process::mp_max(c->count_) - nb_steps_elapsed_)>0 ? Process::mp_max(c->count_) - nb_steps_elapsed_ : 0;
  double total_untracked_time_ts=c->time_alone_.count();
  c = access_std_counter(STD_COUNTERS::total_execution_time);
  double total_time = Process::mp_max(c->total_time_.count());
  double total_untracked_time=c->time_alone_.count();
  double total_comm_time=0.;

  auto write_globalTU_line = [&] (Counter*c_ptr_,std::stringstream & line)
  {
    if (time_loop_ && c_ptr_->count_>0)
      {
        double t_c = c_ptr_->total_time_.count();
        int count = c_ptr_->count_;
        if (Process::je_suis_maitre())
          {
            t_c = Process::mp_max(t_c);
            count = Process::mp_max(count);
          }
        line << "Whose: " << c_ptr_->description_ << " ";
        double t = nb_ts>0 ? t_c/nb_ts : t_c;
        line << t << t_c/total_time*100 ;
        if (nb_ts>0)
          {
            double n = count/nb_ts;
            line << " (" << n << " " << (n==1?"call":"calls") << "/time steps)";
          }
        line << std::endl;
      }
  };

  for (Counter * c_com : std_counters_)
    {
      if (c_com!=nullptr && c_com->count_>0)
        {
          if (c_com->is_comm_)
            {
              if (c_com->family_=="MPI_allreduce")
                {
                  comm_allreduce_q += c_com->quantity_;
                  comm_allreduce_t += c_com->total_time_.count();
                  comm_allreduce_c += c_com->count_;
                }
              if (c_com->family_=="MPI_sendrecv")
                {
                  comm_sendrecv_q += c_com->quantity_;
                  comm_sendrecv_t += c_com->total_time_.count();
                  comm_sendrecv_c += c_com->count_;
                }
              total_comm_time += c_com->time_alone_.count();
            }
        }
    }
  Counter* c_com = nullptr;
  for (auto map_it = custom_counter_map_str_to_counter_.begin(); map_it != custom_counter_map_str_to_counter_.end(); ++map_it)
    {
      c_com = map_it->second;
      if (c_com!=nullptr)
        {
          if (c_com->is_comm_)
            {
              if (c_com->family_=="MPI_allreduce")
                {
                  comm_allreduce_q += c_com->quantity_;
                  comm_allreduce_t += c_com->total_time_.count();
                  comm_allreduce_c += c_com->count_;
                }
              if (c_com->family_=="MPI_sendrecv")
                {
                  comm_sendrecv_q += c_com->quantity_;
                  comm_sendrecv_t += c_com->total_time_.count();
                  comm_sendrecv_c += c_com->count_;
                }
              total_comm_time += c_com->time_alone_.count();
            }
        }
    }
  if(Process::is_parallel())
    {
      min_max_avg_sd_t_q_c_allreduce_comm = compute_min_max_avg_sd(comm_allreduce_t,comm_allreduce_q,comm_allreduce_c);
      min_max_avg_sd_t_q_c_sendrecv_comm =  compute_min_max_avg_sd(comm_sendrecv_t,comm_sendrecv_q,comm_sendrecv_c);
      c = access_std_counter(STD_COUNTERS::mpi_sendrecv);
      min_max_avg_sd_t_q_c_sendrecv = c->compute_min_max_avg_sd_();
      c = access_std_counter(STD_COUNTERS::virtual_swap);
      min_max_avg_sd_t_q_c_echange_espace_virtuel = c->compute_min_max_avg_sd_();
    }
  if (Process::je_suis_maitre())
    {
      if (message == "Computation start-up statistics")
        {
          File_header << "                                                  # Global performance file #"<< std::endl;
          File_header <<  std::endl;
          File_header << "------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
          File_header << "Date:     " << get_date() << std::endl;
          File_header << "OS:     " << get_os() << std::endl;
          File_header << "CPU:     " << get_cpu() << std::endl;
          File_header << "GPU:     " << get_gpu() << std::endl;
          File_header << "Number of processor used = " << nb_procs << std::endl ;
          File_header << "------------------------------------------------------------------------------------------------------------------------------------------" << std::endl << std::endl << std::endl;
          File_header << message << std::endl<< std::endl;
          c = access_std_counter(STD_COUNTERS::total_execution_time);
          File_header << "Total time of the start-up (averaged by proc): " <<  c->total_time_.count() << std::endl;
          File_header << "Percent of untracked time (averaged per proc) during computation start-up: " << Process::mp_sum(total_untracked_time)/Process::mp_sum(total_time) << std::endl<< std::endl;
        }
      else if (message == "Time loop statistics")
        {
          File_header << "------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
          File_header << message << std::endl<< std::endl;
          nb_ts = Process::mp_max(nb_ts);
          if (nb_ts <= 0)
            Process::exit("No time step after cache filling was computed");
          File_header << "The " <<  nb_steps_elapsed_<< " first time steps are not accounted for the computation of the statistics"<< std::endl;
          c = access_std_counter(STD_COUNTERS::total_execution_time);
          File_header << "Total time averaged per proc: " <<  Process::mp_sum(c->total_time_.count())/nb_procs << std::endl;
          File_header << "Number of time steps: " << Process::mp_max(nb_ts) << std::endl;
          c = access_std_counter(STD_COUNTERS::timeloop);
          File_header << "Average time per time steps per proc: " << Process::mp_sum(c->time_ts_.count())/nb_procs << " ; Standard deviation between time steps: " << c->sd_time_per_step_ << std::endl;
          File_header << "Time of cache :" << Process::mp_sum(time_cache_.count())/nb_procs <<" s/proc" <<std::endl << std::endl;
          File_header << "Percent of total time (averaged per proc) tracked by communication counters:" << 100* Process::mp_sum(total_comm_time) / Process::mp_sum(total_time) << std::endl;
          File_header << "Percent of untracked time (averaged per proc) outside of the time loop: " << Process::mp_sum(total_untracked_time)/Process::mp_sum(total_time) << std::endl;
          File_header << "Percent of untracked time (averaged per proc) inside the time loop: " << Process::mp_sum(total_untracked_time_ts)/Process::mp_sum(c->total_time_.count()) << std::endl ;
        }
      else if (message == "Post-treatment statistics")
        {
          File_header << "------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
          File_header << message << std::endl<< std::endl;
          c = access_std_counter(STD_COUNTERS::total_execution_time);
          File_header << "Average time per proc of the post-treatment: " <<  Process::mp_sum(c->total_time_.count())/nb_procs << std::endl;
          File_header << "Percent of untracked time (averaged per proc) during post-treatment: " << Process::mp_sum(total_untracked_time)/nb_procs << std::endl<< std::endl;
          captions <<  std::endl;
          captions << "Max waiting time big    => probably due to a bad partitioning" << std::endl;
          captions << "Communications > 30%    => too many processors or network too slow" << std::endl;
          captions << std::endl;
        }
      else
        Process::exit("You are trying to get stats of an unknown computation step");

      for (Counter *c_ptr :std_counters_)
        {
          if (c_ptr!=nullptr && c_ptr->to_print_in_global_TU_)
            write_globalTU_line(c_ptr,perfs_TU);
        }
      // Loop on the custom counters
      Counter * c_custom_ptr;
      for (auto map_it = custom_counter_map_str_to_counter_.begin(); map_it != custom_counter_map_str_to_counter_.end(); ++map_it)
        {
          c_custom_ptr = map_it->second;
          if (c_custom_ptr->count_ > 0 && c_custom_ptr->to_print_in_global_TU_)
            write_globalTU_line(c_custom_ptr, perfs_TU);
        }
      c = access_std_counter(STD_COUNTERS::virtual_swap);
      if (Process::mp_max(c->count_)>0)
        {
          perfs_TU << "Maximum number of virtual exchanges per time steps :" <<  Process::mp_max(c->count_) << std::endl;
        }
      if (min_max_avg_sd_t_q_c_allreduce_comm[2][1]>0)
        {
          double allreduce_per_ts = (double) min_max_avg_sd_t_q_c_allreduce_comm[2][1]/nb_ts;
          perfs_TU << "Maximum number of MPI allreduce per time step" << allreduce_per_ts << std::endl;
          if (allreduce_per_ts > 30.0)
            {
              perfs_TU << "------------------------------------------------------------------------------------------------------------------------------------------" << std::endl;
              perfs_TU << " Warning: The number of MPI_allreduce calls per time step is high. Contact TRUST support if you plan to run massive parallel calculation. " << std::endl;
              perfs_TU << "------------------------------------------------------------------------------------------------------------------------------------------"<< std::endl;
            }
        }
      c = access_std_counter(STD_COUNTERS::system_solver);
      int tmp = Process::mp_max(c->quantity_);
      if (tmp > 0)
        {
          perfs_TU << "Number of call of the solver per time step: "<< 1.0*tmp / nb_ts << std::endl;
          double avg_time = Process::mp_max(c->total_time_.count()) / tmp;
          perfs_TU << "Average time of the resolution per call: " << avg_time << std::endl;
          perfs_TU << "Average number of iteration per call: "<< Process::mp_max(c->quantity_) / tmp << std::endl;
        }

      c = access_std_counter(STD_COUNTERS::backup_file);
      double total_quantity = static_cast<double>(Process::mp_sum(c->quantity_));
      c = access_std_counter(STD_COUNTERS::backup_file);
      tmp = Process::mp_max(c->count_);
      if (tmp>0)
        {
          perfs_TU << "Number of back-up: " << tmp << std::endl;
          perfs_TU << "Average amount of data per back-up (Mo): " << total_quantity / (tmp *1024*1024) << std::endl;
        }
      // GPU part of the TU :
      auto compute_percent_and_write = [&] (const std::string str)
      {
        double max_time = Process::mp_max(c->total_time_.count());
        double percent = 100*max_time/nb_ts;
        double calls = Process::mp_max(c->count_)/nb_ts;
        double t_ts = max_time/nb_ts;
        double bw = c->quantity_*(1021.*1024.*Process::mp_max(c->total_time_.count()));
        perfs_GPU << str << " time per time step:"<< t_ts <<" ; number of calls per time steps: " << calls << " ; percent of total time: " << percent << " ; bandwidth: " << bw <<  std::endl;
        return percent;
      };
      c = access_std_counter(STD_COUNTERS::gpu_copytodevice);
      tmp = Process::mp_max(c->count_);
      if (tmp>0 && nb_ts >0)
        {
          perfs_GPU << std::endl << "-------------------------------------------------------------------GPU--------------------------------------------------------------------" << std::endl;
          perfs_GPU << "GPU statistics per time step (experimental):" << std::endl;
          c = access_std_counter(STD_COUNTERS::gpu_library);
          double ratio_gpu_library = compute_percent_and_write("Libraries: ");
          c = access_std_counter(STD_COUNTERS::gpu_kernel);
          double ratio_gpu_kernel = compute_percent_and_write("Kernels: ");
          double ratio_gpu = ratio_gpu_kernel+ratio_gpu_library;
          c = access_std_counter(STD_COUNTERS::gpu_copytodevice);
          double ratio_copy = compute_percent_and_write("Copy host to device: ");
          c = access_std_counter(STD_COUNTERS::gpu_copyfromdevice);
          ratio_copy += compute_percent_and_write("Copy device to host: ");
          c = access_std_counter(STD_COUNTERS::timeloop);
          double ratio_comm = 100.0 * (comm_sendrecv_t+ comm_allreduce_t)/c->max_time_per_step_;
          double ratio_cpu = 100 - ratio_copy - ratio_gpu - ratio_comm;
          perfs_GPU << std::setprecision(3) << "GPU: " << ratio_gpu << "% Copy H<->D: " << ratio_copy << "Comm: "<< ratio_comm << " ; CPU & other: " << ratio_cpu << std::setprecision(6)<<std::endl;
          if (ratio_gpu<50)
            {
              Cerr << "==============================================================================================" << std::endl;
              Cerr << "[GPU] Warning: Only " << 0.1*int(10*ratio_gpu) << " % of the time calculation is spent on GPU." << std::endl;
              if (ratio_gpu_library==0) Cerr << "[GPU] First add a GPU solver !" << std::endl;
              else
                Cerr << "[GPU] Probably some algorithms used are not ported yet on GPU. Contact TRUST team." << std::endl;
              Cerr << "==============================================================================================" << std::endl;
            }
        }
      // IO part

      // Estimates latency of MPI allreduce
      Process::barrier();
      std::chrono::time_point t1 = std::chrono::high_resolution_clock::now();
      int i = 0;
      for (i = 0; i < 100; i++)
        Process::mp_sum((int)1);
      std::chrono::time_point t2 =  std::chrono::high_resolution_clock::now();
      std::chrono::duration<double> time = t2-t1;
      double allreduce_peak_perf = time.count();
      allreduce_peak_perf = Process::mp_min(allreduce_peak_perf)/100.0;

      // Estimates bandwidth
      double bandwidth = 1.1e30;
      c = access_std_counter(STD_COUNTERS::mpi_sendrecv);
      if (comm_sendrecv_t > 0)
        bandwidth = c->quantity_/ (c->total_time_.count() + DMINFLOAT);

      double max_bandwidth = Process::mp_max(bandwidth);

      // Calcul du temps d'attente du aux synchronisations
      // On prend le temps total de communication et on retranche le temps
      // theorique calcule a partir de allreduce_peak_perf et de la bande passante maxi
      double theoric_comm_time = comm_allreduce_c * allreduce_peak_perf + comm_sendrecv_c / (max_bandwidth + DMINFLOAT);
      // Je suppose que le temps minimum pour realiser les communications sur un proc
      //  depend du processeur qui a le plus de donnees a envoyer:
      theoric_comm_time = Process::mp_max(theoric_comm_time);

      double total_time_avg, total_time_max;
      if(Process::mp_min(nb_ts) >0)
        {
          c = access_std_counter(STD_COUNTERS::timeloop);
        }
      else
        {
          c = access_std_counter(STD_COUNTERS::total_execution_time);
        }
      total_time_avg = Process::mp_sum(c->total_time_.count())/nb_procs;
      total_time_max = Process::mp_max(c->total_time_.count());
      double wait_time = (comm_sendrecv_t+ comm_allreduce_t)- theoric_comm_time;
      double wait_fraction;
      if (total_time_avg == 0)
        wait_fraction = 0.;
      else
        wait_fraction = wait_time / (total_time_avg + DMINFLOAT);
      wait_fraction = 0.1 * floor(wait_fraction * 1000);
      if (wait_fraction < 0.)
        wait_fraction = 0.;
      if (wait_fraction > 100.)
        wait_fraction = 100.;

      double max_wait_fraction = Process::mp_max(wait_fraction);
      double min_wait_fraction = Process::mp_min(wait_fraction);
      double avg_wait_fraction = Process::mp_sum(wait_fraction)/ Process::nproc();

      c = access_std_counter(STD_COUNTERS::IO_EcrireFicPartageBin);
      int debit_seq = c->total_time_.count()>0 ? static_cast<int>(Process::mp_sum(static_cast<double>(c->quantity_)) / (1024 * 1024) / (c->total_time_.count())) : 0;
      c = access_std_counter(STD_COUNTERS::IO_EcrireFicPartageMPIIO);
      int debit_par = c->total_time_.count()>0 ? static_cast<int>(Process::mp_sum(static_cast<double>(c->quantity_)) / (1024 * 1024) / (c->total_time_.count())) : 0;

      if (debit_seq>0 || debit_par>0)
        perfs_IO << std::endl << "------------------------------------------------------------------IO----------------------------------------------------------------------" << std::endl;
      if (debit_seq>0)
        perfs_IO << "Output write sequential (Mo/s) : " << debit_seq << std::endl;
      if (debit_par>0)
        perfs_IO << "Output write parallel (Mo/s) : " << debit_par << std::endl;
      if(min_max_avg_sd_t_q_c_sendrecv_comm[2][1] > 0)
        {
          c=access_std_counter(STD_COUNTERS::petsc_solver);
          if(Process::mp_max(c->count_)>0)
            {
              perfs_IO<< "---------------------------------------------------------------------------------------------------------"<< std::endl;
              perfs_IO<< "Warning: One or several PETSc solvers are used and thus the communication time below are under-estimated."<< std::endl;
              perfs_IO<< "---------------------------------------------------------------------------------------------------------"<< std::endl;
            }
          double fraction = 0.0;
          fraction = (comm_sendrecv_t + comm_allreduce_t)/ (total_time + DMINFLOAT);
          fraction = 0.1 * floor(fraction * 1000);
          if (fraction > 100.)
            fraction = 100.;
          perfs_IO << "Average of the fraction of the time spent in communications between processors: " << fraction  << std::endl;
          fraction = (min_max_avg_sd_t_q_c_sendrecv_comm[0][1] + min_max_avg_sd_t_q_c_allreduce_comm[0][1])/ (total_time_max + DMINFLOAT);
          fraction = 0.1 * floor(fraction * 1000);
          if (fraction > 100.)
            fraction = 100.;
          perfs_IO << "Max of the fraction of the time spent in communications between processors: " << fraction  << std::endl;
          fraction = (min_max_avg_sd_t_q_c_sendrecv_comm[0][0] + min_max_avg_sd_t_q_c_allreduce_comm[0][0])/ (total_time_max + DMINFLOAT);
          fraction = 0.1 * floor(fraction * 1000);
          perfs_IO << "Min of the fraction of the time spent in communications between processors: " << fraction  << std::endl;
          perfs_IO  << "Time of one mpsum measured by an internal bench over 0.1s (network latency): ";
          if (allreduce_peak_perf == 0.)
            perfs_IO << "not measured (total running time too short <10s)" << std::endl;
          else
            perfs_IO << allreduce_peak_perf << " s" << std::endl;
          perfs_IO << "Network maximum bandwidth on all processors: " << max_bandwidth * 1.e-6 << " MB/s"  << std::endl ;
          perfs_IO << "Total network traffic: "<< comm_sendrecv_q * Process::nproc() / nb_ts * 1e-6 << " MB / time step"  << std::endl;
          perfs_IO << "Average message size: "<< comm_sendrecv_q / comm_sendrecv_c* 1e-3 << " kB" << std::endl;
          perfs_IO << "Min waiting time: " << min_wait_fraction << " % of total time"<< std::endl;;
          perfs_IO << "Max waiting time: " << max_wait_fraction<< " % of total time"<< std::endl;;
          perfs_IO << "Avg waiting time: " << avg_wait_fraction<< " % of total time"<< std::endl;;
        }
    }
  // Concatenate stringtreams in order to print the .TU file
  Nom globalTU(Objet_U::nom_du_cas());
  globalTU +="_new.TU";
  std::string root=Sortie_Fichier_base::root;
  Sortie_Fichier_base::root = "";
  EcrFicPartage file(globalTU, mode_append ? (ios::out | ios::app) : (ios::out));
  Sortie_Fichier_base::root = root;
  file << File_header.str();
  file << perfs_TU.str();
  file << perfs_IO.str();
  file << perfs_GPU.str();
  file << captions.str();
  file.syncfile();

  restart_counters();
}





