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
#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>
#include <time.h>
#include <EntreeSortie.h>
#include "Perf_counters.h"
#include <iomanip>
#include <EcrFicPartage.h>
#include <SChaine.h>
#include <communications.h>
#include <TRUSTArray.h>
#include <Comm_Group_MPI.h>


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

  void begin_count_(int counter_level, std::chrono::time_point<std::chrono::high_resolution_clock> t);

  void end_count_(int count_increment, double quantity_increment, std::chrono::time_point<std::chrono::high_resolution_clock> t_stop);
  //check std::chrono

  inline void restart_t_alone(std::chrono::time_point<std::chrono::high_resolution_clock> t) { last_open_time_alone_ = t; }

  inline void set_parent(const Counter * c) { parent_ = c;}

  bool operator==(const Counter& c) const { return (this==&c); }

  bool operator!=(const Counter& c) const  { return !(*this == c);  }

  /*! @brief update variables : avg_time_per_step_ , min_time_per_step_ , max_time_per_step_ , sd_time_per_step_
   *
   */
  void compute_avg_min_max_var_per_step();

  std::array< std::array<double,4> ,2> compute_min_max_avg_sd_();

protected:
  std::string description_;
  int level_;
  std::string family_ ;
  bool is_comm_;
  bool start_at_the_beginning_of_the_time_step_;
  int count_;
  double quantity_;
  Counter* parent_;
  std::chrono::duration<double> total_time_;
  std::chrono::duration<double> time_alone_;   // time when the counter is open minus the time where an counter of lower lvl was open
  std::chrono::duration<double> time_timestep_;
  std::chrono::time_point<std::chrono::high_resolution_clock> last_open_time_;
  std::chrono::time_point<std::chrono::high_resolution_clock> last_open_time_alone_;
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
  start_at_the_beginning_of_the_time_step_ = false;
  count_ = 0;
  quantity_ = 0.0;
  avg_time_per_step_ = 0.0 ;
  min_time_per_step_ = 0.0 ;
  max_time_per_step_ = 0.0 ;
  sd_time_per_step_ = 0.0 ;
  time_timestep_ = std::chrono::duration<double>::zero();
  total_time_ = std::chrono::duration<double>::zero() ;
  time_alone_ = std::chrono::duration<double>::zero() ;
  last_open_time_alone_= std::chrono::time_point<std::chrono::high_resolution_clock>();
  last_open_time_ = std::chrono::time_point<std::chrono::high_resolution_clock>();
  parent_ = nullptr;
}

Counter::Counter()
{
  description_ = "";
  level_ = -10;
  family_ = "";
  is_comm_ = false;
  start_at_the_beginning_of_the_time_step_ = false;
  count_ = 0;
  quantity_ = 0.0;
  avg_time_per_step_ = 0.0 ;
  min_time_per_step_ = 0.0 ;
  max_time_per_step_ = 0.0 ;
  sd_time_per_step_ = 0.0 ;
  time_timestep_ = std::chrono::duration<double>::zero();
  total_time_ = std::chrono::duration<double>::zero() ;
  last_open_time_alone_= std::chrono::time_point<std::chrono::high_resolution_clock>();
  last_open_time_ = std::chrono::time_point<std::chrono::high_resolution_clock>();
  parent_ = nullptr;
}

void Counter::begin_count_(int counter_level, std::chrono::time_point<std::chrono::high_resolution_clock> t)
{
  if (counter_level != level_)
    {
      Cerr<< "You did not specified the expected counter lvl of : :" << description_ <<std::endl;
      level_ = counter_level;
    }
  if (start_at_the_beginning_of_the_time_step_ == false)
    {
      start_at_the_beginning_of_the_time_step_ = true;
    }

  last_open_time_ = t;
  last_open_time_alone_ = t;
  if (parent_!= nullptr)
    {
      if (*parent_->level_ >0)
        *parent_->time_alone_ += std::chrono::duration<double> (t - last_open_time_alone_);
    }
}

void Counter::end_count_(int count_increment, double quantity_increment, std::chrono::time_point<std::chrono::high_resolution_clock> t_stop)
{
  std::chrono::duration<double> t_tot = t_stop-last_open_time_;
  std::chrono::duration<double> t_alone = t_stop - last_open_time_alone_;
  quantity_ += quantity_increment;
  total_time_ += t_tot;
  time_alone_ += t_alone;
  count_ += count_increment;
  if (parent_!= nullptr)
    {
      if (*parent_->level_ >0)
        {
          *parent_-> last_open_time_alone_ = t_stop;
          parent_ = nullptr;
        }
    }
}

std::array< std::array<double,4> ,2> Counter::compute_min_max_avg_sd_()
{
  assert(Process::is_parallel());
  double min,max,avg,sd ;

  min = Process::mp_min(total_time_.count());
  max = Process::mp_max(total_time_.count());
  avg = (double)Process::mp_sum(total_time_.count())/Process::nproc();
  sd = sqrt(Process::mp_sum((total_time_.count()-avg)*(total_time_.count()-avg))/Process::nproc());

  std::array<double,4> min_max_avg_sd_time = {min,max,avg,sd};

  min = Process::mp_min(quantity_);
  max = Process::mp_max(quantity_);
  avg = (double)Process::mp_sum(quantity_)/Process::nproc();
  sd = sqrt(Process::mp_sum((quantity_-avg)*(quantity_-avg))/Process::nproc());

  std::array<double,4> min_max_avg_sd_quantity = {min,max,avg,sd};

  min = Process::mp_min(count_);
  max = Process::mp_max(count_);
  avg = (double)Process::mp_sum(count_)/Process::nproc();
  sd = sqrt(Process::mp_sum((quantity_-avg)*(quantity_-avg))/Process::nproc());

  std::array<double,4> min_max_avg_sd_count = {min,max,avg,sd};

  return {min_max_avg_sd_time,min_max_avg_sd_quantity,min_max_avg_sd_count};
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
{
  Perf_counters::declare_base_counters();
  two_first_steps_elapsed_ = true;
  end_of_cache_=false;
  max_counter_lvl_to_print_ = 10;
  counters_stop_=false;
  last_opened_counter_ = nullptr;
}

void Perf_counters::declare_base_counters()
{
  // Macro counters
  std_counters_[STD_COUNTERS::temps_total_execution_counter_]= Counter::Counter(0, "Total time");
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


void Perf_counters::check_begin(Counter& c, unsigned int counter_lvl, std::chrono::time_point<std::chrono::high_resolution_clock> t)
{
  if (last_opened_counter_ != nullptr)
    {
      Counter & c_parent = *last_opened_counter_;

      if (counter_lvl != *last_opened_counter_->level_ +1)
        Process::exit("The counter you are trying to start does not have the expected level");
      while (c_parent.parent_ !=nullptr)
        {
          if (c == c_parent)
            Process::exit("You are trying to start a counter that is already running");
          &c_parent = *c_parent.parent_;
        }

      c.set_parent(last_opened_counter_);
    }
  last_opened_counter_ = &c;
}

void Perf_counters::check_end(Counter& c, std::chrono::time_point<std::chrono::high_resolution_clock> t)
{
  if (*last_opened_counter_ != c)
    Process::exit("You are not trying to close the last opened counter");
  last_opened_counter_ = c.parent_;
}

/*!
 *
 * @param std_cnt name in the enumerate STD_COUNTERS that corresponds to the counter you try to start
 * @param counter_lvl wanted lvl of the counter you try to start. It has to be equal to the lvl of the last called counter +1. counter_lvl become the new level of the counter you try to start.
 */
void Perf_counters::begin_count(STD_COUNTERS std_cnt, unsigned int counter_lvl)
{
  Counter & c = std_counters_[std_cnt];
  if (end_of_cache_ || !two_first_steps_elapsed_ || c.level_ == 0)
    {
      std::chrono::time_point<std::chrono::high_resolution_clock> t = std::chrono::high_resolution_clock::now();
      Perf_counters::check_begin(c, counter_lvl,t);
      c.begin_count_(counter_lvl,t);
    }
}

/*!
 *
 * @param custom_count_name key of the map (custom_counter_map_str_to_counter_) of the custom counter you try to close
 * @param counter_lvl
 */
void Perf_counters::begin_count(const std::string& custom_count_name, unsigned int counter_lvl)
{
  Counter & c = custom_counter_map_str_to_counter_.at(custom_count_name);
  if (end_of_cache_ || !two_first_steps_elapsed_ || c.level_ == 0)
    {
      std::chrono::time_point<std::chrono::high_resolution_clock> t = std::chrono::high_resolution_clock::now();
      Perf_counters::check_begin(c, counter_lvl,t);
      c.begin_count_(counter_lvl,t);
    }
}

/*! @brief End the count of a counter and update the counter values
 *
 * @param c is the counter to end the count
 * @param count_increment is the count increment. If not specified, then it is equal to 1
 * @param quantity_increment is the increment of custom variable quantity. If not specified, it is set to 0.
 */
void Perf_counters::end_count(const STD_COUNTERS std_cnt, int count_increment, double quantity_increment)
{
  Counter & c = std_counters_[std_cnt];
  if (end_of_cache_ || !two_first_steps_elapsed_ || c.level_ == 0)
    {
      std::chrono::time_point<std::chrono::high_resolution_clock> t = std::chrono::high_resolution_clock::now();
      Perf_counters::check_end(c, t);
      c.end_count_(count_increment, quantity_increment,t);
    }
}


/*! @brief End the count of a counter and update the counter values
 *
 * @param c is the custom counter to end the count
 * @param count_increment is the count increment. If not specified, then it is equal to 1
 * @param quantity_increment is the increment of custom variable quantity. If not specified, it is set to 0.
 */
void Perf_counters::end_count(const std::string& custom_count_name, int count_increment, double quantity_increment)
{
  Counter & c = custom_counter_map_str_to_counter_.at(custom_count_name);
  if (end_of_cache_ || !two_first_steps_elapsed_ || c.level_ == 0)
    {
      std::chrono::time_point<std::chrono::high_resolution_clock> t = std::chrono::high_resolution_clock::now();
      assert(custom_counter_map_str_to_counter_.count(custom_count_name) > 0);
      Perf_counters::check_end(c, t);
      c.end_count_(count_increment, quantity_increment,t);
    }
}


/*! @brief Compute for each counter avg_time_per_step_, min_time_per_step_, max_time_per_step_ and sd_time_per_step_
 *
 * Called at the end of each time step, and only then
 * @param tstep number of time steps elapsed since the start of the computation
 *
 * The value of the
 */
void Perf_counters::compute_avg_min_max_var_per_step(int tstep)
{
  Perf_counters::stop_counters();
  if (!end_of_cache_ && two_first_steps_elapsed_)
    {
      end_of_cache_ = tstep > 2;
    }
  int step = two_first_steps_elapsed_ ? tstep - 2 : tstep;
  auto compute = [](Counter& c)
                            {
    if (c.start_at_the_beginning_of_the_time_step_ == true)
      {
        c.min_time_per_step_ = (c.min_time_per_step_ < c.time_timestep_) ? c.min_time_per_step_ : c.time_timestep_;
        c.max_time_per_step_ = (c.min_time_per_step_ > c.time_timestep_) ? c.min_time_per_step_ : c.time_timestep_;
        c.avg_time_per_step_ = ((step-1)*c.avg_time_per_step_ + c.time_timestep_)/step;
        c.sd_time_per_step_ += (c.time_timestep_* c.time_timestep_ -2*c.time_timestep_* c.avg_time_per_step_ +  c.avg_time_per_step_ *  c.avg_time_per_step_)/step;
        c.sd_time_per_step_ = sqrt(c.sd_time_per_step_);

        if (c.sd_time_per_step_ < 0)
          c.sd_time_per_step_ = 0;

        c.start_at_the_beginning_of_the_time_step_ = false ;
      }
                            };

  for (Counter &c: std_counters_)
    compute(c);
  if (!custom_counter_map_str_to_counter_.empty())
    {
      for (std::pair<const std::string,Counter> & p: custom_counter_map_str_to_counter_)
        compute(p.second);
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
std::string Perf_counters::get_os()
{
  std::string result;
  struct utsname buffer;
  result += buffer.nodename + "__";
  result += buffer.sysname + "__";
  result += buffer.machine + "__";
  result += buffer.release + "__";
  result += buffer.version;
  return result;
}

/*!
 *
 * @return string that contains cpu model and number of proc
 */

std::string Perf_counters::get_cpu()
{
  system("lscpu 2>/dev/null | grep 'Model name' > cpu_detail.txt");
  system("lscpu 2>/dev/null | grep 'Proc' >> cpu_detail.txt");
  std::stringstream cpu_desc;
  cpu_desc << std::ifstream("cpu_detail.txt").rdbuf();
  system("rm cpu_detail.txt");
  return (cpu_desc.str());
}
/*!
 *
 * @return string with gpu model name
 */
std::string Perf_counters::get_gpu()
{
  std::string gpu_description = "No GPU";

#ifdef TRUST_USE_CUDA
  system("nvidia-smi 2>/dev/null | grep NVIDIA > gpu_detail.txt");
#elif TRUST_USE_HIP
  system("rocminfo 2>/dev/null | grep Marketing > gpu_detail.txt");

  std::stringstream gpu_desc;
  gpu_desc << std::ifstream("gpu_detail.txt").rdbuf();
  system("rm gpu_detail.txt");
  gpu_description = gpu_desc.str();
#endif
  return gpu_description;
}
/*!
 *
 * @return string with complete date and time : DD-MM-YYYY -- hour:minute:second
 */
std::string Perf_counters::get_date()
{
  time_t     now = time(0);
  struct tm  tstruct;
  char       date[80];
  tstruct = *localtime(&now);
  strftime(date, sizeof(date), "%d-%m-%Y -- %X", &tstruct);
  return date;
}


/*! @brief Stop all counters, has to be called on every processor simultaneously
 *
 */
void Perf_counters::stop_counters()
{
  std::chrono::time_point<std::chrono::high_resolution_clock> t_stop = std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time_elapsed_before_stop;

  if (last_opened_counter_ != nullptr)
    {
      Counter & c = *last_opened_counter_;
      c.time_alone_ = t_stop -c.last_open_time_alone_;
      while (c.parent_ !=nullptr)
        {
          time_elapsed_before_stop= t_stop -c.last_open_time_;
          c.time_timestep_ += time_elapsed_before_stop;
          c.total_time_ += time_elapsed_before_stop;
          &c = *c.parent_;
        }

    }
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
      Counter & c = *last_opened_counter_;
      c.last_open_time_alone_ = t_restart;
      while (c.parent_ !=nullptr)
        {
          c.last_open_time_ = t_restart;
          &c = *c.parent_;
        }

    }
  counters_stop_=false;
}

void Perf_counters::set_max_counter_lvl_to_print(unsigned int new_max_counter_lvl_to_print)
{
  assert(new_max_counter_lvl_to_print>1);
  max_counter_lvl_to_print_ = new_max_counter_lvl_to_print;
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

void Perf_counters::print_performance_to_csv(std::string message, bool mode_append)
{
  stop_counters();
  assert(message);
  std::stringstream perfs;   ///< Stringstream that contains stats for each processor
  std::stringstream perfs_globales;   ///< Stringstream that contains stats average on the processors : processor number = -1
  std::stringstream File_header;      ///< Stringstream that contains the lines at the start of the file

  long unsigned int length_line = 23; ///< number of item of a line of the _csv.Tu file
  std::vector<int> item_size(length_line,20);   ///< Contains the the width of the printed string, 20 for numbers by default
  std::vector<std::string> line_items(length_line,"");   ///< Contains the data of a line that we want to print in the _csv.TU file.

  std::stringstream tmp_item; ///< Create a temporary stringstream for converting wanted line items in string to construct the line_items vector and therefore
  int nb_procs =  Process::nproc();
  Counter & c = std_counters_[STD_COUNTERS::temps_total_execution_counter_];
  double time_cache = c.total_time_.count();
  c = std_counters_[STD_COUNTERS::timestep_counter_];
  time_cache -= c.total_time_.count();

  /// We specify the width of large items of lines of the _csv.Tu file for making it readable by human
  item_size[0] = 50;
  item_size[2] = 30;
  item_size[3] = 45;

  if ( (Process::je_suis_maitre()) && (message == "Statistiques d'initialisation du calcul") )
    {
      File_header << "# Detailed performance log file. See the associated validation form for an example of data analysis"<< std::endl;
      File_header << "# Date of the computation :" << Perf_counters::get_date() << std::endl;
      File_header << "# OS used :" << Perf_counters::get_os() << std::endl;
      File_header << "# CPU info:" << Perf_counters::get_cpu() << std::endl;
      File_header << "# GPU info:" << Perf_counters::get_gpu() << std::endl;
      File_header << "# Number of processor used = " << nb_procs << std::endl;
      File_header << "# Filling cache took :" << time_cache << std::endl;
      File_header << "# The time was measured by the following method using std::chrono::high_resolution_clock::now()" << std::endl ;
      File_header << "# By default, only averaged statistics on all processor are printed. For accessing the detail per processor, add 'stat_per_proc_perf_log 1' in the data file"<< std::endl;
      File_header << "# Processor number equal to -1 corresponds to the performance of the calculation averaged on the processors during the simulation step" << std::endl;
      File_header << "# If a counter does not belong in any particular family, then counter family is set to (null)" << std::endl;
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
      line_items[7] = "time_(s)";
      line_items[8] = "t_min";
      line_items[9] = "t_max";
      line_items[10] = "t_SD";
      line_items[11] = "count";
      line_items[15] = "time_per_step";
      line_items[16] = "tps_min";
      line_items[17] = "tps_max";
      line_items[18] = "tps_SD";
      line_items[19] = "Quantity";
      line_items[20] = "q_min";
      line_items[21] = "q_max";
      line_items[22] = "q_SD";
      assert(item_size.size()==length_line);
      assert(line_items.size()==item_size.size());
      /// After filling line_items and item_size, we use the function build_line_csv to build the line at the expected format
      build_line_csv(File_header,line_items,item_size);
    }

  /// Check if all of the processors see the same number of counter, if not print an error message in perfs_globales
  int skip_globals = Objet_U::disable_TU;
  int total_nb_of_counters = std_counters_.size() + custom_counter_map_str_to_counter_.size;
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
      skip_globals = 1; ///< If min_nb_of_counters != max_nb_of_counters, aggregated stats are not printed
    }
  Counter& c_time = std_counters_[STD_COUNTERS::temps_total_execution_counter_];
  std::chrono::duration<double> total_time = c_time.total_time_;
  int nb_ts = std_counters_[STD_COUNTERS::timestep_counter_].count_;

  if (nb_ts == 0)
    {
      if (Process::je_suis_maitre())
        {
          perfs_globales << "The computation didn't start" << std::endl;
        }
      skip_globals = 1; ///< If min_nb_of_counters != max_nb_of_counters, aggregated stats are not printed
    }


  //*****************************************************************************************************************************
  //                                            Lambda function to build line
  //*****************************************************************************************************************************
  int level; ///< Level of details of the counter
  bool is_comm; ///< Equal to 1 if the counter is a communication counter, 0 otherwise
  int count;  ///< number of time the counter is open and closed
  double quantity, min_quantity=0, max_quantity=0; ///< A custom quantity which depends on the counter. Used for example to compute the bandwidth
  std::chrono::duration<double> time;
  double percent_time, min_time=0, max_time=0; ///< Percent of the total time used in the method tracked by the counter
  double SD_time=0, SD_quantity=0; ///< the standard dev of all the prev vars
  double avg_time_per_step, min_time_per_step, max_time_per_step, sd_time_per_step;

  auto fill_items = [&](int proc_nb, const std::string& desc, const std::string& familly)
        						                                            								                                                    {
    tmp_item << message; ///< Convert into string the item we want to print in the line, here the overall simulation step
    line_items[0] = tmp_item.str(); ///< Add the item to the vector line_itmes, used to construct the line of the _csv.TU file

    tmp_item.str(""); ///< Empties the temporary stringstream

    tmp_item<< proc_nb;
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
    tmp_item<< time.count();
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

    tmp_item<< count;
    line_items[11] = tmp_item.str();  ///< Number of time the counter was called on the overall simulation step
    tmp_item.str("");

    tmp_item<< avg_time_per_step;
    line_items[15] = tmp_item.str(); ///< Averaged time elapsed by time step for the operation tracked by the counter on the overall simulation step
    tmp_item.str("");

    tmp_item<< min_time_per_step;
    line_items[16] = tmp_item.str(); ///< Minimum time elapsed by time step for the operation tracked by the counter on the overall simulation step
    tmp_item.str("");

    tmp_item<< max_time_per_step;
    line_items[17] = tmp_item.str(); ///< Maximum time elapsed by time step for the operation tracked by the counter on the overall simulation step
    tmp_item.str("");

    tmp_item<< sqrt(sd_time_per_step);
    line_items[18] = tmp_item.str(); ///< Standard Deviation of time elapsed by time step for the operation tracked by the counter on the overall simulation step
    tmp_item.str("");

    tmp_item<< quantity;
    line_items[19] = tmp_item.str(); ///< Custom variable that depends on the counter the overall simulation step
    tmp_item.str("");

    tmp_item<< min_quantity;
    line_items[20] = tmp_item.str(); ///< Detail per proc, so the min, max and SD on proc is equal to 0
    tmp_item.str("");

    tmp_item<< max_quantity;
    line_items[21] = tmp_item.str(); ///< Detail per proc, so the min, max and SD on proc is equal to 0
    tmp_item.str("");

    tmp_item<< SD_quantity;
    line_items[22] = tmp_item.str(); ///< Detail per proc, so the min, max and SD on proc is equal to 0
    tmp_item.str("");
        						                                            								                                                    };

  std::string family_flag = "";
  int size = std::size(std_counters_);
  for (Counter &c: std_counters_)
    {
      if (c.family_.std::string::find("None")==std::string::npos)
        continue;
      if (family_flag.empty())
        family_flag = c.family_;
      if (!family_flag.compare(c.family_))
        {
          time += c.total_time_;
          count += c.count_;
          quantity += c.quantity_;
          avg_time_per_step += c.avg_time_per_step_;
          min_time_per_step = (c.min_time_per_step_<= min_time_per_step) ? c.min_time_per_step_ : min_time_per_step;
          max_time_per_step = (c.min_time_per_step_<= min_time_per_step) ? c.min_time_per_step_ : min_time_per_step;
          sd_time_per_step = c.sd_time_per_step_;
        }

    }

  for (Counter &c: std_counters_)
    {
      if (c.count_ == 0)
        continue;
      time = c.total_time_;
      count = c.count_;
      quantity = c.quantity_;
      avg_time_per_step = c.avg_time_per_step_;
      min_time_per_step = c.min_time_per_step_;
      max_time_per_step = c.max_time_per_step_;
      sd_time_per_step = c.sd_time_per_step_;

      if (Objet_U::stat_per_proc_perf_log && Process::is_parallel())
        {
          fill_items(Process::me(), c.description_);
          build_line_csv(perfs,line_items,item_size);  ///< Build the line of the stats associated on the counter i for a single proc
        }

      if (Process::je_suis_maitre())
        {
          std::array< std::array<double,4> ,2> table = c.compute_min_max_avg_sd_();
          time = table[0][2];
          min_time = table[0][0];
          max_time = table[0][1];
          SD_time = table[0][3];
          quantity = table[1][2];
          min_quantity = table[1][0];
          max_quantity = table[1][1];
          SD_quantity = table[1][3];
          if (! skip_globals )
            {
              fill_items(-1, c.description_);
              build_line_csv(perfs_globales,line_items,item_size);
            }
        }
      level = 0; ///< Level of details of the counter
      is_comm = false; ///< Equal to 1 if the counter is a communication counter, 0 otherwise
      count = 0;  ///< number of time the counter is open and closed
      quantity =0, min_quantity=0, max_quantity=0; ///< A custom quantity which depends on the counter. Used for example to compute the bandwidth
      time = std::chrono::duration<double>::zero();
      percent_time=0, min_time=0, max_time=0; ///< Percent of the total time used in the method tracked by the counter
      SD_time=0, SD_quantity=0; ///< the standard dev of all the prev vars
      avg_time_per_step = 0, min_time_per_step = 0, max_time_per_step = 0, sd_time_per_step = 0;
    }

  for (std::pair< const std::string , Counter> & pair : custom_counter_map_str_to_counter_)
    {
      Counter& c = pair.second;
      if (c.count_ == 0)
        continue;
      time = c.total_time_;
      count = c.count_;
      quantity = c.quantity_;
      avg_time_per_step = c.avg_time_per_step_;
      min_time_per_step = c.min_time_per_step_;
      max_time_per_step = c.max_time_per_step_;
      sd_time_per_step = c.sd_time_per_step_;

      if (Objet_U::stat_per_proc_perf_log && Process::is_parallel())
        {
          fill_items(Process::me(), c.description_);
          build_line_csv(perfs,line_items,item_size);  ///< Build the line of the stats associated on the counter i for a single proc
        }

      if (Process::is_parallel() && Process::je_suis_maitre())
        {
          std::array< std::array<double,4> ,2> table = c.compute_min_max_avg_sd_();
          time = table[0][2];
          min_time = table[0][0];
          max_time = table[0][1];
          SD_time = table[0][3];
          quantity = table[1][2];
          min_quantity = table[1][0];
          max_quantity = table[1][1];
          SD_quantity = table[1][3];
          if (! skip_globals )
            {
              fill_items(-1, c.description_);
              build_line_csv(perfs_globales,line_items,item_size);
            }
        }
      level = 0; ///< Level of details of the counter
      is_comm = false; ///< Equal to 1 if the counter is a communication counter, 0 otherwise
      count = 0;  ///< number of time the counter is open and closed
      quantity =0, min_quantity=0, max_quantity=0; ///< A custom quantity which depends on the counter. Used for example to compute the bandwidth
      time = std::chrono::duration<double>::zero();
      percent_time=0, min_time=0, max_time=0; ///< Percent of the total time used in the method tracked by the counter
      SD_time=0, SD_quantity=0; ///< the standard dev of all the prev vars
      avg_time_per_step = 0, min_time_per_step = 0, max_time_per_step = 0, sd_time_per_step = 0;

    }



  if (Process::je_suis_maitre())
    {
      Nom CSV(Objet_U::nom_du_cas());
      CSV+="_csv.TU";
      std::string root=Sortie_Fichier_base::root;
      Sortie_Fichier_base::root = "";
      EcrFicPartage file(CSV, mode_append ? (ios::out | ios::app) : (ios::out));
      Sortie_Fichier_base::root = root;
      file << File_header.str();
      file << perfs_globales.str();
      file << perfs.str();
      file.syncfile();
    }

  Perf_counters::restart_counters();
}

void Perf_counters::print_global_TU(std::string message, bool mode_append)
{
  Perf_counters::stop_counters();
  Nom TU(Objet_U::nom_du_cas());
  TU += ".TU";
  std::stringstream perfs_TU;   ///< Stringstream that contains algomerated stats that will be printed in the .TU
  std::stringstream perfs_GPU;   ///< Stringstream that contains algomerated stats that will be printed in the .TU
  std::stringstream perfs_IO;   ///< Stringstream that contains algomerated stats that will be printed in the .TU
  std::stringstream File_header;      ///< Stringstream that contains the File header

  int nb_procs =  Process::nproc();
  double comm_allreduce_q = 0.0,comm_sendrecv_q = 0.0, comm_allreduce_t = 0.0, comm_sendrecv_t = 0.0;
  int comm_allreduce_c = 0,comm_sendrecv_c = 0;
  std::array< std::array<double,4> ,3> min_max_avg_sd_t_q_c_sendrecv_comm ;
  std::array< std::array<double,4> ,3> min_max_avg_sd_t_q_c_allreduce_comm ;
  std::array< std::array<double,4> ,3> min_max_avg_sd_t_q_c_echange_espace_virtuel ;
  std::array< std::array<double,4> ,3> min_max_avg_sd_t_q_c_sendrecv ;
  Counter &c = std_counters_[STD_COUNTERS::timestep_counter_];
  int nb_ts = two_first_steps_elapsed_ ? Process::mp_max(c.count_) - 2 : Process::mp_max(c.count_) ;
  double time_dt = c.total_time_.count();
  if (nb_ts <= 0)
    two_first_steps_elapsed_ ? Process::exit("No time step computed or at least none after filing the cache") : Process::exit("No time step computed");

  double total_time = Process::mp_max(c.total_time_.count());

  auto write_globalTU_line = [&] (std::string str, std::stringstream & line)
                                    {
    if (c.count_>0)
      {
        double t_c = c.total_time_.count();
        int count = c.count_;
        if (Process::je_suis_maitre())
          {
            t_c = Process::mp_max(t_c);
            count = Process::mp_max(count);
          }
        line << "Whose " << str << " ";
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

  if (Process::je_suis_maitre())
    {
      nb_ts = Process::mp_max(nb_ts);
      total_time = Process::mp_max(total_time);
      if (nb_ts <= 0)
        Process::exit("No time step was computed");
      File_header << "# Global performance file #"<< std::endl;
      if (two_first_steps_elapsed_)
        File_header << "The two first time steps are not accounted for the computation of the statistics"<< std::endl;
      File_header <<  std::endl;
      File_header << "Date :" << Perf_counters::get_date() << std::endl;
      File_header << "OS :" << Perf_counters::get_os() << std::endl;
      File_header << "CPU :" << Perf_counters::get_cpu() << std::endl;
      File_header << "GPU :" << Perf_counters::get_gpu() << std::endl;
      File_header << "Number of processor used = " << nb_procs << std::endl << std::endl << std::endl ;
      perfs_TU << message << std::endl<< std::endl;
      c = std_counters_[STD_COUNTERS::temps_total_execution_counter_];
      perfs_TU << "Temps total                       " <<  c.total_time_.count() << std::endl << std::endl;

    }
  for (Counter & c_com : std_counters_)
    {
      if (c_com.is_comm_)
        {
          if (c_com.family_=="MPI_allreduce")
            {
              comm_allreduce_q += c_com.quantity_;
              comm_allreduce_t += c_com.total_time_.count();
              comm_allreduce_c += c_com.count_;
            }
          if (c_com.family_=="MPI_sendrecv")
            {
              comm_sendrecv_q += c_com.quantity_;
              comm_sendrecv_t += c_com.total_time_.count();
              comm_sendrecv_c += c_com.count_;
            }
        }

    }
  for (std::pair< const std::string , Counter> & pair : custom_counter_map_str_to_counter_)
    {
      Counter& c_com = pair.second;
      if (c_com.is_comm_)
        {
          if (c_com.family_=="MPI_allreduce")
            {
              comm_allreduce_q += c_com.quantity_;
              comm_allreduce_t += c_com.total_time_.count();
              comm_allreduce_c += c_com.count_;
            }
          if (c_com.family_=="MPI_sendrecv")
            {
              comm_sendrecv_q += c_com.quantity_;
              comm_sendrecv_t += c_com.total_time_.count();
              comm_sendrecv_c += c_com.count_;
            }
        }

    }


  if (GET_COMM_DETAILS && Process::is_parallel() && Process::je_suis_maitre())
    {
      min_max_avg_sd_t_q_c_allreduce_comm = compute_min_max_avg_sd(comm_allreduce_t,comm_allreduce_q,comm_allreduce_c);
      min_max_avg_sd_t_q_c_sendrecv_comm =  compute_min_max_avg_sd(comm_sendrecv_t,comm_sendrecv_q,comm_sendrecv_c);
      c = std_counters_[STD_COUNTERS::mpi_sendrecv_counter_];
      min_max_avg_sd_t_q_c_sendrecv = c.compute_min_max_avg_sd_();
      c = std_counters_[STD_COUNTERS::echange_vect_counter_];
      min_max_avg_sd_t_q_c_echange_espace_virtuel = c.compute_min_max_avg_sd_();
      comm_allreduce_t = min_max_avg_sd_t_q_c_allreduce_comm[0][2];
      comm_allreduce_q = min_max_avg_sd_t_q_c_allreduce_comm[1][2];
      comm_allreduce_c = min_max_avg_sd_t_q_c_allreduce_comm[2][2];
      comm_sendrecv_t = min_max_avg_sd_t_q_c_sendrecv_comm[0][2];
      comm_sendrecv_q = min_max_avg_sd_t_q_c_sendrecv_comm[1][2];
      comm_sendrecv_c = min_max_avg_sd_t_q_c_sendrecv_comm[2][2];
    }


  // Estimate latency of MPI allreduce
  Process::barrier();
  std::chrono::time_point t1 = std::chrono::high_resolution_clock::now();
  int i = 0;
  for (i = 0; i < 100; i++)
    Process::mp_sum((int)1);
  std::chrono::time_point t2 =  std::chrono::high_resolution_clock::now();
  std::chrono::duration<double> time = t2-t1;
  double allreduce_peak_perf = time.count();
  allreduce_peak_perf = Process::mp_min(allreduce_peak_perf)/100.0;


  // Estimate bandwidth
  double bandwidth = 1.1e30;
  c = std_counters_[STD_COUNTERS::mpi_sendrecv_counter_];
  if (comm_sendrecv_t > 0)
    bandwidth = c.quantity_/ (c.total_time_.count() + DMINFLOAT);

  double max_bandwidth = Process::mp_max(bandwidth);

  // Calcul du temps d'attente du aux synchronisations
  // On prend le temps total de communication et on retranche le temps
  // theorique calcule a partir de allreduce_peak_perf et de la bande passante maxi
  double theoric_comm_time = comm_allreduce_c * allreduce_peak_perf
      + comm_sendrecv_c / (max_bandwidth + DMINFLOAT);
  // Je suppose que le temps minimum pour realiser les communications sur un proc
  //  depend du processeur qui a le plus de donnees a envoyer:
  theoric_comm_time = Process::mp_max(theoric_comm_time);



  double total_time_avg, total_time_max;
  if(two_first_steps_elapsed_)
    {
      c = std_counters_[STD_COUNTERS::timestep_counter_];
    }
  else
    {
      c = std_counters_[STD_COUNTERS::temps_total_execution_counter_];
    }
  total_time_avg = Process::mp_sum(c.total_time_.count())/nb_procs;
  total_time_max = Process::mp_max(c.total_time_.count());
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

  c = std_counters_[STD_COUNTERS::sauvegarde_counter_];
  double total_quantity = Process::mp_sum(c.quantity_);
  c = std_counters_[STD_COUNTERS::IO_EcrireFicPartageBin_counter_];
  int debit_seq = c.total_time_.count()>0 ? (int) (Process::mp_sum(c.quantity_) / (1024 * 1024) / c.total_time_.count()) : 0;
  c = std_counters_[STD_COUNTERS::IO_EcrireFicPartageMPIIO_counter_];
  int debit_par = c.total_time_.count()>0 ? (int) (Process::mp_sum(c.quantity_) / (1024 * 1024) / c.total_time_.count()) : 0;


  if (Process::je_suis_maitre())
    {
      c = std_counters_[STD_COUNTERS::probleme_fluide_];
      write_globalTU_line("probleme thermohydraulique  ", perfs_TU);
      c = std_counters_[STD_COUNTERS::probleme_combustible_];
      write_globalTU_line("probleme combustible        ", perfs_TU);
      perfs_TU << std::endl;
      perfs_TU << "Timesteps                         " << nb_ts << std::endl;
      c = std_counters_[STD_COUNTERS::timestep_counter_];
      perfs_TU << "Secondes / pas de temps           " << c.time_timestep_.count() << "    ;    Standard deviation between time steps: " << c.sd_time_per_step_ << std::endl;
      c = std_counters_[STD_COUNTERS::solv_sys_counter_];
      write_globalTU_line("solveurs Ax=B               ",  perfs_TU);
      c = std_counters_[STD_COUNTERS::diffusion_implicite_counter_];
      write_globalTU_line("solveur diffusion_implicite ", perfs_TU);
      c = std_counters_[STD_COUNTERS::assemblage_sys_counter_];
      write_globalTU_line("assemblage matrice_implicite", perfs_TU);
      c = std_counters_[STD_COUNTERS::mettre_a_jour_counter_];
      write_globalTU_line("mettre_a_jour               ", perfs_TU);
      c = std_counters_[STD_COUNTERS::update_vars_counter_];
      write_globalTU_line("update_vars                 ", perfs_TU);
      c = std_counters_[STD_COUNTERS::update_fields_counter_];
      write_globalTU_line("update_fields               ", perfs_TU);
      c = std_counters_[STD_COUNTERS::convection_counter_];
      write_globalTU_line("operateurs convection       ", perfs_TU);
      c = std_counters_[STD_COUNTERS::diffusion_counter_];
      write_globalTU_line("operateurs diffusion        ", perfs_TU);
      c = std_counters_[STD_COUNTERS::decay_counter_];
      write_globalTU_line("operateurs decroissance     ", perfs_TU);
      c = std_counters_[STD_COUNTERS::gradient_counter_];
      write_globalTU_line("operateurs gradient         ", perfs_TU);
      c = std_counters_[STD_COUNTERS::divergence_counter_];
      write_globalTU_line("operateurs divergence       ", perfs_TU);
      c = std_counters_[STD_COUNTERS::source_counter_];
      write_globalTU_line("operateurs source           ", perfs_TU);
      c = std_counters_[STD_COUNTERS::postraitement_counter_];
      write_globalTU_line("operations postraitement    ", perfs_TU);
      c = std_counters_[STD_COUNTERS::dt_counter_];
      write_globalTU_line("calcul dt                   ", perfs_TU);
      c = std_counters_[STD_COUNTERS::nut_counter_];
      write_globalTU_line("modele turbulence           ", perfs_TU);
      c = std_counters_[STD_COUNTERS::sauvegarde_counter_];
      write_globalTU_line("operations sauvegarde       ", perfs_TU);
      c =  std_counters_[STD_COUNTERS::divers_counter_];
      write_globalTU_line("calcul divers               ", perfs_TU);
    }
  c = std_counters_[STD_COUNTERS::echange_vect_counter_];
  if (Process::mp_max(c.count_)>0)
    {
      perfs_TU << "Maximum number of virtual exchanges per time steps :" <<  Process::mp_max(c.count_) << std::endl;
    }
  if (min_max_avg_sd_t_q_c_allreduce_comm[2][1]>=0)
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
  c = std_counters_[STD_COUNTERS::solv_sys_counter_];
  int tmp = Process::mp_max(c.count_);
  if (tmp > 0)
    {
      perfs_TU << "Number of call of the solver per time step         "<< tmp / nb_ts << std::endl;
      double avg_time = Process::mp_max(c.total_time_.count()) / tmp;
      perfs_TU << "Average time of the resolution per call               " << avg_time << std::endl;
      perfs_TU << "Average number of iteration per call             "<< Process::mp_max(c.quantity_) / tmp << std::endl;
    }

  c = std_counters_[STD_COUNTERS::sauvegarde_counter_];
  tmp = Process::mp_max(c.count_);
  if (tmp>0)
    {
      perfs_TU << "Number of back-up :" << tmp << std::endl;
      perfs_TU << "Average amount of data per back-up (Mo) : " << total_quantity / (tmp *1024*1024) << std::endl;
    }


  // GPU part of the TU :

  auto compute_percent_and_write = [&] (const std::string str)
            {
    double max_time = Process::mp_max(c.total_time_.count());
    double percent = 100*max_time/nb_ts;
    double calls = Process::mp_max(c.count_)/nb_ts;
    double t_ts = max_time/nb_ts;
    double bw = c.quantity_*(1021*1024*Process::mp_max(c.total_time_.count()));
    perfs_GPU << str << " time per time step:"<< t_ts <<" ; number of calls per time steps: " << calls << " ; percent of total time: " << percent << " ; bandwidth: " << bw <<  std::endl;
    return percent;
            };

  c = std_counters_[STD_COUNTERS::gpu_copytodevice_counter_];
  tmp = Process::mp_max(c.count_);
  if (tmp>0)
    {
      perfs_GPU << std::endl << perfs_TU << "-------------------------------------------------------------------GPU--------------------------------------------------------------------" << std::endl;
      perfs_GPU << "This test case was run using GPU. The used GPU is the following one : " << Perf_counters::get_gpu() << std::endl;
      perfs_GPU << "GPU statistics per time step (experimental):" << std::endl;
      c = std_counters_[STD_COUNTERS::gpu_library_counter_];
      double ratio_gpu_library = compute_percent_and_write("Libraries: ");
      c = std_counters_[STD_COUNTERS::gpu_kernel_counter_];
      double ratio_gpu_kernel = compute_percent_and_write("Kernels: ");
      double ratio_gpu = ratio_gpu_kernel+ratio_gpu_library;
      c = std_counters_[STD_COUNTERS::gpu_copytodevice_counter_];
      double ratio_copy = compute_percent_and_write("Copy host to device: ");
      c = std_counters_[STD_COUNTERS::gpu_copyfromdevice_counter_];
      ratio_copy += compute_percent_and_write("Copy device to host: ");
      c = std_counters_[STD_COUNTERS::timestep_counter_];
      double ratio_comm = 100.0 * (comm_sendrecv_t+ comm_allreduce_t)/c.max_time_per_step_;
      double ratio_gpu = 100 - ratio_copy - ratio_gpu - ratio_comm;
      perfs_GPU << "Comm: "<< ratio_comm << " ; CPU & other: " << ratio_gpu << std::endl;
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

#ifdef _OPENMP_TARGET
  else
    {
      Cerr << "============================================================================================" << std::endl;
      Cerr << "[GPU] Warning: Don't use this binary! Slower for your calculation, not at all ported on GPU." << std::endl;
      Cerr << "============================================================================================" << std::endl;
    }
#endif

  // IO part

  perfs_IO << std::endl << "------------------------------------------------------------------IO----------------------------------------------------------------------" << std::endl;
  if (debit_seq>0) perfs_IO << "Output write sequential (Mo/s) : " << debit_seq << std::endl;
  if (debit_par>0) perfs_IO << "Output write parallel (Mo/s) : " << debit_par << std::endl;
  if(comm_sendrecv_c > 0)
    {
      c=std_counters_[STD_COUNTERS::solv_sys_petsc_counter_];
      if(Process::mp_max(c.count_)>0)
        {
          perfs_IO<< "---------------------------------------------------------------------------------------------------------"<< std::endl;
          perfs_IO<< "Warning: One or several PETSc solvers are used and thus the communication time below are under-estimated."<< std::endl;
          perfs_IO<< "---------------------------------------------------------------------------------------------------------"<< std::endl;
        }

    }


  // Concatenate stringtreams in order to print the .TU file
  if(Process::je_suis_maitre())
    {
      Nom globalTU(Objet_U::nom_du_cas());
      globalTU +=".TU";
      std::string root=Sortie_Fichier_base::root;
      Sortie_Fichier_base::root = "";
      EcrFicPartage file(globalTU, mode_append ? (ios::out | ios::app) : (ios::out));
      Sortie_Fichier_base::root = root;
      file << File_header.str();
      file << perfs_TU.str();
      file << perfs_IO.str();
      file << perfs_GPU.str();
      file.syncfile();
    }

  Perf_counters::restart_counters();
}

inline std::array< std::array<double,4> ,3> compute_min_max_avg_sd(double time, double quantity, int count)
{
  assert(Process::is_parallel());
  double min,max,avg,sd ;

  min = Process::mp_min(time);
  max = Process::mp_max(time);
  avg = (double)Process::mp_sum(time)/Process::nproc();
  sd = sqrt(Process::mp_sum((time-avg)*(time-avg))/Process::nproc());

  std::array<double,4> min_max_avg_sd_time = {min,max,avg,sd};

  min = Process::mp_min(quantity);
  max = Process::mp_max(quantity);
  avg = (double)Process::mp_sum(quantity)/Process::nproc();
  sd = sqrt(Process::mp_sum((quantity-avg)*(quantity-avg))/Process::nproc());

  std::array<double,4> min_max_avg_sd_quantity = {min,max,avg,sd};

  min = Process::mp_min(count);
  max = Process::mp_max(count);
  avg = (double)Process::mp_sum(count)/Process::nproc();
  sd = sqrt(Process::mp_sum((quantity-avg)*(quantity-avg))/Process::nproc());

  std::array<double,4> min_max_avg_sd_count = {min,max,avg,sd};

  return {min_max_avg_sd_time,min_max_avg_sd_quantity,min_max_avg_sd_count};
}



