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

#include <Perf_counters.h>

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
#include <iomanip>
#include <EcrFicPartage.h>
#include <memory>
#include <iomanip>
#ifdef TRUST_USE_CUDA
// See https://nvidia.github.io/NVTX/
// See https://stackoverflow.com/questions/23230003/something-between-func-and-pretty-function/29856690#29856690
#include <nvtx3/nvToolsExt.h>
#endif

#define MINFLOAT 1.e-34  // smth small!

/**************************************************************************************************************************
 *
 *      Introduction of the class counter that described the behavior of a single counter in TRUST
 *
 **************************************************************************************************************************/

struct Counter
{

  using clock = std::chrono::high_resolution_clock;
  using time_point = std::chrono::time_point<clock>;
  using duration = std::chrono::duration<double>;
  inline time_point now()
  {
    return clock::now();
  }

  Counter(int counter_level, std::string counter_name, std::string counter_family = "None", bool is_comm = false, bool is_gpu = false);

  void begin_count_(int counter_level, time_point t);

  void end_count_(int count_increment, long int quantity_increment, time_point t_stop);

  inline void set_parent(Counter * parent_counter) { parent_ = parent_counter;}

  inline double get_time_() const {return total_time_.count();}

  inline bool running_() const { return is_running_; }

  /*! @brief update variables : avg_time_per_step_ , min_time_per_step_ , max_time_per_step_ , sd_time_per_step_
   *
   */
  void compute_avg_min_max_var_per_step();

  std::array< std::array<double,4> ,4> compute_min_max_avg_sd_() const;

  void reset();

  const std::string description_;
  int level_;
  const std::string family_ ;
  const bool is_comm_;
  const bool is_gpu_;
  int count_;
  long int quantity_;
  Counter* parent_;
  duration total_time_;
  duration time_alone_;   // time when the counter is open minus the time where an counter of lower lvl was open
  duration time_ts_;  // total time tracked during the current time_steps
  time_point open_time_ts_;
  time_point last_open_time_;
  time_point last_open_time_alone_;
  double avg_time_per_step_;
  double min_time_per_step_;
  double max_time_per_step_;
  double sd_time_per_step_;
  bool is_running_;
}
;

Counter::Counter(int counter_level, std::string counter_name, std::string counter_family , bool is_comm, bool is_gpu)
  :description_(counter_name), level_(counter_level), family_(counter_family), is_comm_(is_comm), is_gpu_(is_gpu),  count_(0),
   quantity_( 0), parent_( nullptr),
   total_time_(duration::zero()),
   time_alone_(duration::zero()),
   time_ts_(duration::zero()),
   open_time_ts_( time_point()),
   last_open_time_( time_point()),
   last_open_time_alone_(time_point()),
   avg_time_per_step_( 0.0 ),
   min_time_per_step_( 0.0 ),
   max_time_per_step_( 0.0 ),
   sd_time_per_step_( 0.0 ),
   is_running_( false)
{
}

void Counter::begin_count_(int counter_level, time_point t)
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
      parent_->time_alone_ +=duration (t - last_open_time_alone_);
      parent_->last_open_time_alone_ =  time_point();
    }
#ifdef TRUST_USE_CUDA
  if (!is_comm_)
    nvtxRangePush(description_.c_str());
#endif

}

void Counter::end_count_(int count_increment, long int quantity_increment, time_point t_stop)
{
  if (!is_running_)
    Process::exit("Last open_time was not properly set"+ description_);
  duration t_tot = t_stop-last_open_time_;
  duration t_alone = t_stop - last_open_time_alone_;
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
  last_open_time_ = time_point();
  last_open_time_alone_ = time_point();
  open_time_ts_ = time_point();
#ifdef TRUST_USE_CUDA
  if (!is_comm_) nvtxRangePop();
#endif
}

std::array< std::array<double,4> ,4> Counter::compute_min_max_avg_sd_() const
{
  assert(Process::is_parallel());
  double qty,cnt,min,max,avg,sd ;
  qty=static_cast<double>(quantity_);
  cnt = static_cast<double>(count_);

  auto l_compute =[&min, &max,&avg,&sd] (double value)
  {
    min = Process::mp_min(value);
    max = Process::mp_max(value);
    avg = Process::mp_sum(value)/Process::nproc();
    sd = sqrt(Process::mp_sum((value-avg)*(value-avg))/Process::nproc());
    std::array<double,4> result = {min,max,avg,sd};
    return result;
  };

  std::array<double,4> min_max_avg_sd_time = l_compute(total_time_.count());

  std::array<double,4> min_max_avg_sd_quantity = l_compute(qty);

  std::array<double,4> min_max_avg_sd_count_ = l_compute(cnt);

  std::array<double,4> min_max_avg_sd_time_alone_ = l_compute(time_alone_.count());

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
  open_time_ts_ = now();
  total_time_ =duration::zero() ;
  last_open_time_alone_= now();
  last_open_time_ = now();
  time_alone_=duration::zero() ;;   // time when the counter is open minus the time where an counter of lower lvl was open
  time_ts_=duration::zero() ;;  // total time tracked during the current time_steps
}



/**************************************************************************************************************************
 *
 * 					Declaration of the class Perf_counters using Pimpl idiom
 *
 **************************************************************************************************************************/
/*! @Brief declare all standard counters of TRUST inside an array
 *
 */
class Perf_counters::Impl
{
public:
  using clock = std::chrono::high_resolution_clock;
  using time_point = std::chrono::time_point<clock>;
  using duration = std::chrono::duration<double>;
  inline time_point now()
  {
    return clock::now();
  }
  Impl();
  ~Impl();
  void create_custom_counter_impl(std::string counter_description, int counter_level, std::string counter_family, bool is_comm, bool is_gpu);
  void begin_count_impl(const STD_COUNTERS& std_cnt, int counter_lvl);
  void begin_count_impl(const std::string& custom_count_name, int counter_lvl);
  void end_count_impl(const STD_COUNTERS& std_cnt, int count_increment, long int quantity_increment);
  void end_count_impl(const std::string& custom_count_name, int count_increment, long int quantity_increment);
  void stop_counters_impl();
  void restart_counters_impl();
  void reset_counters_impl();
  void set_time_steps_elapsed_impl(int time_step_elapsed);
  double get_computation_time_impl();
  double get_total_time_impl(const STD_COUNTERS& name);
  double get_total_time_impl(const std::string& name);
  double get_time_since_last_open_impl(const STD_COUNTERS& name);
  double get_time_since_last_open_impl(const std::string& name);
  void start_timeloop_impl();
  void end_timeloop_impl();
  void start_time_step_impl();
  void end_time_step_impl(unsigned int tstep);
  int get_last_opened_counter_level_impl() const;
  void print_TU_files_impl(const std::string& message);
  void start_gpu_timer_impl();
  void stop_gpu_timer_impl();
  double compute_gpu_time_impl();
  bool is_gpu_verbose_on_impl() const ;
  void set_gpu_verbose_impl(bool on) ;
  bool get_init_device_impl() const ;
  void set_init_device_impl(bool init) ;
  bool get_gpu_timer_impl() const ;
  void add_to_gpu_timer_counter_impl(int to_add=1) ;
  int get_gpu_timer_counter_impl() const ;
  bool get_use_gpu_impl() const {return use_gpu_;}
  bool get_gpu_fence_impl() const {return gpu_fence_;}
  void set_gpu_fence_impl(bool fence) {gpu_fence_=fence;}
  bool running_impl(const STD_COUNTERS name) { return get_counter(name).running_(); }

private:
  Counter& get_counter(const STD_COUNTERS name) ;
  Counter& get_counter(const std::string name);
  void check_begin(Counter& c, int counter_lvl, time_point t);
  void check_end(Counter& c, time_point t);
  double compute_allreduce_peak();
  std::string get_os() const;
  std::string get_cpu() const;
  std::string get_gpu() const;
  std::string get_date() const;
  void print_performance_to_csv(const std::string& message);
  void print_global_TU(const std::string& message);
  std::array<std::unique_ptr<Counter>, static_cast<int>(STD_COUNTERS::NB_OF_STD_COUNTER)> std_counters_; ///< Store standard counters using unique_ptr
  std::map<std::string, std::unique_ptr<Counter>> custom_counter_map_str_to_counter_;  ///< Store custom counters using unique_ptr
  bool end_cache_=false;                 ///< A flag used to know if the two first time steps are over or not
  bool time_loop_=false;                 ///< A flag used to know if we are inside the time loop
  bool counters_stop_=false;             ///< A flag used to know if the counters are paused or not
  //int counter_lvl_to_print_=1;       ///< Counter level that you want to be printed in the global_TU
  duration computation_time_=duration::zero();      ///< Used to compute the total time of the simulation.
  duration time_skipped_ts_=duration::zero();       ///< the duration in seconds of the cache. If cache is too long, use function set_three_first_steps_elapsed in oder to include the stats of the cache in your stats
  Counter* last_opened_counter_=nullptr;   ///< pointer to the last opened counter. Each counter has a parent attribute, which also give the pointer of the counter open before them.
  unsigned int nb_steps_elapsed_=0;  ///< By default, we consider that the two first time steps are used to file the cache, so they are not taken into account in the stats.
  int total_nb_backup_=0;
  double total_data_exchange_per_backup_=0.;
  bool gpu_verbose_ =false;
  bool init_device_ = false;
  bool gpu_timer_ = false;
  bool use_gpu_=false;
  bool gpu_fence_=true;
  time_point gpu_timer_start_;
  int gpu_timer_count_=0;
  int max_str_length_=118;
};
Perf_counters::Impl::~Impl()=default;

Perf_counters::Impl::Impl()
{
  // Initialize standard counters
  std_counters_[static_cast<int>(STD_COUNTERS::total_execution_time)] = std::make_unique<Counter>(-1, "Total time");
  std_counters_[static_cast<int>(STD_COUNTERS::computation_start_up)] = std::make_unique<Counter>(0, "Computation start-up");
  std_counters_[static_cast<int>(STD_COUNTERS::timeloop)] = std::make_unique<Counter>(0, "Time loop");
  std_counters_[static_cast<int>(STD_COUNTERS::backup_file)] = std::make_unique<Counter>(0, "Back-up operations");
  std_counters_[static_cast<int>(STD_COUNTERS::system_solver)] = std::make_unique<Counter>(1, "Linear solver resolutions Ax=B");
  std_counters_[static_cast<int>(STD_COUNTERS::petsc_solver)] = std::make_unique<Counter>(2, "Petsc solver");
  std_counters_[static_cast<int>(STD_COUNTERS::implicit_diffusion)] = std::make_unique<Counter>(1, "Solver for implicit diffusion:");
  std_counters_[static_cast<int>(STD_COUNTERS::compute_dt)] = std::make_unique<Counter>(1, "Computation of the time step dt");
  std_counters_[static_cast<int>(STD_COUNTERS::turbulent_viscosity)] = std::make_unique<Counter>(1, "Turbulence model::update");
  std_counters_[static_cast<int>(STD_COUNTERS::convection)] = std::make_unique<Counter>(1, "Convection operator");
  std_counters_[static_cast<int>(STD_COUNTERS::diffusion)] = std::make_unique<Counter>(1, "Diffusion operator");
  std_counters_[static_cast<int>(STD_COUNTERS::gradient)] = std::make_unique<Counter>(1, "Gradient operator");
  std_counters_[static_cast<int>(STD_COUNTERS::divergence)] = std::make_unique<Counter>(1, "Divergence operator");
  std_counters_[static_cast<int>(STD_COUNTERS::rhs)] = std::make_unique<Counter>(1, "Source_terms");
  std_counters_[static_cast<int>(STD_COUNTERS::postreatment)] = std::make_unique<Counter>(1, "Post-treatment operations");
  std_counters_[static_cast<int>(STD_COUNTERS::restart)] = std::make_unique<Counter>(1, "Read file for restart");
  std_counters_[static_cast<int>(STD_COUNTERS::matrix_assembly)] = std::make_unique<Counter>(1, "Nb matrix assembly for implicit scheme:");
  std_counters_[static_cast<int>(STD_COUNTERS::update_variables)] = std::make_unique<Counter>(1, "Update ::mettre_a_jour");
  // MPI communication counters
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_sendrecv)] = std::make_unique<Counter>(2, "MPI_send_recv", "MPI_sendrecv", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_send)] = std::make_unique<Counter>(2, "MPI_send", "MPI_sendrecv", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_recv)] = std::make_unique<Counter>(2, "MPI_recv", "MPI_sendrecv", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_bcast)] = std::make_unique<Counter>(2, "MPI_broadcast", "MPI_sendrecv", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_alltoall)] = std::make_unique<Counter>(2, "MPI_alltoall", "MPI_sendrecv", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_allgather)] = std::make_unique<Counter>(2, "MPI_allgather", "MPI_sendrecv", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_gather)] = std::make_unique<Counter>(2, "MPI_gather", "MPI_sendrecv", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_partialsum)] = std::make_unique<Counter>(2, "MPI_partialsum", "MPI_allreduce", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_sumdouble)] = std::make_unique<Counter>(2, "MPI_sumdouble", "MPI_allreduce", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_mindouble)] = std::make_unique<Counter>(2, "MPI_mindouble", "MPI_allreduce", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_maxdouble)] = std::make_unique<Counter>(2, "MPI_maxdouble", "MPI_allreduce", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_sumfloat)] = std::make_unique<Counter>(2, "MPI_sumfloat", "MPI_allreduce", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_minfloat)] = std::make_unique<Counter>(2, "MPI_minfloat", "MPI_allreduce", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_maxfloat)] = std::make_unique<Counter>(2, "MPI_maxfloat", "MPI_allreduce", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_sumint)] = std::make_unique<Counter>(2, "MPI_sumint", "MPI_allreduce", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_minint)] = std::make_unique<Counter>(2, "MPI_minint", "MPI_allreduce", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_maxint)] = std::make_unique<Counter>(2, "MPI_maxint", "MPI_allreduce", true);
  std_counters_[static_cast<int>(STD_COUNTERS::mpi_barrier)] = std::make_unique<Counter>(2, "MPI_barrier", "MPI_allreduce", true);

  // GPU counters
  std_counters_[static_cast<int>(STD_COUNTERS::gpu_library)] = std::make_unique<Counter>(2, "GPU_library", "GPU_library", false, true);
  std_counters_[static_cast<int>(STD_COUNTERS::gpu_kernel)] = std::make_unique<Counter>(2, "GPU_kernel", "GPU_kernel", false, true);
  std_counters_[static_cast<int>(STD_COUNTERS::gpu_copytodevice)] = std::make_unique<Counter>(2, "GPU_copyToDevice", "GPU_copy", false, true);
  std_counters_[static_cast<int>(STD_COUNTERS::gpu_copyfromdevice)] = std::make_unique<Counter>(2, "GPU_copyFromDevice","GPU_copy",false, true);
  std_counters_[static_cast<int>(STD_COUNTERS::gpu_malloc_free)] = std::make_unique<Counter>(2, "GPU_allocations"   ,"GPU_alloc",false,true);
  // Scatter
  std_counters_[static_cast<int>(STD_COUNTERS::interprete_scatter)] = std::make_unique<Counter>(2, "Scatter_interprete", "None", true,false);
  std_counters_[static_cast<int>(STD_COUNTERS::virtual_swap)] = std::make_unique<Counter>(2, "DoubleVect/IntVect::virtual_swap", "None", true);
  std_counters_[static_cast<int>(STD_COUNTERS::read_scatter)] = std::make_unique<Counter>(2, "Scatter::read_domaine", "None", true);
  //Parallel meshing
  std_counters_[static_cast<int>(STD_COUNTERS::parallel_meshing)] = std::make_unique<Counter>(0, "Parallel meshing");
  //IO
  std_counters_[static_cast<int>(STD_COUNTERS::IO_EcrireFicPartageBin)] = std::make_unique<Counter>(2, "write", "IO");
  std_counters_[static_cast<int>(STD_COUNTERS::IO_EcrireFicPartageMPIIO)] = std::make_unique<Counter>(2,"MPI_File_write_all", "IO");
  if (nb_steps_elapsed_==0)
    end_cache_=true;
#ifdef TRUST_USE_GPU
  use_gpu_=true;
#endif
}
///////  Private methods of Pimpl

Counter& Perf_counters::Impl::get_counter(const STD_COUNTERS name)
{
  return *std_counters_[static_cast<int>(name)];
}

Counter& Perf_counters::Impl::get_counter(std::string cust_counter_desc)
{
  if (custom_counter_map_str_to_counter_.count(cust_counter_desc)==0)
    Process::exit("You are trying to find a custom counter that does not exists");
  return *custom_counter_map_str_to_counter_.at(cust_counter_desc);
}

void Perf_counters::Impl::check_begin(Counter& c, int counter_lvl, time_point t)
{
  if (last_opened_counter_ != nullptr)
    {
      if (c.is_running_)
        Process::exit("The counter that you are trying to start is already running:" + c.description_);
      int expected_lvl = last_opened_counter_->level_ +1;
      if (c.is_comm_)
        {
          counter_lvl=expected_lvl;
        }
      if (counter_lvl != expected_lvl)
        {
          std::ostringstream error_msg ;
          error_msg << "The counter you are trying to start does not have the expected level, counter running: " << last_opened_counter_->description_ << " counter that you try to open: " << c.description_  << " ; expected level: "<<  expected_lvl << std::endl ;
          Process::exit(error_msg.str());
        }
      if (time_loop_)
        {
          c.open_time_ts_ = t;
        }
      c.set_parent(last_opened_counter_);
    }
  last_opened_counter_ =&c;
}

void Perf_counters::Impl::check_end(Counter& c, time_point t)
{
  if (!c.is_running_ || last_opened_counter_==nullptr)
    Process::exit("You are trying to close a counter that is not running: " + c.description_);
  if (last_opened_counter_ != &c)
    {
      std::string error_msg = "The counter you are trying to close is not the last opened, counter: " + c.description_;
      Process::exit(error_msg);
    }
  if (time_loop_)
    {
      c.time_ts_ += t - c.open_time_ts_;
      c.open_time_ts_ = time_point();
    }
  last_opened_counter_ = c.parent_;
  if (c.is_comm_)
    {
      c.level_=2;
    }
}

double Perf_counters::Impl::compute_allreduce_peak()
{
  Process::barrier();
  time_point t1 = now();
  int i = 0;
  for (i = 0; i < 100; i++)
    Process::mp_sum(static_cast<int>(1));
  time_point t2 =  now();
  duration time = t2-t1;
  double allreduce_peak_perf = time.count();
  return Process::mp_min(allreduce_peak_perf)/100.0;
}


std::string delete_blank_spaces(std::string str)
{
  std::string result;
  int space_count =0;
  for (char ch:str)
    {
      if (ch==' ')
        {
          space_count ++;
          if (space_count <=2)
            result += ch;
        }
      else
        {
          space_count=0;
          result+=ch;
        }
    }
  return(result);
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
std::string Perf_counters::Impl::get_os() const
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
  result = delete_blank_spaces(result);
  return result.substr(0,max_str_length_);
}

/*!
 *
 * @return string that contains cpu model and number of proc
 */

std::string Perf_counters::Impl::get_cpu() const
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
  str= delete_blank_spaces(str);
  return (str.substr(0,max_str_length_));
}
/*!
 *
 * @return string with gpu model name
 */
std::string Perf_counters::Impl::get_gpu() const
{
  std::string gpu_description = "No GPU was used for the computation";
#ifdef TRUST_USE_CUDA
  int result = std::system("nvidia-smi 2>/dev/null | grep NVIDIA > gpu_detail.txt");
  if (result !=0)
    Cerr << "Bash command in Perf_counters::get_gpu failed"<<std::endl;
  std::ostringstream gpu_desc;
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
  std::ostringstream gpu_desc;
  gpu_desc << std::ifstream("gpu_detail.txt").rdbuf();
  result_ = std::system("rm gpu_detail.txt");
  if (result_ !=0)
    Cerr << "Bash command in Perf_counters::get_gpu failed"<<std::endl;
  gpu_description = gpu_desc.str();
#endif
  gpu_description=delete_blank_spaces(gpu_description);
  return gpu_description.substr(0,max_str_length_);
}

/*!
 *
 * @return string with complete date and time : DD-MM-YYYY -- hour:minute:second
 */
std::string Perf_counters::Impl::get_date() const
{
  time_t now = time(0);
  std::ostringstream date;
  struct tm tstruct = *localtime(&now);
  date<< std::put_time(&tstruct, "%d-%m-%Y -- %X");
  std::string result = date.str();
  result = delete_blank_spaces(result);
  return (result.substr(0,max_str_length_));
}

static void build_line_csv(std::ostringstream& lines, const std::array<std::string,24>& line_items, const std::array<int,24>& item_size)
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

void clean_stringstream(std::ostringstream& lines)
{
  lines.str("");
  lines.clear();
}

void Perf_counters::Impl::print_performance_to_csv(const std::string& message)
{
  assert(!message.empty());
  std::ostringstream perfs;   ///< Stringstream that contains stats for each processor
  std::ostringstream perfs_globales;   ///< Stringstream that contains stats average on the processors : processor number = -1
  std::ostringstream file_header;      ///< Stringstream that contains the lines at the start of the file

  const int length_line = 24; ///< number of item of a line of the _csv.Tu file
  std::array<int,length_line> item_size; ///< Contains the the width of the printed string, 20 for numbers by default
  for (int& j:item_size)
    j=20;
  std::array<std::string,length_line> line_items;   ///< Contains the data of a line that we want to print in the _csv.TU file.
  for (std::string& str :line_items)
    str="";
  std::ostringstream tmp_item; ///< Create a temporary ostringstream for converting wanted line items in string to construct the line_items vector and therefore
  int nb_procs =  Process::nproc();

  /// We specify the width of large items of lines of the _csv.Tu file for making it readable by human
  item_size[0] = 50;
  item_size[2] = 40;
  item_size[3] = 45;

  if ( (Process::je_suis_maitre()) && (message == "Computation start-up statistics") )
    {
      file_header << "# Detailed performance log file for case: " << Objet_U::nom_du_cas()<<". See the associated validation form for an example of data analysis"<< std::endl;
      file_header << "# Date of the computation:     " << get_date() << std::endl;
      file_header << "# OS used:     " << get_os() << std::endl;
      file_header << "# CPU info:     " << get_cpu() << std::endl;
      file_header << "# GPU info:     " << get_gpu() << std::endl;
      file_header << "# Number of processor used = " << nb_procs << std::endl;
      file_header << "# The time was measured by the following method using std::chrono::high_resolution_clock::now()" << std::endl ;
      file_header << "# By default, only averaged statistics on all processor are printed. For accessing the detail per processor, add 'stat_per_proc_perf_log 1' in the data file"<< std::endl;
      file_header << "# Processor number equal to -1 corresponds to the performance of the calculation averaged on the processors during the simulation step" << std::endl;
      file_header << "# If a counter does not belong in any particular family, then counter family is set to None" << std::endl;
      file_header << "# Count means the number of time the counter is called during the overall calculation step." << std::endl;
      file_header << "# Min, max and SD accounts respectively for the minimum, maximum and Standard Deviation of the quantity of the previous row." << std::endl;
      file_header << "# Quantity is a custom variable that depends on the counter. It is used to compute bandwidth for communication counters for example. See the table at the end of the introduction on statistics in TRUST form for more details." << std::endl;
      file_header << "#" << std::endl << "#" << std::endl;
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
      build_line_csv(file_header,line_items,item_size);
    }

  /// Check if all of the processors see the same number of counter, if not print an error message in perfs_globales
  bool skip_globals = false;
  int total_nb_of_counters = static_cast<int>(std_counters_.size()) + static_cast<int>(custom_counter_map_str_to_counter_.size());
  int min_total_nb_of_counters = total_nb_of_counters;
  int max_total_nb_of_counters = total_nb_of_counters;

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
  Counter& c_time = get_counter(STD_COUNTERS::timeloop);
  int nb_ts = c_time.count_- nb_steps_elapsed_;

  if (time_loop_ && nb_ts <= 0)
    {
      if (Process::je_suis_maitre())
        {
          perfs_globales << "The computation is shorter than cache" << std::endl;
        }
      skip_globals = true; ///< If min_nb_of_counters != max_nb_of_counters, aggregated stats are not printed
    }

  int level; ///< Level of details of the counter
  bool is_comm; ///< Equal to 1 if the counter is a communication counter, 0 otherwise
  int count;  ///< number of time the counter is open and closed
  long int quantity, min_quantity=0, max_quantity=0; ///< A custom quantity which depends on the counter. Used for example to compute the bandwidth
  double time,time_alone,min_time_alone=0.,max_time_alone=0.,SD_time_alone=0.0;
  double percent_time=0., min_time=0.0, max_time=0.0; ///< Percent of the total time used in the method tracked by the counter
  double SD_time=0.0, SD_quantity=0.0; ///< the standard dev of all the prev vars
  double avg_time_per_step=0., min_time_per_step=0., max_time_per_step=0., sd_time_per_step=0.;

  auto fill_items = [&](int proc_number, const std::string desc, const std::string familly)
  {
    tmp_item << message; ///< Convert into string the item we want to print in the line, here the overall simulation step
    line_items[0] = tmp_item.str(); ///< Add the item to the vector line_itmes, used to construct the line of the _csv.TU file

    tmp_item.str(""); ///< Empties the temporary ostringstream

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
    line_items[12] = tmp_item.str(); ///< Detail per proc, so the min, max, avg and SD on proc is equal to 0
    tmp_item.str("");

    tmp_item<< max_time_alone;
    line_items[13] = tmp_item.str(); ///< Detail per proc, so the min, max, avg and SD on proc is equal to 0
    tmp_item.str("");

    tmp_item<< SD_time_alone;
    line_items[14] = tmp_item.str(); ///< Detail per proc, so the min, max, avg and SD on proc is equal to 0
    tmp_item.str("");

    tmp_item<< count;
    line_items[15] = tmp_item.str();  ///< Number of time the counter was called on the overall simulation step
    tmp_item.str("");

    tmp_item<< avg_time_per_step;
    line_items[16] = tmp_item.str(); ///< Averaged time elapsed by time step for the operation tracked by the counter on the overall simulation step
    tmp_item.str("");

    tmp_item<< min_time_per_step;
    line_items[17] = tmp_item.str(); ///< Minimum time elapsed by time step for the operation tracked by the counter on the overall simulation step
    tmp_item.str("");

    tmp_item<< max_time_per_step;
    line_items[18] = tmp_item.str(); ///< Maximum time elapsed by time step for the operation tracked by the counter on the overall simulation step
    tmp_item.str("");

    tmp_item<< sqrt(sd_time_per_step);
    line_items[19] = tmp_item.str(); ///< Standard Deviation of time elapsed by time step for the operation tracked by the counter on the overall simulation step
    tmp_item.str("");

    tmp_item<< quantity;
    line_items[20] = tmp_item.str(); ///< Custom variable that depends on the counter the overall simulation step
    tmp_item.str("");

    tmp_item<< min_quantity;
    line_items[21] = tmp_item.str(); ///< Detail per proc, so the min, max and SD on proc is equal to 0
    tmp_item.str("");

    tmp_item<< max_quantity;
    line_items[22] = tmp_item.str(); ///< Detail per proc, so the min, max and SD on proc is equal to 0
    tmp_item.str("");

    tmp_item<< SD_quantity;
    line_items[23] = tmp_item.str(); ///< Detail per proc, so the min, max and SD on proc is equal to 0
    tmp_item.str("");
  };

  auto extract_stats = [&](const Counter & c_lambda)
  {
    level = c_lambda.level_; ///< Level of details of the counter
    is_comm = c_lambda.is_comm_; ///< Equal to 1 if the counter is a communication counter, 0 otherwise
    time = c_lambda.total_time_.count();
    time_alone = c_lambda.time_alone_.count();
    count = c_lambda.count_;
    quantity = c_lambda.quantity_;
    avg_time_per_step = c_lambda.avg_time_per_step_;
    min_time_per_step = c_lambda.min_time_per_step_;
    max_time_per_step = c_lambda.max_time_per_step_;
    sd_time_per_step = c_lambda.sd_time_per_step_;
    min_time = 0.;
    max_time = 0.;
    SD_time = 0.;
    min_quantity = 0.;
    max_quantity = 0.;
    SD_quantity = 0.;
    min_time_alone = 0.;
    max_time_alone = 0.;
    SD_time_alone = 0.;
    if (Objet_U::stat_per_proc_perf_log || !Process::is_parallel())
      {
        fill_items(Process::me(),c_lambda.description_, c_lambda.family_);
        build_line_csv(perfs,line_items,item_size);  ///< Build the line of the stats associated on the counter i for a single proc
      }
    if (Process::is_parallel())
      {
        std::array< std::array<double,4> ,4> table = c_lambda.compute_min_max_avg_sd_();
        if (Process::je_suis_maitre())
          {
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
                fill_items(-1,c_lambda.description_,c_lambda.family_);
                build_line_csv(perfs_globales,line_items,item_size);
              }
          }
      }
  };

  for (int i =0 ; i< static_cast<int>(STD_COUNTERS::NB_OF_STD_COUNTER); i++)
    {
      Counter& c_std = *std_counters_[i];
      extract_stats(c_std);
    }

  for (const auto & pair : custom_counter_map_str_to_counter_)
    {
      if (pair.second!=nullptr)
        extract_stats(*pair.second);
    }

  Nom CSV(Objet_U::nom_du_cas());
  CSV+="_csv.TU";
  auto printing_mode = ios::app;
  if (message=="Computation start-up statistics")
    printing_mode= ios::out;
  if (Process::je_suis_maitre() || !Process::is_parallel())
    {
      SFichier file(CSV, printing_mode);
      file << file_header.str();
      file << perfs_globales.str();
      file <<perfs.str();
    }
  else if (Objet_U::stat_per_proc_perf_log && Process::is_parallel())
    {
      EcrFicPartage file(CSV,printing_mode);
      file << perfs.str();
      file.syncfile();
    }
  clean_stringstream(file_header);
  clean_stringstream(perfs);
  clean_stringstream(perfs_globales);
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
inline std::array< std::array<double,4> ,3> compute_min_max_avg_sd(double& time, long int& quantity, int& count)
{
  double qty,cnt,min,max,avg,sd ;
  qty=static_cast<double>(quantity);
  cnt = static_cast<double>(count);
  if( Process::is_parallel())
    {
      auto l_compute =[&min, &max,&avg,&sd] (double value)
      {
        min = Process::mp_min(value);
        max = Process::mp_max(value);
        avg = Process::mp_sum(value)/Process::nproc();
        sd = sqrt(Process::mp_sum((value-avg)*(value-avg))/Process::nproc());
        std::array<double,4> result = {min,max,avg,sd};
        return result;
      };
      std::array<double,4> min_max_avg_sd_time = l_compute(time);
      time =avg;
      std::array<double,4> min_max_avg_sd_quantity = l_compute(qty);
      std::array<double,4> min_max_avg_sd_count = l_compute(cnt);
      count = static_cast<int>(std::floor (avg));
      return {min_max_avg_sd_time,min_max_avg_sd_quantity,min_max_avg_sd_count};
    }
  else
    {
      std::array<double,4> min_max_avg_sd_time = {time,0.,0.,0.};
      std::array<double,4> min_max_avg_sd_quantity = {qty,0.,0.,0.};
      std::array<double,4> min_max_avg_sd_count = {cnt,0.,0.,0.};
      return {min_max_avg_sd_time,min_max_avg_sd_quantity,min_max_avg_sd_count};
    }
}
/*!
 *
 * @param message
 * @param mode_append
 */
void Perf_counters::Impl::print_global_TU(const std::string& message)
{
  assert(!message.empty());
  std::ostringstream perfs_TU;   ///< Stringstream that contains algomerated stats that will be printed in the .TU
  std::ostringstream perfs_GPU;   ///< Stringstream that contains algomerated stats that will be printed in the .TU
  std::ostringstream perfs_IO;   ///< Stringstream that contains algomerated stats that will be printed in the .TU
  std::ostringstream captions;
  std::ostringstream file_header;      ///< Stringstream that contains the File header
  const int counter_description_width = 40;
  const int time_per_step_width= 15;
  const int percent_loop_time_width=11;
  const int count_per_ts_width=15;
  const int level_width=5;
  const int bandwith_width= 10;
  const int tabular_custom_line_width= counter_description_width+3+time_per_step_width+3+percent_loop_time_width+3+count_per_ts_width+3+level_width;
  const int cpu_line_width=counter_description_width+3+time_per_step_width+3+percent_loop_time_width+3+count_per_ts_width;
  const int gpu_line_width=counter_description_width+3+time_per_step_width+3+percent_loop_time_width+3+count_per_ts_width+3+bandwith_width+4;
  const int number_width=15;
  const int text_width =cpu_line_width-count_per_ts_width;
  const int header_txt_width = 10;
  const int message_width = static_cast<int>(message.length());
  const std::string separator = " | ";
  const std::string line_sep_cpu(max_str_length_,'~');
  const std::string line_sep_tabular(cpu_line_width,'-');
  const std::string line_sep_tabular_custom(tabular_custom_line_width,'-');
  const std::string line_sep_gpu(gpu_line_width,'-');
  int nb_procs =  Process::nproc();
  double allreduce_peak_perf = compute_allreduce_peak();
  double comm_allreduce_t = 0.0, comm_sendrecv_t = 0.0;
  long int comm_allreduce_q = 0.0,comm_sendrecv_q = 0.0;
  int comm_allreduce_c = 0,comm_sendrecv_c = 0;
  std::array< std::array<double,4> ,3> min_max_avg_sd_t_q_c_sendrecv_comm ;
  for (std::array<double,4>& arr: min_max_avg_sd_t_q_c_sendrecv_comm)
    for (double & d : arr)
      d=-1;
  std::array< std::array<double,4> ,3> min_max_avg_sd_t_q_c_allreduce_comm = min_max_avg_sd_t_q_c_sendrecv_comm ;
  Counter& c_timeloop = get_counter(STD_COUNTERS::timeloop);
  int nb_ts = c_timeloop.count_;
  nb_ts = std::max(nb_ts,1);
  double time_tl=c_timeloop.total_time_.count();
  double total_untracked_time_ts=c_timeloop.time_alone_.count();
  Counter& c_total_time = get_counter(STD_COUNTERS::total_execution_time);
  Counter& c_mpi_sendrecv = get_counter(STD_COUNTERS::mpi_sendrecv);
  Counter& c_virtual_swap = get_counter(STD_COUNTERS::virtual_swap);
  Counter& c_system_solver= get_counter(STD_COUNTERS::system_solver);
  Counter& c_backup = get_counter(STD_COUNTERS::backup_file);
  Counter& c_todevice = get_counter(STD_COUNTERS::gpu_copytodevice);
  Counter& c_gpu_l = get_counter(STD_COUNTERS::gpu_library);
  Counter& c_gpu_k = get_counter(STD_COUNTERS::gpu_kernel);
  Counter& c_fromdevice = get_counter(STD_COUNTERS::gpu_copyfromdevice);
  Counter& c_io_seq = get_counter(STD_COUNTERS::IO_EcrireFicPartageBin);
  Counter& c_io_par = get_counter(STD_COUNTERS::IO_EcrireFicPartageMPIIO);
  Counter& c_petsc=get_counter(STD_COUNTERS::petsc_solver);
  Counter& c_allocfree=get_counter(STD_COUNTERS::gpu_malloc_free);
  int petcs_count = Process::mp_max(c_petsc.count_);
  int copy_to_device_count = Process::mp_max(c_todevice.count_);
  int max_virtual_swap_c = Process::mp_max(c_virtual_swap.count_);
  double avg_solv_time = Process::mp_max(c_system_solver.total_time_.count());
  double total_time = c_total_time.total_time_.count();
  double total_untracked_time=c_total_time.time_alone_.count();
  double total_quantity = Process::mp_sum(static_cast<double>(c_backup.quantity_));
  int max_nb_backup = Process::mp_max(c_backup.count_);
  double total_comm_time=0.;
  int solver_calls=  Process::mp_max(c_system_solver.count_);
  nb_ts = Process::mp_max(nb_ts);
  double cpu_time=0.;
  double nb_it_per_solver_calls= solver_calls>0 ? static_cast<double>(Process::mp_max(c_system_solver.quantity_)) /solver_calls : static_cast<double>(Process::mp_max(c_system_solver.quantity_)) ;
  if (max_nb_backup>0)
    {
      total_nb_backup_ += c_backup.count_;
      total_data_exchange_per_backup_ += total_quantity / (max_nb_backup *1024*1024);
    }

  auto write_globalTU_line = [&] (const Counter& c_to_print_,std::ostringstream & line)
  {
    if (c_to_print_.count_>0  ) //&& c_to_print_.level_==counter_lvl_to_print_)
      {
        double t_c = c_to_print_.total_time_.count();
        int count = c_to_print_.count_;
        line << std::left <<std::setw(counter_description_width) << c_to_print_.description_ <<separator ;
        double t = nb_ts>0 ? c_to_print_.avg_time_per_step_ : t_c;
        line << std::left << std::setw(time_per_step_width) <<t << separator << std::setprecision(3) << std::setw(percent_loop_time_width) << t_c/total_time*100 ;
        if (nb_ts>0)
          {
            double n = static_cast<double>(count)/nb_ts;
            line << separator <<std::left << std::setw(count_per_ts_width) << std::round(n) << std::setprecision(7);
          }
        line << std::endl;
      }
  };

  auto write_globalTU_line_custom_counters = [&] (const Counter& c_to_print_,std::ostringstream & line)
  {
    if (c_to_print_.count_>0)
      {
        double t_c = c_to_print_.total_time_.count();
        int count = c_to_print_.count_;
        line << std::left <<std::setw(counter_description_width) << "Custom_counter::"+c_to_print_.description_ <<separator ;
        double t = nb_ts>0 ? t_c/nb_ts : t_c;
        line << std::left << std::setw(time_per_step_width) <<t << separator << std::setprecision(3) << std::setw(percent_loop_time_width) << t_c/total_time*100 ;
        if (nb_ts>0)
          {
            double n = static_cast<double>(count)/nb_ts;
            line << separator <<std::left << std::setw(count_per_ts_width) << std::round(n) << std::setprecision(7);
          }
        line << separator <<std::setw(level_width) << c_to_print_.level_ << std::endl;
      }
  };

  for (int i =0 ; i< static_cast<int>(STD_COUNTERS::NB_OF_STD_COUNTER); i++)
    {
      Counter& c_com = *std_counters_[i];
      {
        if (c_com.count_>0)
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
                total_comm_time += c_com.time_alone_.count();
              }
            else if (!c_com.is_gpu_)
              cpu_time+=c_com.time_alone_.count();
          }
      }
    }
  for (const auto & pair : custom_counter_map_str_to_counter_)
    {
      Counter& c_com = *pair.second;
      if (c_com.is_comm_)
        {
          if (c_com.family_=="MPI_allreduce")
            {
              comm_allreduce_q += c_com.quantity_;
              comm_allreduce_t += c_com.total_time_.count();
              comm_allreduce_c += c_com.count_;
            }
          if (pair.second->family_=="MPI_sendrecv")
            {
              comm_sendrecv_q += c_com.quantity_;
              comm_sendrecv_t += c_com.total_time_.count();
              comm_sendrecv_c += c_com.count_;
            }
          total_comm_time += c_com.time_alone_.count();
        }
      else if (!c_com.is_gpu_)
        cpu_time+=c_com.time_alone_.count();
    }
  // IO part
  // Estimates bandwidth
  double bandwidth = 1.1e30;
  if (c_mpi_sendrecv.total_time_.count()>0)
    bandwidth = static_cast<double>(c_mpi_sendrecv.quantity_)/ (c_mpi_sendrecv.total_time_.count() + MINFLOAT);

  double max_bandwidth = Process::mp_max(bandwidth);
  // Compute wait time due to synch
  // We take the total communication time and we substract the theoretical tume computed with allreduce_peak_perf and max bandwidth
  double theoric_comm_time = 0.0;
  if(max_bandwidth)
    theoric_comm_time = static_cast<double>(comm_allreduce_c) * allreduce_peak_perf + static_cast<double>(comm_sendrecv_c) / max_bandwidth;
  // Je suppose que le temps minimum pour realiser les communications sur un proc
  //  depend du processeur qui a le plus de donnees a envoyer:
  theoric_comm_time = Process::mp_max(theoric_comm_time);
  double total_time_avg=0.0, total_time_max=0.0;
  if(nb_ts >0)
    {
      total_time_avg = Process::mp_sum(c_timeloop.total_time_.count())/nb_procs;
      total_time_max = Process::mp_max(c_timeloop.total_time_.count());
    }
  else
    {
      total_time_avg = Process::mp_sum(c_total_time.total_time_.count())/nb_procs;
      total_time_max = Process::mp_max(c_total_time.total_time_.count());
    }
  double wait_time = (comm_sendrecv_t+ comm_allreduce_t)- theoric_comm_time;
  double wait_fraction;
  if (total_time_avg == 0)
    wait_fraction = 0.;
  else
    wait_fraction = wait_time / (total_time_avg + MINFLOAT);
  wait_fraction = 0.1 * floor(wait_fraction * 1000);
  if (wait_fraction < 0.)
    wait_fraction = 0.;
  if (wait_fraction > 100.)
    wait_fraction = 100.;

  double max_wait_fraction = Process::mp_max(wait_fraction);
  double min_wait_fraction = Process::mp_min(wait_fraction);
  double avg_wait_fraction = Process::mp_sum(wait_fraction)/ nb_procs;
  int debit_seq = 0;
  int debit_par = 0;
  double com_time_seq = Process::mp_max(c_io_seq.total_time_.count());
  double com_time_par = Process::mp_max(c_io_par.total_time_.count());
  if (com_time_seq>0)
    debit_seq = static_cast<int>(std::floor(static_cast<double>(Process::mp_sum(c_io_seq.quantity_)) / (1024 * 1024) /com_time_seq));
  if (com_time_par>0)
    debit_par = static_cast<int>(std::floor(static_cast<double>(Process::mp_sum(c_io_par.quantity_)) / (1024 * 1024) /com_time_par));
  if(Process::is_parallel())
    {
      min_max_avg_sd_t_q_c_allreduce_comm = compute_min_max_avg_sd(comm_allreduce_t,comm_allreduce_q,comm_allreduce_c);
      min_max_avg_sd_t_q_c_sendrecv_comm =  compute_min_max_avg_sd(comm_sendrecv_t,comm_sendrecv_q,comm_sendrecv_c);
    }
  if (Process::je_suis_maitre())
    {
      std::string spaces;
      if (message == "Computation start-up statistics")
        {
          spaces.assign((max_str_length_-27)/2,' ');
          file_header <<  spaces <<"# Global performance file #"<< std::endl;
          file_header <<  std::endl;
          file_header << "This is the global file for tracking performance in TRUST. It stores aggregated quantities." <<std::endl;
          file_header << "More detailed statistics can be found in the "<< Objet_U::nom_du_cas() <<"_csv.TU file" <<std::endl;
          file_header << "A jupyter notebook giving detailed information about performance measurement can be found in:" << std::endl << " $TRUST_ROOT/Validation/Rapports_automatiques/Verification/HowTo/" <<std::endl;
          file_header << "For time loop, only level 1 counters statistics are printed in this file by default" << std::endl;
          file_header <<"Time is given in seconds"<< std::endl <<std::endl;
          file_header << line_sep_cpu << std::endl;
          spaces.assign((max_str_length_-26)/2,' ');
          file_header << spaces <<"Context of the computation"<< std::endl;
          file_header << line_sep_cpu << std::endl;
          file_header << std::left << std::setw(header_txt_width)<< "Date:" << get_date() << std::endl;
          file_header << std::left << std::setw(header_txt_width)<< "OS:" << get_os() << std::endl;
          file_header << std::left << std::setw(header_txt_width) << "CPU:" << get_cpu() << std::endl;
          file_header << std::left << std::setw(header_txt_width) << "GPU:" << get_gpu() << std::endl;
          file_header << std::left << std::setw(header_txt_width) << "Nb procs: " << nb_procs << std::endl << std::endl;
          file_header << line_sep_cpu << std::endl;
          spaces.assign((max_str_length_-message_width)/2,' ');
          file_header  << spaces<<message << std::endl;
          file_header << line_sep_cpu << std::endl;
          file_header << std::left << std::setw(text_width)<<"Total time of the start-up: " <<  std::left <<std::setw(number_width) << c_total_time.total_time_.count() << std::endl;
          file_header << std::left << std::setw(text_width)<< "Percent of untracked time during computation start-up: "<<  std::left <<std::setw(number_width) << total_untracked_time/total_time << std::endl;
        }
      else if (message == "Time loop statistics")
        {
          file_header << line_sep_cpu << std::endl;
          spaces.assign((max_str_length_-message_width)/2,' ');
          file_header  << spaces<<message << std::endl;
          file_header << line_sep_cpu << std::endl;
          if (nb_ts <= 0)
            {
              cerr<<"No time step after cache filling was computed"<<std::endl;
              return;
            }
          if (nb_steps_elapsed_>0)
            {
              if (nb_steps_elapsed_>1)
                file_header << "The " <<  nb_steps_elapsed_<< " first time steps are not accounted for the computation of the time loop statistics"<< std::endl;
              else
                file_header << "The first time step is not accounted for the computation of the time loop statistics"<< std::endl;
            }
          file_header <<  std::left <<std::setw(text_width)<< "Total time: "<<  std::left <<std::setw(number_width) << c_total_time.total_time_.count() << std::endl;
          file_header <<  std::left <<std::setw(text_width) <<  "Number of time steps: " <<  std::left <<std::setw(number_width) << nb_ts << std::endl;
          file_header <<  std::left <<std::setw(text_width) << "Average time per time step: " <<  std::left <<std::setw(number_width) << c_total_time.total_time_.count()/nb_ts << endl;
          file_header <<  std::left <<std::setw(text_width) << "Standard deviation between time steps: " <<  std::left <<std::setw(number_width) << c_timeloop.sd_time_per_step_ << std::endl;
          file_header <<  std::left <<std::setw(text_width) << "Time elapsed in the skipped time steps: " <<  std::left <<std::setw(number_width) << time_skipped_ts_.count() <<std::endl << std::endl;
          if (Process::is_parallel())
            file_header <<  std::left <<std::setw(text_width) << "Percent of total time tracked by communication counters:" <<  std::left <<std::setw(number_width) << 100* total_comm_time / total_time << std::endl;
        }
      else if (message == "Post-resolution statistics")
        {
          file_header << line_sep_cpu << std::endl;
          spaces.assign((max_str_length_-message_width)/2,' ');
          file_header  << spaces<<message << std::endl;
          file_header << line_sep_cpu << std::endl;
          file_header <<  std::left <<std::setw(text_width) << "Time of the post-resolution: " <<  std::left <<std::setw(number_width) <<  c_total_time.total_time_.count() << std::endl;
          captions <<  std::endl;
          captions << std::left << std::setw(text_width) << "Total time for the whole computation" << std::left <<std::setw(number_width) << computation_time_.count()<< std::endl<< std::endl;
          if (Process::is_parallel())
            {
              captions << "Max waiting time big    => probably due to a bad partitioning" << std::endl;
              captions << "Communications > 30%    => too many processors or network too slow" << std::endl;
              captions << std::endl;
            }
        }
      else
        Process::exit("You are trying to get stats of an unknown computation step");

      if(message == "Time loop statistics")
        {
          perfs_TU<<std::endl;
          perfs_TU << std::left <<std::setw(counter_description_width) << "Standard counter description" << separator << std::setw(time_per_step_width) << "Time/step" << separator << std::setw(percent_loop_time_width) << "% loop time" << separator << std::setw(count_per_ts_width) << "Call(s)/step"<<std::endl;
          perfs_TU << line_sep_tabular << std::endl;
          for (int i = static_cast<int>(STD_COUNTERS::system_solver); i< static_cast<int>(STD_COUNTERS::postreatment); i++)
            {
              Counter& c_to_print = *std_counters_[i];
              write_globalTU_line(c_to_print,perfs_TU);
            }
          total_time = c_total_time.total_time_.count();
          // Loop on the custom counters
          if (custom_counter_map_str_to_counter_.empty())
            perfs_TU <<std::endl  << std::left <<std::setw(counter_description_width) << "Untracked time" << separator << std::setw(time_per_step_width) << total_untracked_time_ts + total_untracked_time << separator << std::setprecision(3) <<  std::setw(percent_loop_time_width) << 100* (total_untracked_time_ts/time_tl + total_untracked_time/total_time) << separator <<std::endl <<std::endl;
          else
            {
              perfs_TU << std::endl;
              perfs_TU << std::left <<std::setw(counter_description_width) << "Custom counter description" << separator << std::setw(time_per_step_width) << "Time/step" << separator << std::setw(percent_loop_time_width) << "% loop time" << separator << std::setw(count_per_ts_width) << "Call(s)/step"<< separator <<std::setw(level_width) << "Level" <<std::endl;
              perfs_TU << line_sep_tabular_custom<<std::endl;
              for (const auto & pair : custom_counter_map_str_to_counter_)
                write_globalTU_line_custom_counters(*pair.second, perfs_TU);
              perfs_TU <<std::endl   << std::left <<std::setw(counter_description_width) << "Untracked time" << separator << std::setw(time_per_step_width) << total_untracked_time_ts + total_untracked_time << separator << std::setprecision(3) <<  std::setw(percent_loop_time_width) << 100* (total_untracked_time_ts/time_tl + total_untracked_time/total_time) << separator <<std::endl <<std::endl ;
            }
        }
      if (max_virtual_swap_c>0)
        {
          if (message=="Time loop statistics")
            perfs_TU <<  std::left <<std::setw(text_width) << "Total number of virtual exchanges:" <<  std::left <<std::setw(number_width) << max_virtual_swap_c << std::endl;
          else
            perfs_TU <<  std::left <<std::setw(text_width) << "Number of virtual exchanges:" <<  std::left <<std::setw(number_width) << max_virtual_swap_c << std::endl;

        }
      if (min_max_avg_sd_t_q_c_allreduce_comm[2][1]>0 && nb_ts>0)
        {
          double allreduce_per_ts = (double) min_max_avg_sd_t_q_c_allreduce_comm[2][1]/nb_ts;
          perfs_TU <<  std::left <<std::setw(text_width) << "Maximum number of MPI allreduce per time step" <<  std::left <<std::setw(number_width) << allreduce_per_ts << std::endl;
          if (allreduce_per_ts > 30.0 && message=="Time loop statistics")
            {
              perfs_TU << std::endl << line_sep_cpu << std::endl;
              perfs_TU << " Warning: number of MPI_allreduce calls per time step is high. Contact TRUST team to run massive parallel calculations" << std::endl;
              perfs_TU << line_sep_cpu<< std::endl;
            }
        }
      int tmp = c_system_solver.count_;
      if (tmp > 0)
        {
          perfs_TU << std::endl;
          if (!(message=="Time loop statistics"))
            {
              if (nb_ts>0)
                perfs_TU <<  std::left <<std::setw(text_width) << "Number of calls to the linear solver per time step: " <<  std::left <<std::setw(number_width) << static_cast<double>(tmp) / nb_ts << std::endl;
              else
                perfs_TU <<  std::left <<std::setw(text_width) << "Number of call to the linear solver: " <<  std::left <<std::setw(number_width) << tmp << std::endl;
              avg_solv_time = avg_solv_time / tmp;
              perfs_TU <<  std::left <<std::setw(text_width) << "Average time of the resolution of the linear problem per call: " <<  std::left <<std::setw(number_width) << avg_solv_time << std::endl;
            }
          perfs_TU <<  std::left <<std::setw(text_width) << "Average number of iteration of the linear solver per call: " <<  std::left <<std::setw(number_width) <<  nb_it_per_solver_calls << std::endl <<std::endl;
        }
      // GPU part of the TU :
      auto compute_percent_and_write_tabular_line = [& perfs_GPU, & nb_ts, & time_tl, &separator] (const Counter& c_, const std::string str)
      {
        double max_time = c_.time_alone_.count();
        double calls = c_.count_/nb_ts;
        double t_ts = max_time/nb_ts;
        double bw = static_cast<double>(c_.quantity_)/(1024.*1024.*1024*max_time);
        double percent = 100*max_time/time_tl;
        perfs_GPU << std::left << std::setw(counter_description_width) << str <<separator << std::setw(time_per_step_width) << t_ts<<separator  << std::setw(percent_loop_time_width) <<  percent <<separator<< std::setw(count_per_ts_width) << calls <<separator;
        if (bw >1.0e-10)
          perfs_GPU << std::setw(bandwith_width) << bw << "GB/s" << std::endl;
        else
          perfs_GPU << std::endl;
        return percent;
      };
      if (copy_to_device_count>0 && nb_ts >0 && message=="Time loop statistics")
        {
          perfs_GPU << std::endl << line_sep_gpu << std::endl;
          spaces.assign((max_str_length_-14)/2,' ');
          perfs_GPU << spaces <<"GPU statistics" << std::endl;
          perfs_GPU << line_sep_gpu<<std::endl;
          perfs_GPU << std::left <<std::setw(counter_description_width) << "Counter description" << separator <<std::setw(time_per_step_width) << "Time per step" <<separator<< std::setw(percent_loop_time_width) << "% loop time" <<separator<< std::setw(count_per_ts_width) << "Call(s)/step" <<separator<< std::setw(bandwith_width)<< "Bandwidth"<<std::endl;
          perfs_GPU << line_sep_gpu << std::endl;
          double ratio_gpu_library = compute_percent_and_write_tabular_line(c_gpu_l,"Libraries: ");
          double ratio_gpu_kernel = compute_percent_and_write_tabular_line(c_gpu_k,"Kernels: ");
          double ratio_gpu = ratio_gpu_kernel+ratio_gpu_library;
          double ratio_copy = compute_percent_and_write_tabular_line(c_todevice,"Copy host to device: ");
          ratio_copy += compute_percent_and_write_tabular_line(c_fromdevice,"Copy device to host: ");
          double ratio_comm = 100.0 * (total_comm_time)/time_tl;
          double ratio_cpu = 100 * cpu_time/time_tl;
          double ratio_allocfree = compute_percent_and_write_tabular_line(c_allocfree,"Alloc/Free on device: ");
          perfs_GPU << std::setprecision(3) << "GPU: " << ratio_gpu << "% Copy H<->D: " << ratio_copy << "% Alloc/free:" << ratio_allocfree << "% Comm: "<< ratio_comm << "% CPU & other: " << ratio_cpu << std::setprecision(6)<<"%"<<std::endl;
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
      if (message=="Time loop statistics")
        {
          if (debit_seq>0 || debit_par>0)
            {
              perfs_IO << std::endl << line_sep_cpu << std::endl;
              spaces.assign((max_str_length_-message_width-4)/2,' ');
              perfs_IO << spaces << message <<": IO" << std::endl;
              perfs_IO << line_sep_cpu<<std::endl;
            }
          if (debit_seq>0)
            perfs_IO << std::left <<std::setw(text_width) << "Output write sequential: " <<  std::left <<std::setw(number_width) << debit_seq << "MB/s"<< std::endl;
          if (debit_par>0)
            perfs_IO <<  std::left <<std::setw(text_width) << "Output write parallel: " <<  std::left <<std::setw(number_width) << debit_par  << "MB/s" << std::endl;
          if (total_nb_backup_>0)
            {
              perfs_IO <<  std::left <<std::setw(text_width) << "Total number of back-up: " <<  std::left <<std::setw(number_width) << total_nb_backup_ << std::endl;
              perfs_IO <<  std::left <<std::setw(text_width) << "Total amount of data per back-up: " <<  std::left <<std::setw(number_width) << total_data_exchange_per_backup_ << "MB"<< std::endl;
            }
          if(min_max_avg_sd_t_q_c_sendrecv_comm[2][1] > 0)
            {
              if(petcs_count>0)
                {
                  perfs_IO<< std::endl<< "---------------------------------------------------------------------------------------------------------"<< std::endl;
                  perfs_IO<< "Warning: One or several PETSc solvers are used and thus the communication time below are under-estimated."<< std::endl;
                  perfs_IO<< "---------------------------------------------------------------------------------------------------------"<< std::endl<< std::endl;
                }
              double fraction = 0.0;
              fraction = (comm_sendrecv_t + comm_allreduce_t)/ (total_time + MINFLOAT);
              fraction = 0.1 * floor(fraction * 1000);
              if (fraction > 100.)
                fraction = 100.;
              perfs_IO <<  std::left <<std::setw(text_width) << "Average of the fraction of the time spent in communications between processors: " <<  std::left <<std::setw(number_width) << fraction << "%" << std::endl;
              fraction = (min_max_avg_sd_t_q_c_sendrecv_comm[0][1] + min_max_avg_sd_t_q_c_allreduce_comm[0][1])/ (total_time_max + MINFLOAT);
              fraction = 0.1 * floor(fraction * 1000);
              if (fraction > 100.)
                fraction = 100.;
              perfs_IO <<  std::left <<std::setw(text_width) << "Max of the fraction of the time spent in communications between processors: " <<  std::left <<std::setw(number_width) << fraction << "%" << std::endl;
              fraction = (min_max_avg_sd_t_q_c_sendrecv_comm[0][0] + min_max_avg_sd_t_q_c_allreduce_comm[0][0])/ (total_time_max + MINFLOAT);
              fraction = 0.1 * floor(fraction * 1000);
              perfs_IO <<  std::left <<std::setw(text_width) << "Min of the fraction of the time spent in communications between processors: " <<  std::left <<std::setw(number_width) << fraction << "%"  << std::endl;
              perfs_IO  <<  std::left <<std::setw(text_width) << "Time of one mpsum measured by an internal bench over 0.1s (network latency): ";
              if (allreduce_peak_perf == 0.)
                perfs_IO  << "not measured (total running time too short <10s)" << std::endl;
              else
                perfs_IO <<  std::left <<std::setw(number_width) << allreduce_peak_perf << std::endl;
              perfs_IO <<  std::left <<std::setw(text_width) << "Network maximum bandwidth on all processors: "  <<  std::left <<std::setw(number_width) << max_bandwidth * 1.e-6 << "MB/s"  << std::endl ;
              if (nb_ts>0)
                perfs_IO <<  std::left <<std::setw(text_width) << "Total network traffic: " <<  std::left <<std::setw(number_width) << static_cast<double>(comm_sendrecv_q) * Process::nproc() / nb_ts * 1e-6 << "MB/time step"  << std::endl;
              else
                perfs_IO <<  std::left <<std::setw(text_width) << "Total network traffic: " <<  std::left <<std::setw(number_width) << static_cast<double>(comm_sendrecv_q) * Process::nproc()* 1e-6 << "MB"  << std::endl;
              perfs_IO <<  std::left <<std::setw(text_width) << "Average message size: " <<  std::left <<std::setw(number_width) << static_cast<double>(comm_sendrecv_q) / static_cast<double>(comm_sendrecv_c)* 1e-3 << "kB" << std::endl;
              perfs_IO <<  std::left <<std::setw(text_width) << "Min waiting time: "  <<  std::left <<std::setw(number_width) << min_wait_fraction << "% of total time"<< std::endl;;
              perfs_IO <<  std::left <<std::setw(text_width) << "Max waiting time: "  <<  std::left <<std::setw(number_width) << max_wait_fraction<< "% of total time"<< std::endl;;
              perfs_IO <<  std::left <<std::setw(text_width) << "Avg waiting time: " <<  std::left <<std::setw(number_width) << avg_wait_fraction<< "% of total time"<< std::endl;;
            }
        }
      // Concatenate stringtreams in order to print the .TU file
      Nom globalTU(Objet_U::nom_du_cas());
      globalTU +=".TU";
      auto printing_mode = ios::app;
      if (message=="Computation start-up statistics")
        printing_mode= ios::out;
      SFichier file(globalTU, printing_mode);
      file << file_header.str();
      file << perfs_TU.str();
      file << perfs_IO.str();
      file << perfs_GPU.str();
      file << captions.str();
      clean_stringstream(file_header);
      clean_stringstream(perfs_TU);
      clean_stringstream(perfs_GPU);
      clean_stringstream(perfs_IO);
      clean_stringstream(captions);
    }
}

/////////////////////////////////////// Public methods of Pimpl ////////////////////////////////////////////

void Perf_counters::Impl::create_custom_counter_impl(std::string counter_description , int counter_level,  std::string counter_family , bool is_comm, bool is_gpu)
{
  if (counter_level <=0)
    Process::exit("Custom counters should not be set with a zero or negative level value");
  if (custom_counter_map_str_to_counter_.count(counter_description)==0)
    {
      auto result =custom_counter_map_str_to_counter_.emplace(counter_description, std::make_unique<Counter>(counter_level, counter_description, counter_family ,is_comm, is_gpu));
      if (!result.second)
        Process::exit("Failed to insert the new custom counter in the custom counter map");
    }
}

void Perf_counters::Impl::begin_count_impl(const STD_COUNTERS& std_cnt, int counter_lvl)
{
  if (!counters_stop_)
    {
      Counter& c = get_counter(std_cnt);
      if (counter_lvl == -100000)
        counter_lvl = c.level_;
      time_point t = now();
      check_begin(c, counter_lvl,t);
      c.begin_count_(counter_lvl,t);
    }
}

/*!
 *
 * @param custom_count_name key of the map (custom_counter_map_str_to_counter_) of the custom counter you try to close
 * @param counter_lvl
 */
void Perf_counters::Impl::begin_count_impl(const std::string& custom_count_name, int counter_lvl)
{
  if (!counters_stop_)
    {
      Counter& c = get_counter(custom_count_name);
      if (counter_lvl == -100000)
        counter_lvl = c.level_;
      time_point t = now();
      check_begin(c, counter_lvl,t);
      c.begin_count_(counter_lvl,t);
    }
}

/*! @brief End the count of a counter and update the counter values
 *
 * @param c is the counter to end the count
 * @param count_increment is the count increment. If not specified, then it is equal to 1
 * @param quantity_increment is the increment of custom variable quantity. If not specified, it is set to 0.
 */
void Perf_counters::Impl::end_count_impl(const STD_COUNTERS& std_cnt, int count_increment, long int quantity_increment)
{
  if (!counters_stop_)
    {
      Counter& c = get_counter(std_cnt);
      time_point t = now();
      check_end(c, t);
      c.end_count_(count_increment, quantity_increment,t);
      if (c.level_ == -1)
        computation_time_ += c.total_time_;
    }
}


/*! @brief End the count of a counter and update the counter values
 *
 * @param c is the custom counter to end the count
 * @param count_increment is the count increment. If not specified, then it is equal to 1
 * @param quantity_increment is the increment of custom variable quantity. If not specified, it is set to 0.
 */
void Perf_counters::Impl::end_count_impl(const std::string& custom_count_name, int count_increment, long int quantity_increment)
{
  if (!counters_stop_)
    {
      Counter& c = get_counter(custom_count_name);
      time_point t = now();
      assert(custom_counter_map_str_to_counter_.count(custom_count_name) > 0);
      check_end(c, t);
      c.end_count_(count_increment, quantity_increment,t);
    }
}


/*! @brief Stop all counters, has to be called on every processor simultaneously
 *
 */
void Perf_counters::Impl::stop_counters_impl()
{
  time_point t_stop = now();
  duration time_elapsed_before_stop;
  if (counters_stop_)
    Process::exit("The counter are already stop, you can't stop them two times in a row \n");
  if (last_opened_counter_ != nullptr)
    {
      Counter* c = last_opened_counter_;
      c->time_alone_ += t_stop -c->last_open_time_alone_;
      c->last_open_time_alone_ = time_point();
      while (c!=nullptr)
        {
          time_elapsed_before_stop= t_stop -c->last_open_time_;
          if (time_loop_ )
            c->time_ts_ += time_elapsed_before_stop;
          c->total_time_ += time_elapsed_before_stop;
          c = c->parent_;
        }
    }
  counters_stop_=true;
}

/*! @brief Restart all counters, has to be called on every processor simultaneously
 *
 */
void Perf_counters::Impl::restart_counters_impl()
{
  time_point t_restart = now();
  if (counters_stop_==false)
    Process::exit("Try to restart counters but they have never been stopped before");
  if (last_opened_counter_ != nullptr)
    {
      Counter* c = last_opened_counter_;
      c->last_open_time_alone_ = t_restart;
      while (c !=nullptr )
        {
          c->last_open_time_ = t_restart;
          if (time_loop_)
            c->open_time_ts_ = t_restart;
          c = c->parent_;
        }
    }
  counters_stop_=false;
}

void Perf_counters::Impl::reset_counters_impl()
{
  for (int i =0 ; i< static_cast<int>(STD_COUNTERS::NB_OF_STD_COUNTER); i++)
    {
      Counter& c = *std_counters_[i];
      c.reset();
    }
  for (const auto & pair : custom_counter_map_str_to_counter_)
    {
      Counter& c = *pair.second;
      c.reset();
    }
}

void Perf_counters::Impl::set_time_steps_elapsed_impl(int time_step_elapsed)
{
  nb_steps_elapsed_=time_step_elapsed;
}

double Perf_counters::Impl::get_computation_time_impl()
{
  if (counters_stop_)
    Process::exit("The counters are stop, you can't access the total time");
  Counter& c = get_counter(STD_COUNTERS::total_execution_time);
  if (c.is_running_)
    computation_time_+= now() - c.last_open_time_;
  return (computation_time_.count());
}

double Perf_counters::Impl::get_total_time_impl(const STD_COUNTERS& name)
{
  Counter& c = get_counter(name);
  duration t = c.total_time_;
  if (c.is_running_)
    t += now() - c.last_open_time_;
  return (t.count());
}

double Perf_counters::Impl::get_total_time_impl(const std::string& name)
{
  Counter& c = get_counter(name);
  duration t = c.total_time_;
  if (c.is_running_)
    t += now() - c.last_open_time_;
  return (t.count());
}

double Perf_counters::Impl::get_time_since_last_open_impl(const STD_COUNTERS& name)
{
  Counter& c = get_counter(name);
  duration t =duration::zero();
  if (!c.is_running_)
    Process::exit("The counter is not running: " + c.description_);
  t = now() - c.last_open_time_;
  return (t.count());
}

double Perf_counters::Impl::get_time_since_last_open_impl(const std::string& name)
{
  Counter& c = get_counter(name);
  duration t = duration::zero();
  if (!c.is_running_)
    Process::exit("The counter is not running: " + c.description_);
  t = now() - c.last_open_time_;
  return (t.count());
}

void Perf_counters::Impl::start_timeloop_impl()
{
  if (last_opened_counter_==nullptr)
    Process::exit("You are trying to start the time loop before the start-up");
  time_loop_=true;
}

void Perf_counters::Impl::end_timeloop_impl()
{
  if (!time_loop_)
    Process::exit("The time loop has not started, but you are trying to end it");
  time_loop_=false;
}

void Perf_counters::Impl::start_time_step_impl()
{
  assert (last_opened_counter_!=nullptr);
  time_point t = now();
  if (last_opened_counter_ == nullptr)
    Process::exit("You are trying to start a time step outside the time loop");

  for (int i =0 ; i< static_cast<int>(STD_COUNTERS::NB_OF_STD_COUNTER); i++)
    {
      Counter& c_std = *std_counters_[i];
      c_std.time_ts_=duration::zero();
    }
  for (const auto & pair : custom_counter_map_str_to_counter_)
    {
      Counter& c_c = *pair.second;
      c_c.time_ts_=duration::zero();
    }
  Counter* c = last_opened_counter_;
  while (c->parent_ != nullptr && c->level_>=0)
    {
      c->open_time_ts_ = t;
      c = c->parent_;
    }
}

/*! @brief Compute for each counter open during a time step avg_time_per_step_, min_time_per_step_, max_time_per_step_ and sd_time_per_step_
 *
 * Called at the end of each time step, and only then
 * @param tstep number of time steps elapsed since the start of the computation
 *
 * Nota : if the counter is called before the time step loop, it is not accounted for in the computation.
 */
void Perf_counters::Impl::end_time_step_impl(unsigned int tstep)
{
  stop_counters_impl(); ///< stop_counters already updated c->tim_ts_
  if (last_opened_counter_ == nullptr)
    Process::exit("You are trying to compute the statistics of a time steps but have not open any counter");
  if (!time_loop_)
    Process::exit("You are trying to compute time loop statistics outside of the time loop");
  int step = tstep - nb_steps_elapsed_;
  auto compute = [step](Counter& c)
  {
    if (c.level_>=0 && step>0 && c.count_>0)
      {
        c.min_time_per_step_ = (c.min_time_per_step_ < (c.time_ts_).count()) ? c.min_time_per_step_ : (c.time_ts_).count();
        c.max_time_per_step_ = (c.min_time_per_step_ > (c.time_ts_).count()) ? c.min_time_per_step_ : (c.time_ts_).count();
        c.avg_time_per_step_ = ((step-1)*c.avg_time_per_step_ + (c.time_ts_).count())/step;
        c.sd_time_per_step_ += ((c.time_ts_).count()* (c.time_ts_).count() - 2*((c.time_ts_).count())* c.avg_time_per_step_ +  c.avg_time_per_step_ *  c.avg_time_per_step_)/static_cast<double>(step);
        c.sd_time_per_step_ = sqrt(c.sd_time_per_step_);
        if (c.sd_time_per_step_ < 0)
          c.sd_time_per_step_ = 0;
      }
    c.open_time_ts_=time_point();
  };
  if (end_cache_)
    {
      for (int i =0 ; i< static_cast<int>(STD_COUNTERS::NB_OF_STD_COUNTER); i++)
        {
          Counter& c_std = *std_counters_[i];
          compute(c_std);
        }
      if (!custom_counter_map_str_to_counter_.empty())
        {
          for (const auto & pair : custom_counter_map_str_to_counter_)
            {
              Counter& c = *pair.second;
              compute(c);
            }
        }
    }
  else
    {
      end_cache_ = tstep >= nb_steps_elapsed_;
      if (end_cache_)
        {
          time_skipped_ts_ += get_counter(STD_COUNTERS::total_execution_time).total_time_;
          computation_time_ += get_counter(STD_COUNTERS::total_execution_time).total_time_;
          reset_counters_impl();
        }
    }
  restart_counters_impl();
}

int Perf_counters::Impl::get_last_opened_counter_level_impl() const
{
  if (last_opened_counter_==nullptr)
    return -5;
  return last_opened_counter_->level_;
}

void Perf_counters::Impl::print_TU_files_impl(const std::string& message)
{
  if(!Objet_U::disable_TU)
    {
      //Process::barrier();
      stop_counters_impl();
      Counter& c_time = get_counter(STD_COUNTERS::total_execution_time);
      computation_time_ += c_time.total_time_;
      print_global_TU(message);
      print_performance_to_csv(message);
      reset_counters_impl();
      counters_stop_=false;
    }
}

void Perf_counters::Impl::start_gpu_timer_impl()
{
  if (gpu_timer_)
    Process::exit("You try to start the gpu timer and it is already running");
  gpu_timer_start_=now();
  gpu_timer_ = true;
}

void Perf_counters::Impl::stop_gpu_timer_impl()
{
  if(!gpu_timer_)
    Process::exit("You try to stop the GPU timer, but it has not been started yet");
  gpu_timer_=false;
}

double Perf_counters::Impl::compute_gpu_time_impl()
{
  stop_gpu_timer_impl();
  duration d= now() - gpu_timer_start_;
  return d.count();
}

bool Perf_counters::Impl::is_gpu_verbose_on_impl() const
{
  return gpu_verbose_;
}
void Perf_counters::Impl::set_gpu_verbose_impl(bool on)
{
  gpu_verbose_ = on;
}
bool Perf_counters::Impl::get_init_device_impl() const
{
  return init_device_;
}
void Perf_counters::Impl::set_init_device_impl(bool init)
{
  init_device_=init;
}
bool Perf_counters::Impl::get_gpu_timer_impl() const
{
  return gpu_timer_;
}
void Perf_counters::Impl::add_to_gpu_timer_counter_impl(int to_add)
{
  gpu_timer_count_ += to_add;
}
int Perf_counters::Impl::get_gpu_timer_counter_impl() const
{
  return gpu_timer_count_;
}

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
//                                              Methods of Perf_counters
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
Perf_counters::Perf_counters() : pimpl_(std::make_unique<Impl>())
{

}

Perf_counters::~Perf_counters()=default;


double Perf_counters::compute_time(time_point start)
{
  duration d= now() - start;
  {return d.count();}
}

void Perf_counters::create_custom_counter(std::string counter_description , int counter_level,  std::string counter_family , bool is_comm, bool is_gpu)
{
  pimpl_->create_custom_counter_impl(counter_description , counter_level,  counter_family , is_comm, is_gpu);
}

void Perf_counters::begin_count(const STD_COUNTERS& std_cnt, int counter_lvl)
{
  pimpl_->begin_count_impl(std_cnt,counter_lvl);
}

void Perf_counters::begin_count(const std::string& custom_count_name, int counter_lvl)
{
  pimpl_->begin_count_impl(custom_count_name,counter_lvl);
}

void Perf_counters::end_count(const std::string& custom_count_name, int count_increment, long int quantity_increment)
{
  pimpl_->end_count_impl(custom_count_name,count_increment,quantity_increment);
}

void Perf_counters::end_count(const STD_COUNTERS& std_cnt, int count_increment, long int quantity_increment)
{
  pimpl_->end_count_impl(std_cnt,count_increment,quantity_increment);
}

void Perf_counters::stop_counters()
{
  pimpl_->stop_counters_impl();
}

void Perf_counters::restart_counters()
{
  pimpl_->restart_counters_impl();
}

void Perf_counters::reset_counters()
{
  pimpl_->reset_counters_impl();
}

void Perf_counters::print_TU_files(const std::string& message)
{
  pimpl_->print_TU_files_impl(message);
}

double Perf_counters::get_computation_time()
{
  return pimpl_->get_computation_time_impl();
}

double Perf_counters::get_total_time(const STD_COUNTERS& name)
{
  return pimpl_->get_total_time_impl(name);
}

double Perf_counters::get_total_time(const std::string& name)
{
  return pimpl_->get_total_time_impl(name);
}

double Perf_counters::get_time_since_last_open(const STD_COUNTERS& name)
{
  return pimpl_->get_time_since_last_open_impl(name);
}

double Perf_counters::get_time_since_last_open(const std::string& name)
{
  return pimpl_->get_time_since_last_open_impl(name);
}

bool Perf_counters::is_running(const STD_COUNTERS& name)
{
  return pimpl_->running_impl(name);
}

void Perf_counters::start_timeloop()
{
  pimpl_->start_timeloop_impl();
}

void Perf_counters::end_timeloop()
{
  pimpl_->end_timeloop_impl();
}

void Perf_counters::start_time_step()
{
  pimpl_->start_time_step_impl();
}

void Perf_counters::end_time_step(unsigned int tstep)
{
  pimpl_->end_time_step_impl(tstep);
}

void Perf_counters::set_nb_time_steps_elapsed(unsigned int n)
{
  pimpl_->set_time_steps_elapsed_impl(n);
}

int Perf_counters::get_last_opened_counter_level() const
{
  return pimpl_->get_last_opened_counter_level_impl();
}

void Perf_counters::start_gpu_timer()
{
  pimpl_->start_gpu_timer_impl();
}

void Perf_counters::stop_gpu_timer()
{
  pimpl_->stop_gpu_timer_impl();
}

bool Perf_counters::is_gpu_verbose_on() const
{
  return pimpl_->is_gpu_verbose_on_impl();
}

void Perf_counters::set_gpu_verbose(bool on)
{
  pimpl_->set_gpu_verbose_impl(on);
}

bool Perf_counters::get_init_device() const
{
  return pimpl_->get_init_device_impl();
}

void Perf_counters::set_init_device(bool init)
{
  pimpl_->set_init_device_impl(init);
}

bool Perf_counters::get_gpu_timer() const
{
  return pimpl_->get_gpu_timer_impl();
}

void Perf_counters::add_to_gpu_timer_counter(int to_add)
{
  pimpl_->add_to_gpu_timer_counter_impl(to_add);
}

int Perf_counters::get_gpu_timer_counter() const
{
  return pimpl_->get_gpu_timer_counter_impl();
}

double Perf_counters::stop_gpu_timer_and_compute_gpu_time()
{
  return pimpl_->compute_gpu_time_impl();
}

bool Perf_counters::get_use_gpu() const {return pimpl_->get_use_gpu_impl();}

bool Perf_counters::get_gpu_fence() const {return pimpl_->get_gpu_fence_impl();}

void Perf_counters::set_gpu_fence(bool fence) {return pimpl_->set_gpu_fence_impl(fence);}
