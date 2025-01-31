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
#include <mutex>
#include <vector>
#include <map>
#include <tuple>
#include <array>
#include <memory>

#ifndef Perf_counters_included
#define Perf_counters_included

// This file contains all of the needed for the description of the counter associated with the tracking of performance in the TRUST code.
// A class Time is introduced in Perf_counters.cpp for extracting the time
class Counter;

enum class STD_COUNTERS : unsigned int
{
  total_execution_time_ , ///< Lowest level counter that track the total time of the computation
  computation_start_up_ , ///< Track the time before the Resoudre loop
  timestep_ ,   ///< Track time elapsed in the time loop
  system_solver_, ///< Track time elapsed in SolveurSys::resoudre_systeme
  petsc_solver_,  ///< Track the time elapsed using petsc solver
  implicit_diffusion_,  ///< Track time elapsed in Equation_base::conjugue_diff_impl
  compute_dt_ , ///< Track time used to compute the time step dt
  turbulent_viscosity_ ,
  convection_ ,
  diffusion_ ,
  gradient_ ,
  divergence_ ,
  rhs_ ,
  postreatment_ ,
  backup_file_ ,
  restart_ ,
  matrix_assembly_ ,
  update_variables_  ,
  virtual_swap_ ,
  mpi_sendrecv_  ,
  mpi_send_ ,
  mpi_recv_ ,
  mpi_bcast_ ,
  mpi_alltoall_ ,
  mpi_allgather_ ,
  mpi_gather_ ,
  mpi_partialsum_,
  mpi_sumdouble_ ,
  mpi_mindouble_,
  mpi_maxdouble_ ,
  mpi_sumfloat_ ,
  mpi_minfloat_ ,
  mpi_maxfloat_ ,
  mpi_sumint_ ,
  mpi_minint_ ,
  mpi_maxint_ ,
  mpi_barrier_ ,
  gpu_library_ ,
  gpu_kernel_ ,
  gpu_copytodevice_ ,
  gpu_copyfromdevice_ ,
  IO_EcrireFicPartageMPIIO_ ,
  IO_EcrireFicPartageBin_ ,
  interprete_scatter_,
  read_scatter_,
  NB_OF_STD_COUNTER
};


/*!  @brief This class stores and manages counters in TRUST. It is a singleton.
 *
 */
class Perf_counters
{
public:

  /*! @brief The class Perf_counters is based on a singleton pattern. To access to the unique object inside the code, use the getInstance() function
   *
   * @return the unique Perf_counters object
   */
  static Perf_counters& getInstance()
  {
    static Perf_counters counters_stat_ ;
    return counters_stat_;
  }

  /*! @brief Create a new counter and add it to the vector of counters
   *
   * @param counter_level
   * @param counter_description
   * @param counter_family
   * @param is_comm
   * @return create a new counter
   */
  void create_custom_counter(unsigned int counter_level, std::string counter_description, std::string counter_family = "None", bool is_comm=false);

  /*! @brief Start the count of a counter
   *
   */
  //	void begin_count(Counter c);


  /*! Standard counters, start the tracking of the wanted operation
   *
   * @param std_cnt reference to the standard counter
   * @param counter_lvl level of the counter you try to open, warning it changes the value of the counter level associated with counter std_cnt
   */
  void begin_count(const STD_COUNTERS std_cnt, int counter_lvl);

  /*! Custom counters, start the tracking of the wanted operation
   *
   * @param custom_count_name string key in the custom counter map
   * @param counter_lvl level of the counter you try to open, warning it changes the value of the level of the counter associated with custom_count_name
   */
  void begin_count(const std::string& custom_count_name, unsigned int counter_lvl);

  /*! @brief Ensure that the counter you are trying to open is not open yet and that the level is correct and update last_opened_counter_
   *
   * @param c counter you try to open
   * @param counter_lvl
   * @param t time of opening
   */
  void check_begin(Counter& c, unsigned int counter_lvl, std::chrono::time_point<std::chrono::high_resolution_clock> t);

  /*! @brief Used to see if the counter you want to close is indeed the last open and update last_opened_counter_
   *
   * @param c counter you try to close
   * @param t time of closing
   */
  void check_end(Counter& c, std::chrono::time_point<std::chrono::high_resolution_clock> t);


  /*! @brief End the count of a counter and update the counter values
   *
   * @param c is the counter to end the count
   * @param count_increment is the count increment. If not specified, then it is equal to 1
   * @param quantity_increment is the increment of custom variable quantity. If not specified, it is set to 0.
   */
  void end_count(const std::string& custom_count_name, int count_increment=1, int quantity_increment=0);

  /*! @brief End the count of a counter and update the counter values
   *
   * @param c is the counter to end the count
   * @param count_increment is the count increment. If not specified, then it is equal to 1
   * @param quantity_increment is the increment of custom variable quantity. If not specified, it is set to 0.
   */
  void end_count(const STD_COUNTERS std_cnt, int count_increment=1, int quantity_increment=0);

  /*! @brief Stop all counters, has to be called on every processor simultaneously
   *
   */
  void stop_counters();

  /*! @brief Restart all counters, has to be called on every processor simultaneously
   *
   */
  void restart_counters();

  /*! @brief Reset counters to zero, used between the start-up of the computation, the computation itself and the post-processing
   *
   */
  void reset_counters();

  /*!
   *
   * @param is_the_three_first_time_steps_elapsed == True, then the thre first time steps are discarded
   */
  inline void set_three_first_steps_elapsed(bool is_the_three_first_time_steps_elapsed)
  {
    two_first_steps_elapsed_ = is_the_three_first_time_steps_elapsed;
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

  /*!@brief Give as a double the total time (in second) elapsed in the operation tracked by the standard counter call name
   *
   */
  double get_total_time(STD_COUNTERS name);

  /*!@brief Give as a double the total time (in second) elapsed in the operation tracked by the custom counter call name
     *
     */
  double get_total_time(std::string name);

  /*!@brief Give as a double the time (in second) elapsed in the operation tracked by the standard counter call name since the counter was last opened
   *
   */
  double get_time_since_last_open(STD_COUNTERS name);

  /*!@brief Give as a double the time (in second) elapsed in the operation tracked by the standard counter call name since the counter was last opened
   *
   */
  double get_time_since_last_open(std::string name);


  /*!
   * This function aims at starting the tracking of time per time step of counters that have been started before the time loop
   */
  void start_timestep();

  /*!@brief This function compute statistics per time steps of counters used at least once during a time step.
   *
   * @param tstep is the current time step number
   */
  void compute_avg_min_max_var_per_step(int tstep);

  /*! @brief Accessor to the Counter object which pointer is stored in the std_counters_ array
   *
   * @return the reference of a the counter object associated with STD_COUNTERS::name
   */
  inline Counter& access_std_counter(STD_COUNTERS name) {return *std_counters_[static_cast<int>(name)];}

  /*! @brief Accessor to the Counter object which pointer is stored in the std_counters_ array
   *
   * @return the reference of a the counter object associated with custom_counter_map_str_to_counter_[name]
   */
  inline Counter& access_custom_counter(std::string name) {return *custom_counter_map_str_to_counter_.at(name);}

  std::string get_os();

  std::string get_cpu();

  std::string get_gpu();

  std::string get_date();

private:

  Perf_counters();
  ~Perf_counters();
  Perf_counters(const Perf_counters&) = delete;
  Perf_counters& operator=(const Perf_counters&) = delete;
  bool two_first_steps_elapsed_;  ///< By default, we consider that the two first time steps are used to file the cache, so they are not taken into account in the stats.
  bool end_cache_; ///< A flag used to know if the two first time steps are over or not
  bool counters_stop_;  ///< A flag used to know if the counters are paused or not
  std::chrono::duration<double> time_cache_; ///< the duration in seconds of the cache. If cache is too long, use function set_three_first_steps_elapsed in oder to include the stats of the cache in your stats
  Counter * last_opened_counter_; ///< pointer to the last opened counter. Each counter has a parent attribute, which also give the pointer of the counter open before them.
  std::array <Counter*,static_cast<int>(STD_COUNTERS::NB_OF_STD_COUNTER)> std_counters_ ; ///< Array of the pointers to the standard counters of TRUST
  std::map <std::string, Counter* > custom_counter_map_str_to_counter_ ; ///< Map that link the descriptions of the custom counters to their pointers
};

#endif

