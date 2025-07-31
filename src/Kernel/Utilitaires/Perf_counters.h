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
#include <chrono>
#include <memory>
#include <iostream>

#ifndef Perf_counters_included
#define Perf_counters_included

// This file contains all of the needed for the description of the counter associated with the tracking of performance in the TRUST code.

enum class STD_COUNTERS : unsigned int
{
  total_execution_time , ///< Lowest level counter that track the total time of the computation
  computation_start_up , ///< Track the time before the Resoudre loop
  timeloop ,   ///< Track time elapsed in the time loop
  backup_file ,
  system_solver, ///< Track time elapsed in SolveurSys::resoudre_systeme
  convection ,
  diffusion ,
  rhs ,
  gradient ,
  divergence ,
  matrix_assembly ,
  update_variables  ,
  implicit_diffusion,  ///< Track time elapsed in Equation_base::conjugue_diff_impl
  compute_dt , ///< Track time used to compute the time step dt
  turbulent_viscosity ,
  restart ,
  postreatment ,
  petsc_solver,  ///< Track the time elapsed using petsc solver
  mpi_sendrecv  ,
  mpi_send ,
  mpi_recv ,
  mpi_bcast ,
  mpi_alltoall ,
  mpi_allgather ,
  mpi_gather ,
  mpi_partialsum,
  mpi_sumdouble ,
  mpi_mindouble,
  mpi_maxdouble ,
  mpi_sumfloat ,
  mpi_minfloat ,
  mpi_maxfloat ,
  mpi_sumint ,
  mpi_minint ,
  mpi_maxint ,
  mpi_barrier ,
  gpu_library ,
  gpu_kernel ,
  gpu_copytodevice ,
  gpu_copyfromdevice ,
  gpu_malloc_free,
  IO_EcrireFicPartageMPIIO ,
  IO_EcrireFicPartageBin ,
  interprete_scatter,
  virtual_swap ,
  read_scatter,
  parallel_meshing,
  NB_OF_STD_COUNTER
};



/*!  @brief This class stores and manages counters in TRUST. It is a singleton.
 *
 */
class Perf_counters
{
public:
  using clock = std::chrono::high_resolution_clock;
  using time_point = std::chrono::time_point<clock>;
  using duration = std::chrono::duration<double>;
  inline time_point now()
  {
    return clock::now();
  }

////// Time wrapper /////

  inline time_point start_clock() {return now();} ///< Start a clock, return a time_point, not a double

  double compute_time(time_point start); ///< return time since start in seconds

  /*
    static Perf_counters& getInstance()
    {
      static Perf_counters counters_stat_ ;
      return counters_stat_;
    }
  */
  /*! @brief The class Perf_counters is based on a phoenix singleton pattern. To access to the unique object inside the code, use the getInstance() function
   *
   * @return the unique Perf_counters object
   */
  static Perf_counters& getInstance()
  {
    // Static flag for tracking singleton state
    static bool destroyed = false;
    static Perf_counters* instance = nullptr;

    if (instance == nullptr)
      {
        if (destroyed)
          {
            // If the instance has already been destroyed but still asked for somewhere, the instance is reborn
            std::cout<< "[Stats] The singleton pattern had to be reconstructed for avoiding \'Static Initialization Order Fiasco\' " <<std::endl;
            std::cout<< "It does not impact the TU files but be careful with the values of times printed afterwards" <<std::endl;
            static void* memory = ::operator new(sizeof(Perf_counters));
            instance = new (memory) Perf_counters();
            // A cleaning method for avoiding memory leaks
            std::atexit([]()
            {
              instance->~Perf_counters();
            });
          }
        else
          {
            // For the first creation of the instance
            static Perf_counters counters_stat_;
            instance = &counters_stat_;

            // Update destruction flag
            std::atexit([]()
            {
              destroyed = true;
              instance = nullptr;
            });
          }
      }
    return *instance;
  }
  /*! @brief Create a new counter and add it to the map of custom counters
   *
   * @param to_print_in_global_TU : if true, then the statistics associated with the counter will appear in the global_TU file
   * @param counter_level
   * @param counter_description
   * @param counter_family
   * @param is_comm
   * @return create a new counter
   */
  void create_custom_counter(std::string counter_description , int counter_level,  std::string counter_family = "None", bool is_comm=false, bool is_gpu=false);

  /*! Standard counters, start the tracking of the wanted operation
   *
   * @param std_cnt reference to the standard counter
   * @param counter_lvl level of the counter you try to open, warning it changes the value of the counter level associated with counter std_cnt
   */
  void begin_count(const STD_COUNTERS& std_cnt, int counter_lvl = -100000);

  /*! Custom counters, start the tracking of the wanted operation
   *
   * @param custom_count_name string key in the custom counter map
   * @param counter_lvl level of the counter you try to open, warning it changes the value of the level of the counter associated with custom_count_name
   */
  void begin_count(const std::string& custom_count_name, int counter_lvl= -100000);

  /*! @brief End the count of a counter and update the counter values
   *
   * @param c is the counter to end the count
   * @param count_increment is the count increment. If not specified, then it is equal to 1
   * @param quantity_increment is the increment of custom variable quantity. If not specified, it is set to 0.
   */
  void end_count(const std::string& custom_count_name, int count_increment=1, long int quantity_increment=0);

  /*! @brief End the count of a counter and update the counter values
   *
   * @param c is the counter to end the count
   * @param count_increment is the count increment. If not specified, then it is equal to 1
   * @param quantity_increment is the increment of custom variable quantity. If not specified, it is set to 0.
   */
  void end_count(const STD_COUNTERS& std_cnt, int count_increment=1, long int quantity_increment=0);

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

  /*! @brief Function that encapsulate the two functions that writes the TU files
   *
   */
  void print_TU_files(const std::string& message);


  /*!@brief Update computation_time_ and return its value as a double (in seconds)
   *
   */
  double get_computation_time();

  /*!@brief Give as a double the total time (in second) elapsed in the operation tracked by the standard counter call name
   *
   */
  double get_total_time(const STD_COUNTERS& name) ;

  /*!@brief Give as a double the total time (in second) elapsed in the operation tracked by the custom counter call name
     *
     */
  double get_total_time(const std::string& name) ;

  /*!@brief Give as a double the time (in second) elapsed in the operation tracked by the standard counter call name since the counter was last opened
   *
   */
  double get_time_since_last_open(const STD_COUNTERS& name);

  /*!@brief Give as a double the time (in second) elapsed in the operation tracked by the standard counter call name since the counter was last opened
   *
   */
  double get_time_since_last_open(const std::string& name) ;

  /*!@brief Check whether a counter is already running. Should rarely be used!
  *
  */
  bool is_running(const STD_COUNTERS& name);

  /*!
   * @brief Set time_loop_ to true in order to account for cache properly
   */
  void start_timeloop();

  /*!
   * @brief Set time_loop_ to false as we exit the time loop
   */
  void end_timeloop();

  /*!
   * @brief, this function start statistics tracking for a time step. It has to be called at the start of each time step.
   */
  void start_time_step();

  /*!@brief This function compute statistics per time steps of counters used at least once during a time step.
   *
   * @param tstep is the current time step number
   */
  void end_time_step(unsigned int tstep);

  void set_nb_time_steps_elapsed(unsigned int n) ;

  int get_last_opened_counter_level() const ;

  /////// GPU features for a cleaner Device class

  void start_gpu_timer();

  void stop_gpu_timer();

  bool is_gpu_verbose_on() const ;
  void set_gpu_verbose(bool on) ;

  bool get_init_device() const ;

  void set_init_device(bool init) ;

  bool get_gpu_timer() const ;

  void set_gpu_timer(bool timer);

  void add_to_gpu_timer_counter(int to_add) ;

  int get_gpu_timer_counter() const ;

  double stop_gpu_timer_and_compute_gpu_time() ;

  bool get_use_gpu() const;

  bool get_gpu_fence() const;

  void set_gpu_fence(bool fence);

//// end of GPU features

  Perf_counters(const Perf_counters&) = delete;
  Perf_counters& operator=(const Perf_counters&) = delete;

private:

  Perf_counters();
  ~Perf_counters();
  class Impl;
  std::unique_ptr<Impl> pimpl_;
};

/*! @brief An that compact the access to the unique Perf_counters object
 *
 * @return the unique Perf_counters object
 */
inline Perf_counters& statistics() {return Perf_counters::getInstance();}

#endif

