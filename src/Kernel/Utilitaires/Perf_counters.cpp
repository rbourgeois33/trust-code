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

	inline bool is_a_comm_counter()
	{
		return is_comm_;
	}

	void begin_count_(int counter_level, std::chrono::time_point<std::chrono::high_resolution_clock> t);

	void end_count_(int count_increment, double quantity_increment, std::chrono::time_point<std::chrono::high_resolution_clock> t_stop);
	//check std::chrono

	/*! @brief update variables : avg_time_per_step_ , min_time_per_step_ , max_time_per_step_ , sd_time_per_step_
	 *
	 */
	void compute_avg_min_max_var_per_step();

	std::array< std::array<double,4> ,2> compute_min_max_avg_sd();

protected:
	std::string description_;
	int level_;
	std::string family_ ;
	bool is_comm_;
	int is_running_; // 0 if not , 1 if yes and 2 if was running before stop (for Perf_counters::restart_counters)
	bool start_at_the_beginning_of_the_time_step_;
	int count_;
	double quantity_;
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
	is_running_ = 0;
	start_at_the_beginning_of_the_time_step_ = false;
	count_ = 0;
	quantity_ = 0.0;
	avg_time_per_step_ = 0.0 ;
	min_time_per_step_ = 0.0 ;
	max_time_per_step_ = 0.0 ;
	sd_time_per_step_ = 0.0 ;
	time_timestep_ = std::chrono::duration<double>::zero();
	total_time_ = std::chrono::duration<double>::zero() ;
	last_open_time_ = std::chrono::time_point<std::chrono::high_resolution_clock>();
}

Counter::Counter()
{
	description_ = "";
	level_ = -10;
	family_ = "";
	is_comm_ = false;
	is_running_ = 0;
	start_at_the_beginning_of_the_time_step_ = false;
	count_ = 0;
	quantity_ = 0.0;
	avg_time_per_step_ = 0.0 ;
	min_time_per_step_ = 0.0 ;
	max_time_per_step_ = 0.0 ;
	sd_time_per_step_ = 0.0 ;
	time_timestep_ = std::chrono::duration<double>::zero();
	total_time_ = std::chrono::duration<double>::zero() ;
	last_open_time_ = std::chrono::time_point<std::chrono::high_resolution_clock>();
}

void Counter::begin_count_(int counter_level, std::chrono::time_point<std::chrono::high_resolution_clock> t)
{
	if (counter_level != level_)
	{
		Cerr<< "You did not specified the expected counter lvl of : :" <<description_ <<std::endl;
		level_ = counter_level;
	}
	if (is_running_==1)
	{
		Cerr<< "Try to start an already launch counter named :" <<description_ <<std::endl;
		break;
	}
	if (start_at_the_beginning_of_the_time_step_ == false)
	{
		start_at_the_beginning_of_the_time_step_ = true;
	}

	last_open_time_ = t;
	is_running_ = 1;

#ifdef PETSCKSP_H
	Nom info(si.description[id]);
	info+="\n";
	PetscInfo(0,"%s",info.getChar());
#endif

#ifdef TRUST_USE_CUDA
	// Level 1 only to avoid MPI calls
	if (si.counter_level[id]==1)
		nvtxRangePush(si.description[id]);
#endif
}

void Counter::end_count_(int count_increment, double quantity_increment, std::chrono::time_point<std::chrono::high_resolution_clock> t_stop)
{
	assert(is_running_);
	std::chrono::duration<double> t_tot = t_stop-last_open_time_;
	std::chrono::duration<double> t_alone = t_stop - last_open_time_alone_;
	is_running_=false;
	quantity_ += quantity_increment;
	total_time_ += t_tot;
	time_alone_ += t_alone;
	count_ += count_increment;

#ifdef TRUST_USE_CUDA
	// Level 1 only to avoid MPI calls
	if (si.counter_level[id]==1) nvtxRangePop();
#endif
}

std::array< std::array<double,4> ,2> Counter::compute_min_max_avg_sd()
{
	assert(Process::is_parallel());
	double time = total_time_;
	double quantity = quantity_;
	double min_time, max_time, avg_time, sd_time, min_quantity, max_quantity, avg_quantity, sd_quantity ;

	min_time = Process::mp_min(time);
	max_time = Process::mp_max(time);
	avg_time = Process::mp_sum(time)/Process::nproc();
	sd_time = sqrt(Process::mp_sum((time-avg_time)*(time-avg_time))/Process::nproc());

	std::array<double,4> min_max_avg_sd_time = {min_time,max_time,avg_time,sd_time};

	min_quantity = Process::mp_min(quantity);
	max_quantity = Process::mp_max(quantity);
	avg_quantity = Process::mp_sum(quantity)/Process::nproc();
	sd_quantity = sqrt(Process::mp_sum((quantity-avg_quantity)*(quantity-avg_quantity))/Process::nproc());

	std::array<double,4> min_max_avg_sd_quantity = {min_quantity,max_quantity,avg_quantity,sd_quantity};

	return {min_max_avg_sd_time,min_max_avg_sd_quantity};
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
void Perf_counters::begin_count(STD_COUNTERS std_cnt, unsigned int counter_lvl)
{
	std::chrono::time_point<std::chrono::high_resolution_clock> t = std::chrono::high_resolution_clock::now();
	if (counter_lvl != current_highest_counter_lvl_ +1)
	{
		Cerr << "Try to start a counter whose lvl is lower than the one already open"<< std::endl;
		break;
	}
	Counter & c = std_counters_[std_cnt];
	if(!c.is_running_)
	{
		Cerr << "Trying to start an already started counter"<<std::endl;
		break;
	}
	c.begin_count_(counter_lvl,t);
	current_highest_counter_lvl_ = counter_lvl;
	if (counter_lvl>0)
	{
		if (!running_counters_.empty())
		{
			*running_counters_.back()->time_alone_ += t - *running_counters_.back()->last_open_time_alone_ ;
		}
		running_counters_.push_back(&c);
	}
}

/*! Custom counters
 *
 * @param custom_count_name
 */
void Perf_counters::begin_count(const std::string& custom_count_name, unsigned int counter_lvl)
{
	std::chrono::time_point<std::chrono::high_resolution_clock> t = std::chrono::high_resolution_clock::now();
	if (counter_lvl != current_highest_counter_lvl_ +1)
	{
		Cerr << "Try to start a counter whose lvl is lower than the one already open"<< std::endl;
		break;
	}
	Counter & c = custom_counter_map_str_to_counter_.at(custom_count_name);
	if(!c.is_running_)
	{
		Cerr << "Trying to start an already started counter"<<std::endl;
		break;
	}
	c.begin_count_(counter_lvl, t);
	current_highest_counter_lvl_ = counter_lvl;
	if (counter_lvl>0)
	{
		if (!running_counters_.empty())
		{
			*running_counters_.back()->time_alone_ += t - *running_counters_.back()->last_open_time_alone_ ;
		}
		running_counters_.push_back(&c);
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
	std::chrono::time_point<std::chrono::high_resolution_clock> t = std::chrono::high_resolution_clock::now();
	Counter & c = std_counters_[std_cnt];

	if (&c!=running_counters_.back())
	{
		Cerr << "You are trying to close a different counter than the one last opened"<<std::endl;
	}

	if (c.level_ != current_highest_counter_lvl_)
	{
		Cerr << "Try to start a counter whose lvl is lower than the one already open"<< std::endl;
		break;
	}
	if(!c.is_running_)
	{
		Cerr << "Trying to close an already closed counter"<<std::endl;
		break;
	}
	c.end_count_(count_increment, quantity_increment,t);
	current_highest_counter_lvl_ = c.level_ -1;
	if (c.level_>0)
	{
		if (!running_counters_.empty())
		{
			running_counters_.pop_back();
			*running_counters_.back()->last_open_time_alone_ = t;
		}
	}
}


/*! @brief End the count of a counter and update the counter values
 *
 * @param c is the counter to end the count
 * @param count_increment is the count increment. If not specified, then it is equal to 1
 * @param quantity_increment is the increment of custom variable quantity. If not specified, it is set to 0.
 */
void Perf_counters::end_count(const std::string& custom_count_name, int count_increment, double quantity_increment)
{
	std::chrono::time_point<std::chrono::high_resolution_clock> t = std::chrono::high_resolution_clock::now();
	Counter & c = custom_counter_map_str_to_counter_.at(custom_count_name);

	if (&c!=running_counters_.back())
	{
		Cerr << "You are trying to close a different counter than the one last opened"<<std::endl;
		break;
	}

	if (c.level_ != current_highest_counter_lvl_)
	{
		Cerr << "Try to start a counter whose lvl is lower than the one already open"<< std::endl;
		break;
	}
	if(!c.is_running_)
	{
		Cerr << "Trying to close an already closed counter"<<std::endl;
		break;
	}
	c.end_count_(count_increment, quantity_increment,t);
	current_highest_counter_lvl_ = c.level_ -1;
	if (c.level_>0)
	{
		if (!running_counters_.empty())
		{
			running_counters_.pop_back();
			*running_counters_.back()->last_open_time_alone_ = t;
		}
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
	int step = three_first_steps_elapsed_ ? tstep - 3 : tstep;
	for (Counter &c: std_counters_)
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
	for (Counter &c: std_counters_)
	{
		if (c.is_running_)
		{
			time_elapsed_before_stop= t_stop -c.last_open_time_;
			c.time_timestep_ += time_elapsed_before_stop;
			c.total_time_ += time_elapsed_before_stop;
			c.is_running_ = 2;
		}
	}
	for (Counter &c: custom_counter_map_str_to_counter_)
	{
		if (c.is_running_)
		{
			time_elapsed_before_stop = t_stop - c.last_open_time_;
			c.time_timestep_ += time_elapsed_before_stop;
			c.total_time_ += time_elapsed_before_stop;
			c.is_running_ = 2;
		}
	}
}

/*! @brief Restart all counters, has to be called on every processor simultaneously
 *
 */
void Perf_counters::restart_counters()
{
	std::chrono::time_point<std::chrono::high_resolution_clock> t_restart = std::chrono::high_resolution_clock::now();
	for (Counter &c: std_counters_)
	{
		if (c.is_running_ == 2)
		{
			c.last_open_time_ = t_restart;
			c.is_running_ = 1;
		}
	}
	for (Counter &c: custom_counter_map_str_to_counter_)
	{
		if (c.is_running_ == 2)
		{
			c.last_open_time_ = t_restart;
			c.is_running_ = 1;
		}
	}
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
		File_header << "# Number of processor used = " << Process::nproc() << std::endl ;
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

	if (total_time<=std::chrono::duration<double>::zero())
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
			std::array< std::array<double,4> ,2> table = c.compute_min_max_avg_sd();
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

	for (Counter &c: custom_counter_map_str_to_counter_)
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

		if (Process::is_parallel() && Process::je_suis_maitre())
		{
			std::array< std::array<double,4> ,2> table = c.compute_min_max_avg_sd();
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


