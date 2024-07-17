/****************************************************************************
* Copyright (c) 2023, CEA
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
// This file contains all of the needed for the description of the counter associated with the tracking of performance in the TRUST code.

class counter
{
public:
  counter();
  ~counter();
  int get_counter_number(const counter c);
  int get_counter_name(const counter c);
  bool is_a_comm_counter(const counter c);
  void create_new_counter(counter c);
  void change_counter_level(counter c, int new_counter_lvl);
  std::tuple <int, std::string, int, std::string, bool, double, double, double, double, double, double> get_counter_statistics(const counter c);
  std::map <std::string, counter> name_to_counter;

protected:
  int counter_number_;
  std::string counter_name_;
  int counter_level_;
  std::string counter_family_ ;
  bool is_comm_;
  bool is_running;
  double count_;
  double quantity_;
  double current_time_;
  double first_call_time_;
  double total_time_;
  double avg_time_per_step_;
  double min_time_per_step_;
  double max_time_per_step_;
}
  ;


/*!  @brief This class store counters in TRUST
 *
 */

class Perf_counters
{
public:
  Perf_counters();
  ~Perf_counters();

  std::vector<counter> declare_base_counters();

  std::vector<map<std::string,count>> vector_map_name_to_counter;

  counter create_custom_counter();

  /*! @brief Start the count of a counter
     *
     */
  void begin_count(const std::string counter_name);

  /*! @brief End the count of a counter and update the counter values
     *
     */
  void end_count(const std::string counter_name, int count_increment, double quantity_increment);


protected:
  int Number_of_base_counters_;
  int Number_of_counters_;
  std::vector <counter> list_of_counters_ ;
};


