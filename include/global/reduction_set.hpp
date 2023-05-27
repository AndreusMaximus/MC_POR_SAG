#ifndef GLOBAL_REDUCTION_SET
#define GLOBAL_REDUCTION_SET

#include <iostream>
#include <ostream>
#include <cassert>
#include <algorithm>
#include <map>
#include <numeric>
#include <cmath>

#include "jobs.hpp"

#include <chrono>

namespace NP
{

	namespace Global
	{

		template <class Time>
		class LST_container
		{
		public:
			std::vector<Time> eq;
			std::vector<Time> BIW;
			Time HPIW;
			Time LST;
			Time rmax;
			Time ref;
			Time prio;
			std::vector<const Job<Time> *> out_of_range_jobs;
			LST_container(std::vector<Time> base_BIW, Time base_HPIW, Time j_rmax, Time j_prio) : BIW{base_BIW}, HPIW{base_HPIW}, rmax{j_rmax}, prio{j_prio}
			{
				update_eq();
			}
			LST_container() {}

			void update_LST()
			{
				// if the high prio interference is larger than the largest eq
				// Fancy Fill algorithm
				// see whiteboard :*
				LST = rmax + ref + HPIW;
				// std::cout << "HPIW: " << HPIW << " LST " << LST_i << std::endl;

				for (int i = 0; i < BIW.size() - 1; i++)
				{
					if (HPIW >= eq[BIW.size() - 2 - i])
					{
						LST = rmax + ref + BIW[BIW.size() - 1 - i] + floor((HPIW - eq[BIW.size() - 2 - i]) / (BIW.size() - i));
						break;
					}
				}
			}

			void add_high_interference(const Job<Time> *j_i)
			{
				// if its earliest release time is further than our current LST, then we wont add it yet and just add it to our todo list.
				if (j_i->earliest_arrival() > LST)
				{
					insert_sorted(out_of_range_jobs, j_i,
								  [](const Job<Time> *i, const Job<Time> *j) -> bool
								  { return i->earliest_arrival() < j->earliest_arrival(); });
				}
				else
				{
					HPIW += j_i->maximal_cost();
					update_LST();
					for (const Job<Time> *j : out_of_range_jobs)
					{
						if (j->earliest_arrival() <= LST)
						{
							HPIW += j->maximal_cost();
							out_of_range_jobs.erase(out_of_range_jobs.begin());
						}
						else
						{
							update_LST();
							if (j->earliest_arrival() <= LST)
							{
								HPIW += j->maximal_cost();
								out_of_range_jobs.erase(out_of_range_jobs.begin());
							}
							else
							{
								update_LST();
								break;
							}
						}
					}
				}
			}

			void add_interference(const Job<Time> *j_i)
			{
				if (j_i->get_priority() > prio)
				{
					add_low_interference(j_i->maximal_cost());
				}
				else
				{
					add_high_interference(j_i);
				}
				update_LST();
			}

			void add_low_interference(Time low_interference)
			{
				if (low_interference > BIW[0])
				{
					// insert sorted and discard the lowest
					Time swap;
					BIW[0] = low_interference;
					for (int i = 0; i < BIW.size() - 1; i++)
					{
						if (BIW[i] > BIW[i + 1])
						{
							swap = BIW[i + 1];
							BIW[i + 1] = BIW[i];
							BIW[i] = swap;
						}
						else
						{
							break;
						}
					}
					ref = BIW[0];
					update_eq();
				}
			}

			void update_eq()
			{
				for (int i = 1; i < BIW.size(); i++)
				{
					eq.emplace_back(BIW[i] - BIW[i - 1]);
				}
				for (int i = 1; i < BIW.size() - 1; i++)
				{
					eq[i] = eq[i] * (1 + i) + eq[i - 1];
				}
			}
		};

		template <class Time>
		class lst_Job
		{
		public:
			const Job<Time> *J;
			Time LST;
			lst_Job(const Job<Time> *base_J, Time base_LST) : J{base_J}, LST{base_LST} {}
		};

		template <class Time>
		class Reduction_set
		{
		public:
			typedef std::vector<const Job<Time> *> Job_set;
			typedef std::vector<std::size_t> Job_precedence_set;
			typedef std::unordered_map<JobID, Time> Job_map;
			typedef typename Job<Time>::Priority Priority;
			typedef std::unordered_map<JobID, LST_container<Time>> LST_values;

		private:
			std::vector<Interval<Time>> cpu_availability;
			Job_set jobs;
			std::vector<std::size_t> indices;
			Job_precedence_set job_precedence_sets;
			// different sorted sets of jobs
			Job_set jobs_by_latest_arrival;
			Job_set jobs_by_earliest_arrival;

			// Job_set jobs_by_wcet;

			// these ordered sets are required for the reverse fill algorithm
			Job_set jobs_by_rmin_cmin;
			Job_set jobs_by_rmax_cmin;

			Time latest_busy_time;
			Time latest_idle_time;

			// this is the latest- latest start time for all jobs in the reduction set
			Time latest_LST;

			Job_map latest_start_times;
			hash_value_t key;
			Priority max_priority;
			unsigned long num_interfering_jobs_added;
			std::unordered_map<JobID, std::size_t> index_by_job;
			std::map<std::size_t, const Job<Time> *> job_by_index;
			Time total_count = 0;
			LST_values LST_map;

			bool deadline_miss;

		public:
			// CPU availability gaat dus wat anders worden a.d.h.v de meerdere cores
			Reduction_set(std::vector<Interval<Time>> cpu_availability, const Job_set &jobs, std::vector<std::size_t> &indices, const Job_precedence_set job_precedence_sets)
				: cpu_availability{cpu_availability},
				  jobs{jobs},
				  indices{indices},
				  job_precedence_sets{job_precedence_sets},
				  jobs_by_latest_arrival{jobs},
				  jobs_by_earliest_arrival{jobs},
				  // jobs_by_wcet{jobs},
				  jobs_by_rmax_cmin{jobs},
				  jobs_by_rmin_cmin{jobs},
				  key{0},
				  num_interfering_jobs_added{0},
				  index_by_job(),
				  job_by_index(),
				  LST_map()
			{
				// std::cout << "INITIALIZE NEW REDUCITON SET" << std::endl;
				deadline_miss = false;
				std::sort(jobs_by_latest_arrival.begin(), jobs_by_latest_arrival.end(),
						  [](const Job<Time> *i, const Job<Time> *j) -> bool
						  { if(i->latest_arrival() < j->latest_arrival()){
							return true;
						  }else if(i->latest_arrival() > j->latest_arrival()){
							return false;
						  }else{
							return i->get_priority() < j->get_priority();
						  }; });

				std::sort(jobs_by_earliest_arrival.begin(), jobs_by_earliest_arrival.end(),
						  [](const Job<Time> *i, const Job<Time> *j) -> bool
						  { return i->earliest_arrival() < j->earliest_arrival(); });

				// std::sort(jobs_by_wcet.begin(), jobs_by_wcet.end(),
				//		  [](const Job<Time> *i, const Job<Time> *j) -> bool
				//		  { return i->maximal_cost() < j->maximal_cost(); });

				// these sorts are required for the reverse fill algorithm
				std::sort(jobs_by_rmin_cmin.begin(), jobs_by_rmin_cmin.end(),
						  [](const Job<Time> *i, const Job<Time> *j) -> bool
						  { return i->earliest_arrival() + i->minimal_cost() < j->earliest_arrival() + j->minimal_cost(); });

				std::sort(jobs_by_rmax_cmin.begin(), jobs_by_rmax_cmin.end(),
						  [](const Job<Time> *i, const Job<Time> *j) -> bool
						  { return i->earliest_arrival() + i->maximal_cost() < j->earliest_arrival() + j->maximal_cost(); });

				// No clue yet why this exists but i dont want to remove it yet.
				for (int i = 0; i < jobs.size(); i++)
				{
					auto j = jobs[i];
					std::size_t idx = indices[i];

					index_by_job.emplace(j->get_id(), idx);
					job_by_index.emplace(std::make_pair(idx, jobs[i]));
				}

				// latest_busy_time = compute_latest_busy_time();
				latest_idle_time = compute_latest_idle_time();
				// compute_latest_start_times();
				compute_latest_start_time_complex();
				max_priority = compute_max_priority();
				initialize_key();
				// std::cout << "Candidate reduction set contains " << jobs.size() << " jobs" << std::endl;

				// for (const Job<Time>* j : jobs) {
				//	std::cout<<j->get_id()<<" ";
				// }
				// std::cout<<std::endl;
			}

			// For test purposes
			Reduction_set(std::vector<Interval<Time>> cpu_availability, const Job_set &jobs, std::vector<std::size_t> indices)
				: Reduction_set(cpu_availability, jobs, indices, {})
			{
			}

			int get_number_of_jobs()
			{
				return jobs.size();
			}

			void initialize_key()
			{
				for (const Job<Time> *j : jobs)
				{
					key = key ^ j->get_key();
				}
			}

			hash_value_t get_key() const
			{
				return key;
			}

			bool has_potential_deadline_misses() const
			{
				if (deadline_miss)
				{
					return true;
				}
				return false;
			}

			void add_job(const Job<Time> *jx, std::size_t index)
			{
				deadline_miss = false;
				num_interfering_jobs_added++;
				jobs.push_back(jx);

				index_by_job.emplace(jx->get_id(), index);
				job_by_index.emplace(std::make_pair(index, jobs.back()));
				indices.push_back(index);

				insert_sorted(jobs_by_latest_arrival, jx,
							  [](const Job<Time> *i, const Job<Time> *j) -> bool
							  { if(i->latest_arrival() < j->latest_arrival()){
							return true;
						  }else if(i->latest_arrival() > j->latest_arrival()){
							return false;
						  }else{
							return i->get_priority() < j->get_priority();
						  }; });
				insert_sorted(jobs_by_earliest_arrival, jx,
							  [](const Job<Time> *i, const Job<Time> *j) -> bool
							  { return i->earliest_arrival() < j->earliest_arrival(); });
				// insert_sorted(jobs_by_wcet, jx,
				//			  [](const Job<Time> *i, const Job<Time> *j) -> bool
				//			  { return i->maximal_cost() < j->maximal_cost(); });
				insert_sorted(jobs_by_rmax_cmin, jx,
							  [](const Job<Time> *i, const Job<Time> *j) -> bool
							  { return i->latest_arrival() + i->minimal_cost() < j->latest_arrival() + j->minimal_cost(); });
				insert_sorted(jobs_by_rmin_cmin, jx,
							  [](const Job<Time> *i, const Job<Time> *j) -> bool
							  { return i->earliest_arrival() + i->minimal_cost() < j->earliest_arrival() + j->minimal_cost(); });

				// latest_busy_time = compute_latest_busy_time();
				key = key ^ jx->get_key();
				// latest_idle_time = compute_latest_idle_time();
				//  compute_latest_start_times();
				// compute_latest_start_time_complex();

				if (!jx->priority_at_least(max_priority))
				{
					max_priority = jx->get_priority();
				}
			}

			void update_set()
			{
				latest_idle_time = compute_latest_idle_time();
				compute_latest_start_time_complex();
			}

			Job_set get_jobs() const
			{
				return jobs;
			}

			Time get_latest_busy_time() const
			{
				return latest_busy_time;
			}

			Time get_latest_idle_time() const
			{
				return latest_idle_time;
			}

			Time get_latest_LST()
			{
				return latest_LST;
			}

			Job_map get_latest_start_times() const
			{
				return latest_start_times;
			}

			Time get_latest_start_time(const Job<Time> &job) const
			{
				auto iterator = latest_start_times.find(job.get_id());
				return iterator == latest_start_times.end() ? -1 : iterator->second;
			}

			Time earliest_finish_time(const Job<Time> &job) const
			{
				return std::max(cpu_availability[0].min(), job.earliest_arrival()) + job.least_cost();
			}

			Time latest_finish_time(const Job<Time> &job) const
			{
				return get_latest_start_time(job) + job.maximal_cost();
			}

			Priority compute_max_priority() const
			{
				Priority max_prio{};

				for (const Job<Time> *j : jobs)
				{
					if (!j->priority_exceeds(max_prio))
					{
						max_prio = j->get_priority();
					}
				}

				return max_prio;
			}

			/*
			TODO: zorg dat de start_times globally geupdate worden en niet lokaal, hierdoor kunnen we bij de voorheen berekende waardes.
			*/

			// Hierin gaan we voor iedere job de LST berekenen en opslaan dan hoeven we dat niet voor iedere job te doen opnieuw
			void compute_latest_start_times()
			{
				latest_start_times = {};
				for (const Job<Time> *j : jobs_by_latest_arrival)
				{
					Time LST_j = compute_latest_start_time(j);
					latest_start_times.emplace(j->get_id(), LST_j);
				}

				return;
			}

			// Hierin gaan we voor iedere job de LST berekenen en opslaan dan hoeven we dat niet voor iedere job te doen opnieuw
			void compute_latest_start_time_complex()
			{
				// std::cout << "starting complex calculation :0" << std::endl;
				latest_start_times = {};
				// reserve the size for the latest start times preemtively
				latest_start_times.reserve(jobs.size());
				std::vector<lst_Job<Time>> LST_list;
				std::vector<const Job<Time> *> no_LST_list;
				typename std::vector<const Job<Time> *>::iterator lb_LPIW = jobs_by_earliest_arrival.begin();
				for (const Job<Time> *j : jobs_by_latest_arrival)
				{
					Time HPIW = 0;
					std::vector<Time> BIW = complexBIW(LST_list, no_LST_list, lb_LPIW, j, HPIW);
					Time LST_j = complexHPIW(HPIW, lb_LPIW, BIW, j);
					latest_start_times.emplace(j->get_id(), LST_j);
					if (j->exceeds_deadline(LST_j + j->maximal_cost()))
					{
						deadline_miss = true;
						return;
					}
				}

				return;
			}

			Time complexHPIW(Time HPIW, typename std::vector<const Job<Time> *>::iterator it_LPIW, std::vector<Time> BIW, const Job<Time> *j_i)
			{
				// ToDo friday:
				// add the EQ
				std::vector<Time> Ceq;
				for (int i = 1; i < cpu_availability.size(); i++)
				{
					Ceq.emplace_back(BIW[i] - BIW[i - 1]);
				}
				for (int i = 1; i < cpu_availability.size() - 1; i++)
				{
					Ceq[i] = Ceq[i] * (1 + i) + Ceq[i - 1];
				}
				Time ref = BIW[0];
				Time HPIW_i = HPIW;

				// initial LST using the HPIW from previous jobs
				Time LST_i = j_i->latest_arrival() + ref + HPIW;
				for (int i = 0; i < cpu_availability.size() - 1; i++)
				{
					if (HPIW >= Ceq[cpu_availability.size() - 2 - i])
					{
						LST_i = j_i->latest_arrival() + ref + BIW[cpu_availability.size() - 1 - i] + floor((HPIW - Ceq[cpu_availability.size() - 2 - i]) / (cpu_availability.size() - i));
						break;
					}
				}

				// we restart this loop where we ended when calculating the BIW
				// this means that the earliest arrival of all the upcoming jobs are >= than the current job.
				//	meaning that they shoudl never have an lst
				typename std::vector<const Job<Time> *>::iterator it_HPIW = it_LPIW;
				while (it_HPIW != jobs_by_earliest_arrival.end())
				{
					const Job<Time> *j_j = *it_HPIW;
					// only take high/same prio jobs into account
					if (j_j != j_i && j_j->get_priority() <= j_i->get_priority())
					{
						Time current_j_j_lst = j_j->latest_arrival();
						if (j_j->earliest_arrival() <= LST_i)
						{
							if (current_j_j_lst >= j_i->latest_arrival())
							{
								HPIW += j_j->maximal_cost();
								LST_i = j_i->latest_arrival() + ref + HPIW;

								for (int i = 0; i < cpu_availability.size() - 1; i++)
								{
									if (HPIW >= Ceq[cpu_availability.size() - 2 - i])
									{
										LST_i = j_i->latest_arrival() + ref + BIW[cpu_availability.size() - 1 - i] + floor((HPIW - Ceq[cpu_availability.size() - 2 - i]) / (cpu_availability.size() - i));
										break;
									}
								}
								//checking for deadline misses here enables us to discard reduction sets faster
								if (j_i->exceeds_deadline(LST_i + j_i->maximal_cost()))
									return LST_i;
							}
						}
						else
						{
							latest_LST = std::max(latest_LST, LST_i);
							return LST_i;
						}
					}
					it_HPIW++;
				}
				latest_LST = std::max(latest_LST, LST_i);
				return LST_i;
			}

			std::vector<Time> complexBIW(std::vector<lst_Job<Time>> &LST_list, std::vector<const Job<Time> *> &no_LST_list, typename std::vector<const Job<Time> *>::iterator &it_LPIW, const Job<Time> *j_i, Time &pre_calc_high)
			{
				// we will loop over this indefinately
				std::vector<Time> Cmax;
				// number of cores that we have
				int m = cpu_availability.size();
				for (int i = 0; i < m; i++)
				{
					// init this list with only zeroes
					Cmax.emplace_back((Time)0);
				}
				// first, this list with jobs that we have an LST for
				auto lst_it = LST_list.begin();
				while (lst_it != LST_list.end())
				{
					lst_Job<Time> &lst_job = *lst_it;
					// we have to go through the entire list, as we cannot sort i properly since that is reliant on the current j_i->latest_arrival value.
					// but we can remove jobs that we do not need anymore
					// namely jobs that do not cross the current j_i->latest_arrival value anymore.
					if (lst_job.LST + lst_job.J->maximal_cost() <= j_i->latest_arrival())
					{
						// update the iterator and continue
						lst_it = LST_list.erase(lst_it);
						continue;
					}
					// okay so now we know that the job still traverses the r_i^max
					if (lst_job.J->get_priority() > j_i->get_priority())
					{
						// if its a lower prio job, then we add it as follows to the Cmax
						Time C_aug = lst_job.LST + lst_job.J->maximal_cost() - (j_i->latest_arrival() - 1);
						if (C_aug > Cmax[0])
						{
							linear_inserter(Cmax, C_aug);
						}
					}
					else
					{
						// if its a high prio job, then we we can do two things
						if (lst_job.LST >= j_i->latest_arrival())
						{
							// if its larger or equal to the upcoming latest arrival, then add it to the HPIW already
							// we're allowed to do this as we know that the earliest release is before rimax and the LST is after, so it can interfere
							pre_calc_high += lst_job.J->maximal_cost();
						}
						else
						{
							// otherwise we can do the same kind of calculation as performed for the Low prio interference :)
							Time C_aug = lst_job.LST + lst_job.J->maximal_cost() - (j_i->latest_arrival() - 1);
							if (C_aug > Cmax[0])
							{
								linear_inserter(Cmax, C_aug);
							}
						}
					}
					lst_it++;
				}

				// okay now we just have to check all jobs that had a higher rmax than the previous job, and thus did not have an LST yet
				// std::cout << "we have " << LST_list.size() << " jobs without a LST" << std::endl;
				auto no_lst_it = no_LST_list.begin();
				while (no_lst_it != no_LST_list.end())
				{
					const Job<Time> *no_lst_job = *no_lst_it;

					if (no_lst_job == j_i)
					{
						no_lst_it++;
						continue;
					}
					// first check if it might have an LST now
					if (no_lst_job->latest_arrival() >= j_i->latest_arrival())
					{
						// now we can assume that it doesnt have an LST yet (high prio ones will have but does not matter at this point)
						if (no_lst_job->get_priority() <= j_i->get_priority())
						{
							// so if its a higher prio job, then we can add it to HPIW as well

							pre_calc_high += no_lst_job->maximal_cost();
						}
						else
						{
							// otherwise its a lower prio job, and thus can block maximally so lets check if its large enough to block
							Time C = no_lst_job->maximal_cost();
							if (C > Cmax[0])
							{
								linear_inserter(Cmax, C);
							}
						}
					}
					else
					{
						// arriving here means that the job has a rmax lower than the current rmax, and thus MUST have an LST
						Time lst = latest_start_times[no_lst_job->get_id()];
						if (lst + no_lst_job->maximal_cost() <= j_i->latest_arrival())
						{
							// if it will never interfere, then just remove it from the list

							no_lst_it = no_LST_list.erase(no_lst_it);
							continue;
						}
						if (no_lst_job->get_priority() <= j_i->get_priority())
						{
							if (lst >= j_i->latest_arrival())
							{
								// if it has a higher priority then that still means that it can interfere with the current job and therefor we add it to the HPIW and the other list
								pre_calc_high += no_lst_job->maximal_cost();
								// and higher prio jobs also have an LST already so lets move it to the other list
								LST_list.emplace_back(lst_Job<Time>(no_lst_job, lst));
								// and remove it from the current

								no_lst_it = no_LST_list.erase(no_lst_it);
								continue;
							}
							else
							{
								// otherwise we can do the same kind of calculation as performed for the Low prio interference :)
								Time C_aug = lst + no_lst_job->maximal_cost() - (j_i->latest_arrival() - 1);
								if (C_aug > Cmax[0])
								{
									linear_inserter(Cmax, C_aug);
								}
							}
						}
						else
						{
							// so if its a lower prio job then we can
							if (lst >= j_i->latest_arrival())
							{
								Time C = no_lst_job->maximal_cost();
								if (C > Cmax[0])
								{
									linear_inserter(Cmax, C);
								}
							}
							else
							{
								// otherwise we can still do the augmented one
								Time C_aug = lst + no_lst_job->maximal_cost() - (j_i->latest_arrival() - 1);
								if (C_aug > Cmax[0])
								{
									linear_inserter(Cmax, C_aug);
								}
							}
						}
					}
					no_lst_it++;
				}

				while (it_LPIW != jobs_by_earliest_arrival.end())
				{
					const Job<Time> *LPIW_job = *it_LPIW;
					if (LPIW_job == j_i)
					{
						// if we encounter ourselves, then just place it in the no LST list for future
						// std::cout<<"self encounter"<<std::endl;
						no_LST_list.emplace_back(LPIW_job);
						it_LPIW++;
						continue;
					}
					if (LPIW_job->earliest_arrival() >= j_i->latest_arrival())
					{
						// return where we've ended in this exploration so we can continue further at a later stage
						break;
					}
					// all jobs that we encounter here are not present in the lists above since we will start from where we ended last time, so we'll never encounter the same job twice

					if (LPIW_job->get_priority() <= j_i->get_priority())
					{
						if (LPIW_job->latest_arrival() >= j_i->latest_arrival())
						{
							pre_calc_high += LPIW_job->maximal_cost();
							// since its latest arrival is >= rimax, it wont have an LST so add it to the no LST list
							// std::cout<<"\t pushed hi job" << LPIW_job->get_id() << " to no LST list" << std::endl;
							no_LST_list.emplace_back(LPIW_job);
						}
						else
						{
							// since its latest arrival is < rimax it will have an LST, so it might be BIW
							Time high_LST = latest_start_times[LPIW_job->get_id()];
							if (high_LST + LPIW_job->maximal_cost() <= j_i->latest_arrival())
							{
								// we can instanly check if its a waste of a job, then we dont have to calculate further and just continue to the next one
								it_LPIW++;
								continue;
							}
							if (high_LST >= j_i->latest_arrival())
							{
								pre_calc_high += LPIW_job->maximal_cost();
							}
							else
							{
								// if its LST is < rimax, then do the augmented add
								Time C_aug = high_LST + LPIW_job->maximal_cost() - (j_i->latest_arrival() - 1);
								if (C_aug > Cmax[0])
								{
									linear_inserter(Cmax, C_aug);
								}
							}
							// since we have an LST we add it to the LST list
							// std::cout<<"\t pushed hi job " << LPIW_job->get_id() << " to LST list" << std::endl;
							LST_list.emplace_back(lst_Job<Time>(LPIW_job, high_LST));
						}
					}
					else
					{
						// if its a lower prio job then
						if (LPIW_job->latest_arrival() >= j_i->latest_arrival())
						{
							Time C = LPIW_job->maximal_cost();
							if (C > Cmax[0])
							{
								linear_inserter(Cmax, C);
							}
							// std::cout<<"\t pushed lo job" << LPIW_job->get_id() << " to no LST list" << std::endl;
							no_LST_list.emplace_back(LPIW_job);
						}
						else
						{
							// again this means that the latest arrival < rmax so it must ahve an LST
							Time low_LST = latest_start_times[LPIW_job->get_id()];
							if (low_LST + LPIW_job->maximal_cost() <= j_i->latest_arrival())
							{
								// we can instanly check if its a waste of a job, then we dont have to calculate further and just continue to the next one
								it_LPIW++;
								continue;
							}
							if (low_LST >= j_i->latest_arrival())
							{
								// if the LST is more than the current maximum release time then just add the maximum BIW to the list if possible
								Time C = LPIW_job->maximal_cost();
								if (C > Cmax[0])
								{
									linear_inserter(Cmax, C);
								}
							}
							else
							{
								// if its LST is < rimax, then do the augmented add
								Time C_aug = low_LST + LPIW_job->maximal_cost() - (j_i->latest_arrival() - 1);
								if (C_aug > Cmax[0])
								{
									linear_inserter(Cmax, C_aug);
								}
							}
							// add it to the LST list
							// std::cout<<"\t pushed lo job" << LPIW_job->get_id() << " to LST list" << std::endl;
							LST_list.emplace_back(lst_Job<Time>(LPIW_job, low_LST));
						}
					}

					it_LPIW++;
				}
				//  okay now we have the maximum interference that can be caused by jobs that start before j_i
				//  we just need to take the system availabilites into consideration now.
				//  Cmax  is sorted low -> high
				//  A^max is sorted low -> high
				for (int i = 0; i < m; i++)
				{
					// highest value from Cmax
					Time Cmax_i = Cmax[m - 1 - i];
					// compared to the lowest value from A^max
					Time A_max_i = cpu_availability[i].max();
					// compared to the current job offset
					// r_i^max
					Cmax[m - 1 - i] = std::max(j_i->latest_arrival() - 1 + Cmax_i, std::max(A_max_i, j_i->latest_arrival())) - j_i->latest_arrival();
				}
				// make sure its a ascending list
				std::sort(Cmax.begin(), Cmax.end());
				return Cmax;
			}

			void linear_inserter(std::vector<Time> &list, Time object)
			{
				Time swap;
				list[0] = object;
				int m = cpu_availability.size() - 1;
				for (int i = 0; i < m; i++)
				{
					if (list[i] > list[i + 1])
					{
						swap = list[i + 1];
						list[i + 1] = list[i];
						list[i] = swap;
					}
					else
					{
						return;
					}
				}
			}

			Time compute_latest_start_time(const Job<Time> *j_i)
			{

				//  Blocking interfering workload for job j_i

				Time Ceq[cpu_availability.size() - 1];
				Time BIW[cpu_availability.size()];
				std::vector<Time> _Ceq;
				std::vector<Time> _BIW;
				compute_m_blocking_interfering_workload(j_i, BIW);
				std::sort(BIW, BIW + cpu_availability.size());
				for (int i = 0; i < cpu_availability.size(); i++)
				{
					_BIW.emplace_back(BIW[i]);
				}
				//  std::cout<<"managed to do BIW"<<std::endl;
				// Time ref = BIW[cpu_availability.size()-1];
				Time ref = BIW[0];
				// std::cout << "AMAX = [";
				// for(int i = 0; i < cpu_availability.size(); i ++){
				//	std::cout<<cpu_availability[i].max()<<" ";
				// }
				// std::cout<<"]"<<std::endl;
				// std::cout << "REFERENCE POINT = " << ref << std::endl;
				// std::cout << "[";

				// for (int i = 0; i < cpu_availability.size(); i++)
				//{
				//	std::cout << BIW[i] << " ";
				// }
				// std::cout << "]" << std::endl;

				// std::cout << "[";
				for (int i = 1; i < cpu_availability.size(); i++)
				{
					// Ceq[i - 1] = BIW[i] - BIW[i - 1];
					_Ceq.emplace_back(BIW[i] - BIW[i - 1]);
					// std::cout << Ceq[i-1] << " ";
				}
				// std::cout << "]" << std::endl;
				//  last value in Ceq is the largest
				// std::cout << "[";
				for (int i = 1; i < cpu_availability.size() - 1; i++)
				{
					// Ceq[i] = Ceq[i] * (1 + i) + Ceq[i - 1];
					_Ceq[i] = Ceq[i] * (1 + i) + Ceq[i - 1];
					// std::cout << Ceq[i] << " ";
				}
				// std::cout << std::endl;

				Time LST_i = j_i->latest_arrival() + ref;
				// High priority interfering workload for job j_i
				Time HPIW = 0;

				// here we calculate the interference caused by jobs with a higher priority than j_i
				// std::cout << "HPIW: ";
				for (const Job<Time> *j_j : jobs_by_earliest_arrival)
				{
					// std::cout << "use " << j_j->get_id();
					if (j_j != j_i && j_j->get_priority() <= j_i->get_priority())
					{
						// for the improvement we will only count high prio jobs that have a LST >= r_i^max
						Time current_j_j_lst = (j_j->latest_arrival() >= j_i->latest_arrival()) ? j_j->latest_arrival() : latest_start_times[j_j->get_id()];
						// std::cout << " | " << current_j_j_lst << " - " << j_j->latest_arrival() << " : " << j_i->latest_arrival();
						if (j_j->earliest_arrival() <= LST_i)
						{
							if (current_j_j_lst >= j_i->latest_arrival())
							{
								// std::cout << " yes"<<std::endl;;
								HPIW += j_j->maximal_cost();
								// if the high prio interference is larger than the largest eq
								// Fancy Fill algorithm
								// see whiteboard :*
								LST_i = j_i->latest_arrival() + ref + HPIW;
								// std::cout << "HPIW: " << HPIW << " LST " << LST_i << std::endl;

								for (int i = 0; i < cpu_availability.size() - 1; i++)
								{
									if (HPIW >= _Ceq[cpu_availability.size() - 2 - i])
									{
										LST_i = j_i->latest_arrival() + ref + BIW[cpu_availability.size() - 1 - i] + floor((HPIW - _Ceq[cpu_availability.size() - 2 - i]) / (cpu_availability.size() - i));
										break;
									}
								}
							}
						}
						else
						{
							// std::cout << "\t LST for Job " << j_i->get_id() << " = " << LST << std::endl;
							latest_LST = std::max(latest_LST, LST_i);
							// std::cout << "HPIW: "<< HPIW << std::endl;
							// LST_map.emplace(j_i->get_id(), LST_container<Time>(_BIW, HPIW, j_i->latest_arrival(), j_i->get_priority()));
							return LST_i;
						}
					}
				}

				// std::cout << "\t UB on LST for Job " << j_i->get_id() << " = " << LST << std::endl;
				latest_LST = std::max(latest_LST, LST_i);
				// std::cout << "HPIW: "<< HPIW << std::endl;
				// std::cout << std::endl;
				// LST_map.emplace(j_i->get_id(), LST_container<Time>(_BIW, HPIW, j_i->latest_arrival(), j_i->get_priority()));
				return LST_i;
			}

			void show_time_waste()
			{
				std::cout << " current time waste indexing: " << total_count << "ns or " << total_count / 1000000000 << " s" << std::endl;
			}

			// Compute the blocking interfering workload as described in the paper
			void compute_m_blocking_interfering_workload(const Job<Time> *j_i, Time LPIW[])
			{
				// we need to remember the m largest values for the Low Priority Interfering Workload (LPIW)
				// static Time LPIW[cpu_availability.size()];
				for (int i = 0; i < cpu_availability.size(); i++)
				{
					LPIW[i] = 0;
				}

				for (const Job<Time> *j_j : jobs_by_earliest_arrival)
				{
					// Continue if we have the same job
					if (j_j == j_i)
					{
						// std::cout<<"cont."<<std::endl;
						continue;
					}
					// Since the list is sorted by earliest arrival we know that all jobs after this one will also no longer influence the BIW, so we can break here
					if (j_j->earliest_arrival() >= j_i->latest_arrival())
					{
						// std::cout<<"earliest arrival from j_j > latest from j_i "<<std::endl;
						break;
					}
					// If the j_j can arrive before j_i and has a lower priority then we consider it for LPIW
					// we only use the lower prio jobs here as the higher priority jobs are used later
					// we do not only count low prio jobs here anymore as we also will include partial high prio tasks
					// if (j_j->earliest_arrival() < j_i->latest_arrival() && j_j->get_priority() > j_i->get_priority())
					if (j_j->earliest_arrival() < j_i->latest_arrival())
					{
						Time actual_low_prio_interference = 0;
						// At this point i do not properly switch between calculating for low or high jobs
						if (j_j->get_priority() <= j_i->get_priority())
						{
							// high/equal prio jobs
							if (j_j->latest_arrival() >= j_i->latest_arrival())
							{
								// this means that it will be used for the HPIW calculation so the interference is 0
								actual_low_prio_interference = 0;
							}
							else
							{
								// this means that it has its latest release time before the job j_i
								// which in turn means that is MUST have an LST
								// but if its LST is also bigger than the current latest arrival, then we shouldnt count it
								auto s = std::chrono::high_resolution_clock::now();
								Time j_lst = latest_start_times[j_j->get_id()];
								auto e = std::chrono::high_resolution_clock::now();
								auto d = std::chrono::duration_cast<std::chrono::nanoseconds>(e - s);
								total_count += d.count();

								if (j_lst < j_i->latest_arrival())
								{
									actual_low_prio_interference = j_lst + j_j->maximal_cost() - (j_i->latest_arrival() - 1);
								}
							}
						}
						else
						{
							// low prio jobs
							if (j_j->latest_arrival() >= j_i->latest_arrival())
							{
								// this means that it is able to be released JUST before r_i^max so we suffer maximal interference
								actual_low_prio_interference = j_j->maximal_cost();
							}
							else
							{
								// this means that it has its latest release time before the job j_i
								// which in turn means that is MUST have an LST
								auto s = std::chrono::high_resolution_clock::now();
								Time j_lst = latest_start_times[j_j->get_id()];
								auto e = std::chrono::high_resolution_clock::now();
								auto d = std::chrono::duration_cast<std::chrono::nanoseconds>(e - s);
								total_count += d.count();

								if (j_lst > j_i->latest_arrival())
								{
									actual_low_prio_interference = j_j->maximal_cost();
								}
								else
								{
									actual_low_prio_interference = j_lst + j_j->maximal_cost() - (j_i->latest_arrival() - 1);
								}
							}
						}

						// std::cout<<actual_low_prio_interference<<std::endl;
						//  this is the old implementation so comment this out for the improved version
						//  actual_low_prio_interference = j_j->maximal_cost();
						//  std::cout << " actual low prio interference for " << j_j->get_id() << " on " << j_i->get_id() << " = " << actual_low_prio_interference << std::endl;
						//   Linear insert if it is larger than the smallest value in the LPIW array
						if (actual_low_prio_interference > LPIW[cpu_availability.size() - 1])
						{
							Time swap = 0;
							LPIW[cpu_availability.size() - 1] = actual_low_prio_interference;
							// std::cout<<"LPIW: ["<<LPIW[0]<<" ";
							for (int LPIW_i = cpu_availability.size() - 1; LPIW_i > 0; LPIW_i--)
							{
								if (LPIW[LPIW_i] > LPIW[LPIW_i - 1])
								{
									swap = LPIW[LPIW_i - 1];
									LPIW[LPIW_i - 1] = LPIW[LPIW_i];
									LPIW[LPIW_i] = swap;
								}
								/*else
								{
									break;
								}*/

								// std::cout<<LPIW[LPIW_i]<<" ";
							}
							// std::cout<<"]"<<std::endl;
						}
					}
				}

				// iterator integer for the array
				int LPIW_i = 0;
				// Calculate the blocking interfering workload using:
				// BIW_i = max(A_i^max, r_i^max, r_i^max -1 + LP_i)-r_i^max
				// as described in the paper
				for (Interval<Time> it : cpu_availability)
				{
					// it.max is the max of the interval, so A_x^max
					// std::cout<<it.max()<<" - "<< j_i->latest_arrival() << " - " <<j_i->latest_arrival() - 1 + LPIW[LPIW_i];
					LPIW[LPIW_i] = std::max(it.max(), std::max(j_i->latest_arrival(), j_i->latest_arrival() - 1 + LPIW[LPIW_i])) - j_i->latest_arrival();
					// std::cout<<"=>"<<LPIW[LPIW_i]<<std::endl;
					LPIW_i++;
				}
				// return LPIW;
			}

			void update_LSTs(const Job<Time> *j_i)
			{

				latest_idle_time = compute_latest_idle_time();
				Time LST = compute_latest_start_time(j_i);
				latest_start_times.emplace(j_i->get_id(), LST);

				if (j_i->exceeds_deadline(LST + j_i->maximal_cost()))
				{
					// std::cout<<"initial deadline miss for " << j->get_id() << std::endl;
					deadline_miss = true;
					return;
				}
				// okay now we've added the new job properly with an LST for it
				for (const Job<Time> *j : jobs)
				{
					if (j == j_i)
					{
						continue;
					}
					LST_map[j->get_id()].add_interference(j_i);
					Time updated_LST = LST_map[j->get_id()].LST;
					latest_start_times[j->get_id()] = updated_LST;
					latest_LST = std::max(latest_LST, updated_LST);
					if (j->exceeds_deadline(updated_LST + j->maximal_cost()))
					{
						// std::cout<<"initial deadline miss for " << j->get_id() << std::endl;
						deadline_miss = true;
						return;
					}
				}
			}

			Time compute_latest_idle_time()
			{
				// remove
				Time latest_idle_time{-1};
				// this is a counter to check how much jobs will always overlap with the release interval of the current job, if its > m then there will never be an interval at that point
				Time overflow = 0;
				// Save the m-1 largest jobs
				Time Cmax[cpu_availability.size() - 1];
				// Sum to add the remaining workload
				Time Crest = 0;
				// for convenience sake, remember the earliest time a core might become free
				Time Amin = cpu_availability.front().min();

				// std::cout << "Finding idle intervals..." << std::endl;
				// std::cout << jobs_by_latest_arrival.size() << std::endl;
				for (auto it_x = jobs_by_latest_arrival.rbegin(); it_x != jobs_by_latest_arrival.rend(); it_x++)
				{
					const Job<Time> *j_x = *it_x;
					// std::cout << "Looking for idle interval for job " << j_x->get_id() << " possibly ending at " << j_x->latest_arrival() << std::endl;
					//  reset the overflow and Crest values for each job under investigation
					overflow = 0;
					Crest = 0;
					int A_i = 0;
					// Compute the interference from possible previous states
					for (Interval<Time> A : cpu_availability)
					{
						// since we dont use the first value as it is always 0
						if (A_i != 0)
						{
							// Save the interfering workloads as the initial "highest workloads"
							Cmax[A_i - 1] = A.min() - Amin;
							// std::cout << "\t" << Cmax[A_i - 1] << std::endl;
						}
						A_i++;
					}

					// Now we again need to iterate over this list
					// i think i can set it_y to it_x, to limit the search
					for (auto it_y = jobs_by_latest_arrival.rbegin(); it_y != jobs_by_latest_arrival.rend(); it_y++)
					{
						const Job<Time> *j_y = *it_y;
						if (j_y->latest_arrival() < j_x->latest_arrival())
						{
							// this statement checks the overflow condition
							if (j_y->earliest_arrival() + j_y->minimal_cost() >= j_x->latest_arrival())
							{
								overflow++;
							}
							// if the current minimal cost is larger than the smallest saved value then add that value to Crest and insert the new value linearly
							if (j_y->minimal_cost() > Cmax[0])
							{
								Crest += Cmax[0];
								Cmax[0] = j_y->minimal_cost();
								Time swap;
								for (int i = 0; i < cpu_availability.size() - 2; i++)
								{
									if (Cmax[i] > Cmax[i + 1])
									{
										swap = Cmax[i + 1];
										Cmax[i + 1] = Cmax[i];
										Cmax[i] = swap;
									}
									else
									{
										break;
									}
								}
							}
							else
							{
								// if its smaller than all values in Cmax, just add it to Crest
								Crest += j_y->minimal_cost();
							}
						}
					}
					/*
					There is an interval if the workload of all jobs except the m-1 largest jobs spread over all cores
					is smaller than the latest release time of the current job under investigation
					*/
					if (Amin + ceil((double)Crest / (double)cpu_availability.size()) < j_x->latest_arrival() && overflow < cpu_availability.size())
					{

						// std::cout << "\t Latest interval might end at:" << j_x->latest_arrival() << std::endl;
						return j_x->latest_arrival();
					}
				}

				// std::cout << "\t No interval possible in this set" << std::endl;
				//  no jobs can be released at -1 so we can safely return 0 if there is no interval anywhere.
				return 0;
			}

			unsigned long get_num_interfering_jobs_added() const
			{
				return num_interfering_jobs_added;
			}

			bool can_interfere(const Job<Time> &job) const
			{
				// find_if geeft een pointer naar de eerste value die overeenkomt naar de eerset value in die lijst
				auto pos = std::find_if(jobs.begin(), jobs.end(),
										[&job](const Job<Time> *j)
										{ return j->get_id() == job.get_id(); });

				// als die pointer naar de laatste plek wijst, dan betekent dat dat deze dus nog niet in de lijst zit en daarom dus.. buiten de lijst zit
				if (pos != jobs.end())
				{
					// std::cout << "i already use job " << job.get_id() << std::endl;
					return false;
				}

				// rx_min < delta_M
				// latest idle time is het einde van de laatste idle interval in de set, dus die moet steeds herberekend
				// worden voor we nieuwe interfering jobs gaan zoekn
				if (job.earliest_arrival() <= latest_idle_time)
				{
					// std::cout << "Interference due to idle period " << job.get_id() << std::endl;
					return true;
				}

				// min_wcet wordt hier niet gebruikt dus waarom doen we dit?
				// Time min_wcet = min_lower_priority_wcet(job);

				// Zoek de laatste arrival time van deze set om snel te kijken of je daarbuiten valt
				Time max_arrival = jobs_by_latest_arrival.back()->latest_arrival();

				// Als je een lagere dan de minste prio hebben of buiten de set vallen dan hoeven we niet verder te zoeken
				if (!job.priority_exceeds(max_priority) && job.earliest_arrival() >= max_arrival)
				{
					return false;
				}

				// There exists a J_i s.t. rx_min <= LST^hat_i and p_x < p_i
				for (const Job<Time> *j : jobs)
				{
					// Voor iedere job, kijk of de huidige job een hogere prio heeft dan de job en of ie een hogere prio heeft dan de huidige job
					if (job.higher_priority_than(*j) && job.earliest_arrival() <= get_latest_start_time(*j))
					{
						// std::cout << "Job " << job.get_id() << " can be an interfering job for " << j->get_id() << std::endl;
						return true;
					}
				}

				return false;
			}

			void certainly_available_i(Time *Chigh, Time Clow, Time s, std::vector<Time> *CA_values)
			{
				Time tClow = Clow;
				if (cpu_availability.size() == 1)
				{
					CA_values->emplace_back(s + Clow + Chigh[0]);
					return;
				}

				if (Chigh[0] <= 0)
				{
					for (int i = 0; i < cpu_availability.size(); i++)
					{
						if (Chigh[i] > 0)
						{
							CA_values->emplace_back(s + Chigh[i]);
						}
						else
						{
							CA_values->emplace_back(s);
						}
					}
				}
				else
				{
					Time deltaC[cpu_availability.size() - 1];
					for (int delta_i = 0; delta_i < cpu_availability.size() - 1; delta_i++)
					{
						deltaC[delta_i] = Chigh[delta_i + 1] - Chigh[delta_i];
						// std::cout << Chigh[delta_i + 1] <<"-"<< Chigh[delta_i]<< "=" <<Chigh[delta_i + 1] - Chigh[delta_i] << deltaC[delta_i] << std::endl;;
					}
					for (int i = 1; i <= cpu_availability.size(); i++)
					{
						Time eq = Clow;
						for (int d = i - 1; d < cpu_availability.size(); d++)
						{
							if (eq > deltaC[d] * (i + (d - (i - 1))))
							{
								// we have enough to equalize
								eq -= deltaC[d] * (i + (d - (i - 1)));
							}
							else
							{
								// at this point we dont have enough to equalize
								CA_values->emplace_back(s + floor((double)eq / (double)(i + (d - (i - 1)))) + Chigh[d]);
								eq = 0;
								break;
							}
						}

						// std::cout<<"\t seg fault tester"<<std::endl;
						if (eq > 0)
						{
							// we have enough to properly equalize at this point
							CA_values->emplace_back(s + floor((double)eq / (double)cpu_availability.size()) + Chigh[cpu_availability.size() - 1]);
						}
						Clow += Chigh[i - 1];
					}
				}
			}

			std::vector<Time> compute_certainly_available()
			{
				// std::cout<<"\ncomputing CA"<<std::endl;
				std::vector<Time> CA_values;
				Time Chigh[cpu_availability.size()];
				Time Clow = 0;

				Time last_event = -1;
				Time s = std::max(cpu_availability[0].max(), jobs_by_latest_arrival[0]->latest_arrival());

				Time event = s;

				Time prev = 0;
				for (int i = 1; i < cpu_availability.size(); i++)
				{
					prev += cpu_availability[i].max() - cpu_availability[0].max();
				}
				prev = cpu_availability[0].max() + ceil((double)prev / (double)cpu_availability.size());

				// for some reason i need to have this print otherwise i have a seg fault?!?
				// std::cout<<"\t seg fault prevention?!?"<<std::endl;
				//  initialize the chigh values for the first iterative step
				for (int Chigh_i = 0; Chigh_i < cpu_availability.size(); Chigh_i++)
				{
					// std::cout<<cpu_availability[Chigh_i].max()<<" - ";
					Chigh[Chigh_i] = cpu_availability[Chigh_i].max() - s;
					// std::cout<<Chigh[Chigh_i]<<std::endl;
				}
				// std::cout<<"\tChihgh Filled"<<std::endl;
				for (const Job<Time> *j_x : jobs_by_latest_arrival)
				{

					if (j_x->latest_arrival() <= event && j_x->latest_arrival() > last_event)
					{
						// std::cout << j_x->get_id() << std::endl;
						if (j_x->maximal_cost() > Chigh[0])
						{
							if (Chigh[0] > 0)
								Clow += Chigh[0];
							Chigh[0] = j_x->maximal_cost();
							Time swap;
							// its in this point i think,
							for (int i = 0; i < cpu_availability.size() - 1; i++)
							{
								if (Chigh[i] > Chigh[i + 1])
								{
									swap = Chigh[i];
									Chigh[i] = Chigh[i + 1];
									Chigh[i + 1] = swap;
								}
								else
								{
									break;
								}
							}
						}
						else
						{
							Clow += j_x->maximal_cost();
						}

						if (Chigh[0] > 0)
						{
							event = s + ceil((double)(Clow + Chigh[0]) / (double)cpu_availability.size());
						}
						else
						{
							event = s;
						}
					}
					else
					{
						// std::cout << "\t new iteration" << std::endl;
						last_event = event;
						if (Chigh[0] > 0)
						{
							prev = Clow;
							for (int i = 0; i < cpu_availability.size(); i++)
							{
								prev += (Chigh[i] >= 0) ? Chigh[i] : 0;
							}
							prev = s + ceil((double)prev / (double(cpu_availability.size())));
						}
						// std::cout << "prev" << prev << std::endl;
						s = j_x->latest_arrival();

						for (int i = 0; i < cpu_availability.size(); i++)
						{

							if (Chigh[i] > 0)
							{
								Chigh[i] = (0 > ((event + Chigh[i]) - s)) ? prev - s : ((event + Chigh[i]) - s);
							}
							else
							{
								Chigh[i] = prev - s;
							}
							// std::cout << Chigh[i] << " ";
						}
						// std::cout << std::endl;
						//  reset Clow
						Clow = 0;
						if (j_x->maximal_cost() > Chigh[0])
						{
							Clow += Chigh[0];
							Chigh[0] = j_x->maximal_cost();
							Time swap = 0;
							for (int i = 0; i < cpu_availability.size() - 1; i++)
							{
								if (Chigh[i] > Chigh[i + 1])
								{
									swap = Chigh[i];
									Chigh[i] = Chigh[i + 1];
									Chigh[i + 1] = swap;
								}
								else
								{
									break;
								}
							}
						}
						else
						{
							Clow += j_x->maximal_cost();
						}
						event = s;
					}
				}
				certainly_available_i(Chigh, Clow, s, &CA_values);
				return CA_values;
			}

			void possibly_available_pre_filter(int *LFP, int i)
			{
				// create an empty temporary set
				Job_set J_temp{};
				int m = cpu_availability.size();

				// Initialize the LFP array
				for (int j = 0; j < m; j++)
				{
					LFP[j] = 0; // Set a default value here, or you can choose any appropriate initial value
				}
				// std::cout << "\t PA_" << i << " for " << m << " cores i have " << jobs_by_rmax_cmin.size() << " jobs to use" << std::endl;
				for (auto it_x = jobs_by_rmax_cmin.rbegin(); it_x != jobs_by_rmax_cmin.rend(); it_x++)
				{
					const Job<Time> *j_x = *it_x;

					if (J_temp.size() < m - i)
					{
						// remove all elementes from J_temp that exceed the latest finish time of the new job j_x as long as J_temp is not full
						Time newLFP = j_x->latest_arrival() + j_x->minimal_cost();
						// std::cout << "\tnew LFP to add is " << newLFP << std::endl;
						J_temp.erase(std::remove_if(J_temp.begin(),
													J_temp.end(),
													[&newLFP](const Job<Time> *j_y)
													{ return j_y->latest_arrival() > newLFP; }),
									 J_temp.end());
						// afterwards we just add j_x
						J_temp.emplace_back(j_x);
					}
					else
					{
						break;
					}
				}
				// Now that J_temp is filled, we can calculate the m-i latest finish points using only m-i cores
				int LFP_i = 0;
				// std::cout << "\tComputed " << J_temp.size() << " values for the LFP are" << std::endl;
				for (const Job<Time> *j_t : J_temp)
				{
					LFP[LFP_i] = j_t->latest_arrival() + j_t->minimal_cost();
					LFP_i++;
				}
			}

			Time possibly_available_i(int *LFP, int i)
			{
				int m = cpu_availability.size();
				// std::cout << m << "-" << i << "-" << 2 << "=" << m - i - 2 << std::endl;
				int Cmax[m - 1 - i];
				int Crest = 0;
				// just a really small sort, but can be improved upon. so TODO
				std::sort(LFP, LFP + (m - i));
				// This is the latest end point.
				int LEP = LFP[m - i - 1];
				if (m - i - 1 > 0)
				{
					// create the initial Cmax values;
					for (int Cmax_i = 0; Cmax_i < m - i - 1; Cmax_i++)
					{
						Cmax[Cmax_i] = LEP - LFP[Cmax_i];
					}
				}
				// std::cout<<"init CMax"<<std::endl;

				for (auto it_s = jobs_by_rmin_cmin.rbegin(); it_s != jobs_by_rmin_cmin.rend(); it_s++)
				{
					const Job<Time> *j_s = *it_s;
					if (j_s->latest_arrival() < LFP[m - i - 1])
					{
						// we should also check if A^max is less than that value i think.
						if (j_s->earliest_arrival() + j_s->minimal_cost() > LEP - ceil((double)Crest / (double)(m - i)))
						{
							// here we have to use the i'th core
							// std::cout << "\tWe have to use the i^th core at time " << j_s->earliest_arrival() + j_s->minimal_cost() << std::endl;
							return j_s->earliest_arrival() + j_s->minimal_cost();
						}
						// we cannot use Cmax if we have only one core to fill
						if (m - i - 1 > 0)
						{
							if (j_s->minimal_cost() < Cmax[m - i - 2])
							{
								Crest += j_s->minimal_cost();
							}
							else
							{
								// insert it into the Cmax array
								// we know for certain that Cmax[m-i-2] has to be added to Crest
								Crest += Cmax[m - i - 2];
								Cmax[m - i - 2] = j_s->minimal_cost();
								Time swap;
								for (int Cmax_i = m - i - 2; Cmax_i > 0; Cmax_i--)
								{
									if (Cmax[Cmax_i - 1] < Cmax[i])
									{
										swap = Cmax[Cmax_i - 1];
										Cmax[Cmax_i - 1] = Cmax[Cmax_i];
										Cmax[Cmax_i] = swap;
									}
									else
									{
										break;
									}
								}
							}
						}
						else
						{
							Crest += j_s->minimal_cost();
						}
					}
				}
				return cpu_availability[i - 1].min();
			}

			Time possibly_available_m()
			{
				Time Ctot = 0;
				Time maxEnd = 0;
				Time least_start_time = jobs[0]->earliest_arrival();

				for (const Job<Time> *j : jobs)
				{
					maxEnd = (maxEnd < j->earliest_arrival() + j->minimal_cost()) ? j->earliest_arrival() + j->minimal_cost() : maxEnd;
					Ctot += j->minimal_cost();
					least_start_time = std::min(least_start_time, j->earliest_arrival());
				}
				// this returns 0 + that value so the further we go the worse this bound gets
				// We have two options,  or we do it by setting the value to Amin + bound
				//						or we determine the earliest time that A job can start in this set and then take the maximum between that and the Amin
				return std::max(cpu_availability[0].min(), least_start_time) + (ceil((double)Ctot / (double)cpu_availability.size()) > maxEnd) ? ceil((double)Ctot / (double)cpu_availability.size()) : maxEnd;
			}

			std::vector<Time> compute_possibly_available()
			{
				int LFP[cpu_availability.size()];
				std::vector<Time> PA_values{};
				for (int i = 1; i < cpu_availability.size(); i++)
				{
					possibly_available_pre_filter(LFP, i);
					PA_values.push_back(possibly_available_i(LFP, i));
				}
				PA_values.push_back(possibly_available_m());

				return PA_values;
			}

			void compute_new_system_states()
			{
				std::cout << "\tCores are Possibly available at times: ";
				for (Time PA : compute_possibly_available())
				{
					std::cout << PA << " ";
				}
				std::cout << std::endl;
				std::cout << "\tCores are Certainly available at times: ";
				for (Time CA : compute_certainly_available())
				{
					std::cout << CA << " ";
				}
				std::cout << std::endl;
			}

			void created_set()
			{
				std::cout << "I would create the reduction set using the following jobs:" << std::endl;
				std::cout << "\t";
				for (const Job<Time> *j : jobs)
				{
					std::cout << j->get_id() << " ";
				}
				std::cout << std::endl;
				std::cout << "PA old";
				for (Interval<Time> it : cpu_availability)
				{
					std::cout << it.min() << " ";
				}
				std::cout << std::endl;
				std::cout << "CA old";
				for (Interval<Time> it : cpu_availability)
				{
					std::cout << it.max() << " ";
				}
				std::cout << std::endl;

				std::cout << "Where the new system states would be:" << std::endl;
				compute_new_system_states();
			}
		};

		template <typename T, typename Compare>
		typename std::vector<T>::iterator insert_sorted(std::vector<T> &vec, const T &item, Compare comp)
		{
			return vec.insert(std::upper_bound(vec.begin(), vec.end(), item, comp), item);
		}

		template <class Time>
		class Reduction_set_statistics
		{

		public:
			bool reduction_success;

			unsigned long num_jobs, num_interfering_jobs_added;

			std::vector<Time> priorities;

			Reduction_set_statistics(bool reduction_success, Reduction_set<Time> &reduction_set)
				: reduction_success{reduction_success},
				  num_jobs{reduction_set.get_jobs().size()},
				  num_interfering_jobs_added{reduction_set.get_num_interfering_jobs_added()},
				  priorities{}
			{
				for (const Job<Time> *j : reduction_set.get_jobs())
				{
					priorities.push_back(j->get_priority());
				}
			}
		};

	}

};

#endif