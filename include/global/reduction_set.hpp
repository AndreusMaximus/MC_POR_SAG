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

namespace NP
{

	namespace Global
	{

		template <class Time>
		class Reduction_set
		{
		public:
			typedef std::vector<const Job<Time> *> Job_set;
			typedef std::vector<std::size_t> Job_precedence_set;
			typedef std::unordered_map<JobID, Time> Job_map;
			typedef typename Job<Time>::Priority Priority;

		private:
			std::vector<Interval<Time>> cpu_availability;
			Job_set jobs;
			std::vector<std::size_t> indices;
			Job_precedence_set job_precedence_sets;
			// different sorted sets of jobs
			Job_set jobs_by_latest_arrival;
			Job_set jobs_by_earliest_arrival;
			Job_set jobs_by_wcet;

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

		public:
			// CPU availability gaat dus wat anders worden a.d.h.v de meerdere cores
			Reduction_set(std::vector<Interval<Time>> cpu_availability, const Job_set &jobs, std::vector<std::size_t> &indices,const Job_precedence_set job_precedence_sets)
				: cpu_availability{cpu_availability},
				  jobs{jobs},
				  indices{indices},
				  job_precedence_sets{job_precedence_sets},
				  jobs_by_latest_arrival{jobs},
				  jobs_by_earliest_arrival{jobs},
				  jobs_by_wcet{jobs},
				  jobs_by_rmax_cmin{jobs},
				  jobs_by_rmin_cmin{jobs},
				  key{0},
				  num_interfering_jobs_added{0},
				  index_by_job(),
				  job_by_index()
			{
				std::sort(jobs_by_latest_arrival.begin(), jobs_by_latest_arrival.end(),
						  [](const Job<Time> *i, const Job<Time> *j) -> bool
						  { return i->latest_arrival() < j->latest_arrival(); });

				std::sort(jobs_by_earliest_arrival.begin(), jobs_by_earliest_arrival.end(),
						  [](const Job<Time> *i, const Job<Time> *j) -> bool
						  { return i->earliest_arrival() < j->earliest_arrival(); });

				std::sort(jobs_by_wcet.begin(), jobs_by_wcet.end(),
						  [](const Job<Time> *i, const Job<Time> *j) -> bool
						  { return i->maximal_cost() < j->maximal_cost(); });

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
				latest_start_times = compute_latest_start_times();
				max_priority = compute_max_priority();
				// initialize_key();
				std::cout << "Candidate reduction set contains " << jobs.size() << " jobs" << std::endl;
			}

			// For test purposes
			Reduction_set(std::vector<Interval<Time>> cpu_availability, const Job_set &jobs, std::vector<std::size_t> indices)
				: Reduction_set(cpu_availability, jobs, indices, {})
			{
			}


			hash_value_t get_key() const {
				return key;
			}

			bool has_potential_deadline_misses() const {
				for (const Job<Time>* job : jobs) {
					if (job->exceeds_deadline(get_latest_start_time(*job) + job->maximal_cost())) {
						return true;
					}
				}

				return false;
			}

			void add_job(const Job<Time> *jx, std::size_t index)
			{

				num_interfering_jobs_added++;
				jobs.push_back(jx);

				index_by_job.emplace(jx->get_id(), index);
				job_by_index.emplace(std::make_pair(index, jobs.back()));
				indices.push_back(index);

				insert_sorted(jobs_by_latest_arrival, jx,
							  [](const Job<Time> *i, const Job<Time> *j) -> bool
							  { return i->latest_arrival() < j->latest_arrival(); });
				insert_sorted(jobs_by_earliest_arrival, jx,
							  [](const Job<Time> *i, const Job<Time> *j) -> bool
							  { return i->earliest_arrival() < j->earliest_arrival(); });
				insert_sorted(jobs_by_wcet, jx,
							  [](const Job<Time> *i, const Job<Time> *j) -> bool
							  { return i->maximal_cost() < j->maximal_cost(); });
				insert_sorted(jobs_by_rmax_cmin, jx,
							  [](const Job<Time> *i, const Job<Time> *j) -> bool
							  { return i->latest_arrival() + i->minimal_cost() < j->latest_arrival() + j->minimal_cost(); });
				insert_sorted(jobs_by_rmin_cmin, jx,
							  [](const Job<Time> *i, const Job<Time> *j) -> bool
							  { return i->earliest_arrival() + i->minimal_cost() < j->earliest_arrival() + j->minimal_cost(); });

				// latest_busy_time = compute_latest_busy_time();
				latest_idle_time = compute_latest_idle_time();
				latest_start_times = compute_latest_start_times();
				key = key ^ jx->get_key();

				if (!jx->priority_at_least(max_priority))
				{
					max_priority = jx->get_priority();
				}
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


			Time earliest_finish_time(const Job<Time> &job) const {
				return std::max(cpu_availability[0].min(), job.earliest_arrival()) + job.least_cost();
			}

			Time latest_finish_time(const Job<Time> &job) const {
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

			// Hierin gaan we voor iedere job de LST berekenen en opslaan dan hoeven we dat niet voor iedere job te doen opnieuw
			Job_map compute_latest_start_times()
			{
				// std::unordered_map<JobID, Priority> job_prio_map = preprocess_priorities();

				Job_map start_times{};
				for (const Job<Time> *j : jobs)
				{
					start_times.emplace(j->get_id(), compute_latest_start_time(j));
				}

				return start_times;
			}

			Time compute_latest_start_time(const Job<Time> *j_i)
			{
				std::cout << "Computing LST for " << j_i->get_id() << std::endl;
				// Blocking interfering workload for job j_i
				Time BIW = compute_blocking_interfering_workload(j_i);
				Time LST = j_i->latest_arrival() + floor(BIW / cpu_availability.size());
				// High priority interfering workload for job j_i
				Time HPIW = 0;

				// here we calculate the interference caused by jobs with a higher priority than j_i
				for (const Job<Time> *j_j : jobs_by_earliest_arrival)
				{
					if (j_j != j_i && j_j->get_priority() <= j_i->get_priority())
					{
						if (j_j->earliest_arrival() <= LST)
						{
							HPIW += j_j->maximal_cost();
							LST = j_i->latest_arrival() + floor((BIW + HPIW) / cpu_availability.size());
						}
						else
						{
							std::cout << "\t LST for Job " << j_i->get_id() << " = " << LST << std::endl;
							return LST;
						}
					}
				}

				std::cout << "\t UB on LST for Job " << j_i->get_id() << " = " << LST << std::endl;
				latest_LST = std::max(latest_LST, LST);
				return LST;
			}

			// Compute the blocking interfering workload as described in the paper
			Time compute_blocking_interfering_workload(const Job<Time> *j_i)
			{
				// summation for the final value
				Time BIW = 0;
				// we need to remember the m largest values for the Low Priority Interfering Workload (LPIW)
				Time LPIW[cpu_availability.size()];

				for (const Job<Time> *j_j : jobs_by_earliest_arrival)
				{
					// Continue if we have the same job
					if (j_j == j_i)
					{
						continue;
					}
					// Since the list is sorted by earliest arrival we know that all jobs after this one will also no longer influence the BIW, so we can break here
					if (j_j->earliest_arrival() > j_i->latest_arrival())
					{
						break;
					}
					// If the j_j can arrive before j_i and has a lower priority then we consider it for LPIW
					// we only use the lower prio jobs here as the higher priority jobs are used later
					if (j_j->earliest_arrival() < j_i->latest_arrival() && j_j->get_priority() > j_i->get_priority())
					{

						// Linear insert if it is larger than the smallest value in the LPIW array
						if (j_j->maximal_cost() > LPIW[cpu_availability.size() - 1])
						{
							bool updated = false;
							Time swap1 = 0;
							Time swap2 = 0;
							for (int LPIW_i = 0; LPIW_i < cpu_availability.size(); LPIW_i++)
							{
								if (updated)
								{
									swap2 = LPIW[LPIW_i];
									LPIW[LPIW_i] = swap1;
									swap1 = swap2;
								}
								else
								{
									if (LPIW[LPIW_i] < j_j->maximal_cost())
									{
										swap1 = LPIW[LPIW_i];
										LPIW[LPIW_i] = j_j->maximal_cost();
										updated = true;
									}
								}
							}
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
					BIW += std::max(it.max(), std::max(j_i->latest_arrival(), j_i->latest_arrival() - 1 - LPIW[LPIW_i])) - j_i->latest_arrival();
					LPIW_i++;
				}
				return BIW;
			}

			// Returns the smallest wcet among the jobs with a lower priority than job
			Time min_lower_priority_wcet(const Job<Time> &job) const
			{
				auto pos = std::find_if(jobs_by_wcet.begin(), jobs_by_wcet.end(),
										[&job](const Job<Time> *j)
										{ return job.higher_priority_than(*j); });

				if (pos == jobs_by_wcet.end())
				{
					return 0;
				}
				else
				{
					return (*pos)->maximal_cost();
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

				std::cout << "Finding idle intervals..." << std::endl;
				std::cout << jobs_by_latest_arrival.size() << std::endl;
				for (auto it_x = jobs_by_latest_arrival.rbegin(); it_x != jobs_by_latest_arrival.rend(); it_x++)
				{
					const Job<Time> *j_x = *it_x;
					std::cout << "Looking for idle interval for job " << j_x->get_id() << " possibly ending at " << j_x->latest_arrival() << std::endl;
					// reset the overflow and Crest values for each job under investigation
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
							std::cout << "\t" << Cmax[A_i - 1] << std::endl;
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
								for (int i = 0; i < cpu_availability.size() - 1; i++)
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
					if (Amin + ceil(Crest / cpu_availability.size()) < j_x->latest_arrival() && overflow < cpu_availability.size())
					{

						std::cout << "\t Latest interval might end at:" << j_x->latest_arrival() << std::endl;
						return j_x->latest_arrival();
					}
				}

				std::cout << "\t No interval possible in this set" << std::endl;
				// no jobs can be released at -1 so we can safely return 0 if there is no interval anywhere.
				return 0;
			}

			unsigned long get_num_interfering_jobs_added() const {
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
					std::cout << "i already use job " << job.get_id() << std::endl;
					return false;
				}

				// rx_min < delta_M
				// latest idle time is het einde van de laatste idle interval in de set, dus die moet steeds herberekend
				// worden voor we nieuwe interfering jobs gaan zoekn
				if (job.earliest_arrival() <= latest_idle_time)
				{
					std::cout << "Interference due to idle period " << job.get_id() << std::endl;
					return true;
				}

				// min_wcet wordt hier niet gebruikt dus waarom doen we dit?
				Time min_wcet = min_lower_priority_wcet(job);

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
						std::cout << "Job " << job.get_id() << " can be an interfering job for " << j->get_id() << std::endl;
						return true;
					}
				}

				return false;
			}

			void certainly_available_i(Time *Chigh, Time Clow, Time s, std::vector<Time> *CA_values)
			{
				Time tClow = Clow;
				if (Chigh[0] == 0)
				{
					for (int i = 0; i < cpu_availability.size(); i++)
					{
						CA_values->emplace_back(Chigh[i]);
					}
				}
				else
				{
					Time deltaC[cpu_availability.size() - 1];
					for (int delta_i = 0; delta_i < cpu_availability.size() - 1; delta_i++)
					{
						deltaC[delta_i] = Chigh[delta_i + 1] - Chigh[delta_i];
						//std::cout << Chigh[delta_i + 1] <<"-"<< Chigh[delta_i]<< "=" <<Chigh[delta_i + 1] - Chigh[delta_i] << deltaC[delta_i] << std::endl;;
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
								CA_values->emplace_back(s + floor(eq / (i + (d - (i - 1)))) + Chigh[d]);
								eq = 0;
								break;
							}
						}

						if (eq > 0)
						{
							// we have enough to properly equalize at this point
							CA_values->emplace_back(s + floor(eq / cpu_availability.size()) + Chigh[cpu_availability.size() - 1]);
						}
						Clow += Chigh[i - 1];
					}
				}
			}

			std::vector<Time> compute_certainly_available()
			{
				std::vector<Time> CA_values;
				Time Chigh[cpu_availability.size()];
				Time Clow = 0;

				Time event = 0;
				Time last_event = -1;
				Time s = cpu_availability[0].max();

				// initialize the chigh values for the first iterative step
				for (int Chigh_i = 0; Chigh_i < cpu_availability.size(); Chigh_i++)
				{
					Chigh[Chigh_i] = cpu_availability[Chigh_i].max() - cpu_availability[0].max();
				}

				for (const Job<Time> *j_x : jobs_by_latest_arrival)
				{
					if (j_x->latest_arrival() <= event && j_x->latest_arrival() > last_event)
					{
						if (j_x->maximal_cost() > Chigh[0])
						{
							Clow += Chigh[0];
							Chigh[0] = j_x->maximal_cost();
							Time swap;
							for (int i = 0; i < cpu_availability.size(); i++)
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
						event = s + ceil((Clow + Chigh[0]) / cpu_availability.size());
					}
					else
					{
						last_event = event;
						s = j_x->latest_arrival();

						for (int i = 0; i < cpu_availability.size(); i++)
						{
							Chigh[i] = (0 > ((event + Chigh[i]) - s)) ? 0 : ((event + Chigh[i]) - s);
						}
						// reset Clow
						Clow = 0;
						if (j_x->maximal_cost() > Chigh[0])
						{
							Clow += Chigh[0];
							Chigh[0] = j_x->maximal_cost();
							Time swap;
							for (int i = 0; i < cpu_availability.size(); i++)
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

				for (auto it_s = jobs_by_rmin_cmin.rbegin(); it_s != jobs_by_rmin_cmin.rend(); it_s++)
				{
					const Job<Time> *j_s = *it_s;
					if (j_s->latest_arrival() < LFP[m - i - 1])
					{
						// we should also check if A^max is less than that value i think.
						if (j_s->earliest_arrival() + j_s->minimal_cost() > LEP - ceil(Crest / (m - i)))
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
								Time swap;
								for (int i = m - i - 2; i > 0; i++)
								{
									if (Cmax[i - 1] < Cmax[i])
									{
										swap = Cmax[i - 1];
										Cmax[i - 1] = Cmax[i];
										Cmax[i] = swap;
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

				for (const Job<Time> *j : jobs)
				{
					maxEnd = (maxEnd < j->earliest_arrival() + j->minimal_cost()) ? j->earliest_arrival() + j->minimal_cost() : maxEnd;
					Ctot += j->minimal_cost();
				}
				return (ceil(Ctot / cpu_availability.size()) > maxEnd) ? ceil(Ctot / cpu_availability.size()) : maxEnd;
			}

			std::vector<Time> compute_possibly_available()
			{
				int LFP[cpu_availability.size()];
				std::vector<Time> PA_values;
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

			bool has_potential_deadline_misses(){
				return false;
			}

			void created_set(){
				std::cout<<"I would create the reduction set using the following jobs:"<<std::endl;
				std::cout<<"\t";
				for(const Job<Time> * j : jobs){
					std::cout<<j->get_id()<<" ";
				}
				std::cout<<std::endl;
				std::cout<<"Where the new system states would be:"<<std::endl;
				compute_new_system_states();
			}
		};

		template <typename T, typename Compare>
		typename std::vector<T>::iterator insert_sorted(std::vector<T> &vec, const T &item, Compare comp)
		{
			return vec.insert(std::upper_bound(vec.begin(), vec.end(), item, comp), item);
		}

		template<class Time> class Reduction_set_statistics {

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
				for (const Job<Time>* j : reduction_set.get_jobs()) {
					priorities.push_back(j->get_priority());
				}
			}
		};
	}
};

#endif