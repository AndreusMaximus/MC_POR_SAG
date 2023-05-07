#ifndef MERGE_SET_H
#define MERGE_SET_H

#include <iostream>
#include <ostream>
#include <cassert>
#include <algorithm>
#include <map>
#include <numeric>

#include "jobs.hpp"

namespace NP {

	namespace Uniproc {

		template<class Time> class Reduction_set {
		public:

			typedef std::vector<const Job<Time>*> Job_set;
			typedef std::vector<std::size_t> Job_precedence_set;
			typedef std::unordered_map<JobID, Time> Job_map;
			typedef typename Job<Time>::Priority Priority;

		private:

			Interval<Time> cpu_availability;
			Job_set jobs;
			std::vector<std::size_t> indices;
			std::vector<Job_precedence_set> job_precedence_sets;
			Job_set jobs_by_latest_arrival;
			Job_set jobs_by_earliest_arrival;
			Job_set jobs_by_wcet;
			Time latest_busy_time;
			Time latest_idle_time;
			Job_map latest_start_times;
			hash_value_t key;
			Priority max_priority;
			unsigned long num_interfering_jobs_added;
			std::unordered_map<JobID, std::size_t> index_by_job;
			std::map<std::size_t, const Job<Time>*> job_by_index;

		public:

			Reduction_set(Interval<Time> cpu_availability, const Job_set &jobs, std::vector<std::size_t> &indices, const std::vector<Job_precedence_set> &job_precedence_sets)
			: cpu_availability{cpu_availability},
			jobs{jobs},
			indices{indices},
			job_precedence_sets{job_precedence_sets},
			jobs_by_latest_arrival{jobs},
			jobs_by_earliest_arrival{jobs},
			jobs_by_wcet{jobs},
			key{0},
			num_interfering_jobs_added{0},
			index_by_job(),
			job_by_index()
			{
				std::sort(jobs_by_latest_arrival.begin(), jobs_by_latest_arrival.end(),
						  [](const Job<Time>* i, const Job<Time>* j) -> bool { return i->latest_arrival() < j->latest_arrival(); });

				std::sort(jobs_by_earliest_arrival.begin(), jobs_by_earliest_arrival.end(),
						  [](const Job<Time>* i, const Job<Time>* j) -> bool { return i->earliest_arrival() < j->earliest_arrival(); });

				std::sort(jobs_by_wcet.begin(), jobs_by_wcet.end(),
						  [](const Job<Time>* i, const Job<Time>* j) -> bool { return i->maximal_cost() < j->maximal_cost(); });

				for (int i = 0; i < jobs.size(); i++) {
					auto j = jobs[i];
					std::size_t idx = indices[i];

					index_by_job.emplace(j->get_id(), idx);
					job_by_index.emplace(std::make_pair(idx, jobs[i]));
				}

				latest_busy_time = compute_latest_busy_time();
				latest_idle_time = compute_latest_idle_time();
				latest_start_times = compute_latest_start_times();
				max_priority = compute_max_priority();
				initialize_key();

			}

			// For test purposes
            Reduction_set(Interval<Time> cpu_availability, const Job_set &jobs, std::vector<std::size_t> indices)
                    : Reduction_set(cpu_availability, jobs, indices, {})
            {}

			Job_set get_jobs() const {
				return jobs;
			}

			Time get_latest_busy_time() const {
				return latest_busy_time;
			}

			Time get_latest_idle_time() const {
				return latest_idle_time;
			}

			Job_map get_latest_start_times() const {
				return latest_start_times;
			}

			Time get_latest_start_time(const Job<Time> &job) const {
				auto iterator = latest_start_times.find(job.get_id());
				return iterator == latest_start_times.end() ? -1 : iterator->second;
			}

			bool has_potential_deadline_misses() const {
				for (const Job<Time>* job : jobs) {
					if (job->exceeds_deadline(get_latest_start_time(*job) + job->maximal_cost())) {
						return true;
					}
				}

				return false;
			}

			//interference check maar dan met precendence constraints
			bool can_interfere(const Job<Time> &job, const Job_precedence_set &job_precedence_set, const Index_set &scheduled_jobs) {
				if (! job_satisfies_precedence_constraints(job_precedence_set, scheduled_jobs)) {
					return false;
				}

				return can_interfere(job);
			}

			void add_job(const Job<Time>* jx, std::size_t index) {
				num_interfering_jobs_added++;
				jobs.push_back(jx);

				index_by_job.emplace(jx->get_id(), index);
				job_by_index.emplace(std::make_pair(index, jobs.back()));
                indices.push_back(index);

				insert_sorted(jobs_by_latest_arrival, jx,
							  [](const Job<Time>* i, const Job<Time>* j) -> bool { return i->latest_arrival() < j->latest_arrival(); });
				insert_sorted(jobs_by_earliest_arrival, jx,
							  [](const Job<Time>* i, const Job<Time>* j) -> bool { return i->earliest_arrival() < j->earliest_arrival(); });
				insert_sorted(jobs_by_wcet, jx,
							  [](const Job<Time>* i, const Job<Time>* j) -> bool { return i->maximal_cost() < j->maximal_cost(); });

				latest_busy_time = compute_latest_busy_time();
				latest_idle_time = compute_latest_idle_time();
				latest_start_times = compute_latest_start_times();
				key = key ^ jx->get_key();

				if (!jx->priority_at_least(max_priority)) {
					max_priority = jx->get_priority();
				}
			}

			Time earliest_start_time() const {
				return std::max(cpu_availability.min(), (*jobs_by_earliest_arrival.begin())->earliest_arrival());
			}

			Time earliest_finish_time() const {
				Time t = cpu_availability.min();

				for(const Job<Time>* j : jobs_by_earliest_arrival) {
					t = std::max(t, j->earliest_arrival()) + j->least_cost();
				}

				return t;
			}

			Time earliest_finish_time(const Job<Time> &job) const {
				return std::max(cpu_availability.min(), job.earliest_arrival()) + job.least_cost();
			}

			Time latest_start_time() const {
				return std::max(cpu_availability.max(), (*jobs_by_latest_arrival.begin())->latest_arrival());
			}

			Time latest_finish_time(const Job<Time> &job) const {
				return get_latest_start_time(job) + job.maximal_cost();
			}

			Time get_min_wcet() const {
				return jobs_by_wcet.front()->maximal_cost();
			}

			hash_value_t get_key() const {
				return key;
			}

			unsigned long get_num_interfering_jobs_added() const {
				return num_interfering_jobs_added;
			}

		private:

			bool job_satisfies_precedence_constraints(const Job_precedence_set &job_precedence_set, const Index_set &scheduled_jobs) {
				if (job_precedence_set.empty()) {
					return true;
				}

				Index_set scheduled_union_reduction_set{scheduled_jobs};

				for (auto idx : indices) {
					scheduled_union_reduction_set.add(idx);
				}

				Index_set predecessor_indices{};

				for (auto idx : job_precedence_set) {
					predecessor_indices.add(idx);
				}
                // ?????
                auto condition1 = scheduled_union_reduction_set.includes(job_precedence_set);
                auto condition2 = !predecessor_indices.includes(indices);
//				return scheduled_union_reduction_set.includes(job_precedence_set) && !predecessor_indices.includes(indices);
                return condition1 && condition2;
			}

			bool can_interfere(const Job<Time> &job) const {
				auto pos = std::find_if(jobs.begin(), jobs.end(),
										[&job] (const Job<Time>* j) { return j->get_id() == job.get_id(); } );

				// A job can't interfere with itself
				if(pos != jobs.end()) {
					return false;
				}

				// rx_min < delta_M
				if (job.earliest_arrival() <= latest_idle_time) {
					return true;
				}

				Time min_wcet = min_lower_priority_wcet(job);
				Time max_arrival = jobs_by_latest_arrival.back()->latest_arrival();

				if (!job.priority_exceeds(max_priority) && job.earliest_arrival() >= max_arrival) {
					return false;
				}

				// There exists a J_i s.t. rx_min <= LST^hat_i and p_x < p_i
				for (const Job<Time>* j : jobs) {
					if (job.earliest_arrival() <= get_latest_start_time(*j) && job.higher_priority_than(*j)) {
						return true;
					}
				}

				return false;
			}

			Time compute_latest_busy_time() {
				Time t = cpu_availability.max();

				for(const Job<Time>* j : jobs_by_latest_arrival) {
					t = std::max(t, j->latest_arrival()) + j->maximal_cost();
				}

				return t;
			}

			Job_map compute_latest_start_times() {
				std::unordered_map<JobID, Priority> job_prio_map = preprocess_priorities();

				Job_map start_times{};
				for (const Job<Time>* j : jobs) {
					start_times.emplace(j->get_id(), compute_latest_start_time(*j, job_prio_map));
				}

				return start_times;
			}

			// Preprocess priorities for s_i by setting priority of each job to the lowest priority of its predecessors
			std::unordered_map<JobID, Priority> preprocess_priorities() {
				std::unordered_map<JobID, Priority> job_prio_map{};
//				???

                if (!job_precedence_sets.empty()) {
                    auto topo_sorted_jobs = topological_sort<Time>(job_precedence_sets, jobs);
                    for (auto j: topo_sorted_jobs) {
                        const Job_precedence_set &preds = job_precedence_sets[index_by_job.find(j.get_id())->second];
                        Priority max_pred_prio = 0;

                        for (auto pred_idx: preds) {
                            auto iterator = job_by_index.find(pred_idx);

                            // We ignore all predecessors outside the reduction set
                            if (iterator == job_by_index.end()) {
                                continue;
                            }

                            auto pred = iterator->second;
                            max_pred_prio = std::max(max_pred_prio, pred->get_priority());
                        }

                        Priority p = std::max(j.get_priority(), max_pred_prio);
                        job_prio_map.emplace(j.get_id(), p);
                    }
                }else{
                //in this case we don't have precedence constraints, so we just set the priority to the job's own priority
                    for (auto j: jobs) {
                        job_prio_map.emplace(j->get_id(), j->get_priority());
                    }
                }
				return job_prio_map;
			}

			Time compute_latest_start_time(const Job<Time> &i, const std::unordered_map<JobID, Priority> &job_prio_map) {
				Time s_i = compute_si(i, job_prio_map);

				return std::min(s_i, compute_second_lst_bound(i));
			}

			// Upper bound on latest start time (LFT^bar - sum(C_j^max) - C_i^max)
			Time compute_second_lst_bound(const Job<Time> &i) {
				Job_set descendants = get_descendants(i);

				return latest_busy_time - i.maximal_cost() - std::accumulate(descendants.begin(), descendants.end(), 0,
																			 [](int x, auto &y) {
																				 return x + y->maximal_cost();
																			 });
			}

			// Gets all descendants of J_i in J^M
			Job_set get_descendants(const Job <Time> &i) {
				Job_set remaining_jobs{jobs};
				Job_set descendants{};

				std::deque<Job<Time>> queue{};
				queue.push_back(i);

				while (not queue.empty()) {
					Job<Time> &j = queue.front();
					queue.pop_front();
                    if (!job_precedence_sets.empty()) {
                        // if the job set has precedence constraints
                        size_t index_j = index_by_job.find(j.get_id())->second;

                        for (auto k: remaining_jobs) {
                            const Job_precedence_set &preds = job_precedence_sets[index_by_job.find(
                                    k->get_id())->second];
                            // k is a successor of j
                            if (std::find(preds.begin(), preds.end(), index_j) != preds.end()) {
                                descendants.push_back(k);
                                queue.push_back(*k);
                            }
                        }
                    }

					std::remove_if(remaining_jobs.begin(), remaining_jobs.end(), [descendants](auto &x) {
						return std::find(descendants.begin(), descendants.end(), x) != descendants.end();
					});
				}
				return descendants;
			}

			// Upper bound on latest start time (s_i)
			Time compute_si(const Job <Time> &i, const std::unordered_map<JobID, Priority> &job_prio_map) {
				const Job<Time>* blocking_job = nullptr;

				for (const Job<Time>* j : jobs_by_earliest_arrival) {
					if (j->get_id() == i.get_id()) {
						continue;
					}

					// use preprocessed prio level
					if (i.priority_exceeds(job_prio_map.find(j->get_id())->second) && (blocking_job == nullptr || blocking_job->maximal_cost() < j->maximal_cost())) {
						blocking_job = j;
					}
				}

				Time blocking_time = blocking_job == nullptr ? 0 : std::max<Time>(0, blocking_job->maximal_cost() - Time_model::constants<Time>::epsilon());
				Time latest_start_time = std::max(cpu_availability.max(), i.latest_arrival() + blocking_time);

				for (const Job<Time>* j : jobs_by_earliest_arrival) {
					if (j->get_id() == i.get_id()) {
						continue;
					}

					if (j->earliest_arrival() <= latest_start_time && !i.priority_exceeds(job_prio_map.find(j->get_id())->second)) {
						latest_start_time += j->maximal_cost();
					} else if (j->earliest_arrival() > latest_start_time) {
						break;
					}
				}

				return latest_start_time;
			}

			Time compute_latest_idle_time() {
				Time idle_time{-1};
				const Job<Time>* idle_job = nullptr;

				auto first_job_after_eft = std::find_if(jobs_by_latest_arrival.begin(), jobs_by_latest_arrival.end(),
														[this] (const Job<Time>* j) { return j->latest_arrival() > cpu_availability.min(); } );

				// No j with r_max > A1_min, so no idle time before the last job in job_set
				if (first_job_after_eft == jobs_by_latest_arrival.end()) {
					return idle_time;
				}

				for (const Job<Time>* i : jobs_by_latest_arrival) {
					// compute the earliest time at which the set of all jobs with r^max < r_i^max can complete
					Time t = cpu_availability.min();

					for (const Job<Time>* j : jobs_by_earliest_arrival) {
						if ((j->latest_arrival() < i->latest_arrival())) {
							t = std::max(t, j->earliest_arrival()) + j->least_cost();
						}

						if (t >= i->latest_arrival()) {
							break;
						}
					}

					if (t < i->latest_arrival()) {
						if (idle_job == nullptr || i->latest_arrival() > idle_job->latest_arrival()) {
							idle_job = i;
						}
					}
				}

				if (idle_job == nullptr) {
					return idle_time;
				}

				if (idle_job->latest_arrival() == (*jobs_by_latest_arrival.begin())->latest_arrival()) {
					return idle_time;
				} else {
					return idle_job->latest_arrival() - Time_model::constants<Time>::epsilon();
				}
			}

			Priority compute_max_priority() const {
				Priority max_prio {};

				for (const Job<Time>* j : jobs) {
					if (!j->priority_exceeds(max_prio)) {
						max_prio = j->get_priority();
					}
				}

				return max_prio;
			}

			// Returns the smallest wcet among the jobs with a lower priority than job
			Time min_lower_priority_wcet(const Job<Time> &job) const {
				auto pos = std::find_if(jobs_by_wcet.begin(), jobs_by_wcet.end(),
										[&job] (const Job<Time>* j) { return job.higher_priority_than(*j); } );

				if (pos == jobs_by_wcet.end()) {
					return 0;
				} else {
					return (*pos)->maximal_cost();
				}
			}

			void initialize_key() {
				for (const Job<Time>* j : jobs) {
					key = key ^ j->get_key();
				}
			}


		};

		template<typename T, typename Compare>
		typename std::vector<T>::iterator insert_sorted(std::vector<T> &vec, const T &item, Compare comp) {
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
}

#endif