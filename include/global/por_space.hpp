#ifndef POR_GLOB_SCHEDULE_SPACE_H
#define POR_GLOB_SCHEDULE_SPACE_H

#include <unordered_map>
#include <map>
#include <vector>
#include <deque>
#include <list>
#include <algorithm>

#include <iostream>
#include <ostream>
#include <cassert>
#include <queue>

#include "config.h"
#include "problem.hpp"
#include "jobs.hpp"
#include "precedence.hpp"
#include "clock.hpp"

#include "global/state.hpp"
#include "global/space.hpp"
#include "global/por_criterion.hpp"
#include "global/reduction_set.hpp"

#include <chrono>

namespace NP
{

	namespace Global
	{

		template <class Time, class IIP = Null_IIP<Time>, class POR_criterion = POR_criterion<Time>>
		class Por_state_space : public State_space<Time, IIP>
		{

		public:
			typedef typename State_space<Time, IIP>::Problem Problem;
			typedef typename State_space<Time, IIP>::Workload Workload;
			// typedef typename State_space<Time, IIP>::Abort_actions Abort_actions;
			typedef typename State_space<Time, IIP>::State State;
			typedef typename State_space<Time, IIP>::Job_precedence_set Job_precedence_set;

			static Por_state_space explore(
				const Problem &prob,
				const Analysis_options &opts)
			{
				DM("!! MULTI PROCESSOR POR SAG !!" << std::endl);
				// this is a multiprocessor analysis
				assert(prob.num_processors > 1);

				// Preprocess the job such that they release at or after their predecessors release
				auto jobs = preprocess_jobs<Time>(prob.dag, prob.jobs);

				Por_state_space s = Por_state_space(jobs, prob.dag, prob.num_processors,
													opts.timeout, opts.max_depth,
													opts.num_buckets);
				s.group_add = opts.group_add;
				s.limit_fail = opts.limit_fail;
				s.cpu_time.start();
				if (opts.be_naive)
				{
					DM("\texploring naively" << std::endl);
					s.explore_naively();
				}
				else
					s.explore();
				s.cpu_time.stop();
				return s;
			}

			// convenience interface for tests
			static Por_state_space explore_naively(const Workload &jobs)
			{
				Problem p{jobs};
				Analysis_options o;
				o.be_naive = true;
				return explore(p, o);
			}

			// convenience interface for tests
			static Por_state_space explore(const Workload &jobs)
			{
				Problem p{jobs};
				Analysis_options o;
				return explore(p, o);
			}

			unsigned long number_of_por_successes() const
			{
				return reduction_successes;
			}

			unsigned long number_of_por_failures() const
			{
				return reduction_failures;
			}

			unsigned long number_of_jobs_in_por() const
			{
				unsigned long jobs_in_por = 0;
				unsigned long largest_fail = 0;
				unsigned long total_fail = 0;
				for (Reduction_set_statistics<Time> rss : reduction_set_statistics)
				{
					if (rss.reduction_success)
					{
						jobs_in_por += rss.num_jobs;
					}else{
						largest_fail = (rss.num_jobs > largest_fail) ? rss.num_jobs : largest_fail;
						total_fail +=  rss.num_jobs;
					}
				}
				//std::cout << "largest fail set was : " << largest_fail << " AVG fail size is : " << total_fail/reduction_failures <<std::endl;
				return jobs_in_por;
			}

			std::vector<Reduction_set_statistics<Time>> get_reduction_set_statistics() const
			{
				return reduction_set_statistics;
			}

		protected:
			using State_space<Time, IIP>::explore;
			using State_space<Time, IIP>::explore_naively;
			// using State_space<Time, IIP>::dispatch;

			POR_criterion por_criterion;
			unsigned long reduction_successes, reduction_failures, jobs_in_por;
			std::vector<Reduction_set_statistics<Time>> reduction_set_statistics;
			std::vector<Job_precedence_set> job_precedence_sets;

			// Normal state space: jobs, dag_edges, num_cpus, max_cpu_time, max_depth, num_buckets
			Por_state_space(const Workload &jobs,
							const Precedence_constraints &dag_edges,
							unsigned int num_cpus,
							double max_cpu_time = 0,
							unsigned int max_depth = 0,
							std::size_t num_buckets = 1000,
							bool early_exit = true)
				: State_space<Time, IIP>(jobs, dag_edges, num_cpus, max_cpu_time, max_depth, num_buckets), por_criterion(), reduction_successes(0), reduction_failures(0), reduction_set_statistics(), job_precedence_sets(jobs.size())
			{
				for (auto e : dag_edges)
				{
					const Job<Time> &from = lookup<Time>(jobs, e.first);
					const Job<Time> &to = lookup<Time>(jobs, e.second);
					job_precedence_sets[index_of(to, jobs)].push_back(index_of(from, jobs));
				}
			}

			// entrance: this creates the reduction set
			Reduction_set<Time> create_reduction_set(const State &s, typename Reduction_set<Time>::Job_set &eligible_successors)
			{
				std::vector<std::size_t> indices{};
				for (const Job<Time> *j : eligible_successors)
				{
					indices.push_back(this->index_of(*j));
				}
				// okie nu heb ik een basis reductie set
				//========================Creating reduction set===========================//
				Reduction_set<Time> reduction_set = Reduction_set<Time>(s.get_all_core_availabilities(), eligible_successors, indices);

				const Job<Time> *j;
				Time lb = s.core_availability().min();
				while (true)
				{
					if (reduction_set.has_potential_deadline_misses())
					{
						reduction_failures++;
						reduction_set_statistics.push_back(Reduction_set_statistics<Time>{false, reduction_set});

						return reduction_set;
					}
					//========================interfering jobs=================================//
					// an interfering job can only interfere if it can start before the latest latest start time of all jobs in the current candidate set.
					// so lets iterate over the list of jobs sorted by the earliest arrival time from Amin to LST_max

					// i aint using the macros as they require a lot of additions to the states etc
					std::vector<const Job<Time> *> interfering_jobs{};
					bool updated = false;
					// adding the lower bound here does work, but it was not really that much of a time save
					for (auto it = this->jobs_by_earliest_arrival.lower_bound(lb); it != this->jobs_by_earliest_arrival.upper_bound(reduction_set.get_latest_LST()); it++)
					{
						j = it->second;
						const Job<Time> &j_i = *j;
						const Job_precedence_set &preds = this->job_precedence_sets[this->index_of(j_i)];
						if (reduction_set.can_interfere(*j) && s.job_incomplete(this->index_of(j_i)))
						{
							if (!updated)
							{
								lb = j->earliest_arrival();
								updated = true;
							}
							interfering_jobs.push_back(j);
						}
					}

					// Now we have a (possible) set of interfering jobs
					if (!interfering_jobs.empty())
					{
						//Group add is optimization 2 where we add one or all interfering jobs
						//!ALWAYS RUN reduction_set.update_set() AFTER ADDING JOBS!
						if (!this->group_add)
						{
							const Job<Time> *jx = *std::min_element(interfering_jobs.begin(), interfering_jobs.end(),
																	[](const Job<Time> *i, const Job<Time> *j) -> bool
																	{
																		return i->higher_priority_than(*j);
																	});
							reduction_set.add_job(jx, this->index_of(*jx));
							reduction_set.update_set();
						}
						else
						{
							for (const Job<Time> *j_i : interfering_jobs)
							{
								reduction_set.add_job(j_i, this->index_of(*j_i));
							}
							//So only call update_set after adding ALL jobs otherwise there is no time save
							reduction_set.update_set();
						}
					}
					else
					{
						reduction_set_statistics.push_back(Reduction_set_statistics<Time>{true, reduction_set});
						reduction_successes++;
						jobs_in_por += reduction_set.get_number_of_jobs();
						// no more interfering jobs so we can return the reduction set.
						// std::cout << " interference took in total " << interfering_total << " ns " << std::endl;
						// reduction_set.get_timings();
						return reduction_set;
					}
				}
			}

			void process_new_edge(
				const State &from,
				const State &to,
				const Reduction_set<Time> &reduction_set,
				const Interval<Time> &finish_range)
			{
				// update response times
				for (const Job<Time> *j : reduction_set.get_jobs())
				{
					this->update_finish_times(*j, Interval<Time>{reduction_set.earliest_finish_time(*j), reduction_set.latest_finish_time(*j)});
				}
				// update statistics
				this->num_edges++;
#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
				this->edges.push_back(std::make_unique<Reduced_edge>(reduction_set, &from, &to, finish_range));
#endif
			}

			void dispatch_reduction_set_merge(const State &current_state, Reduction_set<Time> &reduction_set)
			{
				std::vector<std::size_t> indices{};
				for (const Job<Time> *j : reduction_set.get_jobs())
				{
					indices.push_back(this->index_of(*j));
				}

				std::vector<Time> PA = reduction_set.compute_possibly_available();
				std::vector<Time> CA = reduction_set.compute_certainly_available();

				const State &next = this->new_or_merged_state(current_state, indices, PA, CA, reduction_set.get_key());

				process_new_edge(current_state, next, reduction_set, {0, 0});
			}

			void dispatch_reduction_set_naive(const State &current_state, Reduction_set<Time> &reduction_set)
			{

				// Interval<Time> finish_range = next_finish_times(reduction_set);

				std::vector<std::size_t> indices{};
				for (const Job<Time> *j : reduction_set.get_jobs())
				{
					indices.push_back(this->index_of(*j));
				}
				std::vector<Time> PA = reduction_set.compute_possibly_available();
				std::vector<Time> CA = reduction_set.compute_certainly_available();

				const State &next = this->new_state(current_state, indices, PA, CA, reduction_set.get_key());
				//	DM("      -----> S" << (states.end() - states.begin()) << std::endl);
				// At this point its an edge with 0,0 but there needs to be a new
				process_new_edge(current_state, next, reduction_set, {0, 0});
			}


			/*
				Here comes optimization 1 actually into play as we want to give the failed set to the next state as well
			*/
			bool dispatch(const State &s, const Job<Time> &j, Time t_wc, std::vector<std::size_t> failed_set, int retry_depth)
			{
				// check if this job has a feasible start-time interval
				auto _st = this->start_times(s, j, t_wc);
				if (_st.first > _st.second)
					return false; // nope

				Interval<Time> st{_st};

				// yep, job j is a feasible successor in state s

				// compute range of possible finish times
				Interval<Time> ftimes = st + j.get_cost();

				// update finish-time estimates
				this->update_finish_times(j, ftimes);
				if (this->aborted)
				{
					// std::cout<<"fail trace is " << s << std::endl;
				}
				//  expand the graph, merging if possible
				//  met be_naive wordt bedoelt dat als ie false is dat ie niet gaat mergen
				//  dus in de toekomst
				const State &next = this->be_naive ? this->new_state(s, this->index_of(j), this->predecessors_of(j),
																	 st, ftimes, j.get_key(), failed_set,retry_depth)
												   : this->new_or_merged_state(s, this->index_of(j), this->predecessors_of(j),
																			   st, ftimes, j.get_key(), failed_set,retry_depth);

				// make sure we didn't skip any jobs
				this->check_for_deadline_misses(s, next);

#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH
				this->edges.emplace_back(&j, &s, &next, ftimes);
#endif
				this->count_edge();
				this->current_job_count++;

				return true;
			}

			// Here we explore our current state where we also try to make our reduction set, this is why this is placed in this file
			void explore(const State &s) override
			{
				bool found_one = false;

				DM("----" << std::endl);

				// (0) define the window of interest

				// earliest time a core is possibly available
				auto t_min = s.core_availability().min();
				// latest time some unfinished job is certainly ready
				auto t_job = this->next_job_ready(s, t_min);
				// latest time some core is certainly available
				auto t_core = s.core_availability().max();
				// latest time by which a work-conserving scheduler
				// certainly schedules some job
				auto t_wc = std::max(t_core, t_job);

				DM(s << std::endl
					 << "t_min: " << t_min << std::endl
					 << "t_job: " << t_job << std::endl
					 << "t_core: " << t_core << std::endl
					 << "t_wc: " << t_wc << std::endl);

				
				// We have to find all elligible succesors here first.

				typename Reduction_set<Time>::Job_set eligible_successors{};

				/*
				eligible successors are:
				1. all jobs that can start before t_min, aka. A_min
				2. all jobs that can start before the first job in the first set is CERTAINLY released
				*/
				// so first all jobs that might be pending before A_min

				// (1) first check jobs that may be already pending
				DM("==== [1] ====" << std::endl);
				bool skip_set = true;
				for (const Job<Time> &j : this->jobs_by_win.lookup(t_min))
				{
					if (j.earliest_arrival() <= t_min && this->ready(s, j))
					{
						eligible_successors.push_back(&j);
						//This is optimization 1: where we look if there exists a job which was not in the previous failed set
						if (!s.job_in_failed_set(this->index_of(j)))
						{
							skip_set = false;
						}
					}
				}

				DM("==== [2] ====" << std::endl);
				// (2) check jobs that are released only later in the interval
				for (auto it = this->jobs_by_earliest_arrival.upper_bound(t_min);
					 it != this->jobs_by_earliest_arrival.end();
					 it++)
				{
					const Job<Time> &j = *it->second;
					DM(j << " (" << index_of(j) << ")" << std::endl);
					// stop looking once we've left the window of interest
					if (j.earliest_arrival() > t_wc)
						break;

					// Job could be not ready due to precedence constraints
					if (!this->ready(s, j))
						continue;

					// Since this job is released in the future, it better
					// be incomplete...
					assert(this->unfinished(s, j));

					eligible_successors.push_back(it->second);
					if (this->limit_fail)
					{
						if (!s.job_in_failed_set(this->index_of(j)))
						{
							skip_set = false;
						}
					}
				}
				// now we can finially create a proper reduction set

				// look if all eligible successors are present in the previously failed set
				// if so, lets not make a reduction

				// if the size is 1 then just do a normal schedule

				std::vector<std::size_t> failed_reduction;
				// if we skip making this set, then keep looking at the previous set in the next state
				bool limit_failures = this->limit_fail;
				//Optimization 1:
				if (skip_set && limit_failures)
				{
					failed_reduction = s.get_previous_failed_set();
				}
				if (!limit_failures || s.get_retry_reduction() == 0)
				{
					skip_set = false;
				}
				//retry reduction is WIP such that we might dynamically state when we want to start new reductions again
				int retry_depth = s.get_retry_reduction() - 1;
				if (eligible_successors.size() > 1 && !skip_set)
				{
					Reduction_set<Time> reduction_set = create_reduction_set(s, eligible_successors);

					if (!reduction_set.has_potential_deadline_misses())
					{
						DM("\n---\nPartial-order reduction is safe" << std::endl);
						// uncomment to print the CA and PA values
						// reduction_set.created_set();
						//  now we must create something to properly schedule the set
						// reduction_set.show_time_waste();
						if (this->be_naive)
						{
							dispatch_reduction_set_naive(s, reduction_set);
						}
						else
						{
							dispatch_reduction_set_merge(s, reduction_set);
						}
						// this->current_job_count += reduction_set.get_jobs().size();
						found_one = true;
						// if there were no deadline misses and we were able to dispatch it normally, then we can return here.
						return;
					}
					else
					{
						retry_depth = reduction_set.get_jobs().size()/1;
						DM("\tPartial order reduction is not safe" << std::endl);
						//Copy all jobs from failed set into a smaller list
						if (limit_failures)
						{
							for (const Job<Time> *j : reduction_set.get_jobs())
							{
								failed_reduction.push_back(this->index_of(*j));
							}
						}
					}
				}

				DM("\n---\nPartial-order reduction is not safe, or just one job" << std::endl);
				for (const Job<Time> *j : eligible_successors)
				{
					// we can use the normal dispatch now so no worries
					// maybe we can give some information to the next state if we have a failed reduction, in such a way that it wont try a next reduction
					//  if that next job was in the failed set.
					found_one |= dispatch(s, *j, t_wc, failed_reduction, retry_depth);
				}

				// check for a dead end
				if (!found_one && !this->all_jobs_scheduled(s))
				{
					// out of options and we didn't schedule all jobs
					DM("POR dead end abortion" << std::endl);
					this->aborted = true;
				}
			}

			/*
			Time earliest_possible_job_release(const State &s, const Reduction_set<Time> &ignored_set)
			{
				DM("      - looking for earliest possible job release starting from: "
				   << s.earliest_job_release() << std::endl);
				const Job<Time> *jp;
				foreach_possibly_pending_job(s, jp)
				{
					const Job<Time> &j = *jp;

					DM("         * looking at " << j << std::endl);

					// skip if it is the one we're ignoring
					bool ignored = false;
					for (const Job<Time> *ignored_job : ignored_set.get_jobs())
					{
						if (&j == ignored_job)
						{
							ignored = true;
							break;
						}
					}

					if (!ignored)
					{
						DM("         * found it: " << j.earliest_arrival() << std::endl);
						// it's incomplete and not ignored => found the earliest
						return j.earliest_arrival();
					}
				}

				DM("         * No more future releases" << std::endl);
				return Time_model::constants<Time>::infinity();
			}*/

			Interval<Time> next_finish_times(const Reduction_set<Time> &reduction_set)
			{

				// standard case -- this job is never aborted or skipped
				return Interval<Time>{
					reduction_set.earliest_finish_time(),
					reduction_set.get_latest_busy_time()};
			}

#ifdef CONFIG_COLLECT_SCHEDULE_GRAPH

			struct Reduced_edge : State_space<Time, IIP>::Edge
			{
				Reduction_set<Time> reduction_set;

				Reduced_edge(const Reduction_set<Time> &reduction_set, const State *src, const State *tgt,
							 const Interval<Time> &fr)
					: State_space<Time, IIP>::Edge(reduction_set.get_jobs()[0], src, tgt, fr), reduction_set{reduction_set}
				{
				}

				bool deadline_miss_possible() const override
				{
					return reduction_set.has_potential_deadline_misses();
				}

				Time earliest_finish_time() const override
				{
					return reduction_set.earliest_finish_time();
				}

				Time latest_finish_time() const override
				{
					return reduction_set.get_latest_busy_time();
				}

				Time earliest_start_time() const override
				{
					return reduction_set.earliest_start_time();
				}

				Time latest_start_time() const override
				{
					return reduction_set.latest_start_time();
				}
			};

			void print_edge(std::ostream &out, const std::unique_ptr<typename State_space<Time, IIP>::Edge> &e, unsigned int source_id, unsigned int target_id) const override
			{
				out << "\tS" << source_id
					<< " -> "
					<< "S" << target_id
					<< "[label=\"";

				auto r = dynamic_cast<Reduced_edge *>(e.get());
				if (r)
				{
					for (auto j : r->reduction_set.get_jobs())
					{
						out << "T" << j->get_task_id()
							<< " J" << j->get_job_id()
							<< "\\nDL=" << j->get_deadline()
							<< "\\n";
					}
				}
				else
				{
					out << "T" << e->scheduled->get_task_id()
						<< " J" << e->scheduled->get_job_id()
						<< "\\nDL=" << e->scheduled->get_deadline();
				}

				out << "\\nES=" << e->earliest_start_time()
					<< "\\nLS=" << e->latest_start_time()
					<< "\\nEF=" << e->earliest_finish_time()
					<< "\\nLF=" << e->latest_finish_time()
					<< "\"";
				if (e->deadline_miss_possible())
				{
					out << ",color=Red,fontcolor=Red";
				}
				out << ",fontsize=8"
					<< "]"
					<< ";"
					<< std::endl;
				if (e->deadline_miss_possible())
				{
					out << "S" << target_id
						<< "[color=Red];"
						<< std::endl;
				}
			}
#endif
		};

	}
}

#endif
