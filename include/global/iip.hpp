#ifndef Global_IIP_HPP
#define Global_IIP_HPP

namespace NP {

	namespace Global {

		template<class Time> class Null_IIP
		{
			public:

			typedef Schedule_state<Time> State;
			typedef State_space<Time, Null_IIP> Space;
			typedef typename State_space<Time, Null_IIP>::Workload Jobs;

			typedef Job_set Scheduled;

			static const bool can_block = false;

			Null_IIP(const Space &space, const Jobs &jobs) {}

			Time latest_start(const Job<Time>& j, Time t, const State& s)
			{
				return Time_model::constants<Time>::infinity();
			}
		};

	}
}

#endif
