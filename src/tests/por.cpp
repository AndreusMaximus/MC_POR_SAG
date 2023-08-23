#include "doctest.h"

#include "uni/por_space.hpp"
#include "uni/por_criterion.hpp"
#include "uni/reduction_set.hpp"

using namespace NP::Uniproc;
using namespace NP;

typedef std::vector<std::size_t> Job_precedence_set;

TEST_CASE("[Partial-order reduction] Example") {
    Job<dtime_t> j1{1, I(6, 8), I(2, 2), 100, 1};
    Job<dtime_t> j2{2, I(7, 10), I(2, 2), 100, 2};
    Job<dtime_t> j3{3, I(9, 11), I(2, 2), 100, 3};
    Job<dtime_t> j4{4, I(5, 9), I(2, 2), 100, 4};
    std::vector<const Job<dtime_t>*> jobs{&j1, &j2, &j3, &j4};
    std::vector<const Job<dtime_t>*> eligible_successors{&j1, &j2, &j4};


	Reduction_set<dtime_t> reduction_set{I(7, 13), eligible_successors, {0, 1, 3}};

    SUBCASE("Latest busy time") {
        CHECK(reduction_set.get_latest_busy_time() == 19);
    }

    SUBCASE("Latest idle time") {
        CHECK(reduction_set.get_latest_idle_time() == -1);
    }

    SUBCASE("Latest start time") {
        CHECK(reduction_set.get_latest_start_time(j1) == 13);
        CHECK(reduction_set.get_latest_start_time(j2) == 15);
        CHECK(reduction_set.get_latest_start_time(j4) == 17);
    }

    SUBCASE("No potential deadline miss") {
        CHECK(!reduction_set.has_potential_deadline_misses());
    }

    SUBCASE("Add job to reduction set") {
        reduction_set.add_job(&j3, 2);

        CHECK(reduction_set.get_latest_busy_time() == 21);
        CHECK(reduction_set.get_latest_idle_time() == -1);
        CHECK(reduction_set.get_latest_start_time(j1) == 13);
        CHECK(reduction_set.get_latest_start_time(j2) == 15);
        CHECK(reduction_set.get_latest_start_time(j3) == 17);
        CHECK(reduction_set.get_latest_start_time(j4) == 19);
        CHECK(!reduction_set.has_potential_deadline_misses());
    }

}

TEST_CASE("[Partial-order reduction] Example 2") {
    Job<dtime_t> j1{1, I(6, 8), I(2, 2), 100, 1};
    Job<dtime_t> j3{3, I(7, 10), I(2, 2), 100, 3};
    Job<dtime_t> j4{4, I(9, 18), I(2, 2), 100, 4};
    Job<dtime_t> j5{5, I(5, 21), I(2, 2), 22, 5};
    Job<dtime_t> j6{6, I(19, 20), I(5, 5), 100, 6};
    std::vector<const Job<dtime_t>*> jobs{&j1, &j3, &j4, &j5, &j6};
    std::vector<const Job<dtime_t>*> eligible_successors{&j1, &j3, &j5};

    Reduction_set<dtime_t> reduction_set{I(7, 13), eligible_successors, {0, 1, 2}};

    SUBCASE("Latest busy time") {
        CHECK(reduction_set.get_latest_busy_time() == 23);
    }

    SUBCASE("Latest idle time") {
        CHECK(reduction_set.get_latest_idle_time() == 20);
    }

    SUBCASE("Latest start time") {
        CHECK(reduction_set.get_latest_start_time(j1) == 13);
        CHECK(reduction_set.get_latest_start_time(j3) == 15);
        CHECK(reduction_set.get_latest_start_time(j5) == 21);
    }

    SUBCASE("Potential deadline miss") {
        CHECK(reduction_set.has_potential_deadline_misses());
    }

    SUBCASE("Jx cannot interfere because r_x > te - min_wcet") {
        Job<dtime_t> jx{2, I(22, 23), I(2, 2), 100, 2};
        CHECK(!reduction_set.can_interfere(jx, {}, {}));
    }

    SUBCASE("Jx cannot interfere because J_x has a lower priority and r_x >= max(r_max)") {
        Job<dtime_t> jx{6, I(22, 23), I(2, 2), 100, 6};
        CHECK(!reduction_set.can_interfere(jx, {}, {}));
    }

    SUBCASE("Jx can interfere because there exists a J_i such that r_x <= s_i and p_x < p_i") {
        Job<dtime_t> jx{2, I(14, 23), I(2, 2), 100, 2};
        CHECK(reduction_set.can_interfere(jx, {}, {}));
    }

    SUBCASE("Jx can interfere because r_x <= l_M") {
        Job<dtime_t> jx{6, I(19, 23), I(2, 2), 100, 6};
        CHECK(reduction_set.can_interfere(jx, {}, {}));
    }

    SUBCASE("Adding job to reduction set doesn't change latest busy time") {
        reduction_set.add_job(&j4, 2);

        CHECK(reduction_set.get_latest_busy_time() == 23);
        CHECK(reduction_set.get_latest_idle_time() == 20);
        CHECK(reduction_set.get_latest_start_time(j1) == 13);
        CHECK(reduction_set.get_latest_start_time(j3) == 15);
        CHECK(reduction_set.get_latest_start_time(j4) == 21);
        CHECK(reduction_set.get_latest_start_time(j5) == 21);
        CHECK(reduction_set.has_potential_deadline_misses());
    }

    SUBCASE("Adding lower priority job to reduction set changes all start times after the idle time") {
        reduction_set.add_job(&j6, 4);

        CHECK(reduction_set.get_latest_busy_time() == 27);
        CHECK(reduction_set.get_latest_idle_time() == 19);
        CHECK(reduction_set.get_latest_start_time(j1) == 13);
        CHECK(reduction_set.get_latest_start_time(j3) == 16);
        CHECK(reduction_set.get_latest_start_time(j5) == 25);
        CHECK(reduction_set.get_latest_start_time(j6) == 22);
        CHECK(reduction_set.has_potential_deadline_misses());
    }


}

TEST_CASE("[Partial-order reduction] Example 3") {
    Job<dtime_t> j1{1, I(6, 8), I(2, 2), 100, 1};
    Job<dtime_t> j2{2, I(7, 11), I(2, 2), 100, 2};
    Job<dtime_t> j4{4, I(5, 10), I(2, 2), 100, 4};
    Job<dtime_t> j5{5, I(9, 11), I(7, 7), 100, 5};
    std::vector<const Job<dtime_t>*> jobs{&j1, &j2, &j4, &j5};
    std::vector<const Job<dtime_t>*> eligible_successors{&j1, &j2, &j4};

    Reduction_set<dtime_t> reduction_set{I(7, 13), eligible_successors, {0, 1, 2}};

    SUBCASE("Latest busy time") {
        CHECK(reduction_set.get_latest_busy_time() == 19);
    }

    SUBCASE("Latest idle time") {
        CHECK(reduction_set.get_latest_idle_time() == 9);
    }

    SUBCASE("Latest start times") {
        CHECK(reduction_set.get_latest_start_time(j1) == 13);
        CHECK(reduction_set.get_latest_start_time(j2) == 15);
        CHECK(reduction_set.get_latest_start_time(j4) == 17);
    }

    SUBCASE("Adding lower priority job to reduction set changes all start times after the idle time") {
        reduction_set.add_job(&j5, 3);

        CHECK(reduction_set.get_latest_start_time(j1) == 14);
        CHECK(reduction_set.get_latest_start_time(j2) == 19);
        CHECK(reduction_set.get_latest_start_time(j4) == 20);
        CHECK(reduction_set.get_latest_start_time(j5) == 19);
    }
}

TEST_CASE("[Partial-order reduction] Example 4") {
    Job<dtime_t> j1{1, I(6, 8), I(2, 2), 100, 1};
    Job<dtime_t> j2{2, I(7, 11), I(2, 2), 100, 2};
    Job<dtime_t> j4{4, I(5, 8), I(2, 2), 100, 4};
    std::vector<const Job<dtime_t>*> eligible_successors{&j1, &j2, &j4};

    Reduction_set<dtime_t> reduction_set{I(7, 13), eligible_successors, {0, 1, 2}};

    SUBCASE("Latest busy time") {
        CHECK(reduction_set.get_latest_busy_time() == 19);
    }

    SUBCASE("Latest idle time") {
        CHECK(reduction_set.get_latest_idle_time() == -1);
    }
}

TEST_CASE("[Partial-order reduction] Example 5") {
    Job<dtime_t> j1{1, I(6, 8), I(2, 2), 100, 1};
    Job<dtime_t> j2{2, I(7, 12), I(2, 2), 100, 2};
    Job<dtime_t> j4{4, I(5, 10), I(2, 2), 100, 4};
    std::vector<const Job<dtime_t>*> eligible_successors{&j1, &j2, &j4};

    Reduction_set<dtime_t> reduction_set{I(7, 13), eligible_successors, {0, 1, 2}};

    SUBCASE("Latest busy time") {
        CHECK(reduction_set.get_latest_busy_time() == 19);
    }

    SUBCASE("Latest idle time") {
        CHECK(reduction_set.get_latest_idle_time() == 11);
    }
}

TEST_CASE("[Partial-order reduction] Example 6") {
    Job<dtime_t> j1{1, I(6, 7), I(2, 2), 100, 1};
    Job<dtime_t> j2{2, I(7, 7), I(2, 2), 100, 2};
    Job<dtime_t> j4{4, I(5, 7), I(2, 2), 100, 4};
    std::vector<const Job<dtime_t>*> eligible_successors{&j1, &j2, &j4};

    Reduction_set<dtime_t> reduction_set{I(7, 13), eligible_successors, {0, 1, 2}};

    SUBCASE("Latest busy time") {
        CHECK(reduction_set.get_latest_busy_time() == 19);
    }

    SUBCASE("No latest idle time") {
        CHECK(reduction_set.get_latest_idle_time() == -1);
    }
}

TEST_CASE("[Partial-order reduction] Example 7") {
    Job<dtime_t> j0{0, I(9, 11), I(1, 1), 21, 0};
    Job<dtime_t> j1{1, I(6, 8), I(2, 2), 20, 1};
    Job<dtime_t> j2{2, I(7, 12), I(2, 2), 26, 2};
    Job<dtime_t> j4{4, I(5, 10), I(10, 10), 28, 4};
    std::vector<const Job<dtime_t>*> eligible_successors{&j0, &j1, &j2, &j4};

    Reduction_set<dtime_t> reduction_set{I(7, 13), eligible_successors, {0, 1, 2, 3}};

    SUBCASE("Latest start time") {
        CHECK(reduction_set.get_latest_start_time(j0) == 20);
        CHECK(reduction_set.get_latest_start_time(j1) == 18);
        CHECK(reduction_set.get_latest_start_time(j2) == 24);
        CHECK(reduction_set.get_latest_start_time(j4) == 18);
    }

    SUBCASE("No potential deadline miss") {
        CHECK(!reduction_set.has_potential_deadline_misses());
    }

}

TEST_CASE("[Partial-order reduction] Example 8") {
    Job<dtime_t> j0{0, I(9, 9), I(1, 1), 100, 0};
    Job<dtime_t> j1{1, I(6, 8), I(1, 1), 100, 1};
    Job<dtime_t> j2{2, I(7, 12), I(2, 2), 24, 2};
    Job<dtime_t> j4{4, I(5, 10), I(10, 10), 100, 4};
    std::vector<const Job<dtime_t>*> eligible_successors{&j0, &j1, &j2, &j4};

    Reduction_set<dtime_t> reduction_set{I(7, 13), eligible_successors, {0, 1, 2, 3}};

    SUBCASE("Latest start time") {
        CHECK(reduction_set.get_latest_start_time(j0) == 18);
        CHECK(reduction_set.get_latest_start_time(j1) == 18);
        CHECK(reduction_set.get_latest_start_time(j2) == 23);
        CHECK(reduction_set.get_latest_start_time(j4) == 17);
    }

    SUBCASE("Potential deadline miss") {
        CHECK(reduction_set.has_potential_deadline_misses());
    }

}

TEST_CASE("[Partial-order reduction] Example 9") {
    Job<dtime_t> j0{0, I(9, 11), I(2, 2), 21, 0};
    Job<dtime_t> j1{1, I(6, 8), I(2, 2), 21, 1};
    Job<dtime_t> j2{2, I(7, 12), I(2, 2), 27, 2};
    Job<dtime_t> j4{4, I(5, 10), I(10, 10), 29, 4};
    std::vector<const Job<dtime_t>*> eligible_successors{&j0, &j1, &j2, &j4};

    Reduction_set<dtime_t> reduction_set{I(7, 13), eligible_successors, {0, 1, 2, 3}};

    SUBCASE("Latest start time") {
        CHECK(reduction_set.get_latest_start_time(j0) == 20);
        CHECK(reduction_set.get_latest_start_time(j1) == 19);
        CHECK(reduction_set.get_latest_start_time(j2) == 25);
        CHECK(reduction_set.get_latest_start_time(j4) == 19);
    }

    SUBCASE("Potential deadline miss") {
        CHECK(reduction_set.has_potential_deadline_misses());
    }

}

TEST_CASE("[Partial-order reduction] Example 10") {
    Job<dtime_t> j0{0, I(9, 10), I(2, 2), 100, 0};
    Job<dtime_t> j1{1, I(6, 8), I(2, 2), 100, 1};
    Job<dtime_t> j2{2, I(7, 13), I(2, 2), 100, 2};
    Job<dtime_t> j4{4, I(5, 10), I(10, 10), 100, 4};
    std::vector<const Job<dtime_t>*> eligible_successors{&j0, &j1, &j2, &j4};

    Reduction_set<dtime_t> reduction_set{I(7, 13), eligible_successors, {0, 1, 2, 3}};

    SUBCASE("Latest start time") {
        CHECK(reduction_set.get_latest_start_time(j0) == 19);
        CHECK(reduction_set.get_latest_start_time(j1) == 19);
        CHECK(reduction_set.get_latest_start_time(j2) == 26);
        CHECK(reduction_set.get_latest_start_time(j4) == 19);
    }
}

TEST_CASE("[Partial-order reduction] Example 11") {
    Job<dtime_t> j1{1, I(6, 10), I(3, 3), 100, 1};
    Job<dtime_t> j2{2, I(7, 13), I(2, 2), 100, 2};
    Job<dtime_t> j4{4, I(5, 14), I(10, 10), 100, 4};
    std::vector<const Job<dtime_t>*> eligible_successors{&j1, &j2, &j4};

    Reduction_set<dtime_t> reduction_set{I(7, 9), eligible_successors, {0, 1, 2}};

    SUBCASE("Latest busy time") {
        CHECK(reduction_set.get_latest_busy_time() == 25);
    }

    SUBCASE("Latest start time") {
        CHECK(reduction_set.get_latest_start_time(j1) == 19);
        CHECK(reduction_set.get_latest_start_time(j2) == 23);
        CHECK(reduction_set.get_latest_start_time(j4) == 15);
    }
}

TEST_CASE("[Partial-order reduction] Example 12") {
    Job<dtime_t> j0{0, I(8, 10), I(2, 2), 100, 0};
    Job<dtime_t> j1{1, I(7, 11), I(2, 2), 100, 1};
    Job<dtime_t> j2{2, I(9, 12), I(2, 2), 100, 2};
    std::vector<const Job<dtime_t>* > eligible_successors{&j0, &j1, &j2};

    Reduction_set<dtime_t> reduction_set{I(7, 13), eligible_successors, {0, 1, 2}};

    SUBCASE("Latest idle time when r_min has different order than r_max") {
        CHECK(reduction_set.get_latest_idle_time() == 11);
    }
}

TEST_CASE("[Partial-order reduction] Latest idle time") {
    Job<dtime_t> t1j125{1, I(124000, 124034), I(158, 198), 125000, 1, 125};
    Job<dtime_t> t1j126{1, I(125000, 125034), I(158, 198), 126000, 1, 126};
    Job<dtime_t> t2j63{2, I(124000, 124188), I(4, 6), 126000, 2, 63};
    Job<dtime_t> t3j26{3, I(125000, 125144), I(394, 493), 130000, 3, 26};
    Job<dtime_t> t7j7{7, I(120000, 120490), I(529, 662), 140000, 7, 7};
    std::vector<const Job<dtime_t>* > jobs{&t1j125, &t1j126, &t2j63, &t3j26, &t7j7};
    std::vector<const Job<dtime_t>* > eligible_successors{&t1j125, &t1j126, &t2j63, &t7j7};

    Reduction_set<dtime_t> reduction_set{I(123389, 124025), eligible_successors, {0, 1, 2, 4}};

    SUBCASE("Latest idle time when r_min has different order than r_max") {
        CHECK(reduction_set.get_latest_idle_time() == 125033);
    }

    SUBCASE("T3J26 can interfere because r_min before end of idle interval") {
        CHECK(reduction_set.can_interfere(t3j26, {}, {}));
    }
}

TEST_CASE("[Partial-order reduction] State space for partial-order reduction") {
    Job<dtime_t> j0{0, I(0, 0), I(7, 13), 100, 0};
    Job<dtime_t> j1{1, I(6, 8), I(2, 2), 100, 1};
    Job<dtime_t> j2{2, I(7, 9), I(2, 2), 100, 2};
    Job<dtime_t> j4{4, I(5, 10), I(2, 2), 100, 4};
    Job<dtime_t> j5{5, I(18, 20), I(1, 3), 100, 5};
    Job<dtime_t> j3{3, I(18, 29), I(4, 4), 100, 3};
    Job<dtime_t> j6{6, I(28, 30), I(3, 4), 100, 6};
    std::vector<Job<dtime_t>> jobs{j0, j1, j2, j3, j4, j5, j6};

    Analysis_options options{};

    auto space = Por_state_space<dtime_t, Null_IIP<dtime_t>, POR_priority_order<dtime_t>>::explore(jobs, options);
    CHECK(1==1);
}

TEST_CASE("[Partial-order reduction] State space for partial-order reduction 2") {
    Job<dtime_t> j0{0, I(0, 0), I(7, 13), 100, 0};
    Job<dtime_t> j1{1, I(6, 8), I(4, 4), 100, 1};
    Job<dtime_t> j2{2, I(7, 10), I(4, 4), 100, 2};
    Job<dtime_t> j3{3, I(5, 21), I(2, 2), 100, 3};
    Job<dtime_t> j4{4, I(14, 16), I(6, 6), 100, 4};
    Job<dtime_t> j5{5, I(14, 16), I(1, 1), 100, 5};
    std::vector<Job<dtime_t>> jobs{j0, j1, j2, j3, j4, j5};

    Analysis_options options{};

    auto space = Por_state_space<dtime_t, Null_IIP<dtime_t>, POR_priority_order<dtime_t>>::explore(jobs, options);
    // j4 should be added before j5
    CHECK(1==1);
}

TEST_CASE("[Partial-order reduction] Keys") {
    Job<dtime_t> j0{0, I(0, 0), I(7, 13), 100, 0};
    Job<dtime_t> j1{1, I(6, 8), I(2, 2), 100, 1};
    Job<dtime_t> j2{2, I(7, 13), I(2, 2), 100, 2};
    Job<dtime_t> j4{4, I(5, 10), I(10, 10), 100, 4};

    std::vector<const Job<dtime_t>*> jobs{&j1, &j2, &j4};
    Reduction_set<dtime_t> reduction_set{I(7, 13), jobs, {1, 2, 3}};

    SUBCASE("Individual jobs") {
        Schedule_state<dtime_t> s0{};
        Schedule_state<dtime_t> s1{s0, j0, 0, I(7, 13), 5};
        Schedule_state<dtime_t> s2{s1, j1, 1, I(9, 15), 5};
        Schedule_state<dtime_t> s3{s2, j2, 2, I(11, 17), 5};
        Schedule_state<dtime_t> s4{s3, j4, 3, I(21, 27), Time_model::constants<dtime_t>::infinity()};

        CHECK(s1.get_key() == j0.get_key());
        CHECK(s2.get_key() == (j0.get_key() ^ j1.get_key()));
        CHECK(s3.get_key() == (j0.get_key() ^ j1.get_key() ^ j2.get_key()));
        CHECK(s4.get_key() == (j0.get_key() ^ j1.get_key() ^ j2.get_key() ^ j4.get_key()));
    }

    SUBCASE("Reduction set") {
        Schedule_state<dtime_t> s0{};
        Schedule_state<dtime_t> s1{s0, j0, 0, I(7, 13), 5};
        Schedule_state<dtime_t> s2{s1, reduction_set, std::vector<std::size_t>{1, 2, 3}, I(21, 27), Time_model::constants<dtime_t>::infinity()};

        CHECK(s1.get_key() == j0.get_key());
        CHECK(s2.get_key() == (j0.get_key() ^ j1.get_key() ^ j2.get_key() ^ j4.get_key()));
        CHECK(s2.get_key() == (j0.get_key() ^ reduction_set.get_key()));
    }

}

TEST_CASE("[Partial-order reduction] Precedence constraints") {
	Job<dtime_t> j0{0, I(0, 0), I(7, 13), 100, 0};
	Job<dtime_t> j1{1, I(6, 8), I(1, 2), 100, 1};
	Job<dtime_t> j2{2, I(7, 10), I(1, 2), 100, 2};
	Job<dtime_t> j3{3, I(12, 13), I(1, 3), 100, 3};
	Job<dtime_t> j4{4, I(5, 9), I(1, 2), 100, 4};
	Job<dtime_t> j5{5, I(25, 29), I(1, 2), 100, 5};

	std::vector<Job<dtime_t>> jobs{j0, j1, j2, j3, j4, j5};
	std::vector<const Job<dtime_t>*> reduction_jobs{&j1, &j2, &j4};
	Reduction_set<dtime_t> reduction_set{I(7, 13), reduction_jobs, {1, 2, 4}};
	Index_set scheduled_jobs{Index_set(), 0};

	SUBCASE("Job can possibly interfere if its predecessors have completed") {
		Job_precedence_set job_precedence_set;
		job_precedence_set.push_back(0); // J0 precedes J3

		CHECK(reduction_set.can_interfere(j3, job_precedence_set, scheduled_jobs));
	}

	SUBCASE("Job can possibly interfere if it has no predecessors") {
		Job_precedence_set job_precedence_set;

		CHECK(reduction_set.can_interfere(j3, job_precedence_set, scheduled_jobs));
	}

	SUBCASE("Job can possibly interfere if a job in the reduction set is its predecessor (case 1)") {
		Job_precedence_set job_precedence_set;
		job_precedence_set.push_back(1); // J1 precedes J3

		CHECK(reduction_set.can_interfere(j3, job_precedence_set, scheduled_jobs));
	}

	SUBCASE("Job can possibly interfere if a job in the reduction set is its predecessor (case 2)") {
		Job_precedence_set job_precedence_set;
		job_precedence_set.push_back(2); // J2 precedes J3

		CHECK(reduction_set.can_interfere(j3, job_precedence_set, scheduled_jobs));
	}

	SUBCASE("Job can possibly interfere if some jobs in the reduction set are its predecessors (case 1)") {
		Job_precedence_set job_precedence_set;
		job_precedence_set.push_back(1); // J1 precedes J3
		job_precedence_set.push_back(2); // J2 precedes J3


		CHECK(reduction_set.can_interfere(j3, job_precedence_set, scheduled_jobs));
	}

	SUBCASE("Job can possibly interfere if some jobs in the reduction set are its predecessors (case 2)") {
		Job_precedence_set job_precedence_set;
		job_precedence_set.push_back(1); // J1 precedes J3
		job_precedence_set.push_back(4); // J4 precedes J3

		CHECK(reduction_set.can_interfere(j3, job_precedence_set, scheduled_jobs));
	}

	SUBCASE("Job can possibly interfere if some jobs in the reduction set are its predecessors (case 3)") {
		Job_precedence_set job_precedence_set;
		job_precedence_set.push_back(2); // J2 precedes J3
		job_precedence_set.push_back(4); // J4 precedes J3

		CHECK(reduction_set.can_interfere(j3, job_precedence_set, scheduled_jobs));
	}

	SUBCASE("Job cannot interfere if all jobs in the reduction set are its predecessors (case 1)") {
		Job_precedence_set job_precedence_set;
		job_precedence_set.push_back(1); // J1 precedes J3
		job_precedence_set.push_back(2); // J2 precedes J3
		job_precedence_set.push_back(4); // J4 precedes J3

		CHECK(! reduction_set.can_interfere(j3, job_precedence_set, scheduled_jobs));
	}

	SUBCASE("Job cannot interfere if all jobs in the reduction set are its predecessors (case 2)") {
		Job_precedence_set job_precedence_set;
		job_precedence_set.push_back(0); // J0 precedes J3
		job_precedence_set.push_back(1); // J1 precedes J3
		job_precedence_set.push_back(2); // J2 precedes J3
		job_precedence_set.push_back(4); // J4 precedes J3

		CHECK(! reduction_set.can_interfere(j3, job_precedence_set, scheduled_jobs));
	}

	SUBCASE("Job cannot interfere if its predecessors have not completed yet") {
		Job_precedence_set job_precedence_set;
		job_precedence_set.push_back(5); // J5 precedes J3

		CHECK(! reduction_set.can_interfere(j3, job_precedence_set, scheduled_jobs));
	}

}

TEST_CASE("Topological sort precedence constraints") {
	Job<dtime_t> j0{0, I(0, 0), I(7, 13), 100, 0};
	Job<dtime_t> j1{1, I(6, 8), I(1, 2), 100, 1};
	Job<dtime_t> j2{2, I(7, 10), I(1, 2), 100, 2};
	Job<dtime_t> j3{3, I(12, 13), I(1, 3), 100, 3};
	Job<dtime_t> j4{4, I(5, 9), I(1, 2), 100, 4};
	Job<dtime_t> j5{5, I(25, 29), I(1, 2), 100, 5};

	std::vector<Job<dtime_t>> jobs{j0, j1, j2, j3, j4, j5};

	SUBCASE("Topological sort 1") {
		std::vector<Job_precedence_set> job_precedence_sets(jobs.size());
		job_precedence_sets[0].push_back(1); // J0 precedes J1
		job_precedence_sets[1].push_back(2);
		job_precedence_sets[2].push_back(3);
		job_precedence_sets[3].push_back(4);
		job_precedence_sets[4].push_back(5);

		auto topo_sorted = topological_sort<dtime_t>(job_precedence_sets, jobs);

		CHECK(topo_sorted[0].get_id() == JobID(5, 0));
		CHECK(topo_sorted[1].get_id() == JobID(4, 0));
		CHECK(topo_sorted[2].get_id() == JobID(3, 0));
		CHECK(topo_sorted[3].get_id() == JobID(2, 0));
		CHECK(topo_sorted[4].get_id() == JobID(1, 0));
		CHECK(topo_sorted[5].get_id() == JobID(0, 0));
	}

	SUBCASE("Topological sort 2") {
		std::vector<Job_precedence_set> job_precedence_sets(jobs.size());
		job_precedence_sets[1].push_back(0);
		job_precedence_sets[2].push_back(1);
		job_precedence_sets[3].push_back(2);
		job_precedence_sets[3].push_back(4);
		job_precedence_sets[4].push_back(1);
		job_precedence_sets[5].push_back(3);

		auto topo_sorted = topological_sort<dtime_t>(job_precedence_sets, jobs);

		CHECK(topo_sorted[0].get_id() == JobID(0, 0));
		CHECK(topo_sorted[1].get_id() == JobID(1, 0));
		CHECK((topo_sorted[2].get_id() == JobID(2, 0) || topo_sorted[2].get_id() == JobID(4, 0)));
		CHECK((topo_sorted[3].get_id() == JobID(2, 0) || topo_sorted[3].get_id() == JobID(4, 0)));
		CHECK(topo_sorted[4].get_id() == JobID(3, 0));
		CHECK(topo_sorted[5].get_id() == JobID(5, 0));
	}

}
