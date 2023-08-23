#include "doctest.h"

#include "global/por_space.hpp"
#include "global/por_criterion.hpp"
#include "global/reduction_set.hpp"

using namespace NP::Global;
using namespace NP;


TEST_CASE("[Partial-order reduction] MC Example 1") {
    //First make some jobs
    Job<dtime_t> j1{1, I(0, 10), I(2, 2), 100, 1};
    Job<dtime_t> j2{2, I(0, 10), I(2, 2), 100, 2};
    Job<dtime_t> j3{3, I(0, 10), I(2, 2), 100, 3};
    Job<dtime_t> j4{4, I(0, 10), I(2, 2), 100, 4};
    //Then make the input job set
    std::vector<const Job<dtime_t>*> jobs{&j1, &j2, &j3, &j4};

    //make the reduction set
	Reduction_set<dtime_t> reduction_set{{I(0, 0), I(0,0)}, jobs, {0, 1, 2, 3}};
    //update the set
    reduction_set.update_set();

    SUBCASE("Latest idle time") {
        CHECK(reduction_set.get_latest_idle_time() == 10);
    }

    SUBCASE("Latest start time") {
        CHECK(reduction_set.get_latest_start_time(j1) >= 11);
        CHECK(reduction_set.get_latest_start_time(j2) >= 11);
        CHECK(reduction_set.get_latest_start_time(j3) >= 12);
        CHECK(reduction_set.get_latest_start_time(j4) >= 12);
    }

    SUBCASE("No potential deadline miss") {
        CHECK(!reduction_set.has_potential_deadline_misses());
    }
    
}


TEST_CASE("[Partial-order reduction] MC Example 2") {
    //First make some jobs
    Job<dtime_t> j1{1, I(0, 10), I(2, 2), 100, 1};
    Job<dtime_t> j2{2, I(0, 10), I(2, 2), 100, 2};
    Job<dtime_t> j3{3, I(0, 10), I(2, 2), 100, 3};
    Job<dtime_t> j4{4, I(0, 10), I(2, 2), 100, 4};
    Job<dtime_t> j5{5, I(10, 10), I(1,5), 16, 5 };
    //Then make the input job set
    std::vector<const Job<dtime_t>*> jobs{&j1, &j2, &j3, &j4};

    //make the reduction set
	Reduction_set<dtime_t> reduction_set{{I(0, 0), I(0,0),I(0, 0), I(0,0)}, jobs, {0, 1, 2, 3}};
    //update the set
    reduction_set.update_set();

    SUBCASE("Latest idle time") {
        CHECK(reduction_set.get_latest_idle_time() == 10);
    }

    SUBCASE("Latest start time") {
        CHECK(reduction_set.get_latest_start_time(j1) >= 10);
        CHECK(reduction_set.get_latest_start_time(j2) >= 10);
        CHECK(reduction_set.get_latest_start_time(j3) >= 10);
        CHECK(reduction_set.get_latest_start_time(j4) >= 10);
    }

    SUBCASE("No potential deadline miss") {
        CHECK(!reduction_set.has_potential_deadline_misses());
    }

    SUBCASE("Add deadline missing job") {
        reduction_set.add_job(&j5, 5);
        reduction_set.update_set();
        CHECK(reduction_set.get_latest_start_time(j5) == 12);
        CHECK(reduction_set.has_potential_deadline_misses());
    }
}

//test for interference based on prio
TEST_CASE("[Partial-order reduction] MC Example prio interfering") {
    //First make some jobs
    Job<dtime_t> j1{1, I(10, 10), I(2, 2), 100, 1};
    Job<dtime_t> j2{2, I(10, 10), I(2, 2), 100, 2};
    Job<dtime_t> j3{3, I(10, 10), I(2, 2), 100, 3};
    Job<dtime_t> j4{4, I(10, 10), I(2, 2), 100, 4};
    Job<dtime_t> j5{5, I(10, 10), I(1,5), 16,  0};
    //Then make the input job set
    std::vector<const Job<dtime_t>*> jobs{&j1, &j2, &j3, &j4};

    //make the reduction set
	Reduction_set<dtime_t> reduction_set{{I(0, 0), I(0,0),I(0, 0), I(0,0)}, jobs, {0, 1, 2, 3}};
    //update the set
    reduction_set.update_set();

    SUBCASE("Latest idle time") {
        CHECK(reduction_set.get_latest_idle_time() == 10);
    }

    SUBCASE("Latest start time") {
        CHECK(reduction_set.get_latest_start_time(j1) >= 10);
        CHECK(reduction_set.get_latest_start_time(j2) >= 10);
        CHECK(reduction_set.get_latest_start_time(j3) >= 10);
        CHECK(reduction_set.get_latest_start_time(j4) >= 10);
    }

    SUBCASE("No potential deadline miss") {
        CHECK(!reduction_set.has_potential_deadline_misses());
    }

    SUBCASE("Check for prio interference") {
        CHECK(reduction_set.can_interfere(j5));
    }
}
//test for interference based on idle interval
TEST_CASE("[Partial-order reduction] MC Example idle interfering") {
    //First make some jobs
    Job<dtime_t> j1{1, I(0, 0), I(2, 2), 100, 1};
    Job<dtime_t> j2{2, I(0, 0), I(2, 2), 100, 2};
    Job<dtime_t> j3{3, I(2, 2), I(2, 2), 100, 3};
    Job<dtime_t> j4{4, I(0, 3), I(2, 2), 100, 4};
    Job<dtime_t> j5{5, I(0, 3), I(2, 2), 100, 100};
    //Then make the input job set
    std::vector<const Job<dtime_t>*> jobs{&j1, &j2, &j3, &j4};

    //make the reduction set
	Reduction_set<dtime_t> reduction_set{{I(0, 0), I(0,0)}, jobs, {0, 1, 2, 3}};
    //update the set
    reduction_set.update_set();

    SUBCASE("Latest idle time") {
        CHECK(reduction_set.get_latest_idle_time() == 3);
    }

    SUBCASE("Latest start time") {
        CHECK(reduction_set.get_latest_start_time(j1) >= 0);
        CHECK(reduction_set.get_latest_start_time(j2) >= 0);
        CHECK(reduction_set.get_latest_start_time(j3) >= 2);
        CHECK(reduction_set.get_latest_start_time(j4) >= 3);
    }

    SUBCASE("No potential deadline miss") {
        CHECK(!reduction_set.has_potential_deadline_misses());
    }
    
    SUBCASE("Check for prio interference") {
        CHECK(reduction_set.can_interfere(j5));
    }
    
}
//re-test the odd edge case, which is supposed to be unschedulable!
TEST_CASE("[Partial-order reduction] State space for partial-order reduction") {
    Job<dtime_t> j0{0, I(0, 100), I(643, 1607), 10000, 0};
    Job<dtime_t> j1{1, I(0, 100), I(153, 382),  10000, 1};
    Job<dtime_t> j2{2, I(0, 100), I(715, 1786), 10000, 2};
    Job<dtime_t> j3{5, I(0, 100), I(530, 1324), 10000, 3};
    Job<dtime_t> j4{3, I(0, 100), I(576, 1440), 15000, 4};
    Job<dtime_t> j5{4, I(0, 100), I(229, 572),  15000, 6};
    Job<dtime_t> j6{6, I(0, 100), I(404, 1009), 15000, 7};
    Job<dtime_t> j7{7, I(0, 100), I(176, 440),  15000, 8};
    Job<dtime_t> j8{8, I(0, 100), I(440, 1098), 25000, 13};
    Job<dtime_t> j9{9, I(0, 100), I(0, 32167),  35000, 24};
    Job<dtime_t> j10{10, I(0, 100), I(71, 177), 50000, 40};
    std::vector<Job<dtime_t>> jobs{j0, j1, j2, j3, j4, j5, j6, j7, j8, j9, j10};

    Analysis_options options{};

    auto space = Por_state_space<dtime_t, Null_IIP<dtime_t>, POR_priority_order<dtime_t>>::explore({jobs,4}, options);
    //This set is schedulable as the POR but not for non-POR with merge
    CHECK(space.is_schedulable());
}
