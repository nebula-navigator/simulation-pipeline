//
// Created by lex on 21/12/18.
//

#ifndef TSUNAMI_SIMPROF_H
#define TSUNAMI_SIMPROF_H

#include <vector>
#include <ctime>
#include <iostream>

class SimProf {

public:
	SimProf() = default;

	~SimProf() = default;

	timespec tt1, tt2;
    double avg;
    double sum = 0.0;
    size_t Nmean = 0;

	void start() {
		clock_gettime(CLOCK_REALTIME, &tt1);
	}

	void stop_store() {
		clock_gettime(CLOCK_REALTIME, &tt2);
		timespec diff;
		diff.tv_sec = tt2.tv_sec - tt1.tv_sec;
		diff.tv_nsec = tt2.tv_nsec - tt1.tv_nsec;
		double time = diff.tv_sec + diff.tv_nsec / 1.0e9;

        sum += time;

        if (Nmean == 0) {
            avg = time;
            Nmean++;
        } else {
            double new_avg = (avg*static_cast<double>(Nmean) + time) / static_cast<double>(Nmean+1);
            Nmean++;
            avg = new_avg;
        }
	}

	void stop_print() {
		clock_gettime(CLOCK_REALTIME, &tt2);
		timespec diff;
		diff.tv_sec = tt2.tv_sec - tt1.tv_sec;
		diff.tv_nsec = tt2.tv_nsec - tt1.tv_nsec;
		double time = diff.tv_sec + diff.tv_nsec / 1.0e9;
		std::cout << "Ind : " << time << std::endl;
	}

	void print_avg(const string &pref = "") const {
		std::cout << pref << " ";
		std::cout  << avg << std::endl;
	}

	void print_sum(const string &pref = "") const {
		std::cout << pref << " ";
		std::cout << sum << std::endl;
	}

    void reset() {
        sum = 0.0;
        Nmean = 0.0;
    }

};


#endif //TSUNAMI_SIMPROF_H
