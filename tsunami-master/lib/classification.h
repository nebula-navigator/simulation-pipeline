//
// Created by lex on 13/01/19.
//

#ifndef TSUNAMI_CLASSIFICATION_H
#define TSUNAMI_CLASSIFICATION_H

#include <cmath>
#include <map>
#include <vector>
#include <fstream>
#include <iomanip>
#include <sstream>
#include "keplerutils.h"
#include "custom_types.hpp"

class Classification {

public:

	Classification(double Mscale, double Lscale) : KU(Mscale, Lscale) {}

	Classification() : Classification(1, 1) {}

	~Classification() = default;

    void set_units(double Mscale, double Lscale);
	void initialize_pinfo(pinfo *xdata, double *rads, double *mass, size_t Npart);

	struct TideTableElement {
		double taulag;
		double kaps;
		double polyt;
        double gyradius;

		inline TideTableElement() : taulag(-1.0), kaps(-1.0), polyt(-1.0), gyradius(-1.0) {};
	};
	std::map<short int, Classification::TideTableElement> tidetable;

protected:
    KeplerUtils KU;

};


class TripleClass : public Classification {

public:

    TripleClass(double Mscale, double Rscale, bool exitAbsorb) :
    Classification(Mscale, Rscale),
    exitAbsorb(exitAbsorb) {
        initialize_triple_stopcond();
    }
    TripleClass(bool exitAbsorb) : Classification(), exitAbsorb(exitAbsorb) {
        initialize_triple_stopcond();
    }

    ~TripleClass()  {
        if(not stopped) {
            info << "t = "<<ctime*KU.Tscale<<" yr: encounter incomplete"<<std::endl;
        } else {
            info << "t = " << std::fixed << ctime * KU.Tscale << std::setprecision(16) << std::scientific << " yr: encounter complete; Nhomo=" << Nhomo
                 << "; Nhomo_fex=" << Nhomo_first_ex << " (" << Nhomo_first_ex_clean << ")"
                 << "; longest_ex=" << longest_ex << " yr; longest_begin=" << longest_begin_time
                 << " yr; first_ex=" << first_ex << " yr; tex_cumulative=" << cumulative_ex_time
                 << " yr; last_ex=" << latest_ex << " yr; last_ex_begin=" << latest_ex_begin_time << " yr; Nex=" << Nex
                 << "; originalbin=" << isoriginalbin << "; hangle=" << hangle << "; status=" << status
                 << "; ahyp=" << a_hyp * KU.Lscale << "; ehyp=" << e_hyp << "; omegabin=" << omega_bin
                 << "; Nbarak=" << Nbarak_first_ex_clean << "; Nsteve=" << Nsteve;
            if (exitAbsorb) {
                info << "; ejection=" << ejection;
            }
            info << std::endl;
        }
        info.close();
        destruct_triple_stopcond();
    }

    void init_info_IO(std::string const &icname, bool cont);

    void collision_end() {
        info << "t = " << std::fixed << ctime * KU.Tscale << std::setprecision(16) << std::scientific << " yr: collision has occurred; Nhomo=" << Nhomo
             << "; Nhomo_fex=" << Nhomo_first_ex << " (" << Nhomo_first_ex_clean << ")"
             << "; longest_ex=" << longest_ex << " yr; longest_begin=" << longest_begin_time
             << " yr; first_ex=" << first_ex << " yr; tex_cumulative=" << cumulative_ex_time
             << " yr; last_ex=" << latest_ex << " yr; last_ex_begin=" << latest_ex_begin_time << " yr; Nex=" << Nex
             << "; status=" << status << "; Nbarak=" << Nbarak_first_ex_clean << "; Nsteve=" << Nsteve;
        info << std::endl;
    }

    void initialize_angmom(double3 *pos, double3 *vel, double *mass, size_t id1, size_t id2)  {
        double3 dp = pos[id1] - pos[id2];
        double3 dv = vel[id1] - vel[id2];
        double mu = mass[id1]*mass[id2] / (mass[id1] + mass[id2]);
        initial_angmom = mu * dp.cross(dv);
    }

    void initialize_triple_stopcond() {
        size_t N = 3;

        others = new size_t* [N];
        for(size_t i=0; i<N; i++)
            others[i] = new size_t [N-1];
        for(size_t i=0; i<N; i++) {
            size_t noth = 0;
            for(size_t k=0; k<N; k++) {
                if(k != i)
                    others[i][noth++] = k;
            }
        }

        pair_index = new size_t* [N];
        for(size_t i=0; i<N; i++)
            pair_index[i] = new size_t [N];
        size_t idpair = 0;
        for(size_t i=0; i<N; i++)
            for(size_t j=0; j<i; j++) {
                pair_index[i][j] = idpair;
                idpair++;
            }

        notinpair = new size_t* [N];
        for(size_t i=0; i<N; i++)
            notinpair[i] = new size_t [N];
        for(size_t i=0; i<N; i++)
            for(size_t j=0; j<i; j++)
                for(size_t k=0; k<N; k++) {
                    if (k == i or k == j) continue;
                    notinpair[i][j] = k;
                }
    }

    void destruct_triple_stopcond() {
        size_t N = 3;

        for(size_t i=0; i<N; i++)
            delete [] others[i];
        delete [] others;
        for(size_t i=0; i<N; i++)
            delete [] pair_index[i];
        delete [] pair_index;
        for(size_t i=0; i<N; i++)
            delete [] notinpair[i];
        delete [] notinpair;
    }

    bool triple_stopcond(double3 *pos, double3 *vel, double *mass, double current_time);

    double compute_hangle(double3 *pos, double3 *vel, size_t id1, size_t id2);

    struct PairOrbit {
        int i, j;  // Always i > j, if i or j is -1 then it's another pair
        double a;
        double e;
        double vdotr;
        double mcom;
        double3 pcom, vcom;
        double t2peri;
        double t2apo;

        PairOrbit *ipair = nullptr;
        PairOrbit *jpair = nullptr;
        size_t Nparents = 0;
        std::vector<PairOrbit*> parents;

        PairOrbit() = default;
        PairOrbit(int i, int j, double3 &dpos, double3 &dvel, double mu) {
            this->semimajor_eccentricity_vdotr(i, j, dpos, dvel, mu);
        }
        ~PairOrbit() = default;

        void semimajor_eccentricity_vdotr(size_t ii, size_t jj,
                                          const double3 &dpos, const double3 &dvel, double mu);

        void get_pair_com(double3 p1, double3 p2, double3 v1,
                          double3 v2, double m1, double m2);

        void time_to_peri(double3 dpos, double mu);

        bool mardl_stable(const PairOrbit& Outer);

        /*inline const PairOrbit& operator = (const PairOrbit &right){ //we do not have 'const' here because 'this' MUST be modified (assignment =)
            this->i = right.i;
            this->j = right.j;
            this->a = right.a;
            this->e = right.e;
            this->vdotr = right.vdotr;
            this->mcom = right.mcom;
            this->pcom = right.pcom;
            this->vcom = right.vcom;
            this->t2peri = right.t2peri;
            this->t2apo = right.t2apo;

            this->ipair = right.ipair;
            this->jpair = right.jpair;

            this->Nparents = right.Nparents;
            this->parents = right.parents;
            return *this;
        } */// Commented, not necessary
    };

    static double compute_ffactor(PairOrbit &Inner, PairOrbit &Outer, double3 *pos, double *mass);

private:
    bool exitAbsorb = true;
    bool ejection, excursion;

    double ctime = 0;
    bool stopped = false;
    double tform, pinn, pout, tdyn;
    double3 initial_angmom;

    size_t **others, **pair_index, **notinpair;
    int pair_i, pair_j;

    // Viraj's variables
    int Nhomo = 0;
    int Nbarak = 0;
    double homothres = 0.33;
    double ffactor_thres = 1e-5;
    double longest_ex = 0;
    double longest_begin_time = 0.0;
    double cumulative_ex_time = 0.0;
    double dhomo_prev = 0.0;
    double homo_prev = 0.0;
    double barak_prev = 0.0;
    double dbarak_prev = 0.0;
    double binsigle_ex_initial_time = 0.0;
    double first_ex = 0.0;
    double latest_ex = 0.0;
    double latest_ex_begin_time = 0.0;
    double binsigle_ex_vdotr = 0.0;
    double ffactor = 1.0;
    int Nhomo_first_ex = 1;
    int Nhomo_first_ex_clean = 1;
    int Nbarak_first_ex_clean = 1;
    int Nbarak_first_ex_temp, Nhomo_first_ex_temp;
    int Nex = 0;
    bool is_first_approach = true;
    bool isoriginalbin = false;
    double hangle = -1.0;
    double a_hyp = -1.0, e_hyp = -1.0;
    double omega_bin = -1.0;
    std::vector<double> pinn_list;
	bool rapid_exchange = 0;
	std::ostringstream ex_string;
    std::ostringstream message;
    int Nsteve;
    double avedist2_prev = 0;
    double davedist2_prev = 0;
    double davedist2 = 0;
    double avedist2_minimum = 0;
    double avedist2_maximum = 0;

    enum sys {
        NONE = 0,
        ISBROKEN_BINSINGLE,
        ISBROKEN_EXPLODED,
        ISINTERACTING_UNBOUND,
        ISINTERACTING_BOUND_BINSINGLE,
        ISINTERACTING_UNBOUND_BINSINGLE,
        ISHIERARCHICAL,
    } status = NONE;
    sys oldstatus = NONE;

    std::string infoname = "info.dat";
    std::ofstream info;

    double ecrit = 1e-3;
    double rpmult_crit = 10;
    bool critical_orbit(PairOrbit *a_e_vdotr, double3 *pos, size_t Npairs);

    bool check_unbound_system(double3 *pos, double3 *vel, double *mass);

    bool check_bound_binary(PairOrbit &Inner, double3 *pos, double3 *vel, double *mass);

    bool is_stable_hierarchical_triple(PairOrbit &Inner, PairOrbit &Outer, double3 *pos, double3 *vel, double *mass);
    bool is_interacting_binary_single(PairOrbit &Inner, PairOrbit &Outer, double3 *pos, double3 *vel, double *mass);

    void compute_orbit_properties(PairOrbit &Inner, PairOrbit &Outer, double3 *pos, double3 *vel);

    void print_coord(double3 *pos, double3 *vel, double *mass, size_t N);
};



#endif //TSUNAMI_CLASSIFICATION_H
