//
// Created by lex on 13/01/19.
//

#include <iomanip>
#include <algorithm>
#include "classification.h"
#include "Nbodyalgorithms.hpp"
#include "keplerutils.h"


void Classification::set_units(double Mscale, double Lscale) {
    KU = KeplerUtils(Mscale, Lscale);
}

void Classification::initialize_pinfo(pinfo *xdata, double *rads, double *mass, size_t Npart) {
	std::cout << "Tides summary:" << std::endl;
	std::cout << " ID" << std::setw(6) << "type" << std::setw(7) << "PN" << std::setw(9) << "Tides" << std::setw(10)
	          << "PolyInd" << std::setw(10) << "kap" << std::setw(18) << "tau (nbody)" << std::setw(7) << "A"
                << std::setw(13) << "sigma" << std::endl;

    std::cout << std::setprecision(3);

	for(size_t i=0; i < Npart; i++) {

		short int it = xdata[i].stype;

		/// STELLAR TIDE
		if((it >= ptype::LOW_MASS_MS) & (it < 100)) {
			xdata[i].hastide = true;
			if(tidetable.find(it) != tidetable.end()) {
				xdata[i].polyt = tidetable[it].polyt;
				xdata[i].kaps = tidetable[it].kaps;
				xdata[i].taulag = tidetable[it].taulag;
                double gyrad = tidetable[it].gyradius * rads[i];
                xdata[i].inert = gyrad*gyrad * mass[i];
			} else {
				xdata[i].polyt = 3.0;
				xdata[i].kaps = 0.014;
				xdata[i].taulag = 0.15;
                double gyrad = 0.25 * rads[i];
                xdata[i].inert = gyrad*gyrad * mass[i];
			}

		/// PLANET TIDE
		} else if((it >= ptype::ROCKY) & (it <= 1000)) {
			xdata[i].hastide = true;
			if (tidetable.find(it) != tidetable.end()) {
				xdata[i].polyt = tidetable[it].polyt;
				xdata[i].kaps = tidetable[it].kaps;
				xdata[i].taulag = tidetable[it].taulag;
                double gyrad = tidetable[it].gyradius * rads[i];
                xdata[i].inert = gyrad*gyrad * mass[i];
            } else {
				xdata[i].polyt = 3.0;
				xdata[i].kaps = 0.1;
				xdata[i].taulag = 0.3;
                double gyrad = 0.1 * rads[i];
                xdata[i].inert = gyrad*gyrad * mass[i];
			}
		}

		/// Convert taulag to Nbody
		xdata[i].taulag = xdata[i].taulag / (KeplerUtils::yr2sec * KU.Tscale);
        double Qt = 2 * xdata[i].kaps / (1 + 2 * xdata[i].kaps);
        double R5 = pow(rads[i], 5);

        xdata[i].Atide = R5 * Qt / (1 - Qt);
        xdata[i].sigmadiss = 2.0 * R5 / (xdata[i].Atide*xdata[i].Atide) * xdata[i].kaps * xdata[i].taulag / 3;

		if((xdata[i].polyt <= 0) and (xdata[i].taulag <= 0) and (xdata[i].kaps <= 0) and (xdata[i].inert <= 0)) {
		    xdata[i].hastide = false;
		    std::cout<<"Particle "<<i<<" has a tidal parameters less or equal than zero: disabling tides on this particle"<<std::endl;
		}
		std::cout << " P" << i << std::setw(5) << it << std::setw(8) << xdata[i].haspn << std::setw(7)
		          << xdata[i].hastide << std::setw(13) << xdata[i].polyt << std::setw(12) << xdata[i].kaps
		          << std::setw(13) << xdata[i].taulag << std::setw(12) << xdata[i].Atide
                << std::setw(12) << xdata[i].sigmadiss << std::endl;
	}
    std::cout << std::setprecision(15);
}



void TripleClass::init_info_IO(std::string const &icname, bool cont) {

    size_t found = icname.find_last_of('.');
    infoname = icname.substr(0, found) + "_info.dat";

    if (!cont) {
        remove(infoname.c_str()); // Removes info file if any (only when not restarting)
    }
    info << std::setprecision(16) << std::fixed;
    info.open(infoname.c_str(), std::ios::out | std::ios::app);
}

void TripleClass::PairOrbit::semimajor_eccentricity_vdotr(size_t ii, size_t jj,
                                                          const double3 &dpos, const double3 &dvel, double mu) {
    this->i = ii;
    this->j = jj;
    this->mcom = mu;
    double muinv = 1. / mu;
    double r = sqrt(dpos * dpos);
    double v2 = dvel * dvel;
    double vcirc2 = mu / r;
    this->a = -mu / (v2 - 2. * vcirc2); // semimajor axis
    double vdiff = v2 - vcirc2;
    this->vdotr = dpos * dvel;
    double3 e_vec = muinv * (vdiff * dpos - this->vdotr * dvel);
    this->e = sqrt(e_vec * e_vec); // eccentricity
}

void TripleClass::PairOrbit::get_pair_com(double3 p1, double3 p2, double3 v1,
                                          double3 v2, double m1, double m2) {
    this->mcom = m1 + m2;
    this->pcom = (m1 * p1 + m2 * p2) / this->mcom;
    this->vcom = (m1 * v1 + m2 * v2) / this->mcom;
}

void TripleClass::PairOrbit::time_to_peri(double3 dpos, double mu) {
    /**E_star = f_to_E(f_star, e_star)
    M_star = E_to_M(E_star, e_star)
    Mtop_star = 2*np.pi - M_star

    T_star = period(a_star, (self.mstar+self.mplanet+self.mbh))
    t_periapse = Mtop_star / (2*np.pi) * T_star**/

    double r = sqrt(dpos * dpos);
    double E;
    if (fabs(e - 1.0) < KeplerUtils::min_e) {
        double w = 0.5 / fabs(a);
        double wr = w * r;
        if (a > 0) {
            t2peri = (asin(sqrt(wr)) - sqrt(wr * (1 - wr))) / sqrt(2 * mu * w * w * w);
        } else {
            t2peri = (sqrt(wr * wr + wr) - log(sqrt(wr) + sqrt(wr + 1))) / sqrt(2 * mu * w * w * w);
        }
    } else {
        if (e < 1.) {
            E = KeplerUtils::acos2(1. - r / a, e, vdotr);
            E = KeplerUtils::mod2pi(E);
        } else {
            E = acosh((1. - r / a) / e);
            if (vdotr > 0.) {
                // Pericenter already passed
                t2peri = -1 * std::numeric_limits<double>::infinity();
                return;
            }
        }
        double M = KeplerUtils::E_to_M(E, e);
        double P = KeplerUtils::period(fabs(a), mu);
        t2peri = (1 - M / (2 * M_PI)) * P;
        t2apo = (1 - M / M_PI) * P;

        //std::cout << "Period (nbody) " << KeplerUtils::period(fabs(a), mu) << std::endl;
    }
}

bool TripleClass::PairOrbit::mardl_stable(const TripleClass::PairOrbit &Outer) {
    double ome = 1.0 - Outer.e;
    if (ome < KeplerUtils::min_e) return false;
    //std::cout << "out.mcom " << out.mcom << "    inn.mcom " <<  inn.mcom << std::endl;

    double qp1 = Outer.mcom / this->mcom;
    //double base = (q*(1.0+out.e)/(ome*ome*ome));
    //std::cout << "base " << base << "   q " << q << std::endl;
    double mardlstab = this->a * 2.8 * pow((qp1 * (1.0 + Outer.e) / (ome * ome * ome)), 0.4);

    //std::cout<<"oute  "<<out.e<<" outa  "<<out.a<<"  inna "<<inn.a<<"  inne "<<inn.e<<std::endl;
    //std::cout<<"Ma  "<<mardlstab<<"  bool "<<(out.a > mardlstab)<<std::endl;
    // If mardlstab is nan (for example when e ~ 1, it returns false, as it should be
    return (Outer.a > mardlstab);
}

void TripleClass::print_coord(double3 *pos, double3 *vel, double *mass, size_t N) {

    info << "t = " << ctime * KU.Tscale << " yr:   x" << std::setw(13) << "y" << std::setw(13) << "z" << std::setw(13)
         << "vx" << std::setw(13) << "vy" << std::setw(13) << "vz" << std::setw(13) << "m" << std::endl;
    for (size_t i = 0; i < N; i++) {
        info << std::setw(24) << pos[i].x << std::setw(13) << pos[i].y << std::setw(13) << pos[i].z << std::setw(13)
             << vel[i].x << std::setw(13) << vel[i].y << std::setw(13) << vel[i].z << std::setw(13)
             << mass[i] << std::endl;
    }
}

std::vector<double> trim_vector(double frac, std::vector<double> vec)
{
    size_t n_remove = ceil(vec.size()*frac);
    std::vector<double> new_vec;
    for (size_t i=0; i<vec.size(); i++) {
        if ((i > n_remove) and (i <= vec.size()-n_remove)) {
            new_vec.push_back(vec[i]);
        }
    }
    return new_vec;
}


bool TripleClass::triple_stopcond(double3 *pos, double3 *vel, double *mass, double current_time) {
    size_t N = 3;

    this->ctime = current_time;

    size_t Npairs = N * (N - 1) / 2;
    PairOrbit a_e_vdotr[Npairs];

    std::vector<size_t> iBound;
    std::vector<double> interdist;

    double avedist2 = 0;
    size_t kbound = 0;
    for (size_t i = 0; i < N; i++) {
        for (size_t j = 0; j < i; j++) {
            double3 dpos = pos[i] - pos[j];
            double3 dvel = vel[i] - vel[j];
            a_e_vdotr[kbound].semimajor_eccentricity_vdotr(i, j, dpos, dvel, mass[i] + mass[j]);
            double dist2 = dpos * dpos;
            avedist2 += dist2;
            interdist.push_back(sqrt(dist2));
            // Checks
            if (a_e_vdotr[kbound].a > 0) {
                a_e_vdotr[kbound].get_pair_com(pos[i], pos[j], vel[i], vel[j], mass[i], mass[j]);

                iBound.push_back(kbound);
            }
            kbound++;
        }
    }

    if (critical_orbit(a_e_vdotr, pos, Npairs))
        return false;


    // Homology radius
    double maxrad = *std::max_element(interdist.begin(), interdist.end());
    double minrad = *std::min_element(interdist.begin(), interdist.end());
    double homorad = minrad / maxrad;
    double dhomo = homorad - homo_prev;
    double barakrad = minrad*minrad / avedist2;
    double dbarakrad = barakrad - barak_prev;

    // Steve factor
    // Let's count extrema of avedist2
    if (status == sys::NONE) avedist2_maximum = avedist2;

    double davedist2 = avedist2 - avedist2_prev;
    if(davedist2 < 0 and davedist2_prev > 0) {
        // Now decreasing, previously increasing - it's a peak
        // avedist2_prev is the maximum
        avedist2_maximum = avedist2_prev;
    }

    if(davedist2 > 0 and davedist2_prev < 0) {
        // Now increasing, previously decreasing - it's a valley
        // avedist2_prev is the current minimum, avedist2_minimum the previous
        if(avedist2_maximum > 2*avedist2_prev and avedist2_maximum > 2*avedist2_minimum) {
            // If the maximum is greater than twice both minima
            Nsteve+=1;
        }
        avedist2_minimum = avedist2_prev;
    }
    avedist2_prev = avedist2;
    davedist2_prev = davedist2;

    if(dhomo < 0) {
        // Now decreasing
        if(dhomo_prev > 0 and homo_prev > homothres) {
            // If was growing at previous step and was above threshold
            Nhomo+=1;
        }
    } // This counts each peak above homothres
    homo_prev = homorad;
    dhomo_prev = dhomo;

    if(dbarakrad < 0) {
        // Now decreasing
        if(dbarak_prev > 0 and barak_prev > homothres*homothres) {
            // If was growing at previous step and was above threshold
            Nbarak+=1;
        }
    }
    barak_prev = barakrad;
    dbarak_prev = dbarakrad;

    bool was_binsingle_ex = false;
    int bs_i, bs_j;
	double bs_ex_initial_time;
    if (status == sys::ISINTERACTING_BOUND_BINSINGLE) {
        was_binsingle_ex = true;
        bs_i = pair_i;
        bs_j = pair_j;
        bs_ex_initial_time = binsigle_ex_initial_time;
    }

    bool stopret;
    // There is a bound pair
    if (not iBound.empty()) {
        PairOrbit PairInner;

        // More than two pairs are bound, system is interacting
        if (iBound.size() > 1) {
            // Get the Harder binary
            double amin = a_e_vdotr[iBound[0]].a;
            PairInner = a_e_vdotr[iBound[0]];
            for (size_t ib = 1; ib < iBound.size(); ib++) {
                if (a_e_vdotr[iBound[ib]].a < amin) {
                    amin = a_e_vdotr[iBound[ib]].a;
                    PairInner = a_e_vdotr[iBound[ib]];
                }
            }
        } else if (iBound.size() == 1) {
            PairInner = a_e_vdotr[iBound[0]];
        }

        stopret = check_bound_binary(PairInner, pos, vel, mass);
    } else {
        // No bound pairs
        stopret = check_unbound_system(pos, vel, mass);
    }

    if (was_binsingle_ex) {
        if (status != sys::ISINTERACTING_BOUND_BINSINGLE or (pair_i != bs_i and pair_j != bs_j)) {
            double ex_length = (ctime - bs_ex_initial_time) * KU.Tscale;
            if (not is_first_approach and (longest_ex < ex_length)) {
                longest_ex = ex_length;
                longest_begin_time = bs_ex_initial_time * KU.Tscale;
            }
            tdyn = Nbodyalgorithms::dynamical_timescale(pos, vel, mass, N);
            const auto median_it = pinn_list.begin() + pinn_list.size() / 2;
            std::nth_element(pinn_list.begin(), median_it , pinn_list.end());
            double pinn_median = *median_it;

            pinn_list.clear();

            if (ex_length > pinn_median * KU.Tscale) {
                info << ex_string.str();
                /*info << "t = " << time * Tscale
                     << " yr: bound binary-single " + hier + ", a="<< Outer.a << " au, e=" << Outer.e << ", next pericenter passage in "
                     << Outer.t2peri * Tscale << " yr (at " << (Outer.t2peri + time) * Tscale
                     << " yr)" << std::endl;*/
                if (not is_first_approach) {
                    cumulative_ex_time += ex_length;
                    Nex += 1;
                    latest_ex = ex_length;
                    latest_ex_begin_time = bs_ex_initial_time * KU.Tscale;
                    if (first_ex == 0.0) {
                        first_ex = ex_length;
                        Nhomo_first_ex = Nhomo;
                        Nhomo_first_ex_clean = Nhomo_first_ex_temp;
                        Nbarak_first_ex_clean = Nbarak_first_ex_temp;
                    }
                    if (exitAbsorb) {
                        ejection = false;
                        stopped = true;
                        stopret = stopped;
                    }
                }
            }
            if(is_first_approach) is_first_approach = false;
        }
    }
    if (not message.str().empty()) {
        info << message.str();
        message = std::ostringstream();
    }

    return stopret;
}


bool TripleClass::critical_orbit(PairOrbit *a_e_vdotr, double3 *pos, size_t Npairs) {

    for (size_t i=0; i<Npairs; i++) {
        double3 dp = pos[a_e_vdotr[i].i] - pos[a_e_vdotr[i].j];
        double rrr = sqrt(dp*dp);
        double rp = a_e_vdotr[i].a*(1-a_e_vdotr[i].e);
        double edelta = fabs(1 -a_e_vdotr[i].e);
        double rp_mult = sqrt(ecrit / edelta);

        if (edelta < ecrit and rrr < rp_mult * rpmult_crit * rp) {
            return true;
        }

    }
    return false;
}


bool TripleClass::check_unbound_system(double3 *pos, double3 *vel, double *mass) {
    size_t N = 3;

    // Check if bound wrt to the others
    bool is_bound_to_com = false;
    for (size_t i = 0; i < N; i++) {
        size_t j = others[i][0];
        size_t k = others[i][1];

        PairOrbit other2com;
        other2com.get_pair_com(pos[j], pos[k], vel[j], vel[k], mass[j], mass[k]);
        double3 dp = other2com.pcom - pos[i];
        double3 dv = other2com.vcom - vel[i];
        PairOrbit aevdotr_other2(i, 0, dp, dv, other2com.mcom + mass[i]);
        if (aevdotr_other2.a > 0) {
            is_bound_to_com = true;
            break;
        }
    }

    if (is_bound_to_com) {
        if (status != sys::ISINTERACTING_UNBOUND) {
            message << "t = " << ctime * KU.Tscale << " yr: no bound pairs but system is interacting" << std::endl;
            status = sys::ISINTERACTING_UNBOUND;
        }
    } else if (status != sys::ISBROKEN_EXPLODED) {
        message << "t = " << ctime * KU.Tscale << " yr: triple has broken, everything is unbound" << std::endl;
        tform = ctime;
        tdyn = Nbodyalgorithms::dynamical_timescale(pos, vel, mass, N);
        //print_coord(pos, vel, mass, N);
        status = sys::ISBROKEN_EXPLODED;
    } else {
        double elapse = ctime - tform;
        if (elapse > tdyn) {
            message << "t = " << ctime * KU.Tscale << " yr: one dynamical timescale elapsed from breakup ("
                 << tdyn * KU.Tscale << " yr), stopping" << std::endl;
            stopped = true;
            return stopped;
        }
    }
    return false;
}

bool TripleClass::check_bound_binary(PairOrbit &Inner, double3 *pos, double3 *vel, double *mass) {
    size_t k = notinpair[Inner.i][Inner.j];
    double3 dp = Inner.pcom - pos[k];
    double3 dv = Inner.vcom - vel[k];
    PairOrbit PairOuter(k, 0, dp, dv, Inner.mcom + mass[k]);

    rapid_exchange = (pair_i != Inner.i and pair_j != Inner.j);
    if (status == sys::NONE) rapid_exchange = false;
    ffactor = compute_ffactor(Inner, PairOuter, pos, mass);

    if (ffactor > 1) {
        status = sys::ISINTERACTING_UNBOUND;
        return false;
    }

    // Below here is only when ffactor < 1
    if (PairOuter.a > 0) {
        return is_stable_hierarchical_triple(Inner, PairOuter, pos, vel, mass);
    } else {
        return is_interacting_binary_single(Inner, PairOuter, pos, vel, mass);
    }
}


bool TripleClass::is_stable_hierarchical_triple(PairOrbit &Inner, PairOrbit &Outer, double3 *pos, double3 *vel,
                                                double *mass) {
    if (Inner.mardl_stable(Outer)) {
        if (status != sys::ISHIERARCHICAL or rapid_exchange) {
            std::string hier =
                    "(" + std::to_string(Inner.i) + "+" + std::to_string(Inner.j) + "," + std::to_string(Outer.i) + ")";
            message << "t = " << ctime * KU.Tscale << " yr: hierarchical triple " + hier << std::endl;
            //print_coord(pos, vel, mass, N);
            tform = ctime;
            pinn = KeplerUtils::period(Inner.a, Inner.mcom);
            pout = KeplerUtils::period(Outer.a, Outer.mcom);
            pair_i = Inner.i;
            pair_j = Inner.j;
            status = sys::ISHIERARCHICAL;
            //message <<std::endl<<"Pinn "<<pinn<<"   ainn "<< Inner.a<<"   aout "<<Outer.a<<"   eout "<<Outer.e<<std::endl;
        } else {
            double elapse = ctime - tform;
            // Stable hierararchical triple, let's stop after 50 outer periods
            if (elapse > 50 * pout) {
                message << "t = " << ctime * KU.Tscale
                     << " yr: elapsed time since hierarchical triple formation is >50 Pout, stopping"
                     << std::endl;

                if ((Inner.i == 1 and Inner.j == 2) or (Inner.j == 1 and Inner.i == 2))
                    isoriginalbin = true;
                hangle = compute_hangle(pos, vel, Inner.i, Inner.j);

                stopped = true;
                return stopped;
            }
        }
    } else {
        std::string hier = "(" + std::to_string(Inner.i) + "+" + std::to_string(Inner.j) + "," + std::to_string(Outer.i) + ")";
        if (status != sys::ISINTERACTING_BOUND_BINSINGLE or rapid_exchange) {
            pinn = KeplerUtils::period(Inner.a, Inner.mcom);
            pinn_list.push_back(pinn);
            binsigle_ex_initial_time = ctime;
            Nhomo_first_ex_temp = Nhomo;
            Nbarak_first_ex_temp = Nbarak;
            Outer.time_to_peri(Inner.pcom - pos[Outer.i], Inner.mcom + mass[Outer.i]);
            binsigle_ex_vdotr = Outer.vdotr;

            /*message << "*t = " << time * Tscale
                 << " yr: bound binary-single " + hier + ", a="<< Outer.a << " au, e=" << Outer.e << ", next pericenter passage in "
                 << Outer.t2peri * Tscale << " yr (at " << (Outer.t2peri + time) * Tscale
                 << " yr)" << std::endl;*/
            pair_i = Inner.i;
            pair_j = Inner.j;
            status = sys::ISINTERACTING_BOUND_BINSINGLE;
        } else if (status == sys::ISINTERACTING_BOUND_BINSINGLE and not rapid_exchange) {
            double oldvdotr = binsigle_ex_vdotr;
            Outer.time_to_peri(Inner.pcom - pos[Outer.i], Inner.mcom + mass[Outer.i]);
            // Check for apocenter passage
            //double3 r = (Inner.pcom - pos[Outer.i]);
            //info << sqrt(r*r) << std::endl;
            if ((oldvdotr > 0.0) and (Outer.vdotr < 0.0)) {
                ex_string = std::ostringstream();
                ex_string << "t = " << binsigle_ex_initial_time * KU.Tscale
                          << " yr: bound binary-single " + hier + ", a="<< Outer.a << " au, e=" << Outer.e << ", next pericenter passage in "
                          << Outer.t2peri * KU.Tscale << " yr (at " << (Outer.t2peri + ctime) * KU.Tscale
                          << " yr)" << std::endl;
            }
            binsigle_ex_vdotr = Outer.vdotr;

            if (exitAbsorb and ffactor < ffactor_thres) {
                if (ex_string.str().empty()) {
                    ex_string << "t = " << binsigle_ex_initial_time * KU.Tscale
                              << " yr: bound binary-single " + hier + ", a="<< Outer.a << " au, e=" << Outer.e << ", next pericenter passage in "
                              << Outer.t2peri * KU.Tscale << " yr (at " << (Outer.t2peri + ctime) * KU.Tscale
                              << " yr)" << std::endl;
                }

                info << ex_string.str();
                message << "t = " << ctime * KU.Tscale
                     << " yr: bound binary-single " + hier + ", ffactor < " << ffactor_thres << ", next pericenter passage in "
                     << Outer.t2peri * KU.Tscale  << " yt, stopping now" << std::endl;

                // We can set these because if we are here, the interaction was not stopped earlier
                Nhomo_first_ex = Nhomo;
                Nhomo_first_ex_clean = Nhomo_first_ex_temp;
                Nbarak_first_ex_clean = Nbarak_first_ex_temp;

                ejection = false;
                stopped = true;
                return stopped;
            }

	        pinn = KeplerUtils::period(Inner.a, Inner.mcom);
	        pinn_list.push_back(pinn); //TODO: SET LIMIT
        }
    }
    return false;
}

bool TripleClass::is_interacting_binary_single(PairOrbit &Inner, PairOrbit &Outer, double3 *pos,
                                               double3 *vel, double *mass) {

    /*ffactor = compute_ffactor(Inner, Outer, pos, mass);

    if (ffactor > 1) {
        status = sys::ISINTERACTING_UNBOUND;
        return false;
    }*/ //NOT NEEDED

    // Triple is unbound, check if they are approaching pericenter
    if (Outer.vdotr > 0) {
        if (status != sys::ISBROKEN_BINSINGLE or rapid_exchange) {
            Nhomo_first_ex_temp = Nhomo;
            Nbarak_first_ex_temp = Nbarak;

            std::string hier =
                    "(" + std::to_string(Inner.i) + "+" + std::to_string(Inner.j) + "," + std::to_string(Outer.i) + ")";
            message << "t = " << ctime * KU.Tscale << " yr: unbound binary-single " + hier + " on diverging orbits"
                 << std::endl;

            pair_i = Inner.i;
            pair_j = Inner.j;
            status = sys::ISBROKEN_BINSINGLE;
        } else {
            double3 dpos = Inner.pcom - pos[Outer.i];
            double dist2 = dpos * dpos;

            // Stable binary + single, let's stop after single is 100 abin far
            if (dist2 > 1e4 * Inner.a * Inner.a) {
                double vinf = -(Inner.mcom + mass[Outer.i]) / Outer.a;
                vinf = sqrt(vinf) * KU.Vscale;
                message << "t = " << ctime * KU.Tscale
                        << " yr: binary-single separation is >100 abin, stopping (escaper id=" << std::setprecision(16) << std::scientific << Outer.i << ", m="
                        << mass[Outer.i] * KU.Mscale << " MSun, vesc=" << vinf << " km/s, abin="
                        << Inner.a * KU.Lscale << " au, ebin=" << Inner.e << ")"
                        << std::endl;

                if ((Inner.i == 1 and Inner.j == 2) or (Inner.j == 1 and Inner.i == 2))
                    isoriginalbin = true;
                hangle = compute_hangle(pos, vel, Inner.i, Inner.j);
                compute_orbit_properties(Inner, Outer, pos, vel);

                if (first_ex == 0.0) {
                    Nhomo_first_ex = Nhomo;
                    Nhomo_first_ex_clean = Nhomo_first_ex_temp;
                    Nbarak_first_ex_clean = Nbarak_first_ex_temp;
                }

                if (exitAbsorb) {

                    ejection = true;
                    stopped = true;
                    return stopped;
                }

                stopped = true;
                return stopped;
            }
        }
    } else {
        if (status != sys::ISINTERACTING_UNBOUND_BINSINGLE or rapid_exchange) {
            std::string hier =
                    "(" + std::to_string(Inner.i) + "+" + std::to_string(Inner.j) + "," + std::to_string(Outer.i) + ")";
            Outer.time_to_peri(Inner.pcom - pos[Outer.i], Inner.mcom + mass[Outer.i]);
            message << "t = " << ctime * KU.Tscale
                 << " yr: unbound binary-single " + hier + " next pericenter passage in "
                 << Outer.t2peri * KU.Tscale << " yr (at " << (Outer.t2peri + ctime) * KU.Tscale
                 << " yr)" << std::endl;
            pair_i = Inner.i;
            pair_j = Inner.j;
            status = sys::ISINTERACTING_UNBOUND_BINSINGLE;
        }
    }
    return false;
}

double TripleClass::compute_hangle(double3 *pos, double3 *vel, size_t id1, size_t id2) {
    double3 dp = pos[id1] - pos[id2];
    double3 dv = vel[id1] - vel[id2];
    double3 h = dp.cross(dv);
    h = h / sqrt(h*h);
    double3 h0 = initial_angmom / sqrt(initial_angmom*initial_angmom);
    double costh = h*h0;
    return acos(costh);
}

double TripleClass::compute_ffactor(PairOrbit &Inner, PairOrbit &Outer, double3 *pos, double *mass) {
    double apodist = (Inner.a * (1 + Inner.e));
    double frel = (mass[Inner.i] * mass[Inner.j]) / (apodist*apodist);

    double3 p3 = Inner.pcom - pos[Outer.i];
    double outerdist = sqrt(p3*p3);
    double ftid = 2*(Inner.mcom)*mass[Outer.i] * apodist / (outerdist*outerdist*outerdist);

    return ftid/frel;
}


void TripleClass::compute_orbit_properties(PairOrbit &Inner, PairOrbit &Outer, double3 *pos, double3 *vel) {

    double3 outpos = pos[Outer.i] - Inner.pcom;
    double3 outvel = vel[Outer.i] - Inner.vcom;
    double out_posvel[6] = {outpos.x, outpos.y, outpos.z, outvel.x, outvel.y, outvel.z};
    double out_mu = Outer.mcom;
    double out_keplpar[6] = {};

    KeplerUtils::cart_to_kepl(out_keplpar, out_posvel, out_mu);
    a_hyp = out_keplpar[0] * KU.Lscale;
    e_hyp = out_keplpar[1];

    double3 innpos = pos[Inner.i] - pos[Inner.j];
    double3 innvel = vel[Inner.i] - vel[Inner.j];
    double inn_posvel[6] = {innpos.x, innpos.y, innpos.z, innvel.x, innvel.y, innvel.z};
    double inn_mu = Inner.mcom;
    double inn_keplpar[6] = {};

    KeplerUtils::cart_to_kepl(inn_keplpar, inn_posvel, inn_mu);
    omega_bin = inn_keplpar[3];
}
