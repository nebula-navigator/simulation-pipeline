//
// Created by lex on 2/27/23.
//

#ifndef TSUNAMI_TSUNAMI_HPP
#define TSUNAMI_TSUNAMI_HPP

#include <array>
#include <iomanip>
#include "config.hpp"
#include "chain.hpp"
#include "custom_types.hpp"
#include "classification.h"
#include "Nbodyalgorithms.hpp"
#include "leapfrog_stepped.hpp"
#include "bulirsch.hpp"
#include "keplerutils.h"
#include "simprof.hpp"

// Quick stuff, will be removed after revision of gravitational_waveform
inline double vdot(const double a[3], const double b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

inline void vcross(const double a[3], const double b[3], double c[3]) {
    c[0] = a[1] * b[2] - a[2] * b[1];
    c[1] = a[2] * b[0] - a[0] * b[2];
    c[2] = a[0] * b[1] - a[1] * b[0];
}

inline void vnormalize(double a[3]) {
    double norm = std::sqrt(vdot(a, a));
    a[0] = a[0] / norm;
    a[1] = a[1] / norm;
    a[2] = a[2] / norm;
}

template<bool profiling, bool debug>
class TsunamiClass {
public:
    /**
     * Mass and distance unit in MSun and au.
     * Defaults to MSun and a.u.
     */
    TsunamiClass() : TsunamiClass(1, 1) {}

    TsunamiClass(double Mscale, double Lscale) :
                 System(Conf),
                 Leapfrog(System),
                 BSExtra(System, Leapfrog, tolerance) {
        std::cout << std::scientific << std::setprecision(15);

        set_units(Mscale, Lscale);
    }

    void set_units(double Mscale, double Lscale) {
        Conf.Lscale = Lscale;
        Conf.Mscale = Mscale;
        this->Mscale = Mscale;
        this->Lscale = Lscale;

        // Time unit (in yr)
        Tscale = sqrt(Lscale*Lscale*Lscale/(Mscale*KeplerUtils::G_yr_msun_au));

        // Velocity unit (in km/s)
        Vscale = Lscale/Tscale * KeplerUtils::au2km / KeplerUtils::yr2sec;
        speed_of_light = KeplerUtils::c_ms * 1e-3 / Vscale;

        System.init_postnewtonians(Mscale, Lscale);
        Marx.set_units(Mscale, Lscale);
    }

    void allocate_arrays() {
        System.Npart = N;
        System.Nchain = N - 1;
        System.allocate_arrays();

        aux_vec1.resize(N);
        aux_vec2.resize(N);
    }

    void initialize_particle_tides() {
        IO::read_tidal_table(Marx.tidetable, tidefile);
        Marx.initialize_pinfo(System.xdata, System.radius, System.mass, N);
    }

    void initialize_regularization(double alpha, double beta, double gamma) {
        Conf.alpha = alpha;
        Conf.beta = beta;
        Conf.gamma = gamma;
        Conf.TTL = (beta != 0);
        System.init_regularization();
    }

    void initialize_chain() {
        System.init_chain();
        System.initialize_chain_arrays();
    }

    void initialize_integrator() {
        Leapfrog.initialize();
        BSExtra.initialize();
        regularization_switch();
        System.init_tdetracker();
    }

    void do_step_leapfrog() {
        Leapfrog.save_step_zero(System.ch_pos, System.ch_vel, System.ch_mass,
                                System.ch_spin, System.ch_xdata);

        System.time0 = time;
        Leapfrog.integrate(nstep, timestep, dtphysical, System.ch_pos,
                           System.ch_vel, System.ch_spin);

        // B0 and omega0 are not updated yet
    }

    void do_step_bulirsh() {
        System.time0 = time;

        if constexpr (profiling) ProfBS.start();
        BSExtra.bs_iterate(dtphysical, timestep);
        if constexpr (profiling) ProfBS.stop_store();

        // B0 and omega0 are not updated yet
    }

    void regularization_switch() {
        // Using last Leapfrog integration as a proxy
        if constexpr (debug) {
            std::cout << " == REG SWITCH == " << std::endl;
            std::cout << "    B = " << Leapfrog.B << std::endl;
            std::cout << "    T = " << Leapfrog.T << std::endl;
            std::cout << "    U = " << Leapfrog.B + Leapfrog.T << std::endl;
            std::cout << "    U/T = " << (Leapfrog.B + Leapfrog.T)/Leapfrog.T << std::endl;
        }
        if (Leapfrog.B < (Conf.Usafe - 1.0) * Leapfrog.T) {
            if (Conf.gamma == 0.0) {
                if constexpr (debug) {
                    std::cout << "    U < " << Conf.Usafe << " T " << std::endl;
                    std::cout << "Setting gamma = 1 at T = " << time << std::endl;
                }
                switch_list.emplace_back(time, 1.0);
                Conf.gamma = 1.0;
            }
        } else {
            if (Conf.gamma == 1.0) {
                if constexpr (debug) {
                    std::cout << "    U > " << Conf.Usafe << " T " << std::endl;
                    std::cout << "Setting gamma = 0 at T = " << time << std::endl;
                }
                switch_list.emplace_back(time, 0.0);
                Conf.gamma = 0.0;
            }
        }
    }

    void update_time_coordinates_chain() {
        if constexpr (profiling) ProfChain.start();

        // If using the Leapfrog only, it should be Leapfrog.B and Leapfrog.omega instead
        Leapfrog.B0 = BSExtra.FinalTable.B;
        Leapfrog.omega0 = BSExtra.FinalTable.omega;

        time += dtphysical;
        Eintegrator = Leapfrog.B0;
        regularization_switch();

        // Chain check has been done earlier

        if (System.chain_has_changed()) {
            System.to_new_chain();
        }
        if constexpr (profiling) ProfChain.stop_store();
    }

    bool iterate_to_collision_bulirsch(double dcoll_thr = 1e-3) {
        size_t id1 = System.CollInfo.collind[0];
        size_t id2 = System.CollInfo.collind[1];
        Nbodyalgorithms::CollType collflag = System.CollInfo.collflag;
        make_ascending(id1, id2);

        double dcoll_current = System.get_distance_between_collided();
        double sumR = System.radius[id1]+System.radius[id2];
        double dcoll_target = sumR * (1 + dcoll_thr);

        double ds_coll = System.CollInfo.colltime_ds;
        double ds_original = BSExtra.dt_prev;
//        std::cout << "current_ds = " << ds_original << "\nds_coll " << ds_coll << std::endl;
//        std::cout << "dcoll/Rsum = " << dcoll_current/sumR << std::endl;
//        std::cout << "coll     i,j " << id1 << "," << id2 << std::endl;

        size_t nits = 0;

        double ds1 = 0.0;
        double ds2 = ds_original;
        double ds_try = (ds1 + ds2)*0.5;

        // If already in a collision, this is a rare situation where only the final state is in a collision
        // return and issue a warning. not a bug per se
        Leapfrog.revert_step(System.ch_pos, System.ch_vel, System.ch_mass,
                             System.ch_spin, System.ch_xdata);
        System.time0 = time;
        double dcoll_initial = System.get_distance_between_collided();
        if (dcoll_initial < sumR) {
            std::cerr << "WARNING: this collision event cannot be iterated close to the collision radius, because it appears that the"
                         " system is already collided at the beginning of the step.\nPlease send the initial conditions to Alessandro\n";
            return true;
        }

        while (true) {
            Leapfrog.revert_step(System.ch_pos, System.ch_vel, System.ch_mass,
                                 System.ch_spin, System.ch_xdata);
            System.time0 = time;

//            std::cout << " -------------------- ds_try/ds0 " << ds_try/ds_original << std::endl;
            double requested_step = ds_try;
            BSExtra.bs_iterate(dtphysical, ds_try);
            dcoll_current = System.get_distance_between_collided();

            if (BSExtra.dt_prev < requested_step and not ((System.CollInfo.collflag) or (dcoll_current < sumR))) {
                // Understep, because of rejection. Let's advance to this time and reset our starting time
//                std::cout << "  asked           = " << requested_step << std::endl;
//                std::cout << "  made step ds    = " << BSExtra.dt_prev << std::endl;
                update_time_coordinates_chain();
                ds_original = ds_original - BSExtra.dt_prev;
                ds_try = requested_step - BSExtra.dt_prev;
                ds2 = ds2 - BSExtra.dt_prev;
//                std::cout << "  NEW original_ds = " << ds_original  << std::endl;
//                std::cout << "  NEW ds_try      = " << ds_try  << std::endl;
//                std::cout << "  NEW ds2         = " << ds2  << std::endl;
                continue;
            }

            if (System.CollInfo.collflag or dcoll_current < sumR) {
                // If we crossed the collision radius, we need to decrease the timestep
                if (System.CollInfo.collflag) {
                    size_t idd1 = System.CollInfo.collind[0];
                    size_t idd2 = System.CollInfo.collind[1];
                    make_ascending(idd1, idd2);
                    if ((idd1 != id1) or (idd2 != id2)) {
                        throw TsuError("CONGRATULATIONS!! You found an edge case that Alessandro predicted, but was too lazy to implement,\n"
                                       "(it's never going to happen, he thought), so he added this error message.\n"
                                       "Now, send him this message: 'FOUND COLLISION THAT CHANGES IDENTITY'\n"
                                       "and he's going to implement this edge case within one working day.");
                    }
                }
                ds2 = BSExtra.dt_prev;
//                std::cout << "  COLL   -  ds/ds0  " << BSExtra.dt_prev/ds_original << "  (dcoll: " << dcoll_current/sumR << ")" << std::endl;
//                std::cout << " ds1/ds0     " << ds1/ds_original << std::endl;
//                std::cout << " ds2/ds0     " << ds2/ds_original << std::endl;
//                std::cout << " coll    i,j " << id1 << "," << id2 << std::endl;
                ds_try = (ds1 + ds2)*0.5;
//                std::cout << " next ds/d0: " << ds_try/ds_original << std::endl;
            } else {
                if ((sumR < dcoll_current) and (dcoll_current < dcoll_target)) {
                    // No collision, but past target. We let's stop
//                    std::cout << "  NOCOLL REACHED - ds/ds0  " << BSExtra.dt_prev/ds_original << "  (dcoll: " << dcoll_current/sumR << ")" << std::endl;
//                    std::cout << "          dcoll_current  " << dcoll_current << std::endl;
//                    std::cout << "          dcoll_target   " << dcoll_target << std::endl;
                    // Let's set the collision data for the Python interface
                    System.CollInfo.collind[0] = id1;
                    System.CollInfo.collind[1] = id2;
                    System.CollInfo.collflag = collflag;
                    break;
                } else {
                    // Bisect between current timestep and ds2
                    ds1 = BSExtra.dt_prev;
//                    std::cout << "  NOCOLL - ds/ds0  " << BSExtra.dt_prev/ds_original << "  (dcoll: " << dcoll_current/sumR << ")" << std::endl;
//                    std::cout << " ds1/ds0     " << ds1/ds_original << std::endl;
//                    std::cout << " ds2/ds0     " << ds2/ds_original << std::endl;
//                    ds_try = (ds1 + ds2)*0.5;
//                    std::cout << " next ds/d0: " << ds_try/ds_original << std::endl;
                }
            }

            nits++;
            if (fabs(ds1 - ds2)/ds_original < 1e-15) {
                timestep = ds_original * 1e-2;
                System.CollInfo.collind[0] = 0;
                System.CollInfo.collind[1] = 0;
                System.CollInfo.collflag = Nbodyalgorithms::NOCOLLISION;
                Leapfrog.revert_step(System.ch_pos, System.ch_vel, System.ch_mass,
                                     System.ch_spin, System.ch_xdata);
                System.time0 = time;
                BSExtra.bs_iterate(dtphysical, ds1);
                System.update_from_chain_to_com();
                System.find_chain();
//                std::cout << "possible false positive " << std::endl;
                return false;
//                std::cout << " we are here " << nits << std::endl;
//                std::cout << " time now " << time << std::endl;
//                size_t nsteps = 1000;
//                double ds = ds_original / static_cast<double>(nsteps);
//                for (size_t i=0; i<nsteps; i++) {
//                    ds_try = ds;
//                    BSExtra.bs_iterate(dtphysical, ds_try);
//                    if (check_collision()) {
//                        dcoll_current = System.get_distance_between_collided();
//                        std::cout << "FOUND COLL " << std::endl;
//                        std::cout << "dcoll_current " << dcoll_current << std::endl;
//                        break;
//                    }
//                    dcoll_current = System.get_distance_between_collided();
//                    std::cout << "distance " << dcoll_current/sumR << std::endl;
//                    update_time_coordinates_chain();
//                }
//                std::cout << " time now " << time * Tscale << std::endl;
//                return check_collision();
            }

            if (nits > 1000) {
                throw TsuError("Too many iterations, probably a false positive collision\nds_coll = " + n2s(ds_coll)
                               + "\nsumR = " + n2s(sumR) + "\ndcoll_target = " + n2s(dcoll_target));
            }
        }
        System.update_from_chain_to_com();
        System.find_chain();
        update_time_coordinates_chain();
        return true;
    }

    void print_profiling() {
        if constexpr (not profiling) {
            throw TsuError("Profiling has been disabled, compile with cmake option -Dprofile=on");
        }
        ProfBS.print_avg("BS");
        double BSonly = ProfBS.avg - Leapfrog.ProfLeap.avg;
        std::cout << "BS only: " << BSonly << std::endl;
        Leapfrog.ProfLeap.print_avg("Leap");
        ProfChain.print_avg("Chain");
    }

    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
    ///       __  _   _ ___ _  _  __  _  _    ____ _  _ _  _  ___ ___ _  __  _  _  ___                               ///
    ///      |__]  \_/   |  |__| |  | |\ |    |___ |  | |\ | |     |  | |  | |\ | [__                                ///
    ///      |      |    |  |  | |__| | \|    |    |__| | \| |___  |  | |__| | \| ___]                               ///
    ///                                                                                                              ///
    ////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

    void add_particle_set(double *pos_in, size_t npos, size_t pncoord,
                          double *vel_in, size_t nvel, size_t vncoord,
                          double *mass_in, size_t nmass,
                          double *rad_in, size_t nrad,
                          long *stype_in, size_t nstype) {

        if((nmass != npos) or (nmass != nvel) or (nmass != nstype) or (nmass != nrad)) {
            throw TsuError("Provided arrays have inconsistent length");
        }
        if((pncoord != 3) or (vncoord != 3)) {
            throw TsuError("Provide 3 coordinates for each particle position and velocity");
        }

        size_t Npart = nmass;
        if(Npart < 2) {
            throw TsuError("Provide at least two particle");
        }

        for(size_t j = 0; j < Npart; j++) {
            if(std::isnan(mass_in[j]) or std::isnan(rad_in[j]) or std::isnan(stype_in[j])) {
                throw TsuError("NaN values in input (mass or radius or type)");
            }
            for(size_t i = 0; i < 3; i++) {
                if(std::isnan(pos_in[i+j*3]) or std::isnan(vel_in[i+j*3])) {
                    throw TsuError("NaN values in input (position or velocity)");
                }
            }
        }

        reallocate_arrays(Npart);

        if (TsunamiConfig::wSpins) {
            if (w_aps) {
                std::cerr << "TsunamiWarning: code was compiled with spins, but no spins are provided in add_particle_set. "
                         "Defaulting spins to zero" << std::endl;
                w_aps = false;
            }
            for(size_t i = 0; i < N; i++) {
                System.spin[i] = {0.0, 0.0, 0.0};
            }
        }

        for(size_t i = 0; i < N; i++) {
            System.pos[i].x = pos_in[3*i];
            System.pos[i].y = pos_in[3*i+1];
            System.pos[i].z = pos_in[3*i+2];
            System.vel[i].x = vel_in[3*i];
            System.vel[i].y = vel_in[3*i+1];
            System.vel[i].z = vel_in[3*i+2];
            System.mass[i] = mass_in[i];
            System.radius[i] = rad_in[i];
            System.xdata[i].stype = static_cast<ptype>(stype_in[i]);
        }
        Nbodyalgorithms::scale_to_cdm(System.pos, System.vel, System.mass, N);

        initialize_regularization(Conf.alpha, Conf.beta, Conf.gamma);
        initialize_chain();
        initialize_integrator();
    }

    void add_particle_set(double *pos_in, size_t npos, size_t pncoord,
                          double *vel_in, size_t nvel, size_t vncoord,
                          double *mass_in, size_t nmass,
                          double *rad_in, size_t nrad,
                          long *stype_in, size_t nstype,
                          double *spin_in, size_t nspin, size_t sncoord) {

        if (not TsunamiConfig::wSpins) {
            if (w_aps) {
                std::cerr << "TsunamiWarning: code was not compiled with spins, ignoring the extra spin argument"
                          << std::endl;
                w_aps = false;
            }
            add_particle_set(pos_in, npos, pncoord, vel_in, nvel, vncoord, mass_in, nmass, rad_in, nrad, stype_in, nstype);
            return;
        }

        if((nmass != npos) or (nmass != nvel) or (nmass != nstype) or (nmass != nrad) or (nmass != nspin)) {
            throw TsuError("Provided arrays have inconsistent length");
        }
        if((pncoord != 3) or (vncoord != 3) or (sncoord != 3)) {
            throw TsuError("Provide 3 coordinates for each particle position and velocity");
        }

        size_t Npart = nmass;
        if(Npart < 2) {
            throw TsuError("Provide at least two particle");
        }

        for(size_t j = 0; j < Npart; j++) {
            if(std::isnan(mass_in[j]) or std::isnan(rad_in[j]) or std::isnan(stype_in[j])) {
                throw TsuError("NaN values in input (mass or radius or type)");
            }
            for(size_t i = 0; i < 3; i++) {
                if(std::isnan(pos_in[i+j*3]) or std::isnan(vel_in[i+j*3]) or std::isnan(spin_in[i+j*3])) {
                    throw TsuError("NaN values in input (position or velocity or spin)");
                }
            }
        }

        reallocate_arrays(Npart);

        for(size_t i = 0; i < N; i++) {
            System.pos[i].x = pos_in[3*i];
            System.pos[i].y = pos_in[3*i+1];
            System.pos[i].z = pos_in[3*i+2];
            System.vel[i].x = vel_in[3*i];
            System.vel[i].y = vel_in[3*i+1];
            System.vel[i].z = vel_in[3*i+2];
            System.spin[i] = {spin_in[3*i], spin_in[3*i+1], spin_in[3*i+2]};
            System.mass[i] = mass_in[i];
            System.radius[i] = rad_in[i];
            System.xdata[i].stype = static_cast<ptype>(stype_in[i]);
        }
        Nbodyalgorithms::scale_to_cdm(System.pos, System.vel, System.mass, N);

        initialize_regularization(Conf.alpha, Conf.beta, Conf.gamma);
        initialize_chain();
        initialize_integrator();
    }


    void sync_internal_state(double *pos_inout, size_t npos, size_t pncoord,
                             double *vel_inout, size_t nvel, size_t vncoord) {

        if((N != npos) or  (N != nvel)) {
            throw TsuError("Provided arrays have lengths different from particle number");
        }
        if((pncoord != 3) or (vncoord != 3)) {
            throw TsuError("Provide 3 coordinates for each particle");
        };

        for(size_t i = 0; i < N; i++) {
            pos_inout[3*i] = System.pos[i].x;
            pos_inout[3*i+1] = System.pos[i].y;
            pos_inout[3*i+2] = System.pos[i].z;
            vel_inout[3*i] = System.vel[i].x;
            vel_inout[3*i+1] = System.vel[i].y;
            vel_inout[3*i+2] = System.vel[i].z;
        }
    }

    void sync_internal_state(double *pos_inout, size_t npos, size_t pncoord,
                             double *vel_inout, size_t nvel, size_t vncoord,
                             double *spin_inout, size_t nspin, size_t sncoord) {
        if (not TsunamiConfig::wSpins) {
            if (w_sis) {
                std::cerr << "TsunamiWarning: code was not compiled with spins, ignoring the extra spin argument" << std::endl;
                w_sis = false;
            }
            sync_internal_state(pos_inout, npos, pncoord, vel_inout, nvel, vncoord);
            return;
        }

        if((N != npos) or  (N != nvel) or (N != nspin)) {
            throw TsuError("Provided arrays have lengths different from particle number");
        }
        if((pncoord != 3) or (vncoord != 3) or (sncoord != 3)) {
            throw TsuError("Provide 3 coordinates for each particle");
        };

        for(size_t i = 0; i < N; i++) {
            pos_inout[3*i] = System.pos[i].x;
            pos_inout[3*i+1] = System.pos[i].y;
            pos_inout[3*i+2] = System.pos[i].z;
            vel_inout[3*i] = System.vel[i].x;
            vel_inout[3*i+1] = System.vel[i].y;
            vel_inout[3*i+2] = System.vel[i].z;
            spin_inout[3*i] = System.spin[i].x;
            spin_inout[3*i+1] = System.spin[i].y;
            spin_inout[3*i+2] = System.spin[i].z;
        }
    }

    void override_masses(double *mass_in, size_t nmass) {
        if(N != nmass) {
            throw TsuError("Provided arrays have lengths different from particle number");
        }

        for(size_t i = 0; i < N; i++) {
            System.mass[i] = mass_in[i];
        }

        // Different from reset_integrator_sameN because we do not reinitialize the chain
        Nbodyalgorithms::scale_to_cdm(System.pos, System.vel, System.mass, N);
        energy = Nbodyalgorithms::energy_calculation(System.pos, System.vel, System.mass, N, pot, kin);
        System.update_chained_data(); // Update chained data
        BSExtra.reset_bulirsch(); 	// Reset BS and its internal leapfrog
    }

    /**
     * Overwrites the stored positions and velocities with the provided ones, leaving unchanged the number of particles
     * If changing also masses, *update masses first*
     * @param pos_in
     * @param npos
     * @param pncoord
     * @param vel_in
     * @param nvel
     * @param vncoord
     */
    void override_position_and_velocities(double *pos_in, size_t npos, size_t pncoord,
                                          double *vel_in, size_t nvel, size_t vncoord) {
        if((N != npos) or  (N != nvel)) {
            throw TsuError("Provided arrays have lengths different from particle number");
        }
        if((pncoord != 3) or (vncoord != 3)) {
            throw TsuError("Provide 3 coordinates for each particle");
        };

        for(size_t i = 0; i < N; i++) {
            System.pos[i].x = pos_in[3*i];
            System.pos[i].y = pos_in[3*i+1];
            System.pos[i].z = pos_in[3*i+2];
            System.vel[i].x = vel_in[3*i];
            System.vel[i].y = vel_in[3*i+1];
            System.vel[i].z = vel_in[3*i+2];
        }

        Nbodyalgorithms::scale_to_cdm(System.pos, System.vel, System.mass, N);
        energy = Nbodyalgorithms::energy_calculation(System.pos, System.vel, System.mass, N, pot, kin);
        System.find_chain();
        System.initialize_chain_arrays(); // Redo chain
        BSExtra.reset_bulirsch(); 	// Reset BS and its internal leapfrog
    }

    /**
     * Initializes tidal parameters. Time lag is given in N-body units
     * @param kaps_in
     * @param nkaps
     * @param taulag_in
     * @param ntaulag
     * @param polyt_in
     * @param npolyt
     */
    void initialize_tidal_parameters(double *kaps_in, size_t nkaps,
                                     double *taulag_in, size_t ntaulag,
                                     double *polyt_in, size_t npolyt) {
        if (TsunamiConfig::wSpins) {
            throw TsuError("Code was compiled with spins\nadd gyration radii as last argument");
        }

        if((N != nkaps) or (N != ntaulag) or (N != npolyt)) {
            throw TsuError("Provided arrays have lengths different from particle number");
        }

        for(size_t i=0; i < N; i++) {
            if((kaps_in[i] <= 0) or (taulag_in[i] <= 0) or (polyt_in[i] <= 0)) {
                System.xdata[i].hastide = false;
                //std::cerr<<"Particle "<<i<<" has a tidal parameter less or equal than zero: disabling tides on this particle"<<std::endl;
            } else {
                System.xdata[i].hastide = true;
                System.xdata[i].polyt = polyt_in[i];
                System.xdata[i].kaps = kaps_in[i];
                System.xdata[i].taulag = taulag_in[i];

                double Qt = 2 * System.xdata[i].kaps / (1 + 2 * System.xdata[i].kaps);
                double R5 = pow(System.radius[i], 5);

                // True parameters used
                System.xdata[i].Atide = R5 * Qt / (1 - Qt);
                System.xdata[i].sigmadiss = 2.0 * R5 / (System.xdata[i].Atide*System.xdata[i].Atide) * System.xdata[i].kaps * System.xdata[i].taulag / 3;
            }
        }

        System.update_chained_data();
    }

    /**
     * Initializes tidal parameters. Time lag is given in N-body units
     * @param kaps_in
     * @param nkaps
     * @param taulag_in
     * @param ntaulag
     * @param polyt_in
     * @param npolyt
     */
    void initialize_tidal_parameters(double *kaps_in, size_t nkaps,
                                     double *taulag_in, size_t ntaulag,
                                     double *polyt_in, size_t npolyt,
                                     double *gyrad_in, size_t ngyrad) {

        if((N != nkaps) or (N != ntaulag) or (N != npolyt) or (N != ngyrad)) {
            throw TsuError("Provided arrays have lengths different from particle number");
        }

        for(size_t i=0; i < N; i++) {
            if(((kaps_in[i] <= 0) or (taulag_in[i] <= 0)) and (polyt_in[i] <= 0)) {
                System.xdata[i].hastide = false;
                //std::cerr<<"Particle "<<i<<" has a tidal parameter less or equal than zero: disabling tides on this particle"<<std::endl;
            } else {
                System.xdata[i].hastide = true;
                System.xdata[i].polyt = polyt_in[i];
                System.xdata[i].kaps = kaps_in[i];
                System.xdata[i].taulag = taulag_in[i];

                double gyrad = gyrad_in[i] * System.radius[i];
                System.xdata[i].inert = gyrad*gyrad * System.mass[i];

                double Qt = 2 * System.xdata[i].kaps / (1 + 2 * System.xdata[i].kaps);
                double R5 = pow(System.radius[i], 5);

                // True parameters used
                System.xdata[i].Atide = R5 * Qt / (1 - Qt);
                System.xdata[i].sigmadiss = 2.0 * R5 / (System.xdata[i].Atide*System.xdata[i].Atide) * System.xdata[i].kaps * System.xdata[i].taulag / 3;
            }
        }

        System.update_chained_data();
    }


/**
 * Evolves the system to the given final time
 * @param[in] tfin
 */
    void evolve_system(double tfin) {
        // All sort of checks
        if(tfin <= time)
            throw TsuError("Trying to evolve a system that has already reached tfin");
        if (System.CollInfo.collflag)
            throw TsuError("Trying to evolve a system that has already collided\nsolve the collision and resume evolution");

        while (not stopcond) { // Start loop

            // Advance system
            try {
                do_step_bulirsh();
            } catch (TsunamiError &error) {
                std::time_t t = std::time(nullptr);
                std::tm tm = *std::localtime(&t);
                std::ostringstream oss;
                oss << std::put_time(&tm, "%d-%m-%Y-%H-%M-%S");
                string fname = "emergency_crash_" + oss.str() + ".bin";
                BSExtra.dt = BSExtra.dt_prev; // If dt screwed up, save the previous one
                save_restart_file(fname);
                throw TsuError("Saved " + fname + ", now ship the file to Alessandro");
            }

#ifdef STOPCOND
            stopcond = StopLogger->check_stopping_condition(pos, vel, mass, rad, ctime);
#endif
            System.update_from_chain_to_com();
            System.find_chain();

            if (check_collision()) {
                if (iterate_to_collision) {
                    bool truepos = iterate_to_collision_bulirsch();
                    if (truepos) break;
                } else {
                    Leapfrog.revert_step(System.ch_pos, System.ch_vel, System.ch_mass,
                                         System.ch_spin, System.ch_xdata);
                    System.revert_chain();
                    System.update_from_chain_to_com();
                    break;
                }
            }

            update_time_coordinates_chain();

            if(time >= tfin) { // If reached target time, break out
                energy = Nbodyalgorithms::energy_calculation(System.pos, System.vel,
                                                             System.mass, N,
                                                             pot, kin);
                deltaE = fabs(log(fabs((kin + (Eintegrator))/(pot))));
                break;
            }
        } // End loop
    }


/**
 * Evolves the system for a single step, making sure it does not integrate beyond tfin
 * @param[in] tfin
 */
    void evolve_system_dtmax(double tfin) {
        // All sort of checks
        if(tfin <= time)
            throw TsuError("Trying to evolve a system that has already reached tfin");
        if (System.CollInfo.collflag)
            throw TsuError("Trying to evolve a system that has already collided\nsolve the collision and resume evolution");

        double dtime = tfin - time;
        double timestep_dt = estimate_ds_from_dt(dtime);
        timestep = (timestep_dt > timestep) ? timestep : timestep_dt;

        do_step_bulirsh();

#ifdef STOPCOND
            stopcond = StopLogger->check_stopping_condition(pos, vel, mass, rad, ctime);
#endif
        System.update_from_chain_to_com();
        System.find_chain();

        update_time_coordinates_chain();


        energy = Nbodyalgorithms::energy_calculation(System.pos, System.vel,
                                                         System.mass, N,
                                                         pot, kin);
        deltaE = fabs(log(fabs((kin + (Eintegrator))/(pot))));
    }

    void commit_parameters() {
        Leapfrog.update_parameters();
    }

    void reallocate_arrays(size_t Npart) {
        if (N == 0 or N != Npart) {
            N = Npart;
            if (System.Npart > 0) System.deallocate_arrays();
            System.Npart = Npart;
            System.Nchain = Npart - 1;
            System.allocate_arrays();
            System.chain.resize(Npart);

            aux_vec1.resize(N);
            aux_vec2.resize(N);
        }
    }

    void revert_step() {
        Leapfrog.revert_step(System.ch_pos, System.ch_vel, System.ch_mass, System.ch_spin, System.ch_xdata);
    }

    void get_chain_vectors(double **chpos_out, int *npos, int *pncoord,
                           double **chvel_out, int *nvel, int *vncoord) {

        if (System.Npart < 2) {
            throw TsuError("No particles found");
        }

        *vncoord = *pncoord = 3;
        *npos = *nvel = System.Nchain;
        auto *tmp_pos = new double [3*System.Nchain];
        auto *tmp_vel = new double [3*System.Nchain];

        for (size_t i=0; i<System.Nchain; i++) {
            tmp_pos[i*System.Nchain + 0] = System.ch_pos[i].x;
            tmp_pos[i*System.Nchain + 1] = System.ch_pos[i].y;
            tmp_pos[i*System.Nchain + 2] = System.ch_pos[i].z;
            tmp_vel[i*System.Nchain + 0] = System.ch_vel[i].x;
            tmp_vel[i*System.Nchain + 1] = System.ch_vel[i].y;
            tmp_vel[i*System.Nchain + 2] = System.ch_vel[i].z;
        }
        *chpos_out = tmp_pos;
        *chvel_out = tmp_vel;
    }
    /**
     * Acceleration of particle i (acc_i_out) from particle j and vice versa (acc_j_out)
     * @param i
     * @param j
     * @param acc_i_out
     * @param acc_j_out
     */
    void get_accelerations_of_particle_pair(size_t i, size_t j,
                                            double acc_i_out[3], double acc_j_out[3]) {

        if ((i >= System.Npart) or (j >= System.Npart)) {
            throw TsuError("Provided particle index does not exist");
        }

        size_t ch_i = System.invchain[i];
        size_t ch_j = System.invchain[j];
        // k chain member connects k and k+1, get the lowest
        make_ascending(ch_i, ch_j);

        auto unchained_acc = aux_vec1.data();
        for (size_t ip = 0; ip < ch_j+1; ip++) {
            unchained_acc[ip] = 0;
        }

        double dummy;
        auto dummy3 = aux_vec2.data();

        double3 dpos = System.ch_pos[ch_i];
        for (size_t ip = ch_i+1; ip < ch_j; ip++) {
            dpos += System.ch_pos[ip];
        }
        System.cache_ind = 0;
        System.twobody_interaction(ch_j, ch_i, dpos,
                                   System.ch_mass, System.ch_radius,
                                   unchained_acc, dummy3, dummy,dummy);

        if (System.Config.wPNs or System.Config.wEqTides or System.Config.wDynTides) {
            double3 dvel = System.ch_vel[ch_i];
            for (size_t ip = ch_i+1; ip < ch_j; ip++) {
                dvel += System.ch_vel[ip];
            }
            System.cache_ind = 0;
            System.twobody_interaction_vdep<false>(ch_j, ch_i, dpos, dvel, System.ch_mass,
                                                System.ch_radius, System.ch_xdata,
                                                unchained_acc, System.spin, dummy3);
        }

        acc_i_out[0] = unchained_acc[System.invchain[i]].x;
        acc_i_out[1] = unchained_acc[System.invchain[i]].y;
        acc_i_out[2] = unchained_acc[System.invchain[i]].z;
        acc_j_out[0] = unchained_acc[System.invchain[j]].x;
        acc_j_out[1] = unchained_acc[System.invchain[j]].y;
        acc_j_out[2] = unchained_acc[System.invchain[j]].z;
    }

    void gravitational_waveform(double &hplus, double &hcross, size_t ip, size_t jp, double theta, double phi, double R) {
        double a1[3], a2[3];
        get_accelerations_of_particle_pair(ip, jp, a1, a2);

        // Now for every other particle
        for (size_t kp=0; kp < N; kp++) {
            if ((kp != jp) and (kp != ip)) {
                double acc_i_temp[3], acc_j_temp[3], dummy[3];
                get_accelerations_of_particle_pair(ip, kp, acc_i_temp, dummy);
                get_accelerations_of_particle_pair(jp, kp, acc_j_temp, dummy);
                for (size_t idim=0; idim<3; idim++) {
                    a1[idim] += acc_i_temp[idim];
                    a2[idim] += acc_j_temp[idim];
                }
            }
        }

        double v1[3], p1[3], v2[3], p2[3];
        System.pos[ip].to_array(p1);
        System.vel[ip].to_array(v1);
        System.pos[jp].to_array(p2);
        System.vel[jp].to_array(v2);
        double &m1 = System.mass[ip];
        double &m2 = System.mass[jp];
        double pcom[3], vcom[3];
        for (size_t k=0; k<3; k++) {
            pcom[k] = (p1[k]*m1 + p2[k]*m2) / (m1 + m2);
            vcom[k] = (v1[k]*m1 + v2[k]*m2)  / (m1 + m2);
        }
        for (size_t k=0; k<3; k++) {
            p1[k] -= pcom[k];
            p2[k] -= pcom[k];
            v1[k] -= vcom[k];
            v2[k] -= vcom[k];
        }

        std::array<std::array<double, 3>, 3> htensor;
        for (size_t ii=0; ii<3; ii++) {
            for (size_t jj=0; jj<3; jj++) {
                double term1 = m1 * (a1[ii] * p1[jj] + p1[ii] * a1[jj] + 2 * v1[ii] * v1[jj]);
                double term2 = m2 * (a2[ii] * p2[jj] + p2[ii] * a2[jj] + 2 * v2[ii] * v2[jj]);
                double diag = (ii == jj) ? -2 * ( m1 * (vdot(a1, p1) + vdot(v1, v1)) +
                                                  m2 * (vdot(a2, p2) + vdot(v2, v2))) / 3 : 0.0;

                htensor[ii][jj] = (2 * System.coeff2pn / R) * (term1 + term2 + diag);
            }
        }

        double n[3] = {sin(theta) * sin(phi), sin(theta) * cos(phi),  cos(theta)};
        double temp[3] = {0.0, 1.0, 0.0};

        // If it's by chance parallel, let's change it
        if (std::fabs(vdot(n, temp)) > 0.9) {
            temp[0] = 1.0;
            temp[1] = 0.0;
        }

        // p, q, n as an orthonormal basis
        double p[3], q[3];
        vcross(temp, n, p);
        vnormalize(p);
        vcross(n, p, q);
        vnormalize(q);

        std::array<std::array<double, 3>, 3> e_plus;
        std::array<std::array<double, 3>, 3> e_cross;
        for (size_t ii=0; ii<3; ii++) {
            for (size_t jj=0; jj<3; jj++) {
                e_plus[ii][jj] = p[ii] * p[jj] - q[ii] * q[jj];
                e_cross[ii][jj] = p[ii] * q[jj] + q[ii] * p[jj];
            }
        }

        hplus = hcross = 0.0;
        for (size_t ii=0; ii<3; ii++) {
            for (size_t jj=0; jj<3; jj++) {
                hplus += htensor[ii][jj] * e_plus[ii][jj];
                hcross += htensor[ii][jj] * e_cross[ii][jj];
            }
        }
        hplus *= 0.5;
        hcross *= 0.5;
    }

    void sync_masses(double *mass_out, size_t nmass) const {
        if (System.Npart < 2) {
            throw TsuError("No particles found");
        }
        if(System.Npart != nmass) {
            throw TsuError("Provided arrays have lengths different from particle number");
        }

        for (size_t i=0; i<N; i++) {
            mass_out[i] = System.mass[i];
        }
    }

    void sync_radii(double *rad_out, size_t nrad) const {
        if (System.Npart < 2) {
            throw TsuError("No particles found");
        }

        if(System.Npart != nrad) {
            throw TsuError("Provided arrays have lengths different from particle number");
        }

        for (size_t i=0; i<System.Npart; i++) {
            rad_out[i] = System.radius[i];
        }
    }

    void sync_eloss(double *eloss_out, size_t neloss) const {
        if (System.Npart < 2) {
            throw TsuError("No particles found");
        }

        if(System.Npart != neloss) {
            throw TsuError("Provided arrays have lengths different from particle number");
        }

        for (size_t i=0; i<System.Npart; i++) {
            eloss_out[i] = System.xdata[i].eloss;
        }
    }

    double estimate_ds_from_dt(const double dt) {
        double dt_u = Leapfrog.calc_dt_pos(1.0, Leapfrog.T + Leapfrog.B, Leapfrog.omega);
        return dt/dt_u;
    }

    Nbodyalgorithms::PTDELog & get_PTDELog() {
        if constexpr (TsunamiConfig::useTDEtracker) {
        } else {
            std::cerr << "TsunamiWarning: code was not compiled with tdetracker, but user is requesting PDTELog\n\trecompile with -Dtdetracker=on" << std::endl;
        }
        return System.PTDEInfo;
    }

    bool check_collision() const {
        return System.CollInfo.collflag;
    }

    /*[[nodiscard]] std::tuple<size_t,size_t> get_collision_indices() const {
        return {System.CollInfo.collind[0], System.CollInfo.collind[1]};
    }*/// Maybe when SWIG will support tuples

    void get_collision_indices(size_t &id1, size_t &id2) const {
        const size_t &ic1 = System.CollInfo.collind[0];
        const size_t &ic2 = System.CollInfo.collind[1];

        id1 = ic1 < ic2 ? ic1 : ic2;
        id2 = ic2 > ic1 ? ic2 : ic1;
    }

    double get_collision_time() const {
        return System.CollInfo.colltime + System.time0;
    }

    void reset_collision_status() {
        System.CollInfo.reset_collision();
    }

    void save_restart_file(const std::string& fname) {
        if (N < 1) throw TsuError("No particles found");

        std::ofstream outfile(fname, std::ios::binary | std::ios::out);
        if (!outfile) throw TsuError("Cannot open file " + fname);

        outfile.write(reinterpret_cast<const char *>(&N), sizeof(size_t));
        outfile.write(reinterpret_cast<const char *>(&timestep), sizeof(double));
        outfile.write(reinterpret_cast<const char *>(&dtphysical), sizeof(double));  // Not necessary but useful info
        outfile.write(reinterpret_cast<char *>(&Conf.wPNs), sizeof(bool));
        outfile.write(reinterpret_cast<char *>(&Conf.wEqTides), sizeof(bool));
        outfile.write(reinterpret_cast<char *>(&Conf.wDynTides), sizeof(bool));
        outfile.write(reinterpret_cast<char *>(&Conf.wExt), sizeof(bool));
        outfile.write(reinterpret_cast<char *>(&Conf.wExt_vdep), sizeof(bool));
        outfile.write(reinterpret_cast<char *>(&Conf.wMassEvol), sizeof(bool));
        outfile.write(reinterpret_cast<char *>(&Conf.Mscale), sizeof(double));
        outfile.write(reinterpret_cast<char *>(&Conf.Lscale), sizeof(double));
        outfile.write(reinterpret_cast<char *>(&Conf.alpha), sizeof(double));
        outfile.write(reinterpret_cast<char *>(&Conf.beta), sizeof(double));
        outfile.write(reinterpret_cast<char *>(&Conf.gamma), sizeof(double));
        outfile.write(reinterpret_cast<char *>(&Conf.Usafe), sizeof(double));
        outfile.write(reinterpret_cast<char *>(&Conf.dcoll), sizeof(double));
        outfile.write(reinterpret_cast<const char *>(&tolerance), sizeof(double));
        outfile.write(reinterpret_cast<const char *>(&BSExtra.kopt_now), sizeof(size_t));
        outfile.write(reinterpret_cast<const char *>(&Leapfrog.B0), sizeof(double));
        outfile.write(reinterpret_cast<const char *>(&Leapfrog.omega0), sizeof(double));
        outfile.write(reinterpret_cast<const char *>(&deltaE), sizeof(double));
        outfile.write(reinterpret_cast<const char *>(&energy), sizeof(double));
        outfile.write(reinterpret_cast<const char *>(&time), sizeof(double));
        outfile.write(reinterpret_cast<const char *>(System.chain.data()), System.Npart*sizeof(size_t));
        outfile.write(reinterpret_cast<const char *>(System.ch_pos), System.Nchain*sizeof(double3));
        outfile.write(reinterpret_cast<const char *>(System.ch_vel), System.Nchain*sizeof(double3));
        outfile.write(reinterpret_cast<const char *>(System.ch_mass), System.Npart*sizeof(double));
        outfile.write(reinterpret_cast<const char *>(System.ch_radius), System.Npart*sizeof(double));
        outfile.write(reinterpret_cast<const char *>(System.ch_xdata), System.Npart*sizeof(pinfo));

        std::vector<double3> zerospin; // If spin is not included
        double3 *spin_data_pointer;
        if (TsunamiConfig::wSpins) {
            spin_data_pointer = System.ch_spin;
        } else {
            for (size_t i=0; i<N; i++) {
                zerospin.emplace_back(0.0, 0.0, 0.0);
            }
            spin_data_pointer = zerospin.data();
        }
        outfile.write(reinterpret_cast<const char *>(spin_data_pointer), System.Npart*sizeof(double3));

        outfile.close();
        if(!outfile.good()) throw TsuError("Error occurred during writing of " + fname);


        /*std::cout << N << std::endl;
        System.print_chain();
        for(size_t i=0; i<System.Npart; i++) std::cout << System.mass[i] << " ";
        std::cout << std::endl;
        for(size_t i=0; i<System.Npart; i++) std::cout << System.radius[i] << " ";
        std::cout << std::endl;
        for(size_t i=0; i<System.Npart; i++) std::cout << System.xdata[i] << std::endl;
        std::cout << std::endl;
        for(size_t i=0; i<System.Npart; i++) std::cout << System.pos[i] << " ";
        std::cout << std::endl;
        for(size_t i=0; i<System.Npart; i++) std::cout << System.vel[i] << " ";
        std::cout << std::endl;
        for(size_t i=0; i<System.Npart; i++) std::cout << System.ch_pos[i] << " ";
        std::cout << std::endl;
        for(size_t i=0; i<System.Npart; i++) std::cout << System.ch_vel[i] << " ";
        std::cout << std::endl;
        std::cout << "energy " << energy << std::endl;*/
    }


    void load_restart_file(const std::string& fname) {
        std::ifstream infile(fname, std::ios::binary | std::ios::in);
        if (!infile) throw TsuError("Cannot open file " + fname);


        size_t Npart;
        infile.read(reinterpret_cast<char *>(&Npart), sizeof(size_t));
        if(Npart < 2) throw TsuError("Cannot initialize integrator with less than two particles");

        reallocate_arrays(Npart);

        double B0, omega0;
        size_t k_opt;
        infile.read(reinterpret_cast<char *>(&timestep), sizeof(double));
        infile.read(reinterpret_cast<char *>(&dtphysical), sizeof(double));  // Not necessary but useful info
        infile.read(reinterpret_cast<char *>(&Conf.wPNs), sizeof(bool));
        infile.read(reinterpret_cast<char *>(&Conf.wEqTides), sizeof(bool));
        infile.read(reinterpret_cast<char *>(&Conf.wDynTides), sizeof(bool));
        infile.read(reinterpret_cast<char *>(&Conf.wExt), sizeof(bool));
        infile.read(reinterpret_cast<char *>(&Conf.wExt_vdep), sizeof(bool));
        infile.read(reinterpret_cast<char *>(&Conf.wMassEvol), sizeof(bool));
        infile.read(reinterpret_cast<char *>(&Conf.Mscale), sizeof(double));
        infile.read(reinterpret_cast<char *>(&Conf.Lscale), sizeof(double));
        infile.read(reinterpret_cast<char *>(&Conf.alpha), sizeof(double));
        infile.read(reinterpret_cast<char *>(&Conf.beta), sizeof(double));
        infile.read(reinterpret_cast<char *>(&Conf.gamma), sizeof(double));
        infile.read(reinterpret_cast<char *>(&Conf.Usafe), sizeof(double));
        infile.read(reinterpret_cast<char *>(&Conf.dcoll), sizeof(double));
        infile.read(reinterpret_cast<char *>(&tolerance), sizeof(double));
        infile.read(reinterpret_cast<char *>(&k_opt), sizeof(size_t));
        infile.read(reinterpret_cast<char *>(&B0), sizeof(double));
        infile.read(reinterpret_cast<char *>(&omega0), sizeof(double));
        infile.read(reinterpret_cast<char *>(&deltaE), sizeof(double));
        infile.read(reinterpret_cast<char *>(&energy), sizeof(double));
        infile.read(reinterpret_cast<char *>(&time), sizeof(double));
        infile.read(reinterpret_cast<char *>(System.chain.data()), System.Npart*sizeof(size_t));
        infile.read(reinterpret_cast<char *>(System.ch_pos), System.Nchain*sizeof(double3));
        infile.read(reinterpret_cast<char *>(System.ch_vel), System.Nchain*sizeof(double3));
        infile.read(reinterpret_cast<char *>(System.ch_mass), System.Npart*sizeof(double));
        infile.read(reinterpret_cast<char *>(System.ch_radius), System.Npart*sizeof(double));
        infile.read(reinterpret_cast<char *>(System.ch_xdata), System.Npart*sizeof(pinfo));

        if (TsunamiConfig::wSpins) {
            infile.read(reinterpret_cast<char *>(System.ch_spin), System.Npart*sizeof(double3));
        } else {
            std::vector<double3> zerospin(System.Npart); // If spin is not included
            infile.read(reinterpret_cast<char *>(zerospin.data()), System.Npart*sizeof(double3));
        }

        infile.close();
        if(!infile.good()) throw TsuError("Error occurred during the reading of " + fname);

        std::cout << Npart << std::endl;
        System.print_chain();

        set_units(Conf.Mscale, Conf.Lscale);

        System.copy_from_chain_to_com();
        System.init_chain();

        initialize_regularization(Conf.alpha, Conf.beta, Conf.gamma);
        initialize_integrator();

        /*for(size_t i=0; i<System.Npart; i++) std::cout << System.mass[i] << " ";
        std::cout << std::endl;
        for(size_t i=0; i<System.Npart; i++) std::cout << System.radius[i] << " ";
        std::cout << std::endl;
        for(size_t i=0; i<System.Npart; i++) std::cout << System.xdata[i] << std::endl;
        std::cout << std::endl;
        for(size_t i=0; i<System.Npart; i++) std::cout << System.pos[i] << " ";
        std::cout << std::endl;
        for(size_t i=0; i<System.Npart; i++) std::cout << System.vel[i] << " ";
        std::cout << std::endl;
        for(size_t i=0; i<System.Npart; i++) std::cout << System.ch_pos[i] << " ";
        std::cout << std::endl;
        for(size_t i=0; i<System.Npart; i++) std::cout << System.ch_vel[i] << " ";
        std::cout << std::endl;
        energy = Nbodyalgorithms::energy_calculation(System.pos, System.vel, System.mass, N, pot, kin);
        std::cout << "energy " << energy << std::endl;*/

        BSExtra.kopt_now = k_opt;
        BSExtra.FinalTable.B = Leapfrog.B0 = B0;
        BSExtra.FinalTable.omega = Leapfrog.omega0 = omega0;
    }


    std::vector<std::pair<double, double>> switch_list;

    /// Number of particles
    size_t N = 0;

    /// Paths and default names
    // Current path parameter file
    string exepath = IO::get_exe_path();
    string paramfile_locpath = "input/tsunami_parameters.txt";
    string fname_locpath = "input/tsunami_default_input.dat";
    string tidetile_locpath = "input/tsunami_tide_table.dat";

    string paramfile = exepath + paramfile_locpath;
    string fname_in = exepath + fname_locpath;
    string tidefile = exepath + tidetile_locpath;
    string out_name = "output.dat";
    string ene_name = "energy.dat";
    string coll_name = "collision.dat";

    /// Parameters
    double timestep = 1.e-13;
    double tolerance = 1.0e-13;
    TsunamiConfig Conf;

    /// Sub-classes
    ChainSys System;
    Classification Marx;
    LeapfrogStepped Leapfrog;
    BSExtrap BSExtra;
    size_t nstep = 2;

    // Useful unit scales
    double Mscale; ///< Mass scale (in Msun)
    double Lscale; ///< Length scale (in au)
    double Tscale; ///< Time scale (in yr)
    double Vscale; ///< Velocity scale (in km/s)
    double speed_of_light;

    double time = 0.0;
    double energy = 0.0;
    double pot = 0.0;
    double kin = 0.0;
    double eoff = 0.0;

    bool stopcond = false;
    bool iterate_to_collision = false;

    double deltaE = 0.0;
    double Eintegrator = 0.0;
    double dtphysical = 0.0; // Physical timestep

    // Profiling utilities
    SimProf ProfBS = SimProf();
    SimProf ProfChain = SimProf();

    // Warnings
    bool w_sis = true;
    bool w_aps = true;

private:
    std::vector<double3> aux_vec1;
    std::vector<double3> aux_vec2;

};

typedef TsunamiClass<TsunamiConfig::useProfiling, (TsunamiConfig::debug_bs or TsunamiConfig::debug_lf or TsunamiConfig::debug_ch)> TsunamiCode;

#endif //TSUNAMI_TSUNAMI_HPP
