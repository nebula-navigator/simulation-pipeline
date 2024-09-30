#ifndef BSINT_H
#define BSINT_H

#include <vector>
#include <cmath>
#include <algorithm>
#include <ctime>
#include "errhand.hpp"
#include "custom_types.hpp"
#include "leapfrog_stepped.hpp"

template<bool withspin, bool debug>
class BSExtrapolator {
    template<bool profiling, bool debug_tsu> friend class TsunamiClass;

public:
    //time_step here should be are_arrays_allocated outside at a very small value (1.0e-13)
    //then it is substituted with the next time_step obtained during the BS iterations

    BSExtrapolator(ChainSys &System,
                   LeapfrogStepped &LeapStepped,
                   double tolerance) : LeapStepped(LeapStepped), System(System) {

        //BS parameters
        set_parameters(tolerance, tolerance);

        // Do not involve the number of particles
        first_step_init();
        init_steps_seq();
        init_matrices();
    }

    ~BSExtrapolator() = default;

    void set_parameters(double abs_err, double rel_err) {
        abs_err_tolerance = abs_err;
        rel_err_tolerance = rel_err;
    }

    void initialize() {
        allocate_arrays();
    }

    void bs_iterate(double &dtphysical, double &timestep) {

        LeapStepped.save_step_zero(System.ch_pos, System.ch_vel, System.ch_mass,
                                   System.ch_spin, System.ch_xdata);

        do { // Keep iterating until a successful BS iteration
            rejection = bs_step(timestep);
            timestep = dt_next;

        } while (rejection);

        // We are here, step is accepted. Update variables
        dtphysical = FinalTable.time;

        // Load output arrays with BS integrated values
        for (size_t i = 0; i < LeapStepped.Neqs; i++) {
            System.ch_pos[i] = FinalTable.pos[i];
            System.ch_vel[i] = FinalTable.vel[i];

            if (LeapStepped.vdep) {
                System.ch_xdata[i].eloss = FinalTable.eloss[i];
                if constexpr (withspin) System.ch_spin[i] = FinalTable.spin[i];
            }
            if (System.Config.wMassEvol) {
                //TODO ADD MASS EVOLUTION HERE
            }
        }

        // Load also eloss and spin, if enabled
        // Instead of another loop, we just check if we included the CoM and Nchain = Npart
        if (not System.Config.wExt) {
            if (LeapStepped.vdep) {
                System.ch_xdata[System.Nchain].eloss = FinalTable.eloss[System.Nchain];
                if constexpr (withspin) System.ch_spin[System.Nchain] = FinalTable.spin[System.Nchain];
            }
            if (System.Config.wMassEvol) {
                //TODO ADD MASS EVOLUTION HERE
            }
        }

        // Update also the numerical evolution of special functions
        // Moved to tsunami.hpp
        // LeapStepped.B0 = FinalTable.B;
        // LeapStepped.omega0 = FinalTable.omega;
    }

    /**
     * Extrapolates over a timestep
     * The variables integrated are:
     *     - System.ch_pos
     *     - System.ch_vel
     *     - System.ch_spin (if withspin)
     *     - System.xdata.eloss
     *     - B
     *     - omega
     *     - dtphysical
     * @param timestep
     * @return
     */
    bool bs_step(const double timestep) {

        dt = timestep;
        rejection = false;

        //std::cout << "\n STEP dt = " << dt << std::endl;
        //if (dt < 1e-16) exit(EXIT_FAILURE);

        std::vector<double> dt_ideal(kmax + 1);
        std::vector<double> work_min(kmax + 1);

        double dt_new = dt;
        size_t new_kopt = kopt_now;

        if constexpr (debug) {
            std::cout << "\n###### BS time: " << System.time0 << " BS ds = " << dt << " last rej: " << last_step_rejected << std::endl;
            std::cout << "BS k = " << 0 << " (" << steps_sequence[0] << ")" << std::endl;
            //std::cout << std::isnan(dt) << std::endl;
            if (std::isnan(dt)) {
                std::cout << "\ndt is nan, TERMINATING" << std::endl;
                exit(EXIT_FAILURE);
            }
        }

        if (std::isnan(dt) or dt <= 0.0) {
            LeapStepped.revert_step(System.ch_pos, System.ch_vel, System.ch_mass,
                                    System.ch_spin, System.ch_xdata);
            throw TsuError("SOMETHING BAD HAPPENED HERE, dt is " + std::to_string(dt));
        }
        
        // First BS step (k = 0, nstep = 2). No need to extrapolate
        System.lastit = false;
        try {
            LeapStepped.integrate(steps_sequence[0], dt, Dtime,
                                  System.ch_pos, System.ch_vel, System.ch_spin);
        } catch (IntegrationError &IE) {
            mitigate_negative_step();
            return true;
        }

        store_integrated_and_reset(FinalTable); // Store extrapolated variables and reset

        if constexpr (debug) {
            std::cout << " time: " << FinalTable.time << std::endl;
        }
        // Higher order BS steps (k = 1 .. kopt_now), here extrapolate
        for (size_t k = 1; k <= kopt_now + 1; k++) {

            System.lastit = (k == kopt_now);
            try {
                LeapStepped.integrate(steps_sequence[k], dt, Dtime,
                                      System.ch_pos, System.ch_vel, System.ch_spin);
            } catch (IntegrationError &IE) {
                mitigate_negative_step();
                return true;
            }
            if constexpr (debug) {
                std::cout << "BS k = " << k << " (" << steps_sequence[k] << ")" << std::endl;
            }

            store_integrated_and_reset(ExtrapTables[k - 1]); // Store integrated variables and reset

            if constexpr (debug) {
                std::cout << " time: " << ExtrapTables[k - 1].time << std::endl;
            }

            // Try to extrapolate pos, vel, spin, B, omega and time
            // ExtrapTables is loaded with the BS current estimation
            // FinalTable contains the previous extrapolated value!
            extrapolate_table(k, ExtrapTables, FinalTable, coeff_extr);

            //estrapolate_rational(k, table, coeff_ratextr, pv_ch_final);

            //extrapolate_scalar(k, tableB, coeff_extr, &Bfinal);
            //extrapolate_scalar(k, tableT, coeff_extr, &Tfinal);

            //now table and pv_ch_final have been both updated
            //details on http://numerical.recipes/webnotes/nr3web21.pdf

            //error = estimate_error(pv_ch_in, pv_ch_final, table, LeapStepped.B0, Bfinal,
            //                       tableB); //mean squared error over all positions and velocities
            error = estimate_error(FinalTable, ExtrapTables[0]);
            if constexpr (debug) {
                std::cout << " error = " << error << std::endl;
            }


            dt_ideal[k] = estimate_ideal_dt(dt, error, k);
            work_min[k] = cost[k] / dt_ideal[k];

            //follow the conditions/checks of Numerical Recipes pag. 925, from eq. (17.3.14) on
            if (k == kopt_now - 1) { //check convergence at kopt - 1 (before koptimal)

                if (error <= 1.0) {//CONVERGENCE!! Thus, accept the step.

                    rejection = false;

                    if ((work_min[k] < param_K2 * work_min[k - 1]) || (kopt_now <= 2)) {
                        new_kopt = std::min<size_t>(kmax - 1,
                                                    std::max<size_t>(2, k + 1));
                        //new_kopt = k+1, that is use the same kopt; still, it must be less than kmax-1 and greater than 2
                        //so, if kopt_now <= 2 the new k_opt will be, at least, 2;
                        //in this case it is not safe to reduce the order of the BS integrator using kopt_new = kopt - 1
                        dt_new = dt_ideal[new_kopt - 1];
                        dt_new *= cost[new_kopt] / static_cast<double>(cost[new_kopt - 1]);
                    } else if (work_min[k - 1] < param_K1 * work_min[k] && k >= 2) {
                        new_kopt = std::min<size_t>(kmax - 1,
                                                    std::max<size_t>(2, k - 1));
                        // new_kopt = k-1: decrease further the order of the BS integration... high order is not needed
                        dt_new = dt_ideal[new_kopt];
                    } else {
                        new_kopt = std::min<size_t>(kmax - 1,
                                                    std::max<size_t>(2, k));
                        dt_new = dt_ideal[new_kopt];
                    }

                    break; // Stop the BS iterations... I have reached convergence
                }

                    //Criterion NR eq. (17.3.17)
                else if (error > mincoeff_for_kopt_convergence(k, steps_sequence)) {
                    //this conditions tell us that the kopt+1 step will fail, so reject the BS iteration (see NR pag. 926)
                    rejection = true;

                    if ((work_min[k] < param_K2 * work_min[k - 1]) || (kopt_now <= 2)) {
                        new_kopt = std::min<size_t>(kmax - 1,
                                                    std::max<size_t>(2, k + 1));
                        //new_kopt = k+1, that is use the same kopt; still, it must be less than kmax-1 and greater than 2
                        //so, if kopt_now <= 2 the new k_opt will be, at least, 2;
                        //in this case it is not safe to reduce the order of the BS integrator using kopt_new = kopt - 1
                        dt_new = dt_ideal[new_kopt - 1];
                        dt_new *= cost[new_kopt] / static_cast<double>(cost[new_kopt - 1]);
                    } else if (work_min[k - 1] < param_K1 * work_min[k] && k >= 2) {
                        new_kopt = std::min<size_t>(kmax - 1,
                                                    std::max<size_t>(2, k - 1));
                        //new_kopt = k-1: decrease further the order of the BS integration... high order is not needed
                        dt_new = dt_ideal[new_kopt];
                    } else {
                        new_kopt = std::min<size_t>(kmax - 1,
                                                    std::max<size_t>(2, k));
                        dt_new = dt_ideal[new_kopt];
                    }

                    break;
                }
            }


            if (k == kopt_now) { //test convergence at kopt

                if (error < 1.0) {//convergence is reached!
                    rejection = false;

                    new_kopt = std::min<size_t>(kmax - 1,
                                                std::max<size_t>(2, k)); // Default case.. leave unchanged
                    dt_new = dt_ideal[new_kopt];

                    if ((work_min[k - 1] < param_K1 * work_min[k])) {
                        new_kopt = std::min<size_t>(kmax - 1,
                                                    std::max<size_t>(2, k - 1)); // Decrease the BS order
                        dt_new = dt_ideal[new_kopt];
                    } else if ((work_min[k] < param_K2 * work_min[k - 1]) && !last_step_rejected) {
                        new_kopt = std::min<size_t>(kmax - 1,
                                                    std::max<size_t>(2, k + 1)); // Increase the BS order
                        dt_new = dt_ideal[new_kopt - 1];
                        dt_new *= cost[new_kopt] / static_cast<double>(cost[new_kopt - 1]);
                    }

                    break; // Break in any case because I reached convergence.

                } else if (error > mincoeff_for_kopt_convergence(k, steps_sequence)) {

                    rejection = true;

                    new_kopt = std::min<size_t>(kmax - 1,
                                                std::max<size_t>(2, k)); // Default case.. leave unchanged
                    dt_new = dt_ideal[new_kopt];

                    if ((work_min[k - 1] < param_K1 * work_min[k])) {
                        new_kopt = std::min<size_t>(kmax - 1,
                                                    std::max<size_t>(2, k - 1)); // Decrease the BS order
                        dt_new = dt_ideal[new_kopt];
                    } else if ((work_min[k] < param_K2 * work_min[k - 1]) && !last_step_rejected) {
                        new_kopt = std::min<size_t>(kmax - 1,
                                                    std::max<size_t>(2, k + 1)); // Increase the BS order
                        dt_new = dt_ideal[new_kopt - 1];
                        dt_new *= cost[new_kopt] / static_cast<double>(cost[new_kopt - 1]);
                    }

                    break;
                }
            }


            if (k == kopt_now + 1) { // Test convergence at k_opt+1
                if (error < 1.0) {   // Convergence at kopt+1
                    rejection = false;

                    new_kopt = std::min<size_t>(kmax - 1,
                                                std::max<size_t>(2, k - 1)); // Default case.. leave the order as it is.
                    dt_new = dt_ideal[new_kopt];

                    if (work_min[k - 2] < param_K1 * work_min[k - 1]) {
                        new_kopt = std::min<size_t>(kmax - 1,
                                                    std::max<size_t>(2, k - 2));
                        dt_new = dt_ideal[new_kopt];
                    } else if ((work_min[k] < param_K2 * work_min[k - 1]) && !last_step_rejected) {
                        new_kopt = std::min<size_t>(kmax - 1,
                                                    std::max<size_t>(2, k));
                        dt_new = dt_ideal[new_kopt - 1] * cost[new_kopt] / static_cast<double>(cost[new_kopt - 1]);
                    }

                } else {
                    rejection = true;

                    new_kopt = std::min<size_t>(kmax - 1,
                                                std::max<size_t>(2, k - 1)); //default case.. leave the order as it is.
                    dt_new = dt_ideal[new_kopt];

                    if (work_min[k - 2] < param_K1 * work_min[k - 1]) {
                        new_kopt = std::min<size_t>(kmax - 1,
                                                    std::max<size_t>(2, k - 2));
                        dt_new = dt_ideal[new_kopt];
                    } else if ((work_min[k] < param_K2 * work_min[k - 1]) && !last_step_rejected) {
                        new_kopt = std::min<size_t>(kmax - 1,
                                                    std::max<size_t>(2, k));
                        dt_new = dt_ideal[new_kopt - 1] * cost[new_kopt] / static_cast<double>(cost[new_kopt - 1]);
                    }
                }

                break;
            }

        } // End BS iterations

        last_step_rejected = rejection;

        dt_next = dt_new;
        dt_prev = dt;
        kopt_now = new_kopt;

        if constexpr (debug) {
            std::cout << " dt_prev " << dt_prev << std::endl;
            std::cout << " dt_next " << dt_next << std::endl;
            std::cout << "rejection: " << rejection << std::endl;
            std::cout << "new kopt : " << kopt_now << std::endl;
        }
        //exit(EXIT_FAILURE);

        return rejection;
    }

    void extrapolate_table(size_t k, BSTable *TT, BSTable &FinalT,
                           std::vector<std::vector<double> > &CC) {

        //Neville recursive formula to calculate the correct coefficient, needed for the extrapolation
        // Vectors with N=Neqs
        for (size_t i = 0; i < LeapStepped.Neqs; i++) {
            for (size_t j = k - 1; j > 0; --j) {
                TT[j - 1].pos[i] = TT[j].pos[i] + CC[k][j] * (TT[j].pos[i] - TT[j - 1].pos[i]);
                TT[j - 1].vel[i] = TT[j].vel[i] + CC[k][j] * (TT[j].vel[i] - TT[j - 1].vel[i]);
            }
            FinalT.pos[i] = TT[0].pos[i] + CC[k][0] * (TT[0].pos[i] - FinalT.pos[i]);
            FinalT.vel[i] = TT[0].vel[i] + CC[k][0] * (TT[0].vel[i] - FinalT.vel[i]);
        }

        // Vectors with N=Npart
        if (LeapStepped.vdep) {
            for (size_t i = 0; i < System.Npart; i++) {
                for (size_t j = k - 1; j > 0; --j) {
                    TT[j - 1].eloss[i] = TT[j].eloss[i] + CC[k][j] * (TT[j].eloss[i] - TT[j - 1].eloss[i]);
                    if constexpr (withspin)
                        TT[j - 1].spin[i] = TT[j].spin[i] + CC[k][j] * (TT[j].spin[i] - TT[j - 1].spin[i]);
                }
                FinalT.eloss[i] = TT[0].eloss[i] + CC[k][0] * (TT[0].eloss[i] - FinalT.eloss[i]);
                if constexpr (withspin) FinalT.spin[i] = TT[0].spin[i] + CC[k][0] * (TT[0].spin[i] - FinalT.spin[i]);
            }
        }

        // Scalars
        for (size_t j = k - 1; j > 0; --j) {
            TT[j - 1].time = TT[j].time + CC[k][j] * (TT[j].time - TT[j - 1].time);
            TT[j - 1].B = TT[j].B + CC[k][j] * (TT[j].B - TT[j - 1].B);
            if (System.Config.TTL) TT[j - 1].omega = TT[j].omega + CC[k][j] * (TT[j].omega - TT[j - 1].omega);
        }
        FinalT.time = TT[0].time + CC[k][0] * (TT[0].time - FinalT.time);
        FinalT.B = TT[0].B + CC[k][0] * (TT[0].B - FinalT.B);
        if (System.Config.TTL) FinalT.omega = TT[0].omega + CC[k][0] * (TT[0].omega - FinalT.omega);
    }

    void extrapolate_scalar(size_t k, double *TT, std::vector<std::vector<double> > &CC, double *extrapolated) {

        for (size_t j = k - 1; j > 0; --j) {
            TT[j - 1] = TT[j] + CC[k][j] * (TT[j] - TT[j - 1]);
            //Neville recursive formula to calculate the correct table[0][i] coefficient, needed for the extrapolation
        }

        *extrapolated = TT[0] + CC[k][0] * (TT[0] - *extrapolated);
    }


    /// Initialize or reset to the first BS iteration
    void first_step_init() {

        //first_step = true; //I am on the first BS step
        rejection = false; //the previous BS step was not rejected
        last_step_rejected = false; //Was the previous step rejected?
        kopt_now = kmax / 2; //estimation of the optimal iterations (half of maximum iterations)
    }

    void reset_bulirsch() {
        first_step_init();
        LeapStepped.initialize_arrays_step_zero();
    }

    double estimate_error(BSTable const &now, BSTable const &now_km1) {
        // Error on time
        double scale = abs_err_tolerance + rel_err_tolerance * fabs(now.time);
        double err = fabs(now.time - now_km1.time) / scale;
        //std::cout << "] new time " << now.time << "  ] km1 time " << now_km1.time << std::endl;
        //std::cout << "err " << err <<  "] scale " << scale << std::endl;
        //err = 0.0;
        if constexpr (debug) {
            std::cout << " terr  = " << err << std::endl;
        }

        // Error on B
        scale = abs_err_tolerance + rel_err_tolerance * std::max(fabs(LeapStepped.B0), fabs(now.B));
        double err_B = fabs(now.B - now_km1.B) / scale;
        if constexpr (debug) {
            std::cout << " Berr  = " << err_B << std::endl;
        }

        err = (err > err_B) ? err : err_B;

        auto abs_err3 = double3(abs_err_tolerance);

        for (size_t i = 0; i < LeapStepped.Neqs; i++) {
            // Pos
            double3 scale3 = abs_err3 + rel_err_tolerance * max(LeapStepped.pos0[i].fabs(), now.pos[i].fabs());
            double3 erri = (now.pos[i] - now_km1.pos[i]).fabs() / scale3;
            double max3 = erri.max();

            if constexpr (debug) {
                std::cout << " p" << i <<  "err = " << max3 << std::endl;
            }

            err = (err > max3) ? err : max3;

            // Vel
            scale3 = abs_err3 + rel_err_tolerance * max(LeapStepped.vel0[i].fabs(), now.vel[i].fabs());
            erri = (now.vel[i] - now_km1.vel[i]).fabs() / scale3;
            max3 = erri.max();

            if constexpr (debug) {
                std::cout << " v" << i <<  "err = " << max3 << std::endl;
            }

            err = (err > max3) ? err : max3;
        }

        if (System.Config.TTL) {
            scale = abs_err_tolerance + rel_err_tolerance * std::max(fabs(LeapStepped.omega0), fabs(now.omega));
            double err_ome = fabs(now.omega - now_km1.omega) / scale;
            err = (err > err_ome) ? err : err_ome;
        }

        if (LeapStepped.vdep) {
            for (size_t i = 0; i < System.Npart; i++) {
                scale = abs_err_tolerance +
                        rel_err_tolerance * std::max(fabs(LeapStepped.eloss0[i]), fabs(now.eloss[i]));
                double err_eloss = fabs(now.eloss[i] - now_km1.eloss[i]) / scale;
                err = (err > err_eloss) ? err : err_eloss;

                if constexpr (debug) {
                    std::cout << " e" << i <<  "err = " << err_eloss << std::endl;
                }

                if constexpr (withspin) {
                    // Pos
                    double3 scale3 =
                            abs_err3 + rel_err_tolerance * min(LeapStepped.spin0[i].fabs(), now.spin[i].fabs());
                    double3 erri = (now.spin[i] - now_km1.spin[i]).fabs() / scale3;
                    double max3 = erri.max();

                    err = (err > max3) ? err : max3;
                }
            }
        }

        return err;
    }

    void mitigate_negative_step() {
        // No matter what, reject and do not consider for extrapolation, redo with halved timestep
        rejection = true;
        last_step_rejected = rejection;
        dt_next = dt * 0.5;
        dt_prev = dt;

        if constexpr (debug) std::cout << " NEG 0 STEP at ds " << dt << std::endl;
        LeapStepped.revert_step(System.ch_pos, System.ch_vel, System.ch_mass,
                                System.ch_spin, System.ch_xdata);
        //kopt_now = new_kopt; Do not change kopt
    }

    double estimate_ideal_dt(double step, double err, size_t k) {

        double expo = exponents[k];
        double fmin = facmin[k];
        double fac(0);

        if (err == double(0))
            fac = 1.0 / fmin;
        else {

            fac = param_S1 * pow(param_S2 / err, expo);

            if (fac < fmin / param_S4) fac = fmin / param_S4;
            else if (fac > 1.0 / fmin) fac = 1.0 / fmin;
        }
        //cout<<"err "<<err<<endl;

        return static_cast<double>(fabs(step * fac));
    }

private:

    void init_matrices() {
        coeff_extr.resize(kmax + 1);
        coeff_ratextr.resize(kmax + 1);
        cost.resize(kmax + 1);
        facmin.resize(kmax + 1);
        exponents.resize(kmax + 1);

        cost[0] = steps_sequence[0] + 1; //see Numerical recipes eq. (17.3.12)
        //::suspect:: cost[0] = steps_sequence[0];

        for (size_t i = 0; i < kmax + 1; i++) {

            coeff_extr[i].resize(i); //Take into account that it is an upper/triangular matrix
            coeff_ratextr[i].resize(i);

            if (i != 0)
                cost[i] = cost[i - 1] + steps_sequence[i];

            exponents[i] = 1.0 / (2.0 * i + 1.0);
            facmin[i] = pow(param_S3, exponents[i]);

            for (size_t k = 0; k < i; ++k) {
                double denom = (double) steps_sequence[i] / (double) steps_sequence[k];
                coeff_extr[i][k] = 1.0 / (denom * denom - double(1)); // extrapolation coefficients: NR eq. (17.3.8)
                coeff_ratextr[i][k] = denom * denom;
                //cout<<" initcoeff = "<<coeff_ratextr[i][k]<<"   "<<i<<"   "<<k<<endl;
            }
        }
    }

    void store_integrated_and_reset(BSTable &Table) {
        for (size_t i = 0; i < LeapStepped.Neqs; i++) {
            Table.pos[i] = System.ch_pos[i];
            Table.vel[i] = System.ch_vel[i];
        }
        if (LeapStepped.vdep) {
            for (size_t i = 0; i < System.Npart; i++) {
                Table.eloss[i] = System.ch_xdata[i].eloss;
                if constexpr (withspin) Table.spin[i] = System.ch_spin[i];
            }
        }
        Table.time = System.time_phys;
        Table.B = LeapStepped.B;
        Table.omega = LeapStepped.omega;

        // Reset integrated values
        LeapStepped.revert_step(System.ch_pos, System.ch_vel, System.ch_mass,
                                System.ch_spin, System.ch_xdata);
    }

    double mincoeff_for_kopt_convergence(size_t kk, std::vector<int> &nseq) {

        double val(1);

        if (kk == kopt_now - 1)
            val = double(nseq[kk + 1] * nseq[kk + 2]) / double(nseq[0] * nseq[0]);
        else if (kk == kopt_now)
            val = double(nseq[kk + 1]) / double(nseq[0]);
        //::suspect:: double(nseq[kk])/double(nseq[0]);

        return val * val;
    }

    ///Initialize the sequence of substeps per cycle;
    void init_steps_seq() {

        steps_sequence.resize(kmax + 1);

        for (int i = 0; i < static_cast<int>(kmax + 1); i++)
            steps_sequence[i] = 2 * (i + 1);

    }

    void allocate_arrays() {
        FinalTable.allocate_arrays(LeapStepped.Neqs, System.Npart);
        for (size_t i = 0; i < kmax; i++) {
            ExtrapTables[i].allocate_arrays(LeapStepped.Neqs, System.Npart);
        }
    }

    bool rejection; //Should I reject the performed iteration?
    bool last_step_rejected; //Was the previous step rejected?
    size_t kopt_now; //My current optimal number of BS iterations

    static const size_t kmax = 8; //Maximum number of BS iterations allowed
    const double param_S1 = 0.94; //S1 param (see NR eq. (17.3.11)), it was 0.94 //1.2 seems the best value for time steps
    const double param_S2 = 0.65; //S2 param (see NR eq. (17.3.11))
    const double param_S3 = 0.02; //S3 param (see NR eq. (17.3.22))
    const double param_S4 = 4.0;  // S4 param (see NR eq. (17.3.22))
    const double param_K1 = 0.8; //parameter for determining the minimum work for each BS iteration
    const double param_K2 = 0.9; //parameter for determining the minimum work for each BS iteration

    double abs_err_tolerance;//tolerance on the absolute errors of BS extrapolations
    double rel_err_tolerance;//tolerance on the relative errors of BS extrapolations

    double dt_prev = 0.0; // Time step of the previous BS iteration
    double dt; // Current time step
    double dt_next; // Time step of the next BS iteration
    double Dtime; //interval of physical time integrated inside the leapfrog algorithm

    double error; //error estimate of the BS iteration;
    std::vector<int> steps_sequence; //Number of substeps per BS iteration
    std::vector<int> cost;//needed to estimate the optimal work per BS step

    // Tableu logic:
    // FinalTable always holds the last diagonal value (current k extrapolation)
    // ExtrapTables always holds the last row of values at k-1, reversed (element 0 is k,k-1), element k is (k,0)
    BSTable ExtrapTables[kmax];
    BSTable FinalTable;

    std::vector<std::vector<double>> coeff_extr;
    std::vector<std::vector<double>> coeff_ratextr;
    std::vector<double> facmin; //to store F (see Numerical recipes eq. (17.3.22))
    std::vector<double> exponents; //to store 1/(2*k+1)

    LeapfrogStepped &LeapStepped;
    ChainSys &System;

    size_t Nreject;
};

typedef BSExtrapolator<TsunamiConfig::wSpins, TsunamiConfig::debug_bs> BSExtrap;


#endif
