//
// Created by lex on 23/12/18.
//

#ifndef TSUNAMI_LEAPFROG_STEPPED_HPP
#define TSUNAMI_LEAPFROG_STEPPED_HPP


#include <cstdlib>
#include <iostream>
#include <memory>
#include "custom_types.hpp"
#include "chain.hpp"
#include "simprof.hpp"

template<bool withspin, bool debug, bool profiling>
class Leapfrog_stepped {
    template<bool withspin_bs, bool debug_bs> friend class BSExtrapolator;
    template<bool profiling_tsu, bool debug_tsu> friend class TsunamiClass;

public:
    explicit Leapfrog_stepped(ChainSys &System) :
            System(System),
            vdep(System.Config.wPNs or
                 System.Config.wEqTides or
                 System.Config.wDynTides),
            midstep(vdep or
                    System.Config.TTL or
                    System.Config.wExt) {
    }

    ~Leapfrog_stepped() {
        if (Neqs != 0) // Deallocate only if not initialized
            deallocate_arrays();
    }

    void initialize() {
        // Adding one extra particle for the CoM
        auto extraP = static_cast<size_t>(System.Config.wExt);
        size_t Neqs_new = System.Nchain + extraP;

        // Deallocate, but only if necessary
        if (Neqs != Neqs_new and Neqs != 0) {
            deallocate_arrays();
        }
        Neqs = Neqs_new;

        allocate_arrays();
        initialize_arrays_step_zero();

        vdep = System.Config.wPNs or System.Config.wEqTides or System.Config.wDynTides;
        midstep = vdep or System.Config.TTL or System.Config.wExt;
    }

    void initialize_arrays_step_zero() {
        System.calc_accelerations_nvdep(System.ch_pos, unchained_pos, acc, acc_ext,
                                        unchained_acc_ext, U, OMEGA, dOMEGA);
        T = System.calc_T(System.ch_vel, unchained_vel); // Initialize T

        omega = omega0 = OMEGA; // Initialize omega0, omega_fake0
        B = B0 = U - T; // Initialize B0
        dB = dB_ext = 0.0; // B variation due to velocity dependent perturbations
        domega = 0.0; // FTL domega

        /* Recompute each time masses change, if using TTL */
        if (System.Config.TTL and System.Config.wMassEvol) System.compute_m2ttl(System.ch_mass);
    }

    void save_step_zero(double3 *pos, double3 *vel, double *mass, double3 *spin, pinfo *xdata) {
        for (size_t i = 0; i < Neqs; i++) {
            pos0[i] = pos[i];
            vel0[i] = vel[i];
        }

        for (size_t i = 0; i < System.Npart; i++) {
            if (vdep) {
                eloss0[i] = xdata[i].eloss;
                if constexpr (withspin) {
                    spin0[i] = spin[i];
                }
            }
            if (System.Config.wMassEvol) {
                mass0[i] = mass[i];
            }
        }
    }

    void revert_step(double3 *pos, double3 *vel, double *mass, double3 *spin, pinfo *xdata) {
        for (size_t i = 0; i < Neqs; i++) {
            pos[i] = pos0[i];
            vel[i] = vel0[i];
        }

        for (size_t i = 0; i < System.Npart; i++) {
            if (vdep) {
                xdata[i].eloss = eloss0[i];
                if constexpr (withspin) {
                    spin[i] = spin0[i];
                }
            }
            if (System.Config.wMassEvol) {
                mass[i] = mass0[i];
            }
        }
    }

    /**
     * Integrates the system over one timestep macrostep using nstep individual steps.
     * The basic scheme is a 2nd order leapfrog + midpoint step
     * General regularized scheme with velocity dependent perturbations AV(P,V,S),
     * gravitational accelerations A(P), non-velocity dependent perturbations Aext(P)
     * torque T(P,V), spin S, velocity dependent perturbations AVext(P,V)
     *
     * T0-Init: Calculate U = U0, T = T0, B = B0, OMEGA = omega = OMEGA0
     * 			Auxiliary velocity Vaux = V(0)
     *
     * T0-Drift: Calculate DeltaT_P from DeltaS/(alpha*(T+B) + beta*omega + gamma)
     *
     * 			 Advance positions P(1/2) = P(0) + V(0) * DeltaT_P
     *
     * T1/2-Kick: Update U, OMEGA; Calculate dOMEGA, A, Aext
     *            Calculate torque T[P(1/2), V(0), S(0)] and AV[P(1/2), V(0), S(0)]
     *
     * 			  Calculate DeltaT_V from DeltaS/(alpha*U + beta*OMEGA + gamma)
     *
     * 			  Advance auxiliary velocities Vaux(1/2) = Vaux(0) + DeltaT_V/2 * { A + AV[P(1/2), V(0), S(0)] }
     * 			  Advance auxiliary spin Saux(1/2) = Saux(0) + DeltaT_V/2 * T[P(1/2), V(0), S(0)]
     *
     * 			  Advance velocities V(1) = V(0) + DeltaT_V * { A + AV[P(1/2), Vaux(1/2), Saux(1/2)] }
     * 			  Advance spin S(1) = S(0) + DeltaT_V * T[P(1/2), Vaux(1/2), Saux(1/2)]
     *
     * 			  Advance TTL omega = omega(0) + DeltaT_V * SUM( Vaux(1/2) * dOMEGA )
     *
     * 			  Advance binding energy B = B0 + DeltaT_V * SUM( Vaux(1/2) * m *
     * 			  { Aext[P(1/2)] + AV[P(1/2), Vaux(1/2), Saux(1/2)] + AVext[P(1/2), Vaux(1/2)] } , Vaux(1/2) )
     *
     * 			  Advance auxiliary velocities Vaux(1) = Vaux(1/2) + DeltaT_V/2 * [A + a(P(1/2), V(1)) ]
     *  		  Advance auxiliary spin Saux(1) = Saux(1/2) + DeltaT_V/2 * T[P(1/2), V(1), S(1)]
     *
     * T1-Drift: Update T
     *
     * 			 Calculate DeltaT_P from DeltaS/(alpha*(T+B) + beta*omega + gamma)
     *
     * 			 Advance positions P(1) = P(1/2) + V(1) * DeltaT_P
     *
     * Comment: we should not set Vaux = V at the end of T1-Drift/beginning of T0-Drift?
     * 			The midpoint is meant only for the BS extrapolation. Does it have any sense to keep
     * 			using the final value Vaux without?
     *
     * @param nsteps Number of substeps
     * @param macrostep Regularized step
     * @param physical_timestep Physical step
     * @param pos_inout Input/Output position
     * @param vel_inout Input/Output velocity
     */
    void integrate(int nsteps, double macrostep, double &physical_timestep,
                   double3 *pos_inout, double3 *vel_inout, double3 *spin_inout) {

        if constexpr (profiling) ProfLeap.start();

        System.time_phys = System.time_fict = 0.0;
        // Initialize starting variables at the beginning of the leapfrog iterations
        omega = omega0; // Initialize omega
        B = B0; // Initialize B
        System.CollInfo.collflag = Nbodyalgorithms::NOCOLLISION;

        if constexpr (debug) {
            std::cout << "LF Init step " << System.time0 << std::endl;
        }

        this->pos_inout = pos_inout;
        this->vel_inout = vel_inout;
        if constexpr (withspin) {
            if (vdep) {
                this->spin_inout = spin_inout;
                for (size_t i = 0; i < System.Npart; i++) {
                    spin_aux[i] = spin_inout[i];
                }
            }
        }

        // Re-initializing vel_aux
        for (size_t i = 0; i < Neqs; i++) {
            vel_aux[i] = vel_inout[i];
        }

        ///LEAPFROG///////////////////////////////////////////
        double local_step = macrostep / static_cast<double>(nsteps);
        advance_positions(0.5 * local_step);

        for (int iter = 1; iter < nsteps; iter++) {
            advance_velocities(local_step);
            advance_positions(local_step);
        }

        advance_velocities(local_step);
        advance_positions(0.5 * local_step);
        /////////////////////////////////////////////////////

        physical_timestep = System.time_phys;

        if constexpr (profiling) ProfLeap.stop_store();
    }

    void update_parameters() {
        vdep = System.Config.wPNs or System.Config.wEqTides or System.Config.wDynTides;
        midstep = vdep or System.Config.TTL or System.Config.wExt;
    }


private:


    double calc_dt_pos(double ds, double TpB, double omega) {
        if constexpr (debug) {
            double denom = (System.Config.alpha * TpB + System.Config.beta * omega + System.Config.gamma);
            std::cout << " U = " << TpB << "  - dt = " << ds << " / " << denom << std::endl;
        }

        return (ds / (System.Config.alpha * TpB + System.Config.beta * omega + System.Config.gamma));
    }

    double calc_dt_vel(double ds, double U, double OMEGA) {
        if constexpr (debug) {
            double denom = (System.Config.alpha * U + System.Config.beta * OMEGA + System.Config.gamma);
            std::cout << " U = " << U << "  - dt = " << ds << " / " << denom << std::endl;
        }

        return (ds / (System.Config.alpha * U + System.Config.beta * OMEGA + System.Config.gamma));
    }

    /**
     * Drift step
     * @param ds Regularized timestep
     */
    void advance_positions(double ds) { //X(s) method: notation from the Cambridge N-body lectures

        T = System.calc_T(vel_inout, unchained_vel); // Compute T

        double dt = calc_dt_pos(ds, T + B, omega); // Get positions' physical timestep
        if ((dt < 0) or std::isinf(dt)) throw IntegrationError(ErrType::DT_NEG_OR_INF);


        if constexpr (debug) {
            std::cout << " DRIFT dt " << dt;
            if (dt < 0.0) std::cout << "  !!!!!!! ";
            std::cout << std::endl;
        }

        System.time_phys += dt; //update local time (private variable)
        System.time_fict += ds;

        for (size_t i = 0; i < Neqs; i++) {
            pos_inout[i] += dt * vel_inout[i];
        }
    }


    /**
     * Kick step
     * @param ds Regularized timestep
     */
    void advance_velocities(double ds) {

        // Calculate accelerations and omega
        System.calc_accelerations_nvdep(pos_inout, unchained_pos, acc,
                                        acc_ext, unchained_acc_ext, U, OMEGA, dOMEGA);

        System.dt = calc_dt_vel(ds, U, OMEGA); // Get velocities' physical timestep
        if ((System.dt < 0) or std::isinf(System.dt)) throw IntegrationError(ErrType::DT_NEG_OR_INF);

        if constexpr (debug) {
            std::cout << " KICK dt " << System.dt;
            if (System.dt < 0.0) std::cout << "  !!!!!!! ";
            std::cout << std::endl;
        }

        // Modified midpoint method (Eqs 23 24 25, Mikkola & Merrit 2006)
        // Auxiliary step to predict velocities (for velocity dependent perturbations and TTL)
        if (midstep) advance_fake_velocities(0.5 * System.dt);
        advance_real_velocities(System.dt); // Perform the real velocity advance
        if (midstep) advance_fake_velocities(0.5 * System.dt); // Auxiliary step to complete predictor of velocities over dt
    }

    /**
     * Advance auxiliary velocity for midpoint method
     * @param dt Physical timestep
     */
    void advance_fake_velocities(double dt) {

        System.calc_accelerations_vdep<false>(pos_inout, unchained_pos,
                                       vel_inout, unchained_vel,
                                       acc_vdep, unchained_acc_vdep,
                                       acc_ext_vdep, unchained_acc_ext_vdep,
                                       spin_inout, torque);

        for (size_t i = 0; i < Neqs; i++) {
            vel_aux[i] += dt * (acc[i] + acc_vdep[i] + acc_ext[i] + acc_ext_vdep[i]);
        }
        if constexpr (withspin) {
            if (vdep) {
                for (size_t i = 0; i < System.Npart; i++) {
                    spin_aux[i] += torque[i] * dt;
                }
            }
        }
    }

    void advance_fake_velocities_wrong(double dt) {
        for (size_t i = 0; i < Neqs; i++) {
            vel_aux[i] = vel_inout[i];
            if constexpr (withspin) spin_aux[i] = spin_inout[i];
        }
    }

    /**
     * Advance real velocity, B and omega using the midpoint method
     * @param dt Physical timestep
     */
    void advance_real_velocities(double dt) {

        if (midstep) {
            System.calc_accelerations_vdep<true>(pos_inout, unchained_pos,
                                           vel_aux, unchained_vel_aux,
                                           acc_vdep, unchained_acc_vdep,
                                           acc_ext_vdep, unchained_acc_ext_vdep,
                                           spin_aux, torque);
            calc_domega_dB_midpoint();

            B -= dt * (dB + dB_ext);
            omega += dt * domega;
        }

        if constexpr (debug) {
            for (size_t i = 0; i < Neqs; i++) {
                std::cout << " CV" << i << "  = " << vel_inout[i] << std::endl;
                std::cout << " CA" << i << "  = " << acc[i] << std::endl;
                std::cout << " CAV" << i << " = " << acc_vdep[i] << std::endl;
            }
        }

        for (size_t i = 0; i < Neqs; i++) {
            vel_inout[i] += dt * (acc[i] + acc_vdep[i] + acc_ext[i] + acc_ext_vdep[i]);
        }
        if constexpr (withspin) {
            if (vdep) {
                if constexpr (debug) {
                    for (size_t i = 0; i < System.Npart; i++) {
                        std::cout << " S" << i << " = " << spin_inout[System.invchain[i]] << std::endl;
                        std::cout << " T" << i << " = " << torque[System.invchain[i]] << std::endl;
                    }
                }

                for (size_t i = 0; i < System.Npart; i++) {
                    spin_inout[i] += torque[i] * dt;
                }
            }
        }

        /** Alternative dB/domega evaluation and advance
         * if (midstep) {
         *      calc_vmid();
         *      calc_domega_dB_midpoint();
         *      B -= dt * (dB + dB_ext);
         *      omega += dt * domega;
         * }
         */

    }

    /**
     * Computes variations of omega and B using the using v_1/2 predicted during the midpoint step
     */
    void calc_domega_dB_midpoint()  {
        // TTL small omega variations
        if (System.Config.TTL) {
            domega = 0.0;
            for (size_t i = 0; i < System.Npart; i++) {
                domega += dOMEGA[i] * unchained_vel_aux[i];
            }
        }

        // Velocity-dependent forces B variations
        if (vdep) {
            dB = 0.0;
            for (size_t i = 0; i < System.Npart; i++) {
                // Update dB variations
                dB += System.ch_mass[i] * (unchained_acc_vdep[i] * unchained_vel_aux[i]);
            }
        }

        // External forces B variations
        if (System.Config.wExt) {
            dB_ext = 0.0;
            for (size_t i = 0; i < System.Npart; i++) {
                // Update dB variations
                dB_ext += System.ch_mass[i] * ((unchained_acc_ext[i] + unchained_acc_ext_vdep[i]) * unchained_vel_aux[i]);
            }
        }
    }

    /**
     * Computes variations of omega and B using the using (v0 + v1)/2
     * Stores into unchained_vel_aux
     */
    void calc_vmid()  {

        // Overwriting unchained_vel_aux with positions at t1
        System.to_individual_chain_ordered_single_cdm(unchained_vel_aux, vel_inout);

        // Unchained_vel still contains velocities at t0
        for (size_t i = 0; i < System.Npart; i++) {
            unchained_vel_aux[i] = 0.5 * (unchained_vel[i] + unchained_vel_aux[i]);
        }
    }

    void allocate_arrays() {
        acc = new double3[Neqs];
        acc_vdep = new double3[Neqs];
        acc_ext = new double3[Neqs];
        acc_ext_vdep = new double3[Neqs];
        vel_aux = new double3[Neqs];

        dOMEGA = new double3[System.Npart];

        unchained_pos = new double3[System.Npart];
        unchained_vel = new double3[System.Npart];
        unchained_vel_aux = new double3[System.Npart];
        unchained_acc_ext = new double3[System.Npart];
        unchained_acc_ext_vdep = new double3[System.Npart];
        unchained_acc_vdep = new double3[System.Npart];

        pos0 = new double3[Neqs];
        vel0 = new double3[Neqs];
        eloss0 = new double[System.Npart];
        if constexpr (withspin) {
            torque = new double3[System.Npart];
            spin_aux = new double3[System.Npart];
            spin0 = new double3[System.Npart];
        }

        if (System.Config.wMassEvol) {
            mass0 = new double[System.Nchain];
        }
    }

    void deallocate_arrays() {
        delete[] acc;
        delete[] acc_vdep;
        delete[] acc_ext;
        delete[] acc_ext_vdep;
        delete[] vel_aux;

        delete[] dOMEGA;

        delete[] unchained_pos;
        delete[] unchained_vel;
        delete[] unchained_vel_aux;
        delete[] unchained_acc_ext;
        delete[] unchained_acc_ext_vdep;
        delete[] unchained_acc_vdep;

        delete[] vel0;
        delete[] pos0;
        delete[] eloss0;
        if constexpr (withspin) {
            delete[] torque;
            delete[] spin_aux;
            delete[] spin0;
        }

        if (System.Config.wMassEvol) {
            delete[] mass0;
        }
    }

    ChainSys &System; // Reference member
    // Default value Neqs means no initalization

    double T, U, OMEGA;
    double dB_ext, dB, domega;
    double3 *dOMEGA;

    double3 *acc, *acc_vdep, *acc_ext, *acc_ext_vdep;
    double3 *torque, *spin_inout;
    double3 *pos_inout, *vel_inout;
    double3 *spin_aux, *vel_aux;

    // Auxiliary vectors, to be used during integration
    double3 *unchained_pos;
    double3 *unchained_vel;
    double3 *unchained_vel_aux;
    double3 *unchained_acc_ext;
    double3 *unchained_acc_ext_vdep;
    double3 *unchained_acc_vdep;

    // Extrapolated scalars
    double B0; ///< Initial binding energy
    double B; ///< Integrated binding energy

    double omega0; ///< Initial TTL value
    double omega; ///< Integrated TTL value

    size_t Neqs = 0; // This is either Nchain or Nchain+1 if we have the center of mass
    double3 *pos0, *vel0;
    double3 *spin0;
    double *mass0, *eloss0;

    bool vdep;
    bool midstep;

    // Profiling utilities
    SimProf ProfLeap = SimProf();
};

typedef Leapfrog_stepped<TsunamiConfig::wSpins, TsunamiConfig::debug_lf, TsunamiConfig::useProfiling> LeapfrogStepped;

#endif //TSUNAMI_LEAPFROG_STEPPED_HPP
