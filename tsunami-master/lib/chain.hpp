//
// Created by lex on 19/12/18.
//

#ifndef TSUNAMI_CHAIN_HPP
#define TSUNAMI_CHAIN_HPP

#include <vector>
#include "config.hpp"
#include "custom_types.hpp"
#include "Nbodyalgorithms.hpp"
#include "IO.h"
#include "keplerutils.h"
#include <algorithm>
#include <iostream>

/**
 * Chain class
 */
template<bool withspin, bool usetdetracker, bool debug>
class Chain {

public:
    Chain(TsunamiConfig &ConfigS) : Config(ConfigS) {
        init_postnewtonians(Config.Mscale, Config.Lscale);
    }

    Chain(size_t N, TsunamiConfig &ConfigS)
            : Config(ConfigS), Npart(N), Nchain(N - 1) {
        init_postnewtonians(Config.Mscale, Config.Lscale);
        allocate_arrays();
    }

    /**
     * ChainedSystem destructor
     */
    ~Chain() {
        if (Npart != 0)  // Deallocate only if not initialized
            deallocate_arrays();
    }

    void init_postnewtonians(double Mscale, double Lscale) {
        // Post-Newtonian init
        double Tscale = sqrt(Lscale*Lscale*Lscale/(Mscale*KeplerUtils::G_yr_msun_au));
        double Vscale = Lscale/Tscale * KeplerUtils::au2km / KeplerUtils::yr2sec;

        vel_light = KeplerUtils::c_ms * 1e-3 / Vscale;
        vel_light2 = vel_light * vel_light;
        coeff1pn = 1.0 / vel_light2;
        coeff2pn = 1.0 / (vel_light2 * vel_light2);
        coeff25pn = 1.0 / (vel_light2 * vel_light2 * vel_light);
        coeff3pn = coeff1pn * coeff2pn;
        coeff35pn = coeff25pn * coeff1pn;
    }

    void init_regularization() {
        // Regularization init, only relevant for TTL
        if (Config.TTL) {
            ttl_tol = 1e-3;         // TTL tolerance value for small-m smooth behaviour
            if (Config.beta == 1.0) {
                ttl_tol = 1e30;    // If using only TTL, use no tolerance
            }
        }
        compute_m2ttl(mass);
    }

    void init_tdetracker() {
        if constexpr (usetdetracker) {
            for (size_t i = 0; i < Npart; i++) {
                if (xdata[i].stype == ptype::HIGH_MASS_MS) {
                    PTDEInfo.stars_ind.push_back(i);
                    Nbodyalgorithms::PTDE_history newPTDEHistory;
                    newPTDEHistory.ind = i;
                    PTDEInfo.PTDEHistories[i] = newPTDEHistory;
                }
            }
        }
    }

    void compute_m2ttl(double *mass) {

        m2ttl = double(0);

        for (size_t i = 0; i < Npart; i++) {
            for (size_t j = 0; j < i; j++) {
                m2ttl += mass[i] * mass[j];
            }
        }
        m2ttl /= (0.5 * Npart * (Npart - 1.0));
    }

    /**
     * Initialize chain and pair list
     */
    void init_chain() {

        // Reset chain, to be safe
        chain.erase(chain.begin(), chain.end());

        // Reset PairList
        PairList.erase(PairList.begin(), PairList.end());

        ParticlePair p;
        for (size_t i = 0; i < Npart; i++) {
            for (size_t j = 0; j < i; j++) {
                double dx = pos[i].x - pos[j].x;
                double dy = pos[i].y - pos[j].y;
                double dz = pos[i].z - pos[j].z;

                p.dist = dx * dx + dy * dy + dz * dz;
                p.i = i;
                p.j = j;
                p.used = false;

                PairList.push_back(p);
            }
        }

        std::sort(PairList.begin(), PairList.end(), Chain::sort_func); // Use the std::sort function

        chain.push_back(PairList[0].i);
        chain.push_back(PairList[0].j);
        PairList[0].used = true; // First pair is used by default
        //cout<<"adamo ed eva "<<ppchain[0].i<<" "<<ppchain[0].j<<endl;
        //now cycle over all the particles and construct the chain
        // Starting from the first pair
        for (int k = 0; k < (int) PairList.size(); k++) {
            if (!PairList[k].used) { // If I still have to use this pair
                if (chain_connection(PairList[k].i, PairList[k].j, PairList[k].used)) {
                    PairList[k].used = true;
                    if (chain.size() <
                        Npart) { k = -1; } // Couple attached to the chain. Inspect all the other remaining
                    // pairs in the correct order (small to big distance)
                }
            }
        }
        invchain.resize(Npart);
        for (size_t i = 0; i < Npart; i++) {
            size_t i_original = chain[i];
            invchain[i_original] = i;
        }
    }

    /**
     * Finds the new chain and stores it into \link Chain::chain \endlink.
     * The old chain get stored into \link Chain::chain_old \endlink.
     * Needs \link Chain::init_chain() \link to have run once to initialize the pair list \link Chain::ParticlePair \endlink.
     * Find chain also checks for collisions at current step
     */
    void find_chain() {
        for (size_t k = 0; k < PairList.size(); k++) {
            size_t p1 = PairList[k].i;
            size_t p2 = PairList[k].j;
            double dx = pos[p1].x - pos[p2].x;
            double dy = pos[p1].y - pos[p2].y;
            double dz = pos[p1].z - pos[p2].z;
            PairList[k].dist = dx * dx + dy * dy + dz * dz;
            double sumrad = Config.dcoll * (radius[p1] + radius[p2]);
            if (PairList[k].dist <= sumrad*sumrad) {
                CollInfo.collflag = Nbodyalgorithms::COLLISION_REAL;
                CollInfo.collind[0] = p1;
                CollInfo.collind[1] = p2;
                CollInfo.colltime = time_phys;
                CollInfo.colltime_ds = time_fict;
//                std::cout << "found inter-timestep collision " << std::endl;
                // We can interrupt chain building, as we will reset the step
                return;
            }

            PairList[k].used = false;
        }

        std::sort(PairList.begin(), PairList.end(), Chain::sort_func); // Use the std::sort function
        
        max_pair = PairList.back();
        min_pair = PairList.front();
        min_pair.dist = sqrt(min_pair.dist);
        max_pair.dist = sqrt(max_pair.dist);
        max_pair_smallest = max_pair.dist < max_pair_smallest.dist ? max_pair : max_pair_smallest;
        min_pair_smallest = min_pair.dist < min_pair_smallest.dist ? min_pair : min_pair_smallest;

        // Save old chain
        chain_old = chain;
        // Reset chain
        chain.erase(chain.begin(), chain.end());

        chain.push_back(PairList[0].i);
        chain.push_back(PairList[0].j);
        PairList[0].used = true; // First pair is used by default
        //cout<<"adamo ed eva "<<ppchain[0].i<<" "<<ppchain[0].j<<endl;
        //now cycle over all the particles and construct the chain
        // Starting from the first pair
        for (int k = 0; k < (int) PairList.size(); k++) {
            if (!PairList[k].used) { // If I still have to use this pair
                if (chain_connection(PairList[k].i, PairList[k].j, PairList[k].used)) {
                    PairList[k].used = true;
                    if (chain.size() <
                        Npart) { k = -1; } // Couple attached to the chain. Inspect all the other remaining
                    // pairs in the correct order (small to big distance)
                }
            }
        }
    }

    void revert_chain() {
        chain = chain_old;
    }

    /**
     * Checks if particle pair can be attached to the current chain
     * @param[in] p1	Particle 1 index
     * @param[in] p2	Particle 2 index
     * @param[out] used Used in chain flag
     * @return			True if particle pair was attached, otherwise False
     */
    bool chain_connection(size_t p1, size_t p2, bool &used) {

        size_t first = chain[0];
        size_t last = chain[chain.size() - 1];

        //    cout<<"inspect "<<p1<<"  "<<p2<<endl;
        if (first == p1) { // Checks if p2 can be attached to head of chain
            auto found = find(chain.begin(), chain.end(), p2);

            if (found != chain.end()) {
                //    cout << "found p2 begin " << p2 << endl;
                used = true;
                return false;
            } else {
                //    cout << "added p2 begin " << p2 << endl;
                chain.insert(chain.begin(), p2);
                return true;
            }
        } else if (first == p2) { // Checks if p1 can be attached to head of chain
            auto found = find(chain.begin(), chain.end(), p1);

            if (found != chain.end()) {
                //    cout << "found p1 begin " << p1 << endl;
                used = true;
                return false;
            } else {
                //    cout << "added p1 begin " << p1 << endl;
                chain.insert(chain.begin(), p1);
                return true;
            }
        }

        if (last == p1) { // Checks if p2 can be attached to tail of chain
            auto found = find(chain.begin(), chain.end(), p2);

            if (found != chain.end()) {
                //    cout << "found p2 last " << p2 << endl;
                used = true;
                return false;
            } else {
                //    cout << "added p2 last " << p2 << endl;
                chain.push_back(p2);
                return true;
            }

        } else if (last == p2) { // Checks if p1 can be attached to tail of chain
            auto found = find(chain.begin(), chain.end(), p1);

            if (found != chain.end()) {
                //    cout << "found p1 last " << p1 << endl;
                used = true;
                return false;
            } else {
                //    cout << "added p1 last " << p1 << endl;
                chain.push_back(p1);
                return true;
            }
        }
        return false;
    }

    /**
     * Initializes the first chain-ordered arrays
     * @param[in] pos original position array
     * @param[in] vel original velocity array
     * @param[in] mass  original mass array
     * @param[in] soft  original softening array
     * @param[in] radius   original radii array
     * @param[in] xdata original particle extra data array
     * @param[in] pcom  center of mass position
     * @param[in] vcom  center of mass velocity
     */
    void initialize_chain_arrays() {

        for (size_t i = 0; i < Nchain; i++) {
            size_t k1 = chain[i];
            size_t k2 = chain[i + 1];

            ch_pos[i] = pos[k2] - pos[k1];
            ch_vel[i] = vel[k2] - vel[k1];

            ch_mass[i] = mass[k1];
            ch_radius[i] = radius[k1];
            ch_xdata[i] = xdata[k1];
            if constexpr (withspin) ch_spin[i] = spin[k1];
        }

        size_t k1 = chain[Nchain];
        ch_mass[Nchain] = mass[k1];
        ch_radius[Nchain] = radius[k1];
        ch_xdata[Nchain] = xdata[k1];
        if constexpr (withspin) ch_spin[Nchain] = spin[k1];

        if (Config.wExt) {
            ch_pos[Nchain] = pcom;
            ch_vel[Nchain] = vcom;
        }
    }


    /**
     * Updates mass, softening, radius, particle extra data, spin stored into the chain with the unchained arrays
     */
    void update_chained_data() {
        for (size_t i = 0; i < Npart; i++) {
            size_t i_original = chain[i];

            ch_mass[i] = mass[i_original];
            ch_radius[i] = radius[i_original];
            ch_xdata[i] = xdata[i_original];
            if constexpr (withspin) ch_spin[i] = spin[i_original];
        }
    }


    /**
     * Converts given chained positions and velocities to chain-ordered individual positions and velocities
     * It rescales to the center of mass
     * @param[out] loc_p Individual positions
     * @param[out] loc_v Individual velocities
     * @param[in] ch_p Chained positions
     * @param[in] ch_v Chained velocities
     */
    void to_individual_chain_ordered_cdm(double3 *loc_p, double3 *loc_v, double3 *ch_p, double3 *ch_v) {

        // Assign 0 to the first element
        loc_p[0] = 0;
        loc_v[0] = 0;

        // Reconstruct all the elements
        for (size_t i = 0; i < Nchain; i++) {
            loc_p[i + 1] = loc_p[i] + ch_p[i];
            loc_v[i + 1] = loc_v[i] + ch_v[i];
        }

        // Evaluate new cdm
        Nbodyalgorithms::scale_to_cdm(loc_p, loc_v, ch_mass, Npart);
    }

    /**
     * Converts given chained positions and velocities to chain-ordered individual positions and velocities
     * Does not rescale to the center of mass
     * @param[out] loc_p Individual positions
     * @param[out] loc_v Individual velocities
     * @param[in] ch_p Chained positions
     * @param[in] ch_v Chained velocities
     */
    void to_individual_chain_ordered_nocdm(double3 *loc_p, double3 *loc_v, double3 *ch_p, double3 *ch_v) {

        // Assign 0 to the first element
        loc_p[0] = 0;
        loc_v[0] = 0;

        // Reconstruct all the elements
        for (size_t i = 0; i < Nchain; i++) {
            loc_p[i + 1] = loc_p[i] + ch_p[i];
            loc_v[i + 1] = loc_v[i] + ch_v[i];
        }
    }

    /**
     * Converts chained velocities/position to chain-ordered individual positions and velocities
     * Scales to the center of mass
     * @param[out] loc_pv Individual positions or velocities
     * @param[in] ch_pv Chained positions or velocities
     */
    void to_individual_chain_ordered_single_cdm(double3 *loc_pv, double3 *ch_pv) {

        // Assign 0 to the first element
        loc_pv[0] = 0;

        // Reconstruct all the elements
        for (size_t i = 0; i < Nchain; i++) {
            loc_pv[i + 1] = loc_pv[i] + ch_pv[i];
        }
        // Evaluate new cdm
        Nbodyalgorithms::scale_to_cdm_single(loc_pv, ch_mass, Npart);
    }

    /**
     * Converts chained velocities/position to chain-ordered individual positions and velocities
     * Does not rescale to the center of mass
     * @param[out] loc_pv Individual positions or velocities
     * @param[in] ch_pv Chained positions or velocities
     */
    void to_individual_chain_ordered_single_nocdm(double3 *loc_pv, double3 *ch_pv) {

        // Assign 0 to the first element
        loc_pv[0] = 0;

        // Reconstruct all the elements
        for (size_t i = 0; i < Nchain; i++) {
            loc_pv[i + 1] = loc_pv[i] + ch_pv[i];
        }
    }

    /**
     * Converts chained positions and velocities to individual ones with the original order
     */
    void update_from_chain_to_com() {

        // Assign 0 to the first element
        size_t individual_ind_next;
        size_t individual_ind = chain[0];
        pos[individual_ind] = 0;
        vel[individual_ind] = 0;

        // Reconstruct all the elements
        for (size_t i = 0; i < Nchain; i++) {
            individual_ind = chain[i];
            individual_ind_next = chain[i + 1];

            pos[individual_ind_next] = pos[individual_ind] + ch_pos[i];
            vel[individual_ind_next] = vel[individual_ind] + ch_vel[i];
        }

        if (Config.wExt) {
            pcom = ch_pos[Nchain];
            vcom = ch_vel[Nchain];
        }

        for (size_t i = 0; i < Npart; i++) {
            size_t i_original = chain[i];
            xdata[i_original].eloss = ch_xdata[i].eloss;
            if constexpr (withspin) spin[i_original] = ch_spin[i];
            if (Config.wMassEvol) mass[i_original] = ch_mass[i];
        }

        // Evaluate new cdm
        Nbodyalgorithms::scale_to_cdm(pos, vel, mass, Npart);
    }

    /**
     * Converts chained positions and velocities to individual ones with the original order
     * The only difference from update_from_chain_to_com is that it copies the entire xdata vector, not just changing quantities
     */
    void copy_from_chain_to_com() {

        // Assign 0 to the first element
        size_t individual_ind_next;
        size_t individual_ind = chain[0];
        pos[individual_ind] = 0;
        vel[individual_ind] = 0;

        // Reconstruct all the elements
        for (size_t i = 0; i < Nchain; i++) {
            individual_ind = chain[i];
            individual_ind_next = chain[i + 1];

            pos[individual_ind_next] = pos[individual_ind] + ch_pos[i];
            vel[individual_ind_next] = vel[individual_ind] + ch_vel[i];
        }
        if (Config.wExt) {
            pcom = ch_pos[Nchain];
            vcom = ch_vel[Nchain];
        }

        for (size_t i = 0; i < Npart; i++) {
            size_t i_original = chain[i];
            xdata[i_original] = ch_xdata[i];
            if constexpr (withspin) spin[i_original] = ch_spin[i];
            mass[i_original] = ch_mass[i];
            radius[i_original] = ch_radius[i];
        }

        // Evaluate new cdm
        Nbodyalgorithms::scale_to_cdm(pos, vel, mass, Npart);
    }


    /**
     * Compares the current chain \link Chain::chain \endlink with the old one \link Chain::chain_old \endlink
     * @return
     */
    bool chain_has_changed() {

        // Start to inspect the pairs... if I find a different pair then I must switch the chain
        if (chain_old[0] != chain[0] &&
            chain_old[0] != chain[Nchain]) // One of the extremes of the chain has changed... so I need to switch
            return true;
        else {
            if (chain_old[0] == chain[0]) {
                for (size_t i = 1; i <= Nchain; i++) {
                    if (chain_old[i] != chain[i])
                        return true;
                }
            } else if (chain_old[0] == chain[Nchain]) {
                for (size_t i = 1; i <= Nchain; i++) {
                    if (chain_old[i] != chain[Nchain - i])
                        return true;
                }
                // If chain is the same but in inverse order, we preserve the old one.
                chain = chain_old;
            }
        }
        return false;
    }

    /**
     * Transforms chain position and velocity arrays (\link Chain::ch_pos \endlink and \link Chain::ch_vel \endlink)
     * from old chain (\link Chain::chain_old \endlink) to the current one (\link Chain::chain \endlink)
     * @param[in,out] mass Mass array in original order
     * @param[in,out] soft Softening array in original order
     * @param[in,out] rad  Radii array in original order
     * @param[in,out] xdata Extra particle array in orginal order
     * @param[in,out] spin Spin particle array in orginal order
     */
    void to_new_chain() {

        size_t k1new, k2new; //new indexes in the new chain

        // Save old chain arrays
        for (size_t i = 0; i < Nchain; i++) {
            ch_pos_old[i] = ch_pos[i];
            ch_vel_old[i] = ch_vel[i];
        }

        // Cycle all over the chain vectors
        for (size_t i = 0; i < Nchain; i++) {
            //the new pair
            size_t p1 = chain[i];
            size_t p2 = chain[i + 1];

            // Search the new pair in the previous chain
            k1new = (size_t) (find(chain_old.begin(), chain_old.end(), p1) - chain_old.begin());
            k2new = (size_t) (find(chain_old.begin(), chain_old.end(), p2) - chain_old.begin());

            ch_pos[i] = 0;
            ch_vel[i] = 0;

            if (k2new < k1new) {
                for (size_t j = k2new; j < k1new; j++) {
                    ch_pos[i] -= ch_pos_old[j];
                    ch_vel[i] -= ch_vel_old[j];
                }
            } else {
                for (size_t j = k1new; j < k2new; j++) {
                    ch_pos[i] += ch_pos_old[j];
                    ch_vel[i] += ch_vel_old[j];
                }
            }
            // Update other vectors
            ch_mass[i] = mass[p1];
            ch_radius[i] = radius[p1];
            ch_xdata[i] = xdata[p1];
            if constexpr (withspin) ch_spin[i] = spin[p1];

        }
        size_t p1 = chain[Nchain];
        ch_mass[Nchain] = mass[p1];
        ch_radius[Nchain] = radius[p1];
        ch_xdata[Nchain] = xdata[p1];
        if constexpr (withspin) ch_spin[Nchain] = spin[p1];

        for (size_t i = 0; i < Npart; i++) {
            size_t i_original = chain[i];
            invchain[i_original] = i;
        }
    }


    /**
     * Prints to screen the chain indices
     */
    void print_chain() {
        std::cout << "Chain: " << std::endl;
        for (size_t i = 0; i < Npart; i++) {
            std::cout << chain[i] << " ";
        }
        std::cout << std::endl;
        /*std::cout<<"Vorder: "<<std::endl;
        for (size_t i = 0; i < Npart; i++) {
            std::cout<<ch_mass[i]<<" "<<ch_xdata[i].stype<<"   ";
        }
        std::cout<<std::endl;*/
    }

    /**
     * Prints to screen the chained vectors and chain indices
     */
    void print_chained_vectors() {
        for (size_t i = 0; i < Nchain; i++) {
            std::cout << chain[i] << "_" << chain[i + 1] << ":   " << ch_pos[i].x << "   " << ch_pos[i].y << "   "
                      << ch_pos[i].z << "   "
                      << ch_vel[i].x << "   " << ch_vel[i].y << "   " << ch_vel[i].z << std::endl;
        }
        std::cout << std::endl;
    }

    ///////////////////////////////////////////////////////////////
    ////                         CHAIN                          ///
    ///////////////////////////////////////////////////////////////

    /**
     * Compute and returns kinetic energy of the chained system
     * @param[in] ch_vel Chained velocities
     * @return T
     */
    double calc_T(double3 *ch_vel, double3 *unch_vel) {

        auto val = double(0);

        to_individual_chain_ordered_single_cdm(unch_vel, ch_vel);
        //TODO We can avoid this. Compute mj * mi * vij / 2 when computing potential
        // SEE Binney & Tremaine p.795 D34

        for (size_t i = 0; i < Npart; i++) {
            val += 0.5 * ch_mass[i] * (unch_vel[i] * unch_vel[i]);
        }

        return val;
    }

    /**
     *
     * Compute non-velocity dependent accelerations, dOMEGA, U and omega from chained positions
     * @param[in] ch_pos
     * @param[out] unch_pos
     * @param[out] ch_acc
     * @param[out] ch_acc_ext
     * @param[out] unch_acc_ext
     * @param[out] U
     * @param[out] OMEGA
     * @param[out] dOMEGA
     */
    void calc_accelerations_nvdep(double3 *ch_pos, double3 *unch_pos,
                                  double3 *ch_acc, double3 *ch_acc_ext,
                                  double3 *unch_acc_ext,
                                  double &U, double &OMEGA, double3 *dOMEGA) {
        U = OMEGA = 0;

        // Reset CoM accelerations
        for (size_t i = 0; i < Npart; i++) {
            unchained_acc[i] = 0;
            if (Config.TTL) dOMEGA[i] = 0;
            if (Config.wExt) unch_acc_ext[i] = 0;
        }
        cache_ind = 0; // Radius storage

        // Chained contribution
        twobody_interaction(1, 0, ch_pos[0],
                            ch_mass, ch_radius,
                            unchained_acc, dOMEGA, U,
                            OMEGA); // First chained

        for (size_t k = 1; k < Nchain; k++) {
            twobody_interaction(k + 1, k, ch_pos[k],
                                ch_mass, ch_radius,
                                unchained_acc, dOMEGA, U,
                                OMEGA);

            double3 tmppos = ch_pos[k];
            for (size_t j = 1; j < k + 1; j++) {
                tmppos += ch_pos[k - j];
                twobody_interaction(k + 1, k - j, tmppos,
                                    ch_mass, ch_radius,
                                    unchained_acc, dOMEGA, U,
                                    OMEGA);
            }
        }

        /* OLD IMPLEMENTATION //TODO MEASURE ACCURACY LOSS/SPEED GAIN
        for (size_t k = 2; k < Chained->Npart; k++) { // Other chained
            twobody_interaction(k, k - 1, ch_pos[k - 1], ch_vel[k - 1],
                                Chained->ch_mass, Chained->ordered_soft, Chained->ch_radius,
                                unchained_acc, this->dOMEGA, U,
                                OMEGA);
            twobody_interaction(k, k - 2, ch_pos[k - 2] + ch_pos[k - 1], ch_vel[k - 2] + ch_vel[k - 1],
                                Chained->ch_mass, Chained->ordered_soft, Chained->ch_radius,
                                unchained_acc, this->dOMEGA, U,
                                OMEGA);
        }

        Chained->to_individual_chain_ordered_nocdm(unchained_pos, unchained_vel, ch_pos, ch_vel);

        // Non-chained contributions
        for (size_t k = 2; k < Chained->Npart; k++) {
            for (size_t j = 0; j < k - 2; j++) {
                twobody_interaction(k, j, unchained_pos[k] - unchained_pos[j], unchained_vel[k] - unchained_vel[j],
                                    Chained->ch_mass, Chained->ordered_soft, Chained->ch_radius,
                                    unchained_acc, this->dOMEGA, U,
                                    OMEGA);
            }
        }*/

        // Convert to chained accelerations
        for (size_t i = 0; i < Nchain; i++) {
            ch_acc[i] = unchained_acc[i + 1] - unchained_acc[i];
        }

        /// External forces, not velocity dependent
        if (Config.wExt) {
            to_individual_chain_ordered_single_cdm(unch_pos, ch_pos);

            // Acceleration from external force
            for (size_t i = 0; i < Npart; i++) {
                unch_pos[i] = ch_pos[Nchain] + unch_pos[i]; // Adding CoM position

                // Evaluate acceleration on particles
                unch_acc_ext[i] = 0.0;
            }

            // Convert to chained accelerations
            for (size_t i = 0; i < Nchain; i++) {
                ch_acc_ext[i] = unch_acc_ext[i + 1] - unch_acc_ext[i];
            }
            // Evaluate acceleration on CoM
            ch_acc_ext[Nchain] = 0.0;
        }
    }

    /**
     * Evaluate external velocity-dependent forces, storing the current unchained velocities into aux_vel.
     * aux_vel is unchained_vel when advancing the fake velocities (pre-midstep) and unchained_vel_fake when advancing
     * the real velocities (at midstep).
     * unchained_vel_fake is then used to evaluate dB and domega at the midstep by calc_domega_dB_midpoint
     * @param ch_pos
     * @param unch_pos
     * @param ch_vel
     * @param unch_vel
     * @param ch_acc_vdep
     * @param unch_acc_vdep
     * @param ch_acc_ext_vdep
     * @param unch_acc_ext_vdep
     */
    template<bool real_step>
    void calc_accelerations_vdep(double3 *ch_pos, double3 *unch_pos,
                                 double3 *ch_vel, double3 *unch_vel,
                                 double3 *ch_acc_vdep, double3 *unch_acc_vdep,
                                 double3 *ch_acc_ext_vdep, double3 *unch_acc_ext_vdep,
                                 double3 *spin, double3 *torque) {

        for (size_t i = 0; i < Npart; i++) {
            unch_acc_vdep[i] = 0;
            if (Config.wExt and Config.wExt_vdep) unch_acc_ext_vdep[i] = 0;
            if constexpr (withspin) torque[i] = 0;
        }
        cache_ind = 0;  // Radius storage

        // unchained_pos was calculated in evaluate_acc_dOMEGA_U_OMEGA, no need to calculate again
        //Chained->to_individual_chain_ordered_cdm(unchained_pos, unchained_vel, ch_pos, ch_vel);
        // Should be if(dissipative). Only dissipative forces need B advancement
        //if (Config.wPNs) to_individual_chain_ordered_single_cdm(unch_vel, ch_vel); // Need COM velocities for PN forces
        //else if (Config.wEqTides or Config.wDynTides) to_individual_chain_ordered_single_nocdm(unch_vel, ch_vel); // Don't need the COM ones, let's save some cpu cicles
        to_individual_chain_ordered_single_cdm(unch_vel, ch_vel);

        /// Chained contribution
        twobody_interaction_vdep<real_step>(1, 0, ch_pos[0], ch_vel[0], ch_mass,
                                            ch_radius, ch_xdata,
                                            unch_acc_vdep, spin, torque);

        cache_ind++;
        for (size_t k = 1; k < Nchain; k++) {
            twobody_interaction_vdep<real_step>(k + 1, k, ch_pos[k], ch_vel[k],
                                                ch_mass,
                                                ch_radius, ch_xdata,
                                                unch_acc_vdep, spin, torque);
            cache_ind++;

            double3 tmppos = ch_pos[k];
            double3 tmpvel = ch_vel[k];
            for (size_t j = 1; j < k + 1; j++) {
                tmppos += ch_pos[k - j];
                tmpvel += ch_vel[k - j];
                twobody_interaction_vdep<real_step>(k + 1, k - j, tmppos, tmpvel,
                                                    ch_mass,
                                                    ch_radius, ch_xdata,
                                                    unch_acc_vdep, spin, torque);
                cache_ind++;

            }
        }

        /*
        twobody_interaction(1, 0, ch_pos[0], ch_vel[0],
                            Chained->ch_mass, Chained->ordered_soft, Chained->ch_radius,
                            unchained_acc, this->dOMEGA, U,
                            OMEGA); // First chained

        // Chained contributions
        for (size_t k = 2; k < Chained->Npart; k++) { // Other chained
            if (usePN)
                twobody_interaction_pn(k, k - 1, ch_pos[k - 1], ch_vel[k - 1], aux_vel,
                                       Chained->ch_mass, unchained_acc_ext_vdep);
            if (useEqTides)
                twobody_interaction_equilibrium_tide(k, k - 1, ch_pos[k - 1], ch_vel[k - 1],
                                                     Chained->ch_mass,
                                                     Chained->ch_radius, Chained->ch_xdata,
                                                     unchained_acc_ext_vdep);
            ++cache_ind;
            if (usePN)
                twobody_interaction_pn(k, k - 2, ch_pos[k - 2] + ch_pos[k - 1], ch_vel[k - 2] + ch_vel[k - 1],
                                       aux_vel,
                                       Chained->ch_mass, unchained_acc_ext_vdep);

            if (useEqTides)
                twobody_interaction_equilibrium_tide(k, k - 2, ch_pos[k - 2] + ch_pos[k - 1],
                                                     ch_vel[k - 2] + ch_vel[k - 1],
                                                     Chained->ch_mass,
                                                     Chained->ch_radius, Chained->ch_xdata,
                                                     unchained_acc_ext_vdep);
            ++cache_ind;
        }

        // Non-chained contributions
        for (size_t k = 2; k < Chained->Npart; k++) {
            for (size_t j = 0; j < k - 2; j++) {
                if (usePN)
                    twobody_interaction_pn(k, j, unchained_pos[k] - unchained_pos[j],
                                           aux_vel[k] - aux_vel[j], aux_vel,
                                           Chained->ch_mass, unchained_acc_ext_vdep);
                if (useEqTides)
                    twobody_interaction_equilibrium_tide(k, j, unchained_pos[k] - unchained_pos[j],
                                                         aux_vel[k] - aux_vel[j],
                                                         Chained->ch_mass,
                                                         Chained->ch_radius, Chained->ch_xdata,
                                                         unchained_acc_ext_vdep);
                ++cache_ind;
            }
        }
        */

        // Return chained accelerations
        for (size_t i = 0; i < Nchain; i++)
            ch_acc_vdep[i] = unch_acc_vdep[i + 1] - unch_acc_vdep[i];

        /// External forces, velocity dependent
        if (Config.wExt and Config.wExt_vdep) {
            to_individual_chain_ordered_single_cdm(unch_vel, ch_vel);

            // Acceleration from external force
            for (size_t i = 0; i < Npart; i++) {
                unch_vel[i] = ch_vel[Nchain] + unch_vel[i]; // Adding CoM velocity
                unch_pos[i] = unch_pos[i]; // May use unch_pos (which is already wrt to CoM)

                // Evaluate acceleration on particles
                unch_acc_ext_vdep[i] = 0.0;
            }

            // Convert to chained accelerations
            for (size_t i = 0; i < Nchain; i++) {
                ch_acc_ext_vdep[i] = unch_acc_ext_vdep[i + 1] - unch_acc_ext_vdep[i];
            }
            // Evaluate acceleration on CoM
            ch_acc_ext_vdep[Nchain] = 0.0;
        }

    }

    void twobody_interaction(size_t k, size_t j, double3 dpos, double *mass, double *rad,
                             double3 *acc, double3 *dOMEGA,
                             double &U, double &OMEGA) {

        double &radj = rad[j];
        double &radk = rad[k];
        double &mj = mass[j];
        double &mk = mass[k];
        double dr = sqrt(dpos * dpos);
        double inv_dr = 1.0 / dr;
        double inv_dr3 = inv_dr * inv_dr * inv_dr;
        double3 inv3dr = inv_dr3 * dpos;
        double m_prod = mj * mk;
        double3 increment;

        // Radius storage
        cache_r[cache_ind] = dr;
        cache_invr[cache_ind++] = inv_dr;

        // Update U
        U += m_prod * inv_dr;

        /** If I use TTL only, m_prod < tolerance*mean_mp2 might never be true, so OMEGA = 0
         * To avoid this I need to specify that if I am using TTL alone the condition
         * m_prod < ttl_tol*mean_mp2 should be always true, ttl_tol = 1e33
         * For TTL only I also use omega_ij = m2ttl.
         */
        if (Config.TTL) {
            double omega_ij = m2ttl; // Maybe also 1 for "vanilla" Mikkola & Aarseth 2001 TTL
            //TODO omega_ij should be zero if mj+mj = 0

            // Update OMEGA and dOMEGA
            if (m_prod < ttl_tol * m2ttl) {
                OMEGA += omega_ij * inv_dr;

                increment = omega_ij * inv3dr;

                dOMEGA[k] -= increment;
                dOMEGA[j] += increment;
            }
        }
        // Update accelerations
        acc[k] -= mj * inv3dr;
        acc[j] += mk * inv3dr;

        if (lastit and dr <= Config.dcoll * (radj + radk)) {
//            std::cout << "#INCOLL i,j " << chain[k] << "," << chain[j] << std::endl;
//            std::cout << "#INCOLL dr/sumR " << dr / (radj + radk) << std::endl;
//            std::cout << "#INCOLL ds      " << time_fict << std::endl;
            CollInfo.collflag = Nbodyalgorithms::COLLISION_REAL;
            CollInfo.collind[0] = chain[k];
            CollInfo.collind[1] = chain[j];
            CollInfo.colltime = time_phys;
            CollInfo.colltime_ds = time_fict;
            return;
        }
    }

    template<bool real_step>
    void twobody_interaction_vdep(size_t k, size_t j, const double3 &dr, const double3 &dv, double *mass,
                                  double *rad, pinfo *xdata, double3 *acc_vdep,
                                  double3 *spin, double3 *torque) {

        /// Common variables
        double &mj = mass[j];
        double &mk = mass[k];

        double &radj = rad[j];
        double &radk = rad[k];
        
        double3 &sj = spin[j];
        double3 &sk = spin[k];
        double &r_inv = cache_invr[cache_ind];
        double m_prod = mj * mk;
        double mtot = mj + mk;
        double v2 = dv * dv;

        ////////////////////////////////////////////////////////////////////////////
        /////////                      PN TERMS                            /////////
        ////////////////////////////////////////////////////////////////////////////
        if (Config.wPNs) {
            if (v2 > vel_light2) {
                throw IntegrationError(ErrType::V_GTR_C);
            } else {
                double3 rhat = dr * r_inv;
                double rdot = dv * rhat;
                double v4 = v2 * v2;
                double v6 = v4 * v2;
                double r_inv2 = r_inv * r_inv;
                double mu = m_prod / mtot;
                double nu = mu/mtot;
                double nu2 = nu * nu;
                double Ut = mtot * r_inv;
                double Ut2 = Ut * Ut;
                double rdot2 = rdot * rdot;
                double rdot4 = rdot2 * rdot2;
                double rdot6 = rdot4 * rdot2;
                double Utnu = Ut * nu;
                double rdot2v2 = rdot2 * v2;
                double v2nu = v2 * nu;

                double A1, B1, A2, B2, A25, B25, A3, B3, A35, B35;

                B35 = B3 = B25 = B2 = B1 = A35 = A3 = A25 = A2 = A1 = 0.0;

                if (Config.pn1) {
                    A1 = coeff1pn * (v2 * (1 + 3 * nu) - 1.5 * rdot2 * nu - 2 * Ut * (2 + nu));
                    B1 = 2 * coeff1pn * rdot * (nu - 2);
                }

                if (Config.pn2) {
                    A2 = coeff2pn * (
                            nu * (1.875 * rdot4 * (1 - 3 * nu) + 4 * v2 * (1.5 * rdot2 - v2) * (nu - 0.75))
                            + Ut * (2 * (nu * v2 * (nu - 3.25) - rdot2 * (1 + 12.5 * nu + nu2))
                                    + Ut * (9 + 21.75 * nu))
                    );

                    B2 = coeff2pn * rdot * (
                            nu * (3 * rdot2 * (1.5 + nu) - v2 * (7.5 + 2 * nu))
                            + Ut * (2 + 20.5 * nu + 4 * nu2)
                    );
                }

                if (Config.pn25) {
                    double combo25 = coeff25pn * Utnu;
                    A25 = -combo25 * rdot * (4.8 * v2 + c136_15 * Ut);
                    B25 = 1.6 * combo25 * (v2 + 3 * Ut);
                }

                if (Config.pn3) {
                    double A3t1, A3t2, A3t3, A3t4, A3t5, A3t6, A3t7, A3t8, A3t9, A3t10;
                    double B3t1, B3t2, B3t3, B3t4, B3t5, B3t6;

                    A3t1 = 2.1875 * rdot6 * (5 * (nu - nu2) - 1);
                    A3t2 = rdot2 * (1 - 4.5 * nu + 4.25 * nu2);
                    A3t3 = v2 * (-1 + 3.95 * nu - 3 * nu2);
                    A3t4 = v6 * (2.75 - 12.25 * nu + 13 * nu2);
                    A3t5 = rdot4 * (79 - 34.5 * nu - 30 * nu2);
                    A3t6 = rdot2 * v2 * (-121 + 16 * nu + 20 * nu2);
                    A3t7 = v4 * (18.75 + 8 * nu - 10 * nu2);
                    A3t8 = rdot2 * (1 + nu * (pn3const1 + 1.375 * nu - 7 * nu2));
                    A3t9 = v2 * nu * (nu2 - pn3const2);
                    A3t10 = -16 - pn3const3 * nu - 35.5 * nu2;

                    A3 = coeff3pn * (
                            nu * (A3t1 + 7.5 * rdot2 * v2 * (A3t2 + A3t3) + A3t4)
                            + Utnu * (A3t5 + A3t6 + A3t7)
                            + Ut2 * (A3t8 + A3t9 + Ut * A3t10)
                    );

                    B3t1 = 15 * rdot4 * (-0.375 + nu + 0.25 * nu2);
                    B3t2 = rdot2 * v2 * (12 - 27.75 * nu - 12 * nu2);
                    B3t3 = v4 * (-8.125 + 19 * nu + 6 * nu2);
                    B3t4 = rdot2 * (c329_6 + 29.5 * nu + 18 * nu2);
                    B3t5 = -v2 * (15 + 27 * nu + 10 * nu2);
                    B3t6 = -4 + nu * (pn3const4 + 25 * nu + 8 * nu2);


                    B3 = coeff3pn * rdot * (
                            nu * (B3t1 + B3t2 + B3t3)
                            + Utnu * (B3t4 + B3t5)
                            + Ut2 * B3t6
                    );
                }

                if (Config.pn35) {
                    double combo35 = coeff35pn * Utnu;
                    A35 = combo35 * rdot * (
                            v4 * (c366_35 + 12 * nu) - v2 * rdot2 * (114 + 12 * nu) + 112 * rdot4
                            + Ut * (v2 * (c692_35 - c724_15 * nu) + rdot2 * (58.8 + 75.2 * nu))
                            + Ut2 * (c3956_35 + 36.8 * nu)
                    );

                    B35 = combo35 * (
                            -v4 * (c626_35 + 2.4 * nu) + v2 * rdot2 * (135.6 + 2.4 * nu) - 120 * rdot4
                            + Ut * (v2 * (c164_21 + 29.6 * nu) - rdot2 * (c82_3 + c848_15 * nu))
                            - Ut2 * (c1060_21 + 20.8 * nu)

                    );
                }

                double3 acc = ((r_inv2 * (A35 + A3 + A25 + A2 + A1)) * rhat +
                                (r_inv2 * (B35 + B3 + B25 + B2 + B1)) * dv);

                if constexpr (withspin) {
                    // spin equations will be defined in terms of "dimensional" spins
                    // because it requires fewer flops; for reference S = G*m*m*chi/c
                    // where S is the dimensional spin, m is the mass, c is the speed of light, and chi is the dimensionless spin
                    double3 a1ss, a1so, a35ss, a35so;
                    double3 t1ssj, t1ssk, t1soj, t1sok, t35ssj, t35ssk;
                    double r_inv3 = r_inv2 * r_inv;
                    double r_inv4 = r_inv2 * r_inv2;
                    double r_inv5 = r_inv3 * r_inv2;
                    double sjsk = sj * sk;
                    double sjn = sj * rhat;
                    double skn = sk * rhat;
                    double sjnskn = sjn * skn;
                    double sjv = sj * dv;
                    double skv = sk * dv;
                    double sjvskv = sjv * skv;
                    double sjkvn = (sjn * skv + skn * sjv);
                    double mk_mj = mk/mj;
                    double mj_mk = 1/mk_mj;
                    double3 sjkn = (sj * skn + sk * sjn);
                    double3 Ln = dr.cross(dv);
                    double3 stot = sj + sk;
                    double3 xi = (mk_mj) * sj + (mj_mk) * sk;

                    if (Config.pn1 and Config.ss) {
                        a1ss = -3 * r_inv4 / mu * (rhat * sjsk + sjkn - 5 * rhat * sjnskn);
                        double3 term1j = (sk - (3 * skn) * rhat);
                        double3 term1k = (sj - (3 * sjn) * rhat);
                        t1ssj = -r_inv3 * term1j.cross(sj);
                        t1ssk = -r_inv3 * term1k.cross(sk);
                    }
                    if (Config.pn1 and Config.so) {
                        double3 term1so = 4 * stot + 3 * xi;
                        a1so = r_inv3 * (
                                (1.5 * r_inv * (Ln * term1so)) * rhat
                                - dv.cross(term1so)
                                + (1.5 * rdot) * rhat.cross(term1so)
                        );

                        double murinv3 = mu * r_inv3;
                        t1soj = (murinv3 * (2 + 1.5 * mk_mj)) * Ln.cross(sj);
                        t1sok = (murinv3 * (2 + 1.5 * mj_mk)) * Ln.cross(sk);

                    }

                    if (Config.pn35 and Config.ss) {
                        double scalar1 = rdot * ((287 * rdot2 - 99 * v2 + 108.2 * Ut) * sjsk
                                                 - (2646 * rdot2 - 714 * v2 + 392.2 * Ut) * sjnskn
                                                 - 336 * sjvskv)
                                         + (1029 * rdot2 - 123 * v2 + 62.9 * Ut) * sjkvn;

                        double scalar2 = (34.2 * v2 - 195 * rdot2 - 67 * Ut) * sjsk
                                         - (174 * v2 - 1386 * rdot2 - 207.6 * Ut) * sjnskn
                                         - 438 * rdot * sjkvn
                                         + 96 * sjvskv;

                        a35ss = r_inv5 * (rhat * scalar1 + dv * scalar2
                                          + (2.7 * v2 - 37.5 * rdot2 - c509_30 * Ut) * (skv * sj + sjv * sk)
                                          + (7.5 * v2 + 38.5 * rdot2 + 19.9 * Ut) * rdot * sjkn
                        );

                        double mr5 = m_prod * r_inv5;
                        double rdot30 = 30 * rdot;
                        t35ssj = mr5 * (c2_3 * skv + rdot30 * skn) * rhat.cross(sj);
                        t35ssk = mr5 * (c2_3 * sjv + rdot30 * sjn) * rhat.cross(sk);

                    }

                    if (Config.pn35 and Config.so) {
                        double LnS = Ln * stot;
                        double Lnxi = Ln * xi;

                        double scalar1 = rdot * r_inv * (
                                (120 * v2 + 280 * rdot2 + 453 * Ut) * LnS
                                + (285 * v2 - 245 * rdot2 + 211 * Ut) * Lnxi
                        );

                        double scalar2 = r_inv * (
                                (87 * v2 - 675 * rdot2 - c901_3 * Ut) * LnS
                                + 4 * (6 * v2 - 75 * rdot2 - 41 * Ut) * Lnxi
                        );

                        double3 term35so1 = -rdot * (c2_3 * (48 * v2 + 15 * rdot2 + 364 * Ut) * stot
                                                     + c1_3 * (375 * v2 - 195 * rdot2 + 640 * Ut) * xi);

                        double v2Ut = v2 * Ut;
                        double rdot2Ut = rdot2 * Ut;
                        double3 term35so2 = .5 * (
                                stot *
                                (31 * v4 - 260 * rdot2v2 + 245 * rdot4 - c689_3 * v2Ut + 537 * rdot2Ut + c4_3 * Ut2)
                                - xi * (29 * v4 - 40 * rdot2v2 - 245 * rdot4 + 211 * v2Ut - 1019 * rdot2Ut - 80 * Ut2)
                        );

                        a35so = -0.2 * mu * r_inv4 * (
                                rhat * scalar1 + dv * scalar2 + dv.cross(term35so1) + rhat.cross(term35so2)
                        );
                    }
                    // dividing by total mass to make it the same units as the linear PN terms,
                    // this is necessary to transform out of the relative reference frame,
                    // applying the acceleration in the global reference frame to each body

                    double3 acc_spin = (coeff1pn * (a1ss + a1so) + coeff35pn * (a35ss + a35so)) / mtot;

                    acc += acc_spin;
                    double3 tj = coeff1pn*(t1ssj + t1soj) + coeff35pn*t35ssj;
                    double3 tk = coeff1pn*(t1ssk + t1sok) + coeff35pn*t35ssk;

                    torque[j] += tj;
                    torque[k] += tk;
                }
                
                acc_vdep[k] -= mj * acc;
                acc_vdep[j] += mk * acc;
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        /////////                   EQUILIBRIUM TIDES                      /////////
        ////////////////////////////////////////////////////////////////////////////
        if (Config.wEqTides) {
            // continue only if something has tides
            pinfo &xdatak = xdata[k];
            pinfo &xdataj = xdata[j];

            if (xdatak.hastide or xdataj.hastide) {

                // Auxiliary variables
                double3 rhat = dr * r_inv;
                double vrad = dv * rhat;
                double invr7 = r_inv * r_inv * r_inv * r_inv * r_inv * r_inv * r_inv;
                double invr8 = invr7 * r_inv;
                double mj_inv = 1.0 / mj;
                double mk_inv = 1.0 / mk;
                // Useful quantities when we have spin
                double3 dissfact;
                if constexpr (withspin) {
                    dissfact = 2*vrad*rhat;
                } else {
                    dissfact = 3*vrad*rhat;
                }

                double3 totforce = 0.0;

                if (xdatak.hastide) {
                    // Dissipation in body k
                    double &Ak = xdatak.Atide;
                    double &sigmak = xdatak.sigmadiss;
                    double mj2 = mj * mj;
                    double Ak2 = Ak * Ak;
                    double nondissk = -6*mj2*Ak*invr7;
                    double dissk = - 4.5 * sigmak * mj2 * Ak2 * invr8;
                    double3 forcek = dissk * dissfact;
                    totforce += forcek;
                    if constexpr (real_step) {
                        xdatak.eloss += mk * (forcek * dv) * dt;
                    }

                    if constexpr (withspin) {
                        double3 &spink = spin[k];
                        double3 orbcouplk = dissk * (dv - spink.cross(dr));
                        totforce += orbcouplk;
                        double3 torqk = - dr.cross(orbcouplk); //(dissk * r) * (r * spink - (dr * spink) * rhat - rhat_cross_v);
                        torque[k] += torqk / xdatak.inert;
                        if constexpr (real_step) {
                            xdatak.eloss += (orbcouplk * dv + spink * torqk) * dt;
                        }
                    }
                    // Non-diss factors
                    totforce += nondissk * rhat;
                }

                if (xdataj.hastide) {
                    // Dissipation in body j
                    double &Aj = xdataj.Atide;
                    double &sigmaj = xdataj.sigmadiss;
                    double mk2 = mk * mk;
                    double Aj2 = Aj * Aj;
                    double nondissj = -6*mk2*Aj*invr7;
                    double dissj = - 4.5 * sigmaj * mk2 * Aj2 * invr8;
                    double3 forcej = dissj * dissfact;
                    totforce += forcej;
                    if constexpr (real_step) {
                        xdataj.eloss += mj * (forcej * dv) * dt;
                    }

                    if constexpr (withspin) {
                        double3 &spinj = spin[j];
                        double3 orbcouplj = dissj * (dv - spinj.cross(dr));
                        totforce += orbcouplj;
                        double3 torqj = - dr.cross(orbcouplj); //(dissk * r) * (r * spink - (dr * spink) * rhat - rhat_cross_v);
                        torque[j] += torqj / xdataj.inert;
                        if constexpr (real_step) {
                            xdataj.eloss += (orbcouplj * dv + spinj * torqj) * dt;
                        }
                    }
                    // Non-diss factors
                    totforce += nondissj * rhat;
                }

                double3 tideforcej = -mj_inv * totforce;
                double3 tideforcek = mk_inv * totforce;

                acc_vdep[j] += tideforcej;
                acc_vdep[k] += tideforcek;
            }
        }
        ////////////////////////////////////////////////////////////////////////////
        /////////                   DYNAMICAL TIDES                        /////////
        ////////////////////////////////////////////////////////////////////////////
        if (Config.wDynTides) {
            // continue only if something has tides
            if (xdata[k].hastide | xdata[j].hastide) {
                // Orbit
                double vcirc2 = mtot * r_inv;
                double a = -mtot / (v2 - 2. * vcirc2);
                double vdiff = v2 - vcirc2;
                double vdotr = dr * dv;
                double mtot_inv = 1. / mtot;
                double3 e_vec = mtot_inv * (vdiff * dr - vdotr * dv);
                double e2 = e_vec * e_vec;
                double e = sqrt(e2); // eccentricity
                double peri = a * (1. - e);

                // Press & Taukolsky 1977
                double p_over_rj = peri / radj;
                double p_over_rk = peri / radk;
                double rj_over_peri = 1.0 / p_over_rj;
                double rk_over_peri = 1.0 / p_over_rk;
                double etaj = sqrt(mj * mtot_inv * p_over_rj * p_over_rj * p_over_rj);
                double etak = sqrt(mk * mtot_inv * p_over_rk * p_over_rk * p_over_rk);

                // Mardling & Aarseth 2001 (A1-A4)
                double mard = 0.5 * fabs(etaj - 2.);
                mard = (0.5 + 0.25 * sqrt(mard * mard * mard));
                double newetaj = etaj * pow(2.0 / (1.0 + e), mard);

                mard = 0.5 * fabs(etak - 2.);
                mard = (0.5 + 0.25 * sqrt(mard * mard * mard));
                double newetak = etak * pow(2. / (1.0 + e), mard);

                double Etid = 0;

                // TIDE ON j
                if ((newetaj > 0) & (newetaj < 10)) {
                    Etid += Etid_polynomial_fits(mk, rj_over_peri, peri, xdata[j], newetaj);
                }

                // TIDE ON k
                if ((newetak > 0) & (newetak < 10)) {
                    Etid += Etid_polynomial_fits(mj, rk_over_peri, peri, xdata[k], newetak);
                }
                double Estar;
                size_t tideind;

                if (xdata[k].hastide) {
                    Estar = mk * mk / radk;
                    tideind = k;
                }
                if (xdata[j].hastide) {
                    Estar = mj * mj / radj;
                    tideind = j;
                }

                double ratio2 = Etid / Estar;
                double ratio3 = peri / rad[tideind];

                if ((ratio2 > 10) and (ratio3 < 1)) {
                    CollInfo.collflag = Nbodyalgorithms::COLLISION_TDE;
                    CollInfo.collind[0] = chain[k];
                    CollInfo.collind[1] = chain[j];
                }

                //double sumrad = radk+radj;
                //double ratio1 = r12 / sumrad;
                //std::cout << "\nr1 / sumrad  " << ratio1 << std::endl;
                //std::cout << "Etid / Estar " << ratio2 << std::endl;
                //std::cout << "peri / rad   " << peri / rad[tideind] << std::endl;

                // Samsing, Leigh & Trani 2018, using n=4
                double Ien = 0.5 * M_PI * (2. + 7. * e2 + e2 * e2);
                double pterm = a * (1. - e2);
                double epsilon =
                        0.5 * Etid * sqrt(mtot_inv * pterm * pterm * pterm * pterm * pterm * pterm * pterm) / Ien;

                double3 r4inv = (epsilon * r_inv * r_inv * r_inv * r_inv) * dv;


                double3 pn25j = (1.0 / mj) * r4inv;
                double3 pn25k = (-1.0 / mk) * r4inv;

                acc_vdep[j] += pn25j;
                acc_vdep[k] += pn25k;
            }
        }

        ////////////////////////////////////////////////////////////////////////////
        /////////                   COLLISIONS                             /////////
        ////////////////////////////////////////////////////////////////////////////

        double &r = cache_r[cache_ind];
        if constexpr (usetdetracker and real_step) {

            if (lastit) {
                bool istark = xdata[k].stype == ptype::HIGH_MASS_MS;
                bool istarj = xdata[j].stype == ptype::HIGH_MASS_MS;
                bool isbhk = xdata[k].stype == ptype::BH;
                bool isbhj = xdata[j].stype == ptype::BH;

                if ((istarj and isbhk) or (istark and isbhj)) {
                    size_t istar, ibh;

                    if (istarj and isbhk) {
                        istar = j;
                        ibh = k;
                    } else {
                        istar = k;
                        ibh = j;
                    }

                    Nbodyalgorithms::PTDE_history &PTDEHist = PTDEInfo.PTDEHistories[chain[istar]];

                    double rt = pow(mass[ibh] / mass[istar], 0.333333333333333333333) * rad[istar];
                    double Psi = (1.47 + exp((mass[istar] - 0.669) / 0.137)) / (1. + 2.34 * exp((mass[istar] - 0.669) / 0.137));
                    Psi *= (0.80 + 0.26 * sqrt((mass[ibh] * Config.Mscale * 1e-6)));
                    //std::cout << Psi << std::endl;
                    //std::cout << dr << "  " << rt* Psi << std::endl;
                    //std::cout << dr / (rt* Psi) << std::endl;

                    double Rt_full = rt * Psi;
                    if (r <= Rt_full) {
                        CollInfo.collflag = Nbodyalgorithms::COLLISION_TDE;
                        CollInfo.collind[0] = istar;
                        CollInfo.collind[1] = ibh;
                        PTDEHist.ptde_candidate = PTDEInfo.ptde_candidate = false;
                        PTDEHist.last_PTDE = Nbodyalgorithms::PTDE_event();
                    } else {
                        if (r <= 2 * rt) {
                            PTDEHist.ptde_candidate = true;
                            PTDEHist.last_PTDE.rmin = std::min(PTDEHist.last_PTDE.rmin, r);
                            double vcirc2 = mtot * r_inv;
                            PTDEHist.last_PTDE.semi = -mtot / (v2 - 2. * vcirc2); // semimajor axis
                            double3 L = dr.cross(dv);
                            double E = 0.5 * v2 - vcirc2;
                            double e = sqrt(1 + 2*E*L*L / (mtot*mtot));

                            PTDEHist.last_PTDE.ecc = e; // eccentricity
                            PTDEHist.last_PTDE.ibh = chain[ibh];

                        } else if (r > 2 * rt and PTDEHist.ptde_candidate) {
                            PTDEInfo.had_ptde = true;
                            PTDEHist.ptde_candidate = PTDEInfo.ptde_candidate = false;
                            PTDEHist.PTDE_list.push_back(PTDEHist.last_PTDE);
                            /*std::cout << "PTDE:" << std::endl;
                            std::cout << "rmin: " << last_PTDE.rmin << std::endl;
                            std::cout << "semi: " << last_PTDE.semi << std::endl;
                            std::cout << "ecc: " << last_PTDE.ecc << std::endl;
                            std::cout << "mbh: " << last_PTDE.bh_mass << std::endl;*/

                            PTDEHist.last_PTDE = Nbodyalgorithms::PTDE_event();
                        }
                    }
                } else if (r <= (radj + radk)) {
                    CollInfo.collflag = Nbodyalgorithms::COLLISION_REAL;
                    CollInfo.collind[0] = chain[k];
                    CollInfo.collind[1] = chain[j];
                    CollInfo.colltime = time_phys;
                    CollInfo.colltime_ds = time_fict;
                    return;
                }
            }
        }
    }


    double Etid_polynomial_fits(double mpert, double r_peri, double peri, pinfo &xdata, double eta) {
        double fA, fB, fC, fD, fE, fF;
        if (xdata.polyt == 3) {
            fA = fA_n3, fB = fB_n3, fC = fC_n3, fD = fD_n3, fE = fE_n3, fF = fF_n3;
        } else if (xdata.polyt == 1.5) {
            fA = fA_n15, fB = fB_n15, fC = fC_n15, fD = fD_n15, fE = fE_n15, fF = fF_n15;
        } else {
            return 0.0;
        }

        double etid = r_peri * r_peri * r_peri * r_peri * r_peri * mpert * mpert / peri;
        // Now computing Teta
        if (eta < 1) { // most likely a disruptive encounter
#ifdef __linux__
            etid *= exp10(fA);
#else
            etid *= pow(fA, 10);
#endif
        } else if (eta < 10) {
            double let = log10(eta);
#ifdef __linux__
            etid *= exp10(fA + let * (fB + let * (fC + let * (fD + let * (fE + let * fF)))));
#else
            etid *= pow(fA + let * (fB + let * (fC + let * (fD + let * (fE + let * fF)))), 10);
#endif
        }
        return etid;
    }

    TsunamiConfig &Config;
    Nbodyalgorithms::CollisionLog CollInfo;
    Nbodyalgorithms::PTDELog PTDEInfo;

    std::vector<size_t> chain; ///< Chain vector: contains the index of particles in the current chain
    std::vector<size_t> chain_old; ///< Previous chain vector: used to compare the current chain and update
    std::vector<size_t> invchain;

    double3 *pos, *vel; ///< Original positions, velocities
    double3 *spin = nullptr; ///< Original spin
    double *mass, *radius; ///< Mass and radii arrays with the unchained particle index
    pinfo *xdata; ///< Particle info arrays with the unchained particle index
    double3 pcom, vcom; ///< Center of mass position and velocity

    double3 *ch_pos, *ch_vel; ///< Chained positions and velocities
    double3 *ch_pos_old, *ch_vel_old; ///< Previous chained positions and velocities
    double *ch_mass, *ch_radius; ///< Mass and radii arrays with the chained particle index
    double3 *ch_spin = nullptr; ///< Chain-ordered spin
    pinfo *ch_xdata; ///< Particle info arrays with the chained particle index

    /// Auxiliary vector
    double3 *unchained_acc;

    double *cache_r; // Avoid calculating radii again
    double *cache_invr;
    size_t cache_ind = 0;

    // Regularization utilities
    double ttl_tol;
    double m2ttl;
    double time0, time_phys, time_fict;
    double dt;

    // PN utilities
    double vel_light, vel_light2;
    double coeff35pn, coeff3pn, coeff25pn, coeff2pn, coeff1pn;
    static constexpr double pn3const1 = 615.0 / 64 * M_PI * M_PI + 22717.0 / 168;
    static constexpr double pn3const2 = 123.0 / 64 * M_PI * M_PI + 20827.0 / 840;
    static constexpr double pn3const3 = -41.0 / 16 * M_PI * M_PI + 1399.0 / 12;
    static constexpr double pn3const4 = -123.0 / 32 * M_PI * M_PI - 5849.0 / 840;
    static constexpr double c136_15 = 136.0 / 15.0;
    static constexpr double c329_6 = 329.0 / 6.0;
    static constexpr double c366_35 = 366.0 / 35.0;
    static constexpr double c692_35 = 692.0 / 35.0;
    static constexpr double c724_15 = 724.0 / 15.0;
    static constexpr double c3956_35 = 3956.0 / 35.0;
    static constexpr double c626_35 = 626.0 / 35.0;
    static constexpr double c164_21 = 164.0 / 21.0;
    static constexpr double c82_3 = 82.0 / 3.0;
    static constexpr double c848_15 = 848.0 / 15.0;
    static constexpr double c1060_21 = 1060.0 / 21.0;
    static constexpr double c509_30 = 509./30.;
    static constexpr double c2_3 = 2./3.;
    static constexpr double c1_3 = 1./3.;
    static constexpr double c4_3 = 4./3.;
    static constexpr double c901_3 = 901./3.;
    static constexpr double c689_3 = 689./3.;

    // Tidal fits from Portegies Zwart et al. 1993
    double fA_n15 = -0.397;
    double fB_n15 = 1.678;
    double fC_n15 = 1.277;
    double fD_n15 = -12.42;
    double fE_n15 = 9.446;
    double fF_n15 = -5.550;

    double fA_n3 = -1.124;
    double fB_n3 = 0.877;
    double fC_n3 = -13.37;
    double fD_n3 = 21.55;
    double fE_n3 = -16.8;
    double fF_n3 = 4.124;

    /** Auxiliary structure to collect all interparticle information to find the chain */
    struct ParticlePair {
        double dist;
        size_t i;
        size_t j;
        bool used;
    };

    double get_distance_between_collided() {
        size_t chid1 = invchain[CollInfo.collind[0]];
        size_t chid2 = invchain[CollInfo.collind[1]];
        // k chain member connects k and k+1, get the lowest
        make_ascending(chid1, chid2);

        double3 dpos = ch_pos[chid1];
        for (size_t ip = chid1+1; ip < chid2; ip++) {
            dpos += ch_pos[ip];
        }

        return dpos.mod();

        // Check if really adjacent
//        if (chid1 + 1 != chid2) {
//            throw TsuError("CONGRATULATIONS!! You found an edge case that Alessandro predicted, but was too lazy to implement,\n"
//                           "(it's never going to happen, he thought), so he added this error message.\n"
//                           "Now, send him this message: 'FOUND COLLISION NOT CHAIN-ADJACENT'\n"
//                           "and he's going to implement this edge case within one working day.");
//        } else {
//            return ch_pos[chid1].mod();
//        }
    }

    ParticlePair max_pair;
    ParticlePair max_pair_smallest {KeplerUtils::dinf};
    ParticlePair min_pair;
    ParticlePair min_pair_smallest {KeplerUtils::dinf};


    /**
      * Auxiliary function to sort the particles' distances to construct the chain
      * @param a		particle a
       * @param b 	particle b
      * @return 		a.dist < b.dist
      */
    static bool sort_func(ParticlePair a, ParticlePair b) { return (a.dist < b.dist); }

    /** Vector of particle pairs */
    std::vector<ParticlePair> PairList;

    size_t Npart = 0; /**< Number of particles*/
    size_t Nchain; /**< Number of chain vectors(Npart -1)*/
    bool lastit = true;

    /** Allocates the chain arrays and chain-ordered vectors */
    void allocate_arrays() {

        // Adding one extra particle to the chain if we have external potentials
        auto extraP = static_cast<size_t>(Config.wExt);

        // Chained quantities
        ch_pos = new double3[Nchain + extraP];
        ch_vel = new double3[Nchain + extraP];
        ch_pos_old = new double3[Nchain + extraP];
        ch_vel_old = new double3[Nchain + extraP];
        ch_mass = new double[Npart];
        ch_radius = new double[Npart];
        ch_xdata = new pinfo[Npart];
        if constexpr (withspin) ch_spin = new double3[Npart];

        // Original quantities
        pos = new double3[Npart];
        vel = new double3[Npart];
        mass = new double[Npart];
        radius = new double[Npart];
        if constexpr (withspin) spin = new double3[Npart];
        xdata = new pinfo[Npart];

        cache_r = new double[(Npart * (Npart - 1)) / 2];
        cache_invr = new double[(Npart * (Npart - 1)) / 2];

        // Auxiliary vectors
        unchained_acc = new double3[Npart];
    }

    /** De-allocates the chain arrays and chain-ordered vectors */
    void deallocate_arrays() {
        delete[] ch_pos;
        delete[] ch_vel;
        delete[] ch_pos_old;
        delete[] ch_vel_old;
        delete[] ch_mass;
        delete[] ch_radius;
        delete[] ch_xdata;
        if constexpr (withspin) delete[] ch_spin;

        delete[] pos;
        delete[] vel;
        delete[] mass;
        delete[] radius;
        if constexpr (withspin) delete[] spin;
        delete[] xdata;

        delete[] cache_r;
        delete[] cache_invr;

        delete[] unchained_acc;
    }

};

typedef Chain<TsunamiConfig::wSpins, TsunamiConfig::useTDEtracker, TsunamiConfig::debug_ch> ChainSys;

#endif //TSUNAMI_CHAIN_HPP
