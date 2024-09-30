//
// Created by lex on 9/3/21.
//

#include "stopcond.h"


StopCond *StopCond::create(StopType type) {
    switch (type) {
        case NONE: return new StopCond();
        case TRIPLE: return new StopCond_Triple();
        case TDE: return new StopCond_TripleTDEs();
        default: return new StopCond();
    }
}


bool StopCond_Triple::check_stopping_condition(double3 *pos, double3 *vel, double *mass, double *rad, double time) {

    std::vector<size_t> iBound;
    TripleClass::PairOrbit Pairs[3];

    size_t kbound = 0;
    for (size_t i = 0; i < 3; i++) {
        for (size_t j = 0; j < i; j++) {
            Pairs[kbound].semimajor_eccentricity_vdotr(i, j, pos[i] - pos[j], vel[i] - vel[j], mass[i] + mass[j]);

            // Checks
            if (Pairs[kbound].a > 0) {
                Pairs[kbound].get_pair_com(pos[i], pos[j], vel[i], vel[j], mass[i], mass[j]);

                iBound.push_back(kbound);
            }
            kbound++;
        }
    }

    if (not iBound.empty()) {
        TripleClass::PairOrbit PairInner;

        // More than two pairs are bound, system is interacting
        if (iBound.size() > 1) {
            // Get the Harder binary
            PairInner = Pairs[iBound[0]];
            for (size_t ib = 1; ib < iBound.size(); ib++) {
                if (Pairs[iBound[ib]].a < PairInner.a) {
                    PairInner = Pairs[iBound[ib]];
                }
            }
        } else if (iBound.size() == 1) {
            PairInner = Pairs[iBound[0]];
        }

        // Here we have identified the most bound binary, PairInner
        size_t k = id_not_in_pair[PairInner.i][PairInner.j]; // This is the outer object index
        double3 dv = PairInner.vcom - vel[k];
        double3 dp = PairInner.pcom - pos[k];

        // PairOuter is the outer binary
        TripleClass::PairOrbit PairOuter(k, 0, dv,
                                         dv,
                                         PairInner.mcom + mass[k]);

        if (PairOuter.a < 0.0) {
            // Outer object is unbound

            if (PairOuter.vdotr > 0.0) {
                if (not is_escaping) {
                    breakup_time = time;
                    is_escaping = true;
                    //std::cout << "Time = "<< time << " - now escaping" << std::endl;
                } else {
                    double3 dout = PairInner.pcom - pos[k];
                    double dist2 = dout*dout;
                    //std::cout << "Time = "<< time << " - still escaping" << std::endl;
                    //std::cout << " d/a = "<< sqrt(dist2)/PairInner.a << " - still escaping" << std::endl;

                    // When outer body is 20 abin far from the binary, stop
                    if (dist2 > 400*PairInner.a*PairInner.a) {
                        // System has broken
                        breakup_time = time;
                        escape_id = k;
                        FinalBinary = PairInner;
                        return true;
                    }
                }
            }
        } else {
            // Outer object is bound
            is_escaping = false;

            if (instability_time == 0.0) {
                double tidal_factor = TripleClass::compute_ffactor(PairInner, PairOuter, pos, mass);
                if(tidal_factor > 1) {
                    instability_time = time;
                }

            }
            // Check for tidal factor
            // ftid = tidal_force / keplerian_force
            // if ftide < 1 ok, continue
            // if ftide > 1, record instability_time = time
            // Needs to happen only once
        }


    } else {
        // No bound pairs, system is interacting because total E<0
        is_escaping = false;
    }

    // If simulation needs to continue
    return false;
}

bool StopCond_TripleTDEs::check_stopping_condition(double3 *pos, double3 *vel, double *mass, double *rad, double time) {
    for (size_t i=0; i<3; i++) {
        if (i == star_id) continue;

        // Record PTDEs : R_t < R < 2 r_t
        double3 dp = pos[i] - pos[star_id];
        double dist = sqrt(dp * dp);
        std::cout << "dist/ptde " << dist/partial_tde_radius_matrix[i][star_id] << std::endl;
        std::cout << "dist/ftde " << dist/full_tde_radius_matrix[i][star_id] << std::endl;

        if (dist < full_tde_radius_matrix[i][star_id]) {
            double3 dv = vel[i] - vel[star_id];
            double mu = mass[star_id] + mass[i];
            TripleClass::PairOrbit Pair(star_id, i, dp, dv, mu);
            double peri = Pair.a * (1 - Pair.e);

            // Rad here is Schw radius
            if (peri < rad[i]) {
                is_dc = true;
            } else {
                is_ftde = true;
            }
            id_bh = static_cast<int>(i);
            return true;
        } else if (not had_ptde and dist < partial_tde_radius_matrix[i][star_id]) {
            had_ptde = true;
        }
    }
    return false;
}

void StopCond_TripleTDEs::set_tidal_radii(size_t i, size_t j, double full_tde_radius, double partial_tde_radius) {
    full_tde_radius_matrix[i][j] = full_tde_radius_matrix[j][i] = full_tde_radius;
    partial_tde_radius_matrix[i][j] = partial_tde_radius_matrix[j][i] = partial_tde_radius;
}

void StopCond_TripleTDEs::set_star_id(size_t id) {
    star_id = id;
}
