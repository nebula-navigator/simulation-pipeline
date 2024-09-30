//
// Created by lex on 9/3/21.
//

#ifndef TSUNAMI_STOPCOND_H
#define TSUNAMI_STOPCOND_H

#include "custom_types.hpp"
#include <array>
#include <memory>
#include "classification.h"

class StopCond {

public:
    StopCond() = default;

    virtual ~StopCond() = default;

    virtual bool check_stopping_condition(double3 *pos, double3 *vel, double *mass, double *rad, double time) {return false;};

    template <typename T, size_t FIRST, size_t SECOND> using array2d = std::array<std::array<T, SECOND>, FIRST>;

    enum StopType {
        NONE = 0,
        TRIPLE = 1,
        TDE = 2,
    } StopCondType = NONE;
    static StopCond *create (StopType type);
};

/// /// /// /// /// /// /// /// ///
/// HIRARCHICAL TRIPLES STOPCOND
/// /// /// /// /// /// /// /// ///
class StopCond_Triple : public StopCond {

public:
    StopCond_Triple() {
        id_not_in_pair[0][1] = id_not_in_pair[1][0] = 2;
        id_not_in_pair[0][2] = id_not_in_pair[2][0] = 1;
        id_not_in_pair[1][2] = id_not_in_pair[2][1] = 0;
    }
    ~StopCond_Triple() override = default;

    array2d<size_t, 3,3> id_not_in_pair;

    bool check_stopping_condition(double3 *pos, double3 *vel, double *mass, double *rad, double time) override;

    bool is_escaping = false;
    double breakup_time = 0.0;
    double instability_time = 0.0;
    size_t escape_id;
    TripleClass::PairOrbit FinalBinary;

};

/// /// /// /// /// /// /// /// ///
/// Tidal disruptions stopcond
/// /// /// /// /// /// /// /// ///
class StopCond_TripleTDEs : public StopCond {

public:
    StopCond_TripleTDEs() = default;
    ~StopCond_TripleTDEs() override = default;

    array2d<double, 3,3> full_tde_radius_matrix;
    array2d<double, 3,3> partial_tde_radius_matrix;
    size_t star_id;

    bool had_ptde = false;
    bool is_ftde = false;
    bool is_dc = false;

    int id_bh = -1;

    bool check_stopping_condition(double3 *pos, double3 *vel, double *mass, double *rad, double time) override;
    void set_tidal_radii(size_t i, size_t j, double full_tde_radius, double partial_tde_radius);
    void set_star_id(size_t id);
};

#endif //TSUNAMI_STOPCOND_H
