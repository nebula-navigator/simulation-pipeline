//
// Created by lex on 10/8/19.
//

#ifndef TSUNAMI_KOZAILIDOV_H
#define TSUNAMI_KOZAILIDOV_H

#include <array>
#include <cmath>
#include "DOP853.hpp"
#include "keplerutils.h"

class Okinami {

public:
	Okinami() : Okinami(1, 1) {}
	Okinami(double Mscale, double Lscale) {
		this->Mscale = Mscale;
		this->Lscale = Lscale;

        // Time unit (in yr)
        Tscale = sqrt(Lscale*Lscale*Lscale/(Mscale*KeplerUtils::G_yr_msun_au));

        // Velocity unit (in km/s)
        Vscale = Lscale/Tscale * KeplerUtils::au2km / KeplerUtils::yr2sec;
        speed_of_light = KeplerUtils::c_ms * 1e-3 / Vscale;

        speed_of_light2 = speed_of_light * speed_of_light;
        coeff1pn = 1.0 / speed_of_light2;
        coeff2pn = 1.0 / (speed_of_light2 * speed_of_light2);
        coeff25pn = 1.0 / (speed_of_light2 * speed_of_light2 * speed_of_light);

		system = std::bind(
                &Okinami::compute_derivatives, this, std::placeholders::_1, std::placeholders::_2,
                std::placeholders::_3);

        Integrator.bind_system(system);
        set_integrator_tolerance();
    }

	~Okinami() = default;

    /// e1, e2, g1, g2, h1, a1, H, a2, m1, m2, m3, R1, R2, R3
    typedef std::array<double, 14> delaunay_state;
    DOP853<delaunay_state, 14> Integrator;

    static double imutual_delaunay(double m1, double  m2, double m3,
									  double a1, double a2,
									  double e1, double e2,
									  double i1, double i2,
									  double Ome1, double Ome2,
									  double &inc1, double &inc2, double &H);

    static void i1_i2_from_itot(double m1, double m2, double m3,
                                double a1, double a2,
                                double e1, double e2,
                                double i_mut,
                                double &inc1, double &inc2);

    static double compute_H(double m1, double m2, double m3,
                            double a1, double a2,
                            double e1, double e2,
                            double i_mut);
    static void compute_H_in_y(delaunay_state &y, double i_mut);

    static constexpr double yr_in_s = 3.1556925993600003e+07;
	static constexpr double G_yr_msun_pc = 4.4994505613511768e-15; // G in yr^-2, msun^-1, parsec^3

	bool octupole = true;
	bool gr = false;
	bool tides = false;
	bool ext = false;
	bool stopgw = false;
    bool stopsec = false;
    bool check_gw = false;
    bool check_sec = false;
    bool pn1 = true;
    bool pn2 = true;
    bool pn25 = true;
    double fgwstop = 0.0;
	short int check_stability = 0;
	bool check_collisions = true;
	bool isdynstable = true;
    bool hascollided = false;
    double dcoll = 1;

    double ctime = 0.0;
    double dt=0.0;
    double dt_done=0.0;

    enum id {
        e1 = 0,
        e2, g1, g2, h1, a1, H, a2, m1, m2, m3, R1, R2, R3
    };

    enum stability_criterion {
        NONE = 0,
        MARDLING = 1,
        PETROVICH = 2,
        AMD = 3
    };

    void get_derivative(const delaunay_state &y_const, delaunay_state &dydt_dev, const double t);
    void compute_derivatives(const delaunay_state &y_const, delaunay_state &dydt_dev, const double t);
	void initialize_constants(double m1, double m2, double a2, double m3);
	void stop_at_gw_frequency(double gwstop);

    static bool is_mardling_stable(double a1, double m1, double m2,
                                   double a2, double m3, double e2, double i_mut);
    bool is_petrovich_stable(double a1, double m1, double m2,
                             double a2, double m3, double e2, double i_mut);
    static bool is_amd_stable(double a1, double a2, double e1, double e2,
                              double m1, double m2, double m3,
                              double H);
    void set_external_derivatives(const delaunay_state &yderiv);
    double dHdt_from_external_changes(double m1, double m2, double m3, double cositot,
                                      double G1, double G2, double H,
                                      double a1, double e1, double ome21,
                                      double a2, double e2, double ome22,
                                      const Okinami::delaunay_state &yderiv);
    void set_tidal_parameters(double k1, double tau1, double k2, double tau2);

    static double get_octupole_strength(double m1, double m2, double a1, double a2, double e2);
    static double get_zkl_timescale(double m1, double m2, double a1, double m3, double a2, double e2);
    static double itot_from_y(const delaunay_state &y_const);

    void evolve_system(delaunay_state &y, double tfin);
    bool stop_function(const delaunay_state &x, double t);
    void set_integrator_tolerance(double rel_err=1e-7, double abs_err=1e-8, double ax=1.0, double adxdt=1.0);

	std::function<void(const delaunay_state&, delaunay_state&, double)> system;

	double e1max = 0, e1min = 1;
	double i_mut;

    /// Failsafes
	bool zero_inclination = false;
	static constexpr double cosi_min_thr = 1e-15;
    bool zero_e1 = false;
    static constexpr double e1_min_thr = 1e-8;
    bool zero_e2 = false;
    static constexpr double e2_min_thr = 1e-8;

	static constexpr double amd3to43 = 4.3267487109222245;
    double Mscale, Lscale, Tscale, Vscale;
    double speed_of_light;

protected:
    double speed_of_light2;
	double coeff25pn, coeff1pn, coeff2pn;

    void check_dynamical_stability(const delaunay_state &y);
    void check_particle_collisions(const delaunay_state &y);
    void check_gw_frequency(const delaunay_state &y);
    void check_quasi_secular(const delaunay_state &y);

    /// Auxiliary functions and variables
	double L2, L1_til, C2_til, C3_til, PN25_fac_a1, PN25_fac_e1, PN25_fac_G1, PN1_fac_g1, PN2_fac_g1, PN15_fac; // Constants

    delaunay_state yderiv_ext = { 0.0 };
    double k1 = 0, tau1 = 0, k2 = 0, tau2 = 0;

	static inline double G_moment(double mx, double my, double a, double e);
	static inline double L_moment(double m1, double m2, double a);
	static inline double L_tilde(double m1, double m2);
	static inline double C2_tilde(double m1, double m2, double m3, double a2);
	static inline double C3_tilde(double m1, double m2, double m3, double a2);
	inline double PN25_factor_a(double m1, double m2);
	inline double PN25_factor_e(double m1, double m2);
	inline double PN25_factor_l(double m1, double m2);
	inline double PN1_factor(double m1, double m2);
    inline double PN2_factor(double m1, double m2);
    inline static double f1(double e2);
    inline static double f2(double e2);
    inline static double f3(double e2);
    inline static double f4(double e2);
    inline static double f5(double e2);
    inline static double pseudosync_spin(double e2, double omegaorb);
};


#endif //TSUNAMI_KOZAILIDOV_H
