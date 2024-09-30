//
// Created by aatrani on 13/06/19.
//

#ifndef TSUNAMI_KEPLERUTILS_H
#define TSUNAMI_KEPLERUTILS_H

#include <limits>
#include <cmath>

class KeplerUtils {

public:

    //KeplerUtils() = default;
    KeplerUtils() :  KeplerUtils(1, 1) {}
    KeplerUtils(double Mscale, double Lscale) : Mscale(Mscale), Lscale(Lscale) {
        // Time unit (in yr)
        Tscale = sqrt(Lscale*Lscale*Lscale/(Mscale*G_yr_msun_au));

        // Velocity unit (in km/s)
        Vscale = Lscale/Tscale * au2km / yr2sec;
        speed_of_light = c_ms * 1e-3 / Vscale;
    }

    ~KeplerUtils()= default;

    static double period(double a, double mu);
    static double mod2pi(double f);
    static double M_to_E(double M, double e);
    static double E_to_nu(double E, double e);
    static double M_to_nu(double M, double e) {return E_to_nu(M_to_E(M, e), e);};
    static double nu_to_E(double nu, double e);
    static double E_to_M(double E, double e);
    static double nu_to_M(double nu, double e) {return E_to_M(nu_to_E(nu, e), e);};
    template<typename T> static int sgn(T val);

    /// Deprecated
    static void cart_to_kepl(double out_keplpar[6], double in_posvel[6], double mu);
    static void kepl_to_cart(double out_posvel[6], double prim_posvel[6], double mu, double a, double e, double i, double Ome, double ome, double nu);
    ///

    void cart_to_kepl(double out_keplpar[6], double in_pos[3], double in_vel[3], double m1, double m2,
                      bool pn1=false, bool pn2=false, bool pn3=false, bool memmcorr=false);
    void kepl_to_cart(double out_posvel[6], double prim_posvel[6], double m1, double m2, double a, double e, double i,
                      double ome, double Ome, double nu, bool pn1=false, bool pn2=false);

    void kepl_to_cart2(double out_posvel[6], double prim_posvel[6], double m1, double m2, double a, double e, double i,
                      double ome, double Ome, double nu, bool pn1=false, bool pn2=false, bool pn3=false);


    static double acos2(double num, double den, double posflag);

    static constexpr double dtiny = std::numeric_limits<double>::min();
    static constexpr double dhuge = std::numeric_limits<double>::max();
    static constexpr double dinf = std::numeric_limits<double>::infinity();
    static constexpr double min_e = 1e-7;
    static constexpr double min_e2 = 1e-15;
    static constexpr double min_i = 1e-7;
    static constexpr double M_2PI = 2.*M_PI;
    static constexpr double G_s_kg_m = 6.6743e-11; ///< G in s^-2, kg^-1, m^3, CODATA2018
    static constexpr double au2km = 1.495978707e+8; ///< Nominal IAU2012
    static constexpr double yr2sec = 3.15576e+7; ///< Nominal
    static constexpr double c_ms = 2.99792458e+8; ///< Nominal
    static constexpr double GMSun = 1.3271244e+20; ///< Nominal solar mass parameter, m^3 s^-2, IAU2015
    static constexpr double GMEarth = 3.986004e+14; ///< Nominal Earth mass parameter, m^3 s^-2, IAU2015
    static constexpr double GMJup = 1.2668653e+17; ///< Nominal Jupiter mass parameter, m^3 s^-2, IAU2015
    static constexpr double RSun_km = 6.957e+5; ///< Nominal, IAU2015
    static constexpr double REarth_m = 6.3781e+6; ///< Nominal equatorial value, IAU2015
    static constexpr double RJup_m = 7.1492e+7; ///< Nominal equatorial value, IAU2015
    static constexpr double pc2au = (648000/M_PI); ///< Nominal, IAU2015
    static constexpr double pc2km = au2km * pc2au; ///< Nominal, IAU2015

    // Derived quantities
    static constexpr double RSun2au = RSun_km / au2km;
    static constexpr double RJup2au = RJup_m / au2km * 1e-3;
    static constexpr double REarth2au = REarth_m / au2km * 1e-3;
    static constexpr double MSun_kg = GMSun / G_s_kg_m;
    static constexpr double MEarth_kg = GMEarth / G_s_kg_m;
    static constexpr double MJup_kg = GMJup / G_s_kg_m;

    static constexpr double MJup2MSun = MJup_kg / MSun_kg;
    static constexpr double MEarth2MSun = MEarth_kg / MSun_kg;
    static constexpr double G_yr_msun_au = G_s_kg_m * ((yr2sec*yr2sec) * MSun_kg / (au2km*au2km*au2km) * 1e-9); ///< G in yr^-2, msun^-1, au^3

    double Tscale; ///< Tscale in yr
    double Mscale; ///< Mscale in MSun
    double Lscale; ///< Lscale in au
    double Vscale; ///< Vscale in km/s
    double speed_of_light; ///< c in N-body units


private:

    static inline double ecc_to_mean(double E, double M, double e) {
        return (E - e * sin(E) - M);
    }
    static inline double ecc_to_mean_prime(double E, double e) {
        return (1. - e * cos(E));
    }
    static inline double ecc_to_mean_h(double E, double M, double e) {
        return (E - e * sinh(E) + M);
    }
    static inline double ecc_to_mean_h_prime(double E, double e) {
        return (1. - e * cosh(E));
    }

    static constexpr int maxit = 100;
};


#endif //TSUNAMI_KEPLERUTILS_H
