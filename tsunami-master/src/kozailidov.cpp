//
// Created by lex on 10/8/19.
//

#include "kozailidov.h"
#include "keplerutils.h"
#include "errhand.hpp"
#include <array>
#include <functional>

double Okinami::G_moment(double mx, double my, double a, double e) {
	double mtot = mx + my;
	double ome2 = 1 - e*e;
	double L = mx * my * sqrt(a * ome2 / mtot);
	return L;
}
double Okinami::L_moment(double mx, double my, double a) {
	double mtot = mx + my;
	double L = mx * my * sqrt(a/mtot); //*G**0.5
	return L;
}
double Okinami::L_tilde(double m1, double m2) {
	double Lti = m1*m2/sqrt(m1+m2); //*G**0.5
	return Lti;
}
double Okinami::C2_tilde(double m1, double m2, double m3, double a2) {
	double C2ti = m1*m2*m3/16; //*constants.G
	C2ti = C2ti / ( (m1+m2)*a2*a2*a2 );
	return C2ti;
}
double Okinami::C3_tilde(double m1, double m2, double m3, double a2) {
	double C3ti = -15*m1*m2*m3*(m1-m2); //*constants.G
	C3ti = C3ti/( 64*((m1+m2)*(m1+m2))*a2*a2*a2*a2 );
	return C3ti;
}
double Okinami::PN25_factor_a(double m1, double m2) {
	double PN25_fa = -64*m1*m2*(m1+m2)*coeff25pn/5; //*G*G*G
	return PN25_fa;
}
double Okinami::PN25_factor_e(double m1, double m2) {
	double PN25_fa = -304*m1*m2*(m1+m2)*coeff25pn/15; //*G*G*G
	return PN25_fa;
}
double Okinami::PN25_factor_l(double m1, double m2) {
	double PN25_fa = -32*m1*m1*m2*m2*sqrt(m1+m2)*coeff25pn/5; //*G^7/2
	return PN25_fa;
}
double Okinami::PN1_factor(double m1, double m2) {
	double mu_3 = (m1+m2); //*G
	mu_3 = mu_3*mu_3*mu_3;
	double PN1_fa = 3*coeff1pn*sqrt(mu_3);
	return PN1_fa;
}
double Okinami::PN2_factor(double m1, double m2) {
    double mu_5 = (m1+m2); //*G
    mu_5 = mu_5*mu_5*mu_5*mu_5*mu_5;
    double PN2_fa = 3*coeff2pn*sqrt(mu_5)/4;
    return PN2_fa;
}
double Okinami::f1(double e2) {
    return 1+e2*(15.5+e2*(31.875+e2*(11.5625+e2*0.390625)));
}
double Okinami::f2(double e2) {
    return 1.0+e2*(7.5+e2*(5.625+e2*0.3125));
}
double Okinami::f3(double e2) {
    return 1.0+e2*(3.75+e2*(1.875+e2*7.8125e-02));
}
double Okinami::f4(double e2) {
    return 1+e2*(1.5+e2*0.125);
}
double Okinami::f5(double e2) {
    return 1.0+e2*(3.+e2*0.375);
}
double Okinami::pseudosync_spin(double e2, double mmean) {
    double syncspin = f2(e2) / (f5(e2) *  pow(1.0 - e2, 1.5)) * mmean;
    return syncspin;
}

void Okinami::set_tidal_parameters(double k1, double tau1, double k2, double tau2) {
    this->k1 = k1;
    this->k2 = k2;
    this->tau1 = tau1;
    this->tau2 = tau2;
}

double Okinami::imutual_delaunay(double m1, double m2, double m3, double a1, double a2, double e1, double e2,
                                 double i1, double i2, double Ome1, double Ome2,
                                 double &inc1, double &inc2, double &H) {
	double G1 = G_moment(m1, m2, a1, e1);
	double G2 = G_moment(m1+m2, m3, a2, e2);
	double H1 = G1*cos(i1);
	double H2 = G2*cos(i2);

	double K1 = sqrt(G1 * G1 - H1 * H1);
	double K2 = sqrt(G2 * G2 - H2 * H2);

	double Lx = K1 * sin(Ome1) + K2 * sin(Ome2);
	double Ly = -K1 * cos(Ome1) - K2 * cos(Ome2);
	double Lz = H1 + H2;
	double Ltot = sqrt(Lx * Lx + Ly * Ly + Lz * Lz);

	double itot = (Ltot * Ltot - G1 * G1 - G2 * G2) / (2 * G1 * G2);
	itot = acos(itot);

	double H1_x = (Ltot * Ltot + G1 * G1 - G2 * G2) / (2 * Ltot);
	double H2_x = (Ltot * Ltot + G2 * G2 - G1 * G1) / (2 * Ltot);

	inc1 = acos(H1_x / G1);
	inc2 = acos(H2_x / G2);
	H = Ltot;
	return itot;
}

double Okinami::compute_H(double m1, double m2, double m3,
                          double a1, double a2,
                          double e1, double e2,
                          double i_mut) {
    double G1 = G_moment(m1, m2, a1, e1);
    double G2 = G_moment(m1 + m2, m3, a2, e2);
    double H = sqrt(G1*G1 + G2*G2 + 2*G1*G2*cos(i_mut));
    return H;
}


void Okinami::i1_i2_from_itot(double m1, double m2, double m3,
                              double a1, double a2,
                              double e1, double e2,
                              double i_mut,
                              double &inc1, double &inc2) {

    //std::cout <<  cos(i_mut) << std::endl;
    if (1-cos(i_mut) < cosi_min_thr) {
        inc1 = inc2 = 0.0;
        return;
    }

    double G1 = G_moment(m1, m2, a1, e1);
    double G2 = G_moment(m1 + m2, m3, a2, e2);
    double H = sqrt(G1*G1 + G2*G2 + 2*G1*G2*cos(i_mut));

    double cosi1 = (H*H + G1*G1 - G2*G2) / (2*H*G1);
    double cosi2 = (H*H + G2*G2 - G1*G1) / (2*H*G2);

    inc1 = acos(cosi1);
    inc2 = acos(cosi2);
}


void Okinami::compute_H_in_y(delaunay_state &y, double i_mut) {

    y[id::H] = compute_H(y[id::m1], y[id::m2], y[id::m3],
                         y[id::a1], y[id::a2],
                         y[id::e1], y[id::e2],
                         i_mut);
}

void Okinami::initialize_constants(double m1, double m2, double a2, double m3) {
    L1_til = L_tilde(m1, m2);
	L2 = L_moment(m1+m2, m3, a2);
	C2_til = C2_tilde(m1, m2, m3,  a2);
	C3_til = C3_tilde(m1, m2, m3, a2);
	PN25_fac_a1 = PN25_factor_a(m1, m2);
	PN25_fac_e1 = PN25_factor_e(m1, m2);
	PN25_fac_G1 = PN25_factor_l(m1, m2);
	PN1_fac_g1 = PN1_factor(m1, m2);
    PN2_fac_g1 = PN2_factor(m1, m2);
}

void Okinami::set_external_derivatives(const delaunay_state &yderiv) {
    yderiv_ext = yderiv;
    // Do not mess with dHdt
    yderiv_ext[id::H] = 0;
}

double Okinami::dHdt_from_external_changes(double m1, double m2, double m3, double cositot,
                                           double G1, double G2, double H,
                                           double a1, double e1, double ome21,
                                           double a2, double e2, double ome22,
                                           const delaunay_state &yderiv) {
    double minn = m1 + m2;
    double mout = m1 + m2 + m3;

    double dminndt = yderiv[id::m1] + yderiv[id::m2];
    double dmoutdt = dminndt + yderiv[id::m3];

    double dG1dt = yderiv[id::m1]/m1 + yderiv[id::m2]/m2 -0.5*dminndt/minn + 0.5*yderiv[id::a1]/a1 - e1 * yderiv[id::e1] / ome21;
    double dG2dt = yderiv[id::m3]/m3 + dminndt/(m1+m2)   -0.5*dmoutdt/mout + 0.5*yderiv[id::a2]/a2 - e2 * yderiv[id::e2] / ome22;

    dG1dt *= G1;
    dG2dt *= G2;

    double dHdt = (dG1dt * (G1 + G2 * cositot) + dG2dt * (G2 + G1 * cositot))/H;

    return dHdt;
}



double Okinami::itot_from_y(const delaunay_state &y_const) {
	double e1 = y_const[id::e1];
	double e2 = y_const[id::e2];
	double a1 = y_const[id::a1];
	double H  = y_const[id::H];

	double a1_05 = sqrt(a1);
	double e2_2 = e2*e2;
	double ome22 = (1-e2_2);
	double ome22_05 = sqrt(ome22);
	double e1_2 = e1*e1;
	double ome21 = 1-e1_2;
	double ome21_05 = sqrt(ome21);

    double L1_til = L_tilde(y_const[id::m1], y_const[id::m2]);
    double L2 = L_moment(y_const[id::m1] + y_const[id::m2],
                         y_const[id::m3], y_const[id::a2]);

    double G1 = L1_til*a1_05*ome21_05;
	double G2 = L2 * ome22_05;
	double cositot = (H*H - G1*G1 - G2*G2)/(2*G1*G2);
    cositot = (1-cositot > cosi_min_thr) ? cositot : 1;
	return acos(cositot);
}

double Okinami::get_octupole_strength(double m1, double m2, double a1, double a2, double e2) {
    double eps_oct = e2 / (1 - e2*e2) * (a1 / a2);
    eps_oct *= (m1 - m2) / (m1 + m2);
    return eps_oct;
}

double Okinami::get_zkl_timescale(double m1, double m2, double a1, double m3, double a2, double e2) {
    double m_inn = m1 + m2;
    double mtot = m_inn + m3;
    double Pinn = 2 * M_PI * sqrt(a1 * a1 * a1 / m_inn);
    double Pout = 2 * M_PI * sqrt(a2 * a2 * a2 / mtot);
    double tzkl = Pout * Pout / Pinn * mtot / m3 * pow(1 - e2 * e2, 1.5);
    return tzkl;
}

void Okinami::get_derivative(const delaunay_state &y_const,
                             delaunay_state &dydt_dev,
                             const double t) {
    Integrator.return_derivative(y_const, dydt_dev, t);
}

void Okinami::compute_derivatives(const delaunay_state &y_const,
                                  delaunay_state &dydt_dev,
                                  const double t) {
	// e1, e2, g1, g2, h1, a1, H = y_const
	// a2, m1, m2, m3, R1, R2, R3
	double e1 = y_const[id::e1];
	double e2 = y_const[id::e2];
	double g1 = y_const[id::g1];
	double g2 = y_const[id::g2];
	double h1 = y_const[id::h1];
	double a1 = y_const[id::a1];
	double H  = y_const[id::H];
    double a2 = y_const[id::a2];
    double m1 = y_const[id::m1];
    double m2 = y_const[id::m2];
    double m3 = y_const[id::m3];
    double R1 = y_const[id::R1];
    double R2 = y_const[id::R2];
    double R3 = y_const[id::R3];

    /// Initialize non-constant
    if (ext) {
        L1_til = L_tilde(m1, m2);
        L2 = L_moment(m1+m2, m3, a2);
        C2_til = C2_tilde(m1, m2, m3,  a2);
        C3_til = C3_tilde(m1, m2, m3, a2);
    }

	/// Auxiliary variables
	double e2_2 = e2*e2;
	double ome22 = (1-e2_2);
	double ome22_05 = sqrt(ome22);
	//double ome22_3 = ome22*ome22*ome22;
	//double ome22_15 = sqrt(ome22_3);
	double ome22_15 = ome22*ome22_05;
	double C2 = C2_til*a1*a1/ome22_15;

	double a1_05 = sqrt(a1);
	double e1_2 = e1*e1;
	double ome21 = 1-e1_2;
	double ome21_05 = sqrt(ome21);

	/// Moments
	double G1 = L1_til*a1_05*ome21_05;
	double G2 = L2 * ome22_05;

	/// Angles
	double cositot = (H*H - G1*G1 - G2*G2)/(2*G1*G2);

	//cout << "cositot " << cositot << endl;
    /// Check for near-zero inclination
    if((1 - cositot) < cosi_min_thr) {
        cositot = 1;
        zero_inclination = true;
    } else {
        zero_inclination = false;
    }

    /// Check for negative semimajor axis
    if (a1 < 0 and not hascollided) {
        std::cout << "System has collided, a1 < 0" << std::endl;
        hascollided = true;
        for (size_t i=0; i<14; i++) dydt_dev[i] = 0.0;
        return;
    }

    // Stopping integration
    //if (hascollided or isdynstable or stopsec or stopgw) {
    //    for (size_t i=0; i<14; i++) dydt_dev[i] = 0.0;
    //    return;
    //}

    /// Check for zero inner eccentricity
    zero_e1 = e1 < e1_min_thr;
    zero_e2 = e2 < e2_min_thr;

    double cositot_2 = cositot*cositot;
    double sinitot_2 = 1 - cositot_2;
	double sinitot = sqrt(sinitot_2);
	double sin2itot = 2*sinitot*cositot;
	double cosg1 = cos(g1);
	double sing1 = sin(g1);
	double cosg2 = cos(g2);
	double sing2 = sin(g2);
	double cos2g1 = 2*cosg1*cosg1-1; //cos(2 * g1);
	double sin2g1 = 2*sing1*cosg1; //sin(2 * g1);
	double sing1_sing2 = sing1 * sing2;
	double cosphi = -cosg1 * cosg2 - cositot * sing1_sing2; // also theta2

    double C3, B, A; /// Octupole variables
	if(octupole) {
		//double ome22_5 = ome22_3*ome22*ome22;
		//double ome22_25 = sqrt(ome22_5);
		double ome22_25 = ome22*ome22*ome22_05;
		C3 = C3_til*a1*a1*a1/ome22_25;
		B = 2 + e1_2 * (5 - 7 * cos2g1);
		A = 4 + 3 * e1_2 - 2.5 * B * sinitot_2;
	}

	double de1dt = C2 * (1 - e1_2) / G1 * (30 * e1 * sinitot_2 * sin2g1);
	if (octupole) {
		de1dt += C3 * e2 * (1 - e1_2) / G1 * (35 * cosphi * sinitot_2 * e1_2 * sin2g1
		                                      - 10 * cositot * sinitot_2 * cosg1 * sing2 * (1 - e1_2)
		                                      - A * (sing1 * cosg2 - cositot * cosg1 * sing2));
	}

	double de2dt = 0;
	if (octupole) {
		de2dt += -C3 * e1 * (1 - e2_2) / G2 * (10 * cositot * sinitot_2 * (1 - e1_2) * sing1 * cosg2
		                                             + A * (cosg1 * sing2 - cositot * sing1 * cosg2));
	}

	double dg1dt = 6 * C2 * (1.0 / G1 * (4 * cositot_2 + (5 * cos2g1 - 1)
	                                                     * (1 - e1_2 - cositot_2)) +
	                         cositot / G2 * (2 + e1_2 * (3 - 5 * cos2g1)));
	if(octupole and not zero_e1) {
		dg1dt += -C3 * e2 * (e1 * (1.0 / G2 + cositot / G1) *
		                     (sing1_sing2 * (10 * (3 * cositot_2 - 1) * (1 - e1_2) + A)
		                      - 5 * B * cositot * cosphi) - (1 - e1_2) / (e1 * G1)
		                                                                 * (sing1_sing2 * 10 * cositot * sinitot_2 *
		                                                                    (1 - 3 * e1_2)
		                                                                    + cosphi * (3 * A - 10 * cositot_2 + 2)));
	}

	double dg2dt = 3 * C2 * (2 * cositot / G1 * (2 + e1_2 * (3 - 5 * cos2g1))
	                         + 1.0 / G2 * (4 + 6 * e1_2 + (5 * cositot_2 - 3) * (2 + e1_2 * (3 - 5 * cos2g1))));
	if (octupole and not zero_e2) {
		dg2dt += C3 * e1 * (sing1_sing2 * ((4 * e2_2 + 1) / (e2 * G2) * 10 * cositot * sinitot_2 * (1 - e2)
		                                   -
		                                   e2 * (1.0 / G1 + cositot / G2) * (A + 10 * (3 * cositot_2 - 1) * (1 - e1_2)))
		                    + cosphi *
		                      (5 * B * cositot * e2 * (1.0 / G1 + cositot / G2) + (4 * e2_2 + 1) / (e2 * G2) * A));
	}

    double dh1dt = 0;
    if(not zero_inclination) {
        dh1dt += -3 * C2 * H / (G1 * G2 * sinitot) * (2 + 3 * e1_2 - 5 * e1_2 * cos2g1) * sin2itot;
        if (octupole) {
            dh1dt += -C3 * e1 * e2 *
                     (5 * B * cositot * cosphi - A * sing1_sing2 + 10 * (1 - 3 * cositot_2) * (1 - e1_2) * sing1_sing2) *
                     H / (G1 * G2);
        }
    }

	double da1dt = 0;
	double dHdt = 0;
    double da2dt = 0;
    double dm1dt = 0;
    double dm2dt = 0;
    double dm3dt = 0;
    double dR1dt = 0;
    double dR2dt = 0;
    double dR3dt = 0;

	if(gr and not hascollided) {
	    if (ext) {
            PN25_fac_a1 = PN25_factor_a(m1, m2);
            PN25_fac_e1 = PN25_factor_e(m1, m2);
            PN25_fac_G1 = PN25_factor_l(m1, m2);
            PN1_fac_g1 = PN1_factor(m1, m2);
            PN2_fac_g1 = PN2_factor(m1, m2);
	    }
        double mtot = m1+m2;
        double eta = m1*m2/(mtot*mtot);

        double ome21_2 = ome21*ome21;
        double ome21_3 = ome21_2*ome21;
		//double ome21_7 = ome21_3*ome21_3*ome21;
		double ome21_35 = ome21_3*ome21_05;
		double a1_3 = a1 * a1 * a1;
        if (pn25) da1dt += PN25_fac_a1 * (1 + e1_2 * (73.0 / 24 + 37 * e1_2 / 96)) / (a1_3 * ome21_35);
        //double PN25_a1 = PN25_fac_a1 * (1 + e1_2 * (73.0 / 24 + 37 * e1_2 / 96)) / (a1_3 * ome21_35);
        //std::cout << "PN25 da1dt" << PN25_a1 << std::endl;

		//double ome21_5 = ome21_3*ome21*ome21;
		//double ome21_25 = sqrt(ome21_5);
		double ome21_25 = ome21_2*ome21_05;
		//double PN25_e1 = PN25_fac_e1 * e1 * (1 + 121 * e1_2 / 304) / (a1_3 * a1 * ome21_25);
        if (pn25) de1dt += PN25_fac_e1 * e1 * (1 + 121 * e1_2 / 304) / (a1_3 * a1 * ome21_25);;
        //std::cout << "PN25 de1dt" << PN25_e1 << std::endl;

		//double a1_5 = a1_3*a1*a1;
		double a1_2 = a1*a1;
		double a1_25 = a1_2*a1_05;
        if (pn1) dg1dt += PN1_fac_g1 / (ome21 * a1_25);
        //double PN1_g1 = PN1_fac_g1 / (ome21 * a1_25);
        //std::cout << "PN1 dg1dt" << PN1_g1 << std::endl;

        double a1_35 = a1*a1_25;
        if (pn2) dg1dt += - PN2_fac_g1 * (10 + 4*eta - (1+10*eta)*e1_2) / (ome21_2*a1_35);

		//double a1_7 = a1_5*a1*a1;
		//double a1_35 = sqrt(a1_7);
        if (pn25) dHdt += PN25_fac_G1 * (G1+G2*cositot)/H * (1+7*e1_2/8) / (ome21_2*a1_35);
        //double PN25_H = PN25_fac_G1 * (G1+G2*cositot)/H * (1+7*e1_2/8) / (ome21*ome21*a1_35);
        //std::cout << "PN25 dHdt" << PN25_H << std::endl;
    }

    if(tides and not hascollided) {
        double a1_3 = a1 * a1 * a1;
        double a1_12 = a1_3*a1_3*a1_3*a1_3;
        double ome21_15 = ome21*ome21_05;
        double ome21_3 = ome21*ome21*ome21;
        double ome21_65 = ome21_3*ome21_3*ome21_05;
        double ome21_75 = ome21_65*ome21;

        // Fabricky and Tremaine 2007, Eq 29
        dg1dt += 15 * sqrt((m1+m2)/(a1_12*a1)) * (8 + 12*e1_2 + e1_2*e1_2) / (ome21_3*ome21*ome21) *
                (m2/m1 * k1 * R1*R1*R1*R1*R1 + m1/m2 * k2 * R2*R2*R2*R2*R2);

        double meanm = sqrt((m1 + m2) / a1_3);
        double spin = pseudosync_spin(e1_2, meanm);
        double q1 = m2 / m1;
        double R1_over_a = R1/a1;
        double kT1 = m1*tau1*k1 / (R1*R1*R1);
        double q2 = m1 / m2;
        double R2_over_a = R2/a1;
        double kT2 = m2*tau2*k2 / (R2*R2*R2);
        
        double common_tides = kT1*q1*(1+q1) * pow(R1_over_a, 8) + kT2*q2*(1+q2) * pow(R2_over_a, 8);
        double da1dt_tides = common_tides * -6 * a1 / ome21_75 * (f1(e1_2) - ome21_15 * f2(e1_2) * spin/meanm);
        da1dt += da1dt_tides;

        double de1dt_tides = common_tides * -27 * e1 / ome21_65 * (f3(e1_2) - 11./18.* ome21_15 * f4(e1_2) *spin/meanm);
        de1dt += de1dt_tides;
    }

	// e1, e2, g1, g2, h1, a1, H = y_const
	dydt_dev[id::e1] = de1dt;
	dydt_dev[id::e2] = de2dt;
	dydt_dev[id::g1] = dg1dt;
	dydt_dev[id::g2] = dg2dt;
	dydt_dev[id::h1] = dh1dt;
	dydt_dev[id::a1] = da1dt;
	dydt_dev[id::H] = dHdt;
    dydt_dev[id::a2] = da2dt;
    dydt_dev[id::m1] = dm1dt;
    dydt_dev[id::m2] = dm2dt;
    dydt_dev[id::m3] = dm3dt;
    dydt_dev[id::R1] = dR1dt;
    dydt_dev[id::R2] = dR2dt;
    dydt_dev[id::R3] = dR3dt;

    // External derivatives
    if (ext) {
        for (size_t i=0; i<14; i++)
            dydt_dev[i] += yderiv_ext[i];

        delaunay_state yext = yderiv_ext;
        yext[id::a2] = - a2 * (dydt_dev[id::m1] + dydt_dev[id::m2] + dydt_dev[id::m3]) / (m1+m2+m3);
        dydt_dev[id::a2] += yext[id::a2];

        double dHdt_ext = dHdt_from_external_changes(m1, m2, m3, cositot, G1, G2, H, a1, e1, ome21, a2, e2, ome22,
                                                     yext);
        dydt_dev[id::H] += dHdt_ext;
    }
}

void Okinami::stop_at_gw_frequency(double gwstop_frequency) {
    if (gwstop_frequency <= 0.0) {
        check_gw = false;
        fgwstop = 0.0;
    } else {
        check_gw = true;
        fgwstop = gwstop_frequency;
    }
}

void Okinami::check_dynamical_stability(const delaunay_state &y) {
    i_mut = itot_from_y(y);

    switch (check_stability) {
        case stability_criterion::MARDLING :
            isdynstable = is_mardling_stable(y[id::a1], y[id::m1], y[id::m2],
                                             y[id::a2], y[id::m3], y[id::e2], i_mut) and isdynstable;
            break; //optional
        case stability_criterion::PETROVICH :
            //TODO
            break; //optional
        case stability_criterion::AMD :
            isdynstable = is_amd_stable(y[id::a1], y[id::a2], y[id::e1], y[id::e2],
                                        y[id::m1], y[id::m2], y[id::m3],
                                        y[id::H]) and isdynstable;
            break; //optional

    }
}

void Okinami::check_gw_frequency(const delaunay_state &y) {
    double e1 = y[id::e1];
    double a1 = y[id::a1];
    double m1 = y[id::m1];
    double m2 = y[id::m2];

    double p = a1 * (1-e1*e1);
    double f = sqrt(m1+m2) * pow(1 + e1, 1.1954) / (sqrt(p*p*p) * M_PI);
    f = f/(yr_in_s * Tscale); // To Hz

    stopgw = (f > fgwstop) or stopgw;
}


void Okinami::check_quasi_secular(const delaunay_state &y) {
    double e1 = y[id::e1];
    double e2 = y[id::e2];
    double a1 = y[id::a1];
    double a2 = y[id::a2];
    double m1 = y[id::m1];
    double m2 = y[id::m2];
    double m3 = y[id::m3];

    double afac = (a1 / (a2 * (1 - e2)));
    double ome1_crit = 5 * M_PI * m3 / (m1+m2) * afac*afac*afac;
    ome1_crit *= ome1_crit;
    double ome1 = 1 - e1;

    stopsec = (ome1 < ome1_crit) or stopsec;
}

void Okinami::check_particle_collisions(const delaunay_state &y) {
    double e1 = y[id::e1];
    double a1 = y[id::a1];
    double R1 = y[id::R1];
    double R2 = y[id::R2];

    double peri = a1 * (1 - e1);
    double sumrad = R1 + R2;

    hascollided = dcoll*sumrad > peri or hascollided;
}

bool Okinami::is_mardling_stable(double a1, double m1, double m2,
                                 double a2, double m3, double e2, double i_mut) {
    double q = m3 / (m1+m2);
    double ome2 = 1 - e2;
    double aratio_stable = 2.8 / ome2 * (1 - 0.3 * i_mut / M_PI);
    aratio_stable *= pow( (1.0+q)*(1+e2)/sqrt(ome2), 0.4 );
    double aratio = a2 / a1;

    bool isstable = aratio > aratio_stable;
    if (not isstable) {
        std::cout << "MARDLING INSTABILITY" << std::endl;
        std::cout << "a_ratio = " << aratio << std::endl;
        std::cout << "a_stabi = " << aratio_stable << std::endl;
        std::cout << "a1 = " << a1 << "; a2 = " << a2  << "; e2 = " << e2 << "; i_mut = " << i_mut / M_PI * 180 << std::endl;
    }
    return isstable;
}

bool Okinami::is_petrovich_stable(double a1, double m1, double m2, double a2, double m3, double e2, double i_mut) {
    return true;
}

bool Okinami::is_amd_stable(double a1, double a2, double e1, double e2,
                            double m1, double m2, double m3,
                            double H) {
    double H_2 = H*H;
    double G1 = G_moment(m1, m2, a1, e1);
    double G2 = G_moment(m1+m2, m3, a2, e2);
    double G1_2 = G1*G1;
    double G2_2 = G2*G2;

    double cosi1 = (H_2 + G1_2 - G2_2)/(2*H*G1);
    double cosi2 = (H_2 + G2_2 - G1_2)/(2*H*G2);

    double Lambda1 = m2*sqrt(m1*a1);
    double Lambda2 = m3*sqrt(m1*a2);

    double C = Lambda1 * (1-sqrt(1-e1*e1)*cosi1) + Lambda2 * (1-sqrt(1-e2*e2)*cosi2);

    double eps = (m2+m3)/m1;
    double gamma = m2/m3;
    double alpha = a1/a2;
    double gammap1 = (1 + gamma);

    double Chill = gamma * sqrt(alpha) + 1 - pow(1+gamma, 1.5) *
            sqrt(alpha/(alpha+gamma)*(1 + (amd3to43 * gamma * pow(eps, 0.6666666666666666))/(gammap1*gammap1)));

    bool isstable = C < Chill;
    if (not isstable) {
        double i1 = acos(cosi1);
        double i2 = acos(cosi2);
        std::cout << "AMD INSTABILITY" << std::endl;
        std::cout << "C = " << C << std::endl;
        std::cout << "C_thr = " << Chill << std::endl;
        std::cout << "a1 = " << a1 << "; e1 = " << e1  << "; i1 = " << i1 <<std::endl;
        std::cout << "a2 = " << a2 << "; e2 = " << e2  << "; i2 = " << i2 <<std::endl;
    }

    return C < Chill;
}

void Okinami::set_integrator_tolerance(double rel_err, double abs_err, double ax, double adxdt) {
    Integrator.rel_err = rel_err;
    Integrator.abs_err = abs_err;
    Integrator.ax = ax;
    Integrator.adxdt = adxdt;
}

void Okinami::evolve_system(delaunay_state &y, double tfin) {

    if(tfin <= ctime) return;
    if (dt == 0.0) Integrator.estimate_first_timestep(y, ctime, dt);

    do {
        double dtmax = tfin - ctime;
        dt = std::min(dtmax, dt);

        Integrator.iterate_step(y, ctime, dt, dt_done);
        if (Integrator.reject) throw TsuError("Maximum iteration reached");

        ctime += dt_done;

        y[id::g1] = KeplerUtils::mod2pi(y[id::g1]);
        y[id::g2] = KeplerUtils::mod2pi(y[id::g2]);
        y[id::h1] = KeplerUtils::mod2pi(y[id::h1]);

        bool stop = stop_function(y, ctime);
        if (stop) break;

    } while (ctime < tfin);
}

bool Okinami::stop_function(const delaunay_state &x, double t) {
    e1max = std::max(x[0], e1max);
    e1min = std::min(x[0], e1min);

    if (check_stability) check_dynamical_stability(x);
    if (check_collisions) check_particle_collisions(x);
    if (check_gw) check_gw_frequency(x);
    if (check_sec) check_quasi_secular(x);

    return hascollided or not isdynstable or stopgw or stopsec;
}