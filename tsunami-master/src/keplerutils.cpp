//
// Created by aatrani on 13/06/19.
//

#include "keplerutils.h"
#include "errhand.hpp"


template <typename T> int KeplerUtils::sgn(T val) {
    return (T(0) < val) - (val < T(0));
}

double KeplerUtils::mod2pi(double f){
    f = fmod(M_2PI + fmod(f, M_2PI), M_2PI);
    return f;
}

double KeplerUtils::acos2(double num, double den, double posflag){
    double arccos;
    double cose = num/den;
    if(cose > -1.0 and cose < 1.0) {
        arccos = acos(cose);
        if (posflag < 0.0) {
            arccos = -arccos;
        }
    } else if (den==0) {
        arccos = 0.0;
    } else {
        arccos = (cose <= -1.0) ? M_PI : 0.0;
    }
    return arccos;
}

double KeplerUtils::period(double a, double mu) {
    return M_2PI*sqrt(fabs((a*a*a)/mu));
}

double KeplerUtils::M_to_E(double M, double e){
    // Newton-Raphson method
    double E;
    if(e < 1.){
        M = mod2pi(M);
        E = M + sgn(sin(M)) * 0.85 * e;
        double g = ecc_to_mean(E, M, e);
        for(int i=0; i<maxit; i++){
            E = E - g/ecc_to_mean_prime(E, e);
            g = ecc_to_mean(E, M, e);
            if(fabs(g) < 1e-16){
                break;
            }
        }
        E = mod2pi(E);
        return E;
    }
    else{
        E = M/fabs(M)*log(2.*fabs(M)/e + 1.8);

        double g = ecc_to_mean_h(E, M, e);
        for(int i=0; i<maxit; i++){
            E = E - g/ecc_to_mean_h_prime(E, e);
            g = ecc_to_mean_h(E, M, e);
            if(fabs(g) < 1e-16){
                break;
            }
        }
        return E;
    }
}

double KeplerUtils::E_to_nu(double E, double e) {
    if(e > 1.){
        return 2.*atan2(sqrt(1.+e)*sinh(0.5*E), sqrt(1.-e)*cosh(0.5*E));
    }
    else{
        return 2.* atan2(sqrt(1.+e)*sin(0.5*E), (sqrt(1.-e)*cos(0.5*E)));
    }
}

double KeplerUtils::nu_to_E(double nu, double e) {
    if(e > 1.) {
        return 2.0 * atanh( sqrt((e-1.0)/(e+1.0)) * tan(nu*0.5) );
    } else {
        return mod2pi(atan2(sqrt(1.0-e*e) *sin(nu), (e + cos(nu))));
    }
}

double KeplerUtils::E_to_M(double E, double e) {
    if(e > 1.) {
        return e * sinh(E) - E;
    } else {
        return E - e * sin(E);
    }
}

void KeplerUtils::kepl_to_cart(double out_posvel[6], double prim_posvel[6], double mu, double a, double e, double i, double ome, double Ome, double nu){
    if(e == 1.){
        throw TsuError("e == 1");
    }
    if(e < 0.){
        throw TsuError("e < 0");
    }
    if(e > 1.){
        if(a > 0.){
            throw TsuError("e > 1 but a > 0");
        }
    }
    else{
        if(a < 0.){
            throw TsuError("e < 1 but a < 0");
        }
    }
    if(e*cos(nu) < -1.){
        throw TsuError("e*cos(nu) < -1");
    }

    double cf = cos(nu);

    double r = a*(1-e*e)/(1 + e*cf);
    double v0 = sqrt(mu/(a*(1.-e*e)));

    double cO = cos(Ome);
    double sO = sin(Ome);
    double co = cos(ome);
    double so = sin(ome);
    double sf = sin(nu);
    double ci = cos(i);
    double si = sin(i);

    out_posvel[0] = prim_posvel[0] + r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
    out_posvel[1] = prim_posvel[1] + r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
    out_posvel[2] = prim_posvel[2] + r*(so*cf+co*sf)*si;

    out_posvel[3] = prim_posvel[3] + v0*((e+cf)*(-ci*co*sO - cO*so) - sf*(co*cO - ci*so*sO));
    out_posvel[4] = prim_posvel[4] + v0*((e+cf)*(ci*co*cO - sO*so)  - sf*(co*sO + ci*so*cO));
    out_posvel[5] = prim_posvel[5] + v0*((e+cf)*co*si - sf*si*so);
}

void KeplerUtils::kepl_to_cart(double out_posvel[6], double prim_posvel[6], double m1, double m2, double a, double e,
                               double i, double ome, double Ome, double nu, bool pn1, bool pn2) {
    if(e == 1.){
        throw TsuError("e == 1");
    }
    if(e < 0.){
        throw TsuError("e < 0");
    }
    if(e > 1.){
        if (pn1 or pn2) {
            // Disabling PN corrections for hyperbolic orbit
            pn1 = pn2 = false;
        }
        if(a > 0.){
            throw TsuError("e > 1 but a > 0");
        }
    }
    else{
        if(a < 0.){
            throw TsuError("e < 1 but a < 0");
        }
    }
    if(e*cos(nu) < -1.){
        throw TsuError("e*cos(nu) < -1");
    }

    double cf = cos(nu);
    double sf = sin(nu);

    double e2 = e*e;
    double mu = m1 + m2;
    double p = a * (1 - e2);
    double r0 = p/(1 + e*cf);
    double v0 = sqrt(mu / (a * (1. - e * e)));

    double r = r0;

    double vx, vy;

    if (pn1 or pn2) {
        double eta = m1*m2/(mu*mu);
        double e3 = e*e2;
        double e4 = e2*e2;
        double zeta = mu / p;
        double zeta05 = sqrt(zeta);
        double zeta2 = zeta*zeta;
        double zeta3 = zeta2*zeta;
        double c2_inv = 1/(speed_of_light*speed_of_light);
        double c4_inv = c2_inv*c2_inv;
        double eta2 = eta*eta;

        double c2f = cos(2*nu);
        double c3f = cos(3*nu);
        double c4f = cos(4*nu);
        double s2f = sin(2*nu);
        double s3f = sin(3*nu);
        double s4f = sin(4*nu);

        double A = 1 + 1.0/3*(-eta + 1.75*e2*(4-3*eta));
        double B = 1.0/3*e*((9-4*eta) + (1-3*eta)*e2);
        double C =  - 0.25*eta*e2;
        double D = 1 - 65.0/12*eta + e2/24*(356-319*eta+48*eta2) + e4/192*(256-265*eta+459*eta*eta);
        double E = e*(1.0/12*(96-231*eta+8*eta2) + e2/48*(323 - 351*eta + 180*eta2) + eta/12*(5+12*eta)*e4);
        double F = - e2*(1.0/24*(60+159*eta-16*eta2) + eta/48*(29-27*eta)*e2);
        double G = - e3*0.0625*(1 + 19*eta - 4*eta2);
        double H = - e4*0.015625*eta*(1 -3*eta);

        double m_r_corr = 0.0;
        double r2nudot_corr = 1;
        double rdot_corr = e*sf;

        if (pn1) {
            m_r_corr += c2_inv*zeta2*(A + B*cf + C*c2f);
            r2nudot_corr += - zeta*c2_inv*(2.0/3*((3-eta) + (1-3*eta)*e2) + (4-2*eta)*e*cf);
            rdot_corr += c2_inv*zeta*(B*sf + 2*C*s2f);
        }
        if (pn2) {
            m_r_corr += c4_inv*zeta3*(D + E*cf + F*c2f + G*c3f + H*c4f);
            r2nudot_corr +=  zeta2*c4_inv*( 1.0/3*(0.5*(6+53*eta+2*eta2) - 0.125*(28 + 117*eta - 12*eta2)*e2 + 0.5*(2-17*eta+6*eta2)*e4
                                                   - ( 6-53*eta-2*eta2 - 0.125*(32 - 211*eta + 54*eta2)*e2)*e*cf)
                                            + 0.125*(36-13*eta+4*eta2)*e2*c2f - eta*0.125*(3+2*eta)*e3*c3f);
            rdot_corr += c4_inv*zeta2*(E*sf + 2*F*s2f + 3*G*s3f + 4*H*s4f);

        }
        // Position correction
        r /= (1 + m_r_corr*r0/mu);

        double r2nudot = mu / zeta05 * r2nudot_corr;
        double rnudot = r2nudot / r;
        double rdot = r2nudot/p * rdot_corr;

        // Velocity correction
        vx = rdot * cf - rnudot * sf;
        vy = rdot * sf + rnudot * cf;

    } else {
        vx = - v0 * sf;
        vy = v0 * (e + cf);
    }

    //std::cout << "m_r corr fact " << m_r_corr << std::endl;
    //std::cout << "m_r  " << mu/r0 << std::endl;
    //std::cout << "r/r0 " << 1/(1 + m_r_corr*r0/mu) << std::endl;
    //std::cout << "Vcorr " << Vcorr << std::endl;

    double cO = cos(Ome);
    double sO = sin(Ome);
    double co = cos(ome);
    double so = sin(ome);
    double ci = cos(i);
    double si = sin(i);

    out_posvel[0] = prim_posvel[0] + r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
    out_posvel[1] = prim_posvel[1] + r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
    out_posvel[2] = prim_posvel[2] + r*(so*cf+co*sf)*si;

    out_posvel[3] = prim_posvel[3] + (vy * (-ci * co * sO - cO * so) + vx * (co * cO - ci * so * sO));
    out_posvel[4] = prim_posvel[4] + (vy * (ci * co * cO - sO * so) + vx * (co * sO + ci * so * cO));
    out_posvel[5] = prim_posvel[5] + (vy * co * si + vx * si * so);
}


void KeplerUtils::kepl_to_cart2(double out_posvel[6], double prim_posvel[6], double m1, double m2, double a, double e,
                               double i, double ome, double Ome, double nu, bool pn1, bool pn2, bool pn3) {
    if(e == 1.){
        throw TsuError("e == 1");
    }
    if(e < 0.){
        throw TsuError("e < 0");
    }
    if(e > 1.){
        if (pn1 or pn2) {
            // Disabling PN corrections for hyperbolic orbit
            pn1 = pn2 = false;
        }
        if(a > 0.){
            throw TsuError("e > 1 but a > 0");
        }
    }
    else{
        if(a < 0.){
            throw TsuError("e < 1 but a < 0");
        }
    }
    if(e*cos(nu) < -1.){
        throw TsuError("e*cos(nu) < -1");
    }

    double cf = cos(nu);
    double sf = sin(nu);

    double e2 = e*e;
    double mu = m1 + m2;
    double p = a * (1 - e2);
    double r0 = p/(1 + e*cf);
    double v0 = sqrt(mu / (a * (1. - e * e)));

    double r = r0;

    double vx, vy;

    if (pn1 or pn2) {
        double eta = m1*m2/(mu*mu);
        double e3 = e*e2;
        double e4 = e2*e2;
        double zeta = mu / p;
        double zeta05 = sqrt(zeta);
        double zeta2 = zeta*zeta;
        double zeta3 = zeta2*zeta;
        double c2_inv = 1/(speed_of_light*speed_of_light);
        double c4_inv = c2_inv*c2_inv;
        double eta2 = eta*eta;

        double c2f = cos(2*nu);
        double c3f = cos(3*nu);
        double c4f = cos(4*nu);
        double s2f = sin(2*nu);
        double s3f = sin(3*nu);
        double s4f = sin(4*nu);

        double A = 1 + 1.0/3*(-eta + 1.75*e2*(4-3*eta));
        double B = 1.0/3*e*((9-4*eta) + (1-3*eta)*e2);
        double C =  - 0.25*eta*e2;
        double D = 1 - 65.0/12*eta + e2/24*(356-319*eta+48*eta2) + e4/192*(256-265*eta+459*eta*eta);
        double E = e*(1.0/12*(96-231*eta+8*eta2) + e2/48*(323 - 351*eta + 180*eta2) + eta/12*(5+12*eta)*e4);
        double F = - e2*(1.0/24*(60+159*eta-16*eta2) + eta/48*(29-27*eta)*e2);
        double G = - e3*0.0625*(1 + 19*eta - 4*eta2);
        double H = - e4*0.015625*eta*(1 -3*eta);

        double m_r_corr = 0.0;
        double r2nudot_corr = 1;
        double rdot_corr = e*sf;

        if (pn1) {
            m_r_corr += c2_inv*zeta2*(A + B*cf + C*c2f);
            r2nudot_corr += - zeta*c2_inv*(2.0/3*((3-eta) + (1-3*eta)*e2) + (4-2*eta)*e*cf);
            rdot_corr += c2_inv*zeta*(B*sf + 2*C*s2f);
        }
        if (pn2) {
            m_r_corr += c4_inv*zeta3*(D + E*cf + F*c2f + G*c3f + H*c4f);
            r2nudot_corr +=  zeta2*c4_inv*( 1.0/3*(0.5*(6+53*eta+2*eta2) - 0.125*(28 + 117*eta - 12*eta2)*e2 + 0.5*(2-17*eta+6*eta2)*e4
                                                   - ( 6-53*eta-2*eta2 - 0.125*(32 - 211*eta + 54*eta2)*e2)*e*cf)
                                            + 0.125*(36-13*eta+4*eta2)*e2*c2f - eta*0.125*(3+2*eta)*e3*c3f);
            rdot_corr += c4_inv*zeta2*(E*sf + 2*F*s2f + 3*G*s3f + 4*H*s4f);

        }
        // Position correction
        r /= (1 + m_r_corr*r0/mu);

        double r2nudot = mu / zeta05 * r2nudot_corr;
        double rnudot = r2nudot / r;
        double rdot = r2nudot/p * rdot_corr;

        // Velocity correction
        vx = rdot * cf - rnudot * sf;
        vy = rdot * sf + rnudot * cf;

    } else {
        vx = - v0 * sf;
        vy = v0 * (e + cf);
    }

    //std::cout << "m_r corr fact " << m_r_corr << std::endl;
    //std::cout << "m_r  " << mu/r0 << std::endl;
    //std::cout << "r/r0 " << 1/(1 + m_r_corr*r0/mu) << std::endl;
    //std::cout << "Vcorr " << Vcorr << std::endl;

    double cO = cos(Ome);
    double sO = sin(Ome);
    double co = cos(ome);
    double so = sin(ome);
    double ci = cos(i);
    double si = sin(i);

    out_posvel[0] = prim_posvel[0] + r*(cO*(co*cf-so*sf) - sO*(so*cf+co*sf)*ci);
    out_posvel[1] = prim_posvel[1] + r*(sO*(co*cf-so*sf) + cO*(so*cf+co*sf)*ci);
    out_posvel[2] = prim_posvel[2] + r*(so*cf+co*sf)*si;

    out_posvel[3] = prim_posvel[3] + (vy * (-ci * co * sO - cO * so) + vx * (co * cO - ci * so * sO));
    out_posvel[4] = prim_posvel[4] + (vy * (ci * co * cO - sO * so) + vx * (co * sO + ci * so * cO));
    out_posvel[5] = prim_posvel[5] + (vy * co * si + vx * si * so);
}


void KeplerUtils::cart_to_kepl(double out_keplpar[6], double in_pos[3], double in_vel[3], double m1, double m2,
                               bool pn1, bool pn2, bool pn3, bool memmcorr) {

    double dx = in_pos[0];
    double dy = in_pos[1];
    double dz = in_pos[2];

    double r = sqrt ( dx*dx + dy*dy + dz*dz );
    if(r < dtiny) {
        throw TsuError("particle have same position");
    }

    double mu = m1 + m2;
    if(mu < dtiny) {
        throw TsuError("total mass is zero");
    }

    double dvx = in_vel[0];
    double dvy = in_vel[1];
    double dvz = in_vel[2];
    double v2 = dvx*dvx + dvy*dvy + dvz*dvz;

    double rvr = (dx*dvx + dy*dvy + dz*dvz);
    double vr = rvr/r;

    double muinv = 1/mu;
    double mu_r = mu / r;

    // Specific energy
    double E = 0.5 * v2 - mu_r;

    // Specific angular momentum
    double Lx = (dy * dvz - dz * dvy);
    double Ly = (dz * dvx - dx * dvz);
    double Lz = (dx * dvy - dy * dvx);
    double L0 = sqrt (Lx * Lx + Ly * Ly + Lz * Lz );
    double L = L0;

    double a, e, e2;
    bool haspn = pn1 or pn2 or pn3;
    if (haspn) {
        double eta = m1*m2/(mu*mu);
        double inv_c2 = 1/(speed_of_light*speed_of_light);
        double E1 = 0.0;
        double E2 = 0.0;
        double E3 = 0.0;
        double L1c = 0.0;
        double L2c = 0.0;
        double L3c = 0.0;

        double mu_r2 = mu_r*mu_r;
        double v4 = v2*v2;
        double vr2 = vr*vr;
        double eta2 = eta*eta;

        double mu_r3 = mu_r2*mu_r;
        double v6 = v4*v2;
        double vr4 = vr2*vr2;
        double pi2 = M_PI*M_PI;
        double eta3 = eta2*eta;

        if (pn1) {
            E1 = 0.5 * mu_r2 + 0.375 * (1 - 3*eta)*v4 + 0.5 * (3 + eta)*v2*mu_r + 0.5*eta*mu_r*vr2;
            L1c = (3 + eta)*mu_r + 0.5*(1-3*eta)*v2;
        }

        if (pn2) {
            E2 = -0.25*(2 + 15*eta)*mu_r3
                    + 0.3125 * (1 - 7*eta + 13*eta2) * v6
                    + 0.125 * (14 - 55*eta + 4*eta2) * mu_r2 * v2
                 + 0.125 * (4 + 69*eta + 12*eta2) * mu_r2 * vr2
                 + 0.125 * (21 - 23*eta - 27*eta2) * mu_r * v4
                 + 0.25*eta*(1 - 15*eta) * mu_r * v2 * vr2
                 - 0.375*eta*(1 - 3*eta) * mu_r * vr4;
            L2c = 0.25*(14 - 41*eta +4*eta2) * mu_r2
                    + 0.375*(1 - 7*eta +13*eta2) * v4
                    + 0.5*(7 - 10*eta - 9*eta2) * mu_r * v2
                    - 0.5*eta*(2 + 5*eta) * mu_r * vr2;
        }

        if (pn3) {
            double v8 = v4*v4;
            double vr6 = vr4*vr2;
            double mu_r4 = mu_r2*mu_r2;

            E3 = (0.375 + 18469.0/840 * eta) * mu_r4
                    + (1.25 - (6747.0/280.0 - 0.640625 * pi2)*eta - 5.25*eta2 + 0.5*eta3) * mu_r3 * v2
                    + (1.5 + (2321.0/280.0 - 1.921875 * pi2)*eta + 12.75*eta2 + 3.5*eta3) * mu_r3 * vr2
                    + (0.2734375 - 3.2265625*eta + 13.015625*eta2 - 17.6640625*eta3) * v8
                    + (8.4375 - 12.125*eta + 25.375*eta2 - 6.75*eta3) * mu_r2*v4
                    + (0.75 + 15.5*eta - 50.9375*eta2 - 20.25*eta3) * mu_r2*v2*vr2
                    - eta * (731.0/48.0 - 10.25*eta - 6.0*eta2) * mu_r2*vr4
                    + (3.4375 - 13.4375*eta + 10.375*eta2 + 20.3125*eta3) * mu_r*v6
                    + eta*(0.3125 - 1.5625*(eta - eta2)) * mu_r*vr6
                    - eta*(1.3125 + 4.6875*eta - 23.4375*eta2) * mu_r*v4*vr2
                    - eta*(0.5625 - 5.25*eta + 10.3125*eta2) * mu_r*v2*vr4;
            L3c = (2.5 - (5199.0/280.0 - 1.28125 * pi2)*eta - 7*eta2 + eta3) * mu_r3
                    + (0.3125 - 3.6875*eta + 14.875*eta2 - 20.1875*eta3) * v6
                    + (11.25 - 322.0/12.0*eta + 26.25*eta2 - 9.0*eta3) * mu_r2*v2
                    + (0.5 - 287.0/24.0*eta - 39.625*eta2 - 13.5*eta3) * mu_r2*vr2
                    + (4.125 - 17.75*eta + 13.25*eta2 + 24.375*eta3) * mu_r*v4
                    - eta*(3.0 - 1.75*eta - 18.75*eta2) * mu_r*v2*vr2
                    + eta*(0.75*(1.0-eta) - 4.125*eta2) * mu_r*vr4;
        }

        E += inv_c2 * (E1 + inv_c2*(E2 + inv_c2*E3));
        L *= 1 + inv_c2 * (L1c + inv_c2*(L2c + inv_c2*L3c));

        if (memmcorr) {
            // Eq 25a,b from Memmesheimer
            double h = L/mu;
            double Efact = -0.5*E*inv_c2;
            double Lfact = -2*E*h*h;
            double inv_Lfact = 1.0/Lfact;

            a = -mu/(2*E);
            e2 = 1 - Lfact;

            double acorr = 1;
            double ecorr = 0;
            if (pn1) {
                acorr += Efact * (-7 + eta);
                ecorr += Efact * (24 - 4*eta + 5*(-3 + eta)*Lfact);
            }

            if (pn2) {
                acorr += Efact*Efact * (1 + eta2 + 16.0*inv_Lfact*(-4 + 7*eta));
                ecorr += 2*Efact*Efact * (60.0 + 148*eta + 2*eta2 - (80 - 45*eta + 4*eta2)*Lfact
                        - 32*inv_Lfact*(-4 + 7*eta));
            }

            if (pn3) {
                acorr += Efact*Efact*Efact * (1.0 - eta + eta3
                        + inv_Lfact * (256.0 + 41.0*pi2*eta -215408.0*eta/105.0 + 448.0*eta2)
                        - 4*inv_Lfact*inv_Lfact*(512.0 - 176024*eta/105.0 + 41.0*pi2*eta + 144.0*eta2));
                ecorr += Efact*Efact*Efact * (-32.0 + 181264*eta/105.0 + 82.0*pi2*eta - 640.0*eta2
                        + Lfact * (-1488.0 + 1120*eta + 195*eta2 + 4*eta3)
                        - 80.0*inv_Lfact * (9.6 - 4226.0*eta/21.0 + 8.2*pi2*eta + 21.6*eta2)
                        + 16*inv_Lfact*inv_Lfact * (512.0 - 176024*eta/105 + 41.0*pi2*eta + 144.0*eta2));
            }

            a = a * acorr;
            e2 = e2 + ecorr;
            e = sqrt(e2);
        }
    }

    if (not memmcorr or not haspn) {
        double h = L/mu;
        a = -mu/(2*E);
        e2 = 1 + 2*E*h*h;
        if (e2 > min_e2) {
            e = sqrt(e2);
        } else {
            e = 0.0;
        }
    }

    double vdiff2 = v2 - mu_r;

    double ex = muinv*( vdiff2*dx - rvr*dvx);
    double ey = muinv*( vdiff2*dy - rvr*dvy);
    double ez = muinv*( vdiff2*dz - rvr*dvz);
    //double e = sqrt( ex*ex + ey*ey + ez*ez );

    // inclination
    double inc = acos2(Lz, L0, 1.);

    // line of the nodes
    double nx = -Ly;							// vector pointing along the ascending node = zhat cross h
    double ny =  Lx;
    double n = sqrt( nx*nx + ny*ny );

    // Longitude of ascending node
    double Ome = acos2(nx, n, ny);

    // planar case
    double ome, nu, true_long, peri_long, ome_nu;
    if(inc < min_i or inc > M_PI - min_i){
        true_long = acos2(dx, r, dy);		//true longitude (planar)
        peri_long = acos2(ex, e, ey);		//longitude of pericenter (planar)

        // prograde
        if(inc < M_PI/2.){
            ome = peri_long - Ome;   		// argument of pericenter (planar, prograde)
            nu = true_long - peri_long;		// true anomaly (planar, prograde)
        } else{
            ome = Ome - peri_long;			// argument of pericenter (planar, retrograde)
            nu = peri_long - true_long;		// true anomaly (planar, retrograde)
        }

    }
        // non-planar case
    else{
        ome_nu = acos2(nx*dx + ny*dy, n*r, dz);
        ome = acos2(nx*ex + ny*ey, n*e, ez);
        if(inc < M_PI/2.){
            peri_long = Ome + ome;			// longitude of pericenter (inclined, prograde)
            nu = ome_nu - ome;				// true anomaly (inclined, prograde)
            true_long = Ome + ome_nu;		// true longitude (inclided, prograde)
        }
        else{
            peri_long = Ome - ome;			// longitude of pericenter (inclined, retrograde)
            nu = ome_nu - ome;				// true anomaly (inclined, retrograde)
            true_long = Ome - ome_nu;		// true longitude (inclided, retrograde)
        }
    }
    out_keplpar[0] = a;
    out_keplpar[1] = e;
    out_keplpar[2] = inc;
    out_keplpar[3] = mod2pi(ome);
    out_keplpar[4] = mod2pi(Ome);
    out_keplpar[5] = mod2pi(nu);
}


void KeplerUtils::cart_to_kepl(double out_keplpar[6], double in_posvel[6], double mu) {
    double dx,dy,dz,r,dvx,dvy,dvz,v2,vcirc2,vdiff2;
    double lx,ly,lz,l,vr,rvr,muinv,ex,ey,ez,nx,ny,n;
    double true_long,peri_long,ome_nu;

    if(mu < dtiny) {
        throw TsuError("total mass is zero");
    }

    dx = in_posvel[0];
    dy = in_posvel[1];
    dz = in_posvel[2];
    dvx = in_posvel[3];
    dvy = in_posvel[4];
    dvz = in_posvel[5];
    r = sqrt ( dx*dx + dy*dy + dz*dz );
    if(r < dtiny) {
        throw TsuError("particle have same position");
    }

    v2 = dvx*dvx + dvy*dvy + dvz*dvz;

    vcirc2 = mu/r;

    // semimajor axis
    double a = -mu/( v2 - 2.*vcirc2 );

    //angular momentum vector
    lx = (dy*dvz - dz*dvy);
    ly = (dz*dvx - dx*dvz);
    lz = (dx*dvy - dy*dvx);
    l = sqrt ( lx*lx + ly*ly + lz*lz );

    vdiff2 = v2 - vcirc2;
    vr = (dx*dvx + dy*dvy + dz*dvz)/r;
    rvr = r*vr;
    muinv = 1./mu;

    // eccentricity
    ex = muinv*( vdiff2*dx - rvr*dvx );
    ey = muinv*( vdiff2*dy - rvr*dvy );
    ez = muinv*( vdiff2*dz - rvr*dvz );
    double e = sqrt( ex*ex + ey*ey + ez*ez );

    // inclination
    double inc = acos2(lz, l, 1.);

    // line of the nodes
    nx = -ly;							// vector pointing along the ascending node = zhat cross h
    ny =  lx;
    n = sqrt(nx*nx + ny*ny);

    // Longitude of ascending node
    double Ome = acos2(nx, n, ny);

    // planar case
    double ome, nu;
    if(inc < min_i or inc > M_PI - min_i) {
        true_long = acos2(dx, r, dy);		//true longitude (planar)
        peri_long = acos2(ex, e, ey);		//longitude of pericenter (planar)

        // prograde
        if(inc < 0.5*M_PI) {
            ome = peri_long - Ome;   		// argument of pericenter (planar, prograde)
            nu = true_long - peri_long;		// true anomaly (planar, prograde)
        }
        else {
            ome = Ome - peri_long;			// argument of pericenter (planar, retrograde)
            nu = peri_long - true_long;		// true anomaly (planar, retrograde)
        }
    }
        // non-planar case
    else{
        ome_nu = acos2(nx*dx + ny*dy, n*r, dz);
        ome = acos2(nx*ex + ny*ey, n*e, ez);
        if(inc < 0.5*M_PI){
            peri_long = Ome + ome;			// longitude of pericenter (inclined, prograde)
            nu = ome_nu - ome;				// true anomaly (inclined, prograde)
            true_long = Ome + ome_nu;		// true longitude (inclided, prograde)
        }
        else{
            peri_long = Ome - ome;			// longitude of pericenter (inclined, retrograde)
            nu = ome_nu - ome;				// true anomaly (inclined, retrograde)
            true_long = Ome - ome_nu;		// true longitude (inclided, retrograde)
        }
    }
    out_keplpar[0] = a;
    out_keplpar[1] = e;
    out_keplpar[2] = inc;
    out_keplpar[3] = ome;
    out_keplpar[4] = Ome;
    out_keplpar[5] = nu;
}