#ifndef ARTYPES_H
#define ARTYPES_H

//Define some vectors.  They will be used by CPU/ARchain only.
//No GPUs here.

//this is a template. In this way I can change precision
//of the declared variables with just one change in the last line
//of this file: typedef mychainvectors<double> chvec;

//All quantities
#include <iostream>
#include <vector>
#include <cmath>
#include "config.hpp"

template<typename P>
struct chvec {

    //x,  y,  z,  for positions/velocities
    P x, y, z;

    //in the following, 'this' is the object that call the overloaded operators
    constexpr chvec(const chvec &) = default;

    //Default constructor initializes all elements to 0
    //0 here is intended as a type P number
    inline chvec() { x = y = z = P(0); }

    //Overloaded constructor. Explicitly pass elements to use for initialization.
    inline chvec(const P &xx, const P &yy, const P &zz) {
        x = xx;
        y = yy;
        z = zz;
    }

    //Overloaded constructor. Pass a common constant.
    inline chvec(const P &r) {
        x = r;
        y = r;
        z = r;
    }

    //overload common operations (chvec1 + chvec2)
    inline const chvec<P> operator+(const chvec<P> &right) const { //const here means that 'this' must not be modified
        return chvec<P>(this->x + right.x, this->y + right.y, this->z + right.z);
    }

    //overload common operations (chvec1 - chvec2)
    inline const chvec<P> operator-(const chvec<P> &right) const { //const here means that 'this' must not be modified
        return chvec<P>(this->x - right.x, this->y - right.y, this->z - right.z);
    }

    //overload common operations (chvec1 * number)
    //number can be of any type (generic type K)
    template<typename K>
    inline const chvec<P> operator*(const K &r) const { //const here means that 'this' must not be modified
        return chvec<P>(P(r) * this->x, P(r) * this->y, P(r) * this->z);
    }

    //overload common operations (chvec1 * chvec2): DOT PRODUCT
    inline const P operator*(const chvec<P> &right) const { //const here means that 'this' must not be modified
        return (this->x * right.x + this->y * right.y + this->z * right.z);
    }

    //overload common operations (chvec1/chvec2)
    inline const chvec<P> operator/(const chvec<P> &right) const { //const here means that 'this' must not be modified
        return chvec<P>(this->x / right.x, this->y / right.y, this->z / right.z);
    }

    //overload common operations (chvec1 = chvec2)
    inline const chvec<P> &
    operator=(const chvec<P> &right) { //we do not have 'const' here because 'this' MUST be modified (assignment =)
        this->x = right.x;
        this->y = right.y;
        this->z = right.z;
        return *this;
    }

    //overload common operations (chvec1 = real)
    template<typename K>
    inline const chvec<P> &
    operator=(const K &r) { //we do not have 'const' here because 'this' MUST be modified (assignment =)
        this->x = P(r);
        this->y = P(r);
        this->z = P(r);
        return *this;
    }

    //overload common operations (chvec1 += chvec2)
    inline const chvec<P> &
    operator+=(const chvec<P> &right) { //we do not have 'const' here because 'this' MUST be modified (assignment =)
        *this = *this + right; //call overloaded operators + and =
        return *this;
    }

    //overload common operations (chvec1 -= chvec2)
    inline const chvec<P> &
    operator-=(const chvec<P> &right) { //we do not have 'const' here because 'this' MUST be modified (assignment =)
        *this = *this - right; //call overloaded operators + and =
        return *this;
    }

    //overload common operations (chvec1 *= chvec2)
    inline const chvec<P> &
    operator*=(const chvec<P> &right) { //we do not have 'const' here because 'this' MUST be modified (assignment =)
        *this = *this * right; //call overloaded operators + and =
        return *this;
    }

    //overload common operations (chvec1 / number)
    //number can be of any type (generic type K)
    template<typename K>
    inline const chvec<P> operator/(const K &r) const { //const here means that 'this' must not be modified
        P val = P(1) / P(r);
        return chvec<P>(this->x * val, this->y * val, this->z * val);
    }

    //overload common operations (chvec1 /= number)
    template<typename K>
    inline const chvec<P> &
    operator/=(const K &r) { //we do not have 'const' here because 'this' MUST be modified (assignment =)
        *this = *this / P(r); //call overloaded operators + and =
        return *this;
    }

    //obtain absolute values
    inline chvec fabs() const {
        return chvec<P>(std::fabs(this->x), std::fabs(this->y), std::fabs(this->z));
    }

    //obtain absolute value 
    inline const P mod() const {
        return sqrt((*this) * (*this));
    }

    //cross product
   inline chvec<P> cross(const chvec &right) const { //const here means that 'this' must not be modified
        return chvec<P>(this->y * right.z - this->z * right.y, this->z * right.x - this->x * right.z,
                        this->x * right.y - this->y * right.x);
    }

   inline double max() const {
        double maxval = this->x;
        maxval = (maxval > this->y) ? maxval : this->y;
        maxval = (maxval > this->z) ? maxval : this->z;
        return maxval;
    }

    inline void to_array(double a[3]) const {
        a[0] = this->x;
        a[1] = this->y;
        a[2] = this->z;
    }
};

//overload common operations (chvec1 * number)
//number can be of any type (generic type K)
//this must be a non-member function. In this way we do not have the pointer 'this' and the operator * takes just the two arguments r and right
//We could have also declared this function as a FRIEND member function, since friend functions do not have the 'this' pointer.
template<typename P, typename K>
inline const chvec<P> operator*(const K &r, const chvec<P> &right) { //const here means that 'this' must not be modified
    return chvec<P>(right * P(r)); //call overloaded operator * (member function)
}

template<typename P>
inline const chvec<P> operator-(const chvec<P> &right) { //const here means that 'this' must not be modified
    return chvec<P>(-right.x, -right.y, -right.z); //call overloaded operator * (member function)
}

template<typename P>
int sgn(P val) {
    return (P(0) <= val) - (val < P(0));
}

template<typename P>
inline std::ostream &operator<<(std::ostream &output, const chvec<P> &right) {
    output << right.x << "   " << right.y << "   " << right.z;
    return output;
}

template<typename P>
inline const chvec<P> min(const chvec<P> &left, const chvec<P> &right) {

    chvec<P> tmp;

    tmp.x = (left.x < right.x) ? left.x : right.x;
    tmp.y = (left.y < right.y) ? left.y : right.y;
    tmp.z = (left.z < right.z) ? left.z : right.z;

    return tmp;
}

template<typename P>
inline const chvec<P> max(const chvec<P> &left, const chvec<P> &right) {

    chvec<P> tmp;

    tmp.x = (left.x < right.x) ? right.x : left.x;
    tmp.y = (left.y < right.y) ? right.y : left.y;
    tmp.z = (left.z < right.z) ? right.z : left.z;

    return tmp;
}

typedef chvec<double> double3;

template <bool withspin>
struct BStable {
    std::vector<double3> pos;
    std::vector<double3> vel;
    std::vector<double3> spin;
    std::vector<double> eloss;

    double B;
    double omega;
    double time;

    void allocate_arrays(size_t Neqs, size_t Npart) {
        pos.resize(Neqs);
        vel.resize(Neqs);
        eloss.resize(Npart);

        if constexpr (withspin) spin.resize(Npart);
    }
};

typedef BStable<TsunamiConfig::wSpins> BSTable;


enum ptype { // stars
    LOW_MASS_MS = 0,
    HIGH_MASS_MS = 1,
    HGAP = 2,
    GIANT_BRANCH = 3,
    CORE_HE_BURN = 4,
    EARLY_AGBe = 5,
    TP_AGB = 6,
    NAKED_HE = 7,
    NAKED_HE_HGAP = 8,
    NAKED_HE_GIANT = 9,
    HE_WD = 10,
    CO_WD = 11,
    ONE_WD = 12,
    NS = 13,
    BH = 14,
    // planets
    ROCKY = 100,
    GAS_GIANT = 101,
    // null
    UNCLASSIFIED = -1,
};

// particle infos
struct pinfo {
    ptype stype; // stellar type
    double polyt; // polytropic index
    double kaps; // 2nd degree love number (apsidal constant)
    double inert; // inertia
    double taulag; // time lag
    double eloss; //
    bool haspn;
    bool hastide;

    double sigmadiss = 0.0;
    double Atide = 0.0;

    //Default constructor initializes all elements to null
    inline pinfo() {
        stype = ptype::UNCLASSIFIED;
        kaps = 0.0;
        polyt = 0.0;
        inert = 0.0;
        taulag = 0.0;
        eloss = 0.0;
        haspn = false;
        hastide = false;
    }

    //overload common operations (pinfo1 = pinfo2)
    /*inline const pinfo & operator=(const pinfo &right) { //we do not have 'const' here because 'this' MUST be modified (assignment =)
        this->stype = right.stype;
        this->polyt = right.polyt;
        this->kaps = right.kaps;
        this->inert = right.inert;
        this->taulag = right.taulag;
        this->eloss = right.eloss;
        this->haspn = right.haspn;
        this->hastide = right.hastide;
        return *this;
    }*/

    inline bool isplanet() const {
        return (this->stype >= 100) and (this->stype <= 1000);
    }

};

inline std::ostream& operator<< (std::ostream& output, const pinfo &right) {
    output << right.stype << "   " << right.polyt <<"   " << right.kaps <<"   " << right.inert <<"   "
           << right.taulag  << "   " << right.haspn <<"   " << right.hastide;
    return output;
}

// swap function
template <typename T>
void make_ascending(T& t1, T& t2) {
    if (t1 > t2)
        std::swap(t1, t2);
}


// collision info - for now only used for the fewbody interface
struct collinfo {
    short int collind[2]; // indexes of particles
    short int colltype; // collsion type - as BSEXIT

    //Default constructor initializes all elements to null
    inline collinfo() {
        collind[0] = collind[1] = -1;
        colltype = -1;
    }
};


#endif
