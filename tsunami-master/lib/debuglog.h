//
// Created by lex on 8/19/21.
//

#ifndef TSUNAMI_DEBUGLOG_H
#define TSUNAMI_DEBUGLOG_H

#include <string>
#include <fstream>
#include <sstream>
#include <iomanip>

class DebugLog {
public:
    DebugLog(std::string devlog_name = "devlog.dat") {
        devlogname = devlog_name;

#ifdef DEVLOG
        devlog.open(devlogname.c_str(), std::ios::out);
        devlog<<std::scientific<<std::setprecision(16);
#endif
    }

    ~DebugLog() {
#ifdef DEVLOG
        devlog.close();
#endif
    }

    void debug(std::string debugline) {
#ifdef DEVLOG
        devlog << debugline << std::endl;
#endif
    }

    template <typename T> const std::string n2s(T val, const unsigned int precision=6) {

        std::ostringstream stream;
        stream << std::scientific << std::setprecision(precision);
        stream << val;

        if (stream.fail())
            debug("Cannot convert into a string");

        return stream.str();
    }

private:
    std::string devlogname;
    std::ofstream devlog;
};


#endif //TSUNAMI_DEBUGLOG_H
