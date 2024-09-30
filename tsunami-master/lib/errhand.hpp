#ifndef ERRHAND_H
#define ERRHAND_H

#include <iostream>
#include <string>
#include <sstream>
#include <utility>

#ifdef TELEGRAM
#include <bot_alert.h>
#endif

#define s2n(K, a) TsunamiError::_s2n<K>(a, __FILE__, __LINE__)
#define n2s(a) TsunamiError::_n2s(a, __FILE__, __LINE__)
#define TsuError(a) TsunamiError(a, __FILE__, __LINE__)

class TsunamiError : public std::exception {

public:
    TsunamiError() = default; //overloaded constructor

    TsunamiError(std::string err_msg_input, const char *file_input, int line_input) {
        err_msg = std::move(err_msg_input);
        file = to_string_loc(file_input);
        line = line_input;

        print_to_stdout();
        invoke_telegram_bot();
    }

    std::string get_msg() const { return err_msg; }

    std::string get_file() const { return file; }

    int get_line() const { return line; }


    static std::string to_string_loc(const char *obj) {
        std::ostringstream out;
        out << obj;
        return out.str();
    }

    inline void invoke_telegram_bot() {
#ifdef TELEGRAM
        bot = new bot_alert;
        string message = "<b>An error occurred</b>\n<b>Message</b> = " + get_msg() + "\n" + "<b>File</b> = " + get_file() + "\n" + "<b>Line</b> = " + bot->to_string_loc(get_line()) + "\n";

        bot->init();
        bot->send_message(message);
        bot->send_sticker("sad");
#endif
    }

    inline void print_to_stdout() {
        std::cout << " An error occurred: " << std::endl;
        std::cout << "          Message ----> " << get_msg() << std::endl;
        std::cout << "          File    ----> " << get_file() << std::endl;
        std::cout << "          Line    ----> " << get_line() << std::endl;
    }

    template<typename T>
    static const T _s2n(std::string str, const char *file_input, const int line_input) {
        T value;

        std::stringstream stream(str);
        stream >> value;

        if (stream.fail()) throw TsunamiError("Cannot convert string: [" + str + "]", file_input, line_input);

        return value;
    }

    template<typename T>
    static const std::string _n2s(T val, const char *file_input, const int line_input) {
        std::ostringstream stream;
        stream << val;

        if (stream.fail()) throw TsunamiError("Cannot convert number: [" + stream.str() + "]", file_input, line_input);

        return stream.str();
    }

private:
    std::string err_msg, file;
    int line{};
#ifdef TELEGRAM
    bot_alert *bot;
#endif

};

enum ErrType {
    DT_NEG_OR_INF,
    V_GTR_C,
};

class IntegrationError : public std::exception {

public:
    IntegrationError() = default;
    explicit IntegrationError(ErrType err) : error(err) {};

    ErrType error;
    double dt;

};

#endif
