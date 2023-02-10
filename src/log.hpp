#ifndef LOG_H_
#define LOG_H_

#include <fstream>
#include <iomanip> // setw
#include <iostream>

class Logger
{
  public:
    std::ofstream cao;
    bool is_screen = true;

    Logger(std::string filename)
    {
        cao.open(filename.c_str());
        if(!cao) throw std::runtime_error(filename + " : " + strerror(errno));
        cao.precision(3);
        cao.flags(std::ios::fixed | std::ios::right);
    }

    ~Logger()
    {
        cao.close();
    }

    template<class S>
    Logger & operator<<(const S & val)
    {
        cao << val;
        if(is_screen) std::cout << val;
        return *this;
    }

    Logger & operator<<(std::ostream & (*pfun)(std::ostream &))
    {
        pfun(cao);
        if(is_screen) pfun(std::cout);
        return *this;
    };

    template<class S>
    void printSpace(std::ostream & os, const S & val)
    {
        if(std::is_integral_v<std::decay_t<decltype(val)>>)
            os << std::setw(2) << val;
        else if(std::is_floating_point_v<std::decay_t<decltype(val)>>)
            os << std::fixed << val;
        else
            os << val << " ";
    }

    template<typename... Args>
    void print(const Args &... args)
    {
        std::cout.precision(3);
        std::cout.flags(std::ios::fixed | std::ios::right);
        (..., printSpace(std::cout, args));
        std::cout << std::endl;
        (..., printSpace(cao, args));
        cao << std::endl;
    }

    template<typename... Args>
    void warn(const Args &... args)
    {
        std::cout << "\x1B[33m";
        (..., printSpace(std::cout, args));
        std::cout << "\033[0m" << std::endl;
        (..., printSpace(cao, args));
        cao << std::endl;
    }

    template<typename... Args>
    void error(const Args &... args)
    {
        std::cout << "\x1B[31m";
        (..., printSpace(std::cout, args));
        std::cout << "\033[0m" << std::endl;
        (..., printSpace(cao, args));
        cao << std::endl;
    }

    template<typename... Args>
    void done(const Args &... args)
    {
        auto printSpace = [](std::ostream & os, const auto & val) -> void { os << val << " "; };
        std::cout << "\x1B[32m";
        (..., printSpace(std::cout, args));
        std::cout << "\033[0m" << std::endl;
        (..., printSpace(cao, args));
        cao << std::endl;
    }
};

#endif // LOG_H_
