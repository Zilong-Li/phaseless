#ifndef LOG_H_
#define LOG_H_

#include <fstream>
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
    Logger & warn(const S & val)
    {
        cao << val;
        std::cout << "\x1B[33m" << val << "\033[0m";
        return *this;
    }

    template<class S>
    Logger & error(const S & val)
    {
        cao << val;
        std::cout << "\x1B[31m" << val << "\033[0m";
        return *this;
    }

    template<class S>
    Logger & done(const S & val)
    {
        cao << val;
        std::cout << "\x1B[32m" << val << "\033[0m";
        return *this;
    }
};

#endif // LOG_H_
