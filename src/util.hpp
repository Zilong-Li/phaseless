#ifndef UTIL_H_
#define UTIL_H_

#include <cstdint>
#include <string>
#include <vector>

//******************************************************************************
//                               String Utils
//******************************************************************************

inline std::vector<std::string> split_string(const std::string & s, const std::string & separators)
{
    std::vector<std::string> ret;
    bool is_seperator[256] = {false};
    for(auto & ch : separators)
    {
        is_seperator[(unsigned int)ch] = true;
    }
    int begin = 0;
    for(int i = 0; i <= (int)s.size(); i++)
    {
        if(is_seperator[(uint8_t)s[i]] || i == (int)s.size())
        {
            ret.push_back(std::string(s.begin() + begin, s.begin() + i));
            begin = i + 1;
        }
    }
    return ret;
}

inline std::string trim_string(const std::string & s)
{
    int begin = 0, end = (int)s.size();
    while(begin < end && s[begin] == ' ') begin++;
    while(begin < end && s[end - 1] == ' ') end--;
    return std::string(s.begin() + begin, s.begin() + end);
}

inline bool ends_with(std::string const & str, std::string const & ending)
{
    if(ending.size() > str.size())
        return false;
    else
        return std::equal(ending.begin(), ending.end(), str.end() - ending.size());
}

inline bool starts_with(std::string const & str, std::string const & ending)
{
    if(ending.size() > str.size())
        return false;
    else
        return std::equal(ending.begin(), ending.end(), str.begin());
}

#endif // UTIL_H_
