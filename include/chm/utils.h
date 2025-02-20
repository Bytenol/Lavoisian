#ifndef __BYTENOL_CHM_UTILS_H__
#define __BYTENOL_CHM_UTILS_H__

#include <vector>
#include <string>

namespace chm
{

    std::vector<std::string> Split(const std::string& str, char delimeter);

    // template<typename T>
    // std::string GetError(const T& instance);
}

#endif 