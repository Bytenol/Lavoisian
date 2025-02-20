#include <chm/utils.h>

namespace chm
{
    std::vector<std::string> Split(const std::string& str, char delimeter)
    {
        std::vector<std::string> res;
        std::string text = str;
        size_t n;
        while((n = text.find(delimeter)) < text.size()) {
            res.push_back(text.substr(0, n));
            text = text.substr(n + 1, text.size());
        }
        res.push_back(text);
        return res;
    }
}

// template<typename T>
// std::string GetError(const T& instance)
// {
//     return instance.error;
// }