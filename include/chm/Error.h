#ifndef __BYTENOL_CHM_ERROR_H__
#define __BYTENOL_CHM_ERROR_H__

#include <stdexcept>

namespace chm
{

    class DBFailedException: std::logic_error
    {
        public:
            DBFailedException(const char* msg): std::logic_error(msg){}
    };

}

#endif 