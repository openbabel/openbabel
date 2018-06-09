#ifndef COLOR_H
#define COLOR_H

#include <iostream>
#include <cstdlib>
#include <cstring>
extern "C" {
#include <unistd.h>
}

namespace color {

    enum style {
        Reset     = 0,
        bold      = 1,
        dim       = 2,
        italic    = 3,
        underline = 4,
        blink     = 5,
        reversed  = 6,
        conceal   = 7,
        crossed   = 8
    };
    enum fg {
        def     = 39,
        black   = 30,
        red     = 31,
        green   = 32,
        yellow  = 33,
        blue    = 34,
        magenta = 35,
        cyan    = 36,
        gray    = 37
    };
    enum bg {
        def1     = 49,
        black1   = 40,
        red1     = 41,
        green1   = 42,
        yellow1  = 43,
        blue1    = 44,
        magenta1 = 45,
        cyan1    = 46,
        gray1    = 47
    };
}


namespace {
    bool isAllowed = false;
    bool isTerminal()
    {
        return isatty(STDERR_FILENO);
    }
    bool supportsColor()
    {
    	if(const char *env_p = std::getenv("TERM")) {
        	const char *const terms[] = {
       		    "xterm", "xterm-256", "xterm-256color", "vt100",
       		    "color", "ansi",      "cygwin",         "linux"};
		for(int i = 0; i < 8; i++ ){
			if(std::strcmp(env_p, terms[i]) == 0) return true;
	   	}
    	}
        return false;
    }
    std::ostream &operator<<(std::ostream &os, color::style v)
    {
        return isAllowed ? os << "\e[" << static_cast<int>(v) << "m" : os;
    }
    std::ostream &operator<<(std::ostream &os, color::fg v)
    {
        return isAllowed ? os << "\e[" << static_cast<int>(v) << "m" : os;
    }
    std::ostream &operator<<(std::ostream &os, color::bg v)
    {
        return isAllowed ? os << "\e[" << static_cast<int>(v) << "m" : os;
    }
    namespace init {
        void color()
        {
            isAllowed = isTerminal() && supportsColor() ? true : false;
        }
    }
}
#endif /* ifndef COLOR_H*/

