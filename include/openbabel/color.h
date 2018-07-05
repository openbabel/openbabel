#ifndef COLOR_H
#define COLOR_H

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <sstream>
extern "C" {
#include <unistd.h>
}
#include <stdio.h>

#define ANSI_COLOR_RED     "\x1b[31m"
#define ANSI_COLOR_GREEN   "\x1b[32m"
#define ANSI_COLOR_YELLOW  "\x1b[33m"
#define ANSI_COLOR_BLUE    "\x1b[34m"
#define ANSI_COLOR_MAGENTA "\x1b[35m"
#define ANSI_COLOR_CYAN    "\x1b[36m"
#define ANSI_COLOR_RESET   "\x1b[0m"
#define ANSI_COLOR_BOLD    "\x1b[1m"
#define ANSI_COLOR_DIM     "\x1b[2m"

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

struct streamstate {
    std::string coutstate;
    std::string cerrstate;
    std::string clogstate;
};
extern streamstate state;
bool update(std::ostream &os,int v);

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
        return isAllowed && update(os,static_cast<int>(v))? os << "\e[" << static_cast<int>(v) << "m" : os;
    }
    std::ostream &operator<<(std::ostream &os, color::fg v)
    {   
        return isAllowed && update(os,static_cast<int>(v))? os << "\e[" << static_cast<int>(v) << "m" : os;
    }
    std::ostream &operator<<(std::ostream &os, color::bg v)
    {   
        return isAllowed && update(os,static_cast<int>(v))? os << "\e[" << static_cast<int>(v) << "m" : os;
    }
    namespace init {
        void color()
        {
            isAllowed = isTerminal() && supportsColor() ? true : false;
        }
    }
}
#endif /* ifndef COLOR_H*/
