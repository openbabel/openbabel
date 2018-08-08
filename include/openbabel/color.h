#ifndef COLOR_H
#define COLOR_H

// the following snippet of code detects the current OS and
// defines the appropriate macro that is used to wrap some
// platform specific things
#if defined(_WIN32) || defined(_WIN64)
#   define COLOR_OS_WINDOWS
#elif defined(__APPLE__)
#   define COLOR_OS_MACOS
#elif defined(__unix__) || defined(__unix)
#   define COLOR_OS_LINUX
#else
#   error unsupported platform
#endif

#include <iostream>
#include <cstdlib>
#include <cstring>
#include <sstream>
#include <cstdio>

#if defined(COLOR_OS_MACOS) || defined(COLOR_OS_LINUX)
#   include <unistd.h>
#elif defined(COLOR_OS_WINDOWS)
#   include <io.h>
#   include <windows.h>
    #ifndef ENABLE_VIRTUAL_TERMINAL_PROCESSING
    #define ENABLE_VIRTUAL_TERMINAL_PROCESSING 0x0004
    #endif
#endif


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

namespace {
    bool isAllowed = false;    
    bool update(std::ostream &os,int v){
      if(&os == &std::cout){
        std::stringstream s;
        s << "\e[" << v << "m";
        state.coutstate = s.str();
        std::clog << state.clogstate;
        return true;
      }
      else if(&os == &std::clog){
        std::stringstream s;
        s << "\e[" << v << "m";
        state.clogstate = s.str();
        std::cout << state.coutstate;
        return true;
      }
      return false;
    }
    inline
    FILE* get_standard_stream(const std::ostream& stream)
    {
        if (&stream == &std::cout)
            return stdout;
        else if ((&stream == &std::cerr) || (&stream == &std::clog))
            return stderr;

        return 0;
    }
    bool isTerminal(const std::ostream& stream)
    {
        FILE* std_stream = get_standard_stream(stream);

            // Unfortunately, fileno() ends with segmentation fault
            // if invalid file descriptor is passed. So we need to
            // handle this case gracefully and assume it's not a tty
            // if standard stream is not detected, and 0 is returned.
            if (!std_stream)
                return false;

        #if defined(COLOR_OS_MACOS) || defined(COLOR_OS_LINUX)
            return ::isatty(fileno(std_stream));
        #elif defined(COLOR_OS_WINDOWS)
            return ::_isatty(_fileno(std_stream));
        #endif    
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
namespace init {
        int color()
        {
            #if defined(COLOR_OS_WINDOWS)
                HANDLE hOut = GetStdHandle(STD_OUTPUT_HANDLE);
                if (hOut == INVALID_HANDLE_VALUE)
                {
                    return GetLastError();
                }

                DWORD dwMode = 0;
                if (!GetConsoleMode(hOut, &dwMode))
                {
                    return GetLastError();
                }

                dwMode |= ENABLE_VIRTUAL_TERMINAL_PROCESSING;
                if (!SetConsoleMode(hOut, dwMode))
                {
                    return GetLastError();
                }
            #endif
                return 0;
        }
        void checkterm(std::ostream &s){
            isAllowed = isTerminal(s);
        }
    }
    std::ostream &operator<<(std::ostream &os, color::style v)
    {
        init::checkterm(os);   
        return isAllowed && update(os,static_cast<int>(v))? os << "\e[" << static_cast<int>(v) << "m" : os;
    }
    std::ostream &operator<<(std::ostream &os, color::fg v)
    {  
	init::checkterm(os); 
        return isAllowed && update(os,static_cast<int>(v))? os << "\e[" << static_cast<int>(v) << "m" : os;
    }
    std::ostream &operator<<(std::ostream &os, color::bg v)
    {  
	init::checkterm(os); 
        return isAllowed && update(os,static_cast<int>(v))? os << "\e[" << static_cast<int>(v) << "m" : os;
    }
}
#endif /* ifndef COLOR_H*/
