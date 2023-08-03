#include <iostream>
#include <cstdlib>

void report_error(const char* msg, const char* file, int line, const char* func_name, bool require)
{
    std::cout << file << ":" << line << ": " << msg << " (FAIL)" << std::endl;
    if (require)
      exit(-1);
}

const char* ob_expr(const char *expr) { return expr; }
