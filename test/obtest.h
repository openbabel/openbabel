#include <iostream>
#include <cstdlib>

#ifdef _MSC_VER
#define FUNCTION_SIGNATURE __FUNCSIG__
#else
#define FUNCTION_SIGNATURE __PRETTY_FUNCTION__
#endif

/*inline*/ void report_error(const char* msg, const char* file, int line, const char* func_name, bool require = false)
{
    std::cout << file << ":" << line << ": " << msg << " (FAIL)" << std::endl; 
    if (require)
      exit(-1);
}

template <typename T1, typename T2>
void ob_compare(T1 a, T2 b, const char *file, int line, const char *func_name)
{
  if (!(a == b))
    std::cout << file << ":" << line << ": " << a << " == " << b << " (FAIL)" << std::endl; 
}

#define OB_ASSERT(exp) \
  ( (exp) ? static_cast<void>(0) : report_error(#exp, __FILE__, __LINE__, FUNCTION_SIGNATURE, false) ) 

#define OB_REQUIRE(exp) \
  ( (exp) ? static_cast<void>(0) : report_error(#exp, __FILE__, __LINE__, FUNCTION_SIGNATURE, true) ) 

#define OB_COMPARE(a,b) \
  ob_compare(a, b, __FILE__, __LINE__, FUNCTION_SIGNATURE)
