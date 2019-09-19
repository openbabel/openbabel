#ifndef OB_TEST_H
#define OB_TEST_H

#include <iostream>
#include <cstdlib>

#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/shared_ptr.h>

#ifdef _MSC_VER
#define FUNCTION_SIGNATURE __FUNCSIG__
#else
#define FUNCTION_SIGNATURE __PRETTY_FUNCTION__
#endif

void report_error(const char* msg, const char* file, int line, const char* func_name, bool require = false);

template <typename T1, typename T2>
void ob_compare(T1 a, T2 b, const char *expr, const char *file, int line, const char *func_name)
{
  if (!(a == b))
    std::cout << file << ":" << line << ": " << expr << " [" << a << " == " << b << "] (FAIL)" << std::endl;
}

#define OB_ASSERT(exp) \
  ( (exp) ? static_cast<void>(0) : report_error(#exp, __FILE__, __LINE__, FUNCTION_SIGNATURE, false) )

#define OB_REQUIRE(exp) \
  ( (exp) ? static_cast<void>(0) : report_error(#exp, __FILE__, __LINE__, FUNCTION_SIGNATURE, true) )

const char* ob_expr(const char *expr);
#define OB_EXPR(expr) ob_expr(#expr)

#define OB_COMPARE(a,b) \
  ob_compare(a, b, OB_EXPR( a == b ), __FILE__, __LINE__, FUNCTION_SIGNATURE)



// some utility functions
typedef obsharedptr<OpenBabel::OBMol> OBMolPtr;

struct OBTestUtil
{
  static std::string GetFilename(const std::string &filename)
  {
    std::string path = TESTDATADIR + filename;
    return path;
  }

  static OBMolPtr ReadFile(const std::string &filename)
  {
    std::string file = GetFilename(filename);

    std::ifstream ifs;
    ifs.open(file.c_str());
    OB_REQUIRE( ifs );

    OpenBabel::OBConversion conv;
    OpenBabel::OBFormat *format = conv.FormatFromExt(file.c_str());
    OB_REQUIRE(format);
    OB_REQUIRE(conv.SetInFormat(format));

    OBMolPtr mol(new OpenBabel::OBMol);
    OB_REQUIRE(conv.Read(mol.get(), &ifs));

    return mol;
  }

  static std::string ReadFileContent(const std::string &filename)
  {
    std::string fn = GetFilename(filename);

    std::ifstream ifs(fn.c_str());
    std::stringstream buffer;
    buffer << ifs.rdbuf();
    std::string content = buffer.str();

    return content;
  }
};

#endif // OB_TEST_H
