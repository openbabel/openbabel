#ifdef WIN32
  #include <Windows.h> // QueryPerformanceFrequency & QueryPerformanceCounter
#else
  #include <sys/timeb.h> // Linux and MacOSX
#endif

#include <iostream>
#include <sstream>
#include <fstream>
#include <ctime>

#include "obtest.h"

static bool headerWritten = false;

class BenchmarkLoop 
{
  public:
    BenchmarkLoop(const std::string &name = std::string())
    {
#ifdef WIN32
      QueryPerformanceFrequency(&m_frequency); // counts per second
#endif
      m_time = 0;
      m_loops = 0;
      m_start = getMilliCount();
      m_name = name;
    }
    ~BenchmarkLoop()
    {
      unsigned int milliSecsPerIter = m_time / m_loops;
      
      if (!m_name.empty()) {
        std::cout << "Benchmark: " << m_name << std::endl;

        // create header
        std::time_t now;
        std::time(&now);
        std::string dateStr = std::ctime(&now);
        if (dateStr[dateStr.size()-1] == '\n')
          dateStr = dateStr.substr(0, dateStr.size()-1);
        std::string header = "\"" + dateStr + "\"";

        // write time to a file so progress can be monitored
        std::ofstream ofs;
        ofs.open("benchmark.txt.new");
 
        // read the current benchmark.txt, update lines and write to .new file
        std::ifstream ifs;
        ifs.open("benchmark.txt");
        std::string line;
        bool foundBenchmark = false, foundHeader = false;
        if (headerWritten)
          foundHeader = true;
        while (std::getline(ifs, line)) {
          if (line.substr(0, 10) == "__Header__" && !headerWritten) {
            foundHeader = true;
            headerWritten = true;
            std::stringstream ss;
            ss << line << " " << header;
            line = ss.str();          
          }
          if (line.substr(0, m_name.size()) == m_name) {
            foundBenchmark = true;
            std::stringstream ss;
            ss << line << " " << milliSecsPerIter;
            line = ss.str();
          }
          
          ofs << line << std::endl;
        }
        ifs.close();
        ofs.close();

        if (!foundBenchmark) {
          // benchmark not found, append the new line
          ofs.open("benchmark.txt", std::ios::out | std::ios::app);
          if (!foundHeader)
            ofs << "__Header__ " << header << std::endl;
          ofs << m_name << " " << milliSecsPerIter << std::endl;
        } else {
          // benchmark found and updated .new file, copy it now
          ifs.open("benchmark.txt.new");
          ofs.open("benchmark.txt");
          if (!foundHeader)
            ofs << "__Header__" << header << std::endl;
          ofs << ifs.rdbuf();
        }
     }
      
      // print out in Xms, Xs or XmYs
      std::stringstream duration;
      if (milliSecsPerIter > 1000) {
        unsigned int secsPerIter = milliSecsPerIter / 1000;
        if (secsPerIter > 60) {
          unsigned int minutesPerIter = secsPerIter / 60;
          duration << minutesPerIter << "m" << secsPerIter % 60 << "s";        
        } else {
          duration << secsPerIter << "s";
        }
      } else {
        duration << milliSecsPerIter << "ms";
      }
      std::cout << duration.str() << " per iteration (" << m_loops << " iterations total)" << std::endl;
    }
    bool done()
    {
      if (m_time > 500) // iterate untill at least 500ms passed
        return true;
      return false;
    }
    void next()
    {
      m_time = getMilliSpan(m_start);
      ++m_loops;
    }
    unsigned int getMilliCount()
    {
#ifdef WIN32
      LARGE_INTEGER performanceCount;
      QueryPerformanceCounter(&performanceCount); // value in counts
      LARGE_INTEGER nCount = performanceCount / m_frequency * 1000;
      return nCount;
#else
      timeb tb;
      ftime( &tb );
      unsigned int nCount = tb.millitm + (tb.time & 0xfffff) * 1000;
      return nCount;
#endif
    }
    unsigned int getMilliSpan(unsigned int start)
    {
      int nSpan = getMilliCount() - start;
#ifndef WIN32
      if ( nSpan < 0 ) // handle overflow...
        nSpan += 0x100000 * 1000;
#endif
      return nSpan;
    }
  private:
    unsigned int m_time, m_start, m_loops;
    std::string m_name;
#ifdef WIN32
    LARGE_INTEGER m_frequency;
#endif
};

#define OB_BENCHMARK \
  for (BenchmarkLoop loop; !loop.done(); loop.next())

#define OB_NAMED_BENCHMARK(name) \
  for (BenchmarkLoop loop(#name); !loop.done(); loop.next())
