extern "C" {
#include <openbabel/color.h>
}
streamstate state = {"","",""};

// bool update(std::ostream &os,int v){
//   if(&os == &std::cout){
//     std::stringstream s;
//     s << "\e[" << v << "m";
//     state.coutstate = s.str();
//     std::clog << state.clogstate;
//     return true;
//   }
//   else if(&os == &std::clog){
//     std::stringstream s;
//     s << "\e[" << v << "m";
//     state.clogstate = s.str();
//     std::cout << state.coutstate;
//     return true;
//   }
//   return false;
// }
