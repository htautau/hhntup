#include "StrUtil.h"
namespace AnalysisFramework
{
namespace External
{
namespace GRLStrUtil {

  void
  trim (string& input) {
    // trim leading and trailing whitespace
    string::size_type position = input.find_first_not_of(" \t\n");
    if ( position == std::string::npos ) return; // skip, it's all whitespace
    input.erase(0, position);
    position= input.find_last_not_of(" \t\n");
    if ( position != std::string::npos)
      input.erase( position+1 );
  }

  void
  split (const string& input, string& first, string& second) {
    // split input in two
    string::size_type position = input.find_first_of(" \t\n");
    if ( position==std::string::npos ) {
      first  = input;
      second = "";
    } else {
      first  = input.substr(     0, position );
      second = input.substr( position+1, input.size()-position );
      // trim leading whitespace of second
      position= second.find_first_not_of(" \t\n");
      second.erase(0, position);
    }
  }

  vector<string>
  split (string input) {
    trim(input);
    vector<string> splitVec;
    string first, second;
    do {
      split(input,first,second);
      if (!first.empty()) splitVec.push_back(first);
      input = second;
    } while(!input.empty());
    return splitVec;
  }
}
}
}
