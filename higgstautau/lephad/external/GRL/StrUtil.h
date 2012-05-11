#include <string>
#include <vector>

using namespace std;
namespace AnalysisFramework
{
namespace External
{
namespace GRLStrUtil {
  void trim (string& input);
  void split (const string& input, string& first, string& second);
  vector<string> split(string input);
}
}
}
