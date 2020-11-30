#include "gen_anchors.h"
#include <string>

using namespace std;

int main()
{
  string xml_dir = R"(C:\Users\DELL\Desktop\Annotations\)";
  int width = 960;
  int height = 540;
  int num_of_clusters = 25;

  gen_anchors(xml_dir, num_of_clusters, width, height);

  system("pause");
  return 0;
}
