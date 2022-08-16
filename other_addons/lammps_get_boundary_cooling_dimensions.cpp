#include "stringlib.h"

using namespace std;

int main(int argc, char **argv) {

    if(argc != 3) {
      cout << "This script is used in au np in silica problem.\n";
      cout << "Input: <dump file>  \n";
      cout << "       <width of the boundary cooling region> [Å]\n";
      cout << "Output: <minx> <maxx> <miny> <maxy> of the region\n\n";      
      cout << "usage: \n" << argv[0] << " <dump file> <width>\n";
      return 1;
    }

    string in_dump = argv[1];
    double width = stod(argv[2]);
    
    auto boxsizex_arr=split_array_2<double>(exec("grep -A 1 \"ITEM: BOX BOUNDS\" " + in_dump + " | tail -1"));
    auto boxsizey_arr=split_array_2<double>(exec("grep -A 2 \"ITEM: BOX BOUNDS\" " + in_dump + " | tail -1"));
    auto boxsizez_arr=split_array_2<double>(exec("grep -A 3 \"ITEM: BOX BOUNDS\" " + in_dump + " | tail -1"));

    cout << boxsizex_arr[0][0] + width << " ";
    cout << boxsizex_arr[0][1] - width << " ";    

    cout << boxsizey_arr[0][0] + width << " ";
    cout << boxsizey_arr[0][1] - width << "\n";    

    
}
