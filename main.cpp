#define _USE_MATH_DEFINES
#include <iostream>
#include "cloth_4_3_compute\cloth_4_3_compute.h"
using namespace std;


int main(int argc, char* argv[]) {
    cout << "Cloth simulation:" << endl;
    cloth_4_3_compute::main(argc, argv);
}