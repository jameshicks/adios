#include <fstream>
#include <iostream>
#include <string>
#include <string.h>
#include <stdlib.h>

int main(int argc, char** argv) {
    using std::cout; using std::endl;
    if (argc < 3) { std::cerr << "Not enough arguments" << std::endl; return 1; }

    std::string filename = argv[1];
    std::ifstream vcffile(filename);
    if (!vcffile) { std::cerr << "Cant open file: " << filename << std::endl;  return 1; }

    std::string desired = argv[2];

    std::string line; 
    while (std::getline(vcffile, line)) {
        if (line[0] == '#') continue;

        char* origdup = strdup(line.c_str());
        char* dup = origdup;

        int tokidx = -1;
        char* token;
        while (tokidx != 7) {
            token = strsep(&dup, "\t");
            tokidx++;
        }

        while (char* subtok = strsep(&token, ";")) {
            char* flag = strsep(&subtok, "=");
            if (desired.compare(flag) == 0) {
                char* val = strsep(&subtok,"=");
                cout << val << '\n';
                break;
            } 
        }

        free(origdup);
    }
}