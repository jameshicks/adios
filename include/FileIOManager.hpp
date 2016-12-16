#include <string>
#include <stdexcept>
#include <stdio.h>
#include <string.h>

#define READBUF_SIZE 10000

class FileReader {
private:
    FILE* f;
    std::string filename;
public:
    inline bool good(void) { return !(feof(f) || ferror(f)); }
    inline bool eof(void) { return feof(f); }
    FileReader(const std::string& fn);
    ~FileReader(void);
    std::string getline(void);
};
