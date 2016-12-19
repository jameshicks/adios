#include <string>
#include <stdexcept>
#include <stdio.h>
#include <string.h>

#define READBUF_SIZE 10000

class CFileWrapper {
public:
    FILE* f;
    std::string filename;
    ~CFileWrapper(void);    
    void openfile(const std::string& filename, const std::string& mode);
    void closefile(void);
    inline bool good(void) { return !(feof(f) || ferror(f)); }
    inline bool eof(void) { return feof(f); }
};

class FileReader : public CFileWrapper {
public:
    FileReader(const std::string& fn);
    std::string getline(void);
};

