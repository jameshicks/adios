#ifndef FILEIOMANAGER_HPP
#define FILEIOMANAGER_HPP

#include <string>
#include <vector>
#include <stdexcept>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>
#define READBUF_SIZE 10000

class CFileWrapper {
public:
    FILE* f;
    std::string filename;
    ~CFileWrapper(void);    
    void openfile(const std::string& filename, bool write);
    void closefile(void);
    inline bool good(void) { return !(feof(f) || ferror(f)); }
    inline bool eof(void) { return feof(f); }
};

class FileReader : public CFileWrapper {
public:
    FileReader(const std::string& fn);
    std::string getline(void);
};

class DelimitedFileWriter : public CFileWrapper {
private:
    void write(const std::string& s);
    void write(char v);
public:
    char delim;
    inline bool is_stdout(void) { return f == stdout; }
    DelimitedFileWriter(const std::string& fn, char delimiter);
    void writetoks(const std::vector<std::string>& toks);
};

#endif