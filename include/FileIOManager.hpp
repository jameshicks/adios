#ifndef FILEIOMANAGER_HPP
#define FILEIOMANAGER_HPP

#include <string>
#include <vector>
#include <stdexcept>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <errno.h>

#include "config.h"

#ifdef HAVE_ZLIB
#include "zlib.h"
#endif

#define READBUF_SIZE 10000

class CFileWrapper
{
public:
    FILE* f = NULL;
    std::string filename;
    ~CFileWrapper(void);
    virtual void openfile(const std::string& filename, bool write);
    virtual void closefile(void);
    virtual bool good(void);
    virtual bool eof(void);
};

class FileReader : public CFileWrapper
{
public:
    FileReader(void) { return; }
    FileReader(const std::string& fn);
    std::string getline(void);
};

#ifdef HAVE_ZLIB

class GZFileReader : public FileReader
{

private:
    gzFile gzf;
public:
    GZFileReader(const std::string& filename);
    void openfile(const std::string& filename, bool write);
    void closefile(void);
    std::string getline(void);
    
    bool good(void);

    bool eof(void);
};


#endif

class DelimitedFileWriter : public CFileWrapper
{
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