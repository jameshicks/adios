#ifndef FILEIOMANAGER_HPP
#define FILEIOMANAGER_HPP

#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <unistd.h>
#include <fstream>

#include <errno.h>

#include "config.h"

#ifdef HAVE_ZLIB
#include "zlib.h"
#endif

#define READBUF_SIZE 10000

class FileObject
{
public:
    virtual void openfile(const std::string& filename, bool write) = 0;
    virtual void closefile(void) = 0;
    virtual bool good(void) = 0;
    virtual bool eof(void) = 0;
    virtual std::string getline(void) = 0;
    virtual ~FileObject(void) {};

};

class UncompressedFile : public FileObject
{
public:
    FILE* f = NULL;
    std::string filename;
    UncompressedFile(void);
    UncompressedFile(const std::string& filename);
    ~UncompressedFile(void);
    void openfile(const std::string& filename, bool write);
    void closefile(void);
    bool good(void);
    bool eof(void);
    std::string getline(void);
};



#ifdef HAVE_ZLIB

class GZFile : public FileObject
{
public:
    gzFile gzf;
    std::string filename;
    GZFile(void);
    GZFile(const std::string filename);
    void openfile(const std::string& filename, bool write);
    void closefile(void);
    bool good(void);
    bool eof(void);
    std::string getline(void);

};


#endif

class DelimitedFileWriter : public UncompressedFile
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

class Logstream {
public:
    std::ofstream logfile;
    Logstream(const std::string& fn);
};

template <typename T> 
Logstream& operator<<(Logstream& l, const T& rhs ) {
    std::cout << rhs;
    l.logfile << rhs;
    l.logfile.flush();

    return l;
}

#endif