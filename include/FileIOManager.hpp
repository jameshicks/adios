#ifndef FILEIOMANAGER_HPP
#define FILEIOMANAGER_HPP

#include <string>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <stdio.h>
#include <string.h>
#include <unistd.h>

#include <errno.h>

#include "config.h"

#ifdef HAVE_ZLIB
#include "zlib.h"
#endif

#define READBUF_SIZE 10000

class FileWrapper {
public:
    virtual void openfile(const std::string& filename, bool write) =0;
    virtual void closefile(void) =0;
    virtual bool good(void) =0;
    virtual bool eof(void) =0;
    virtual ~FileWrapper(void) =0;

};

class CFileWrapper : virtual public FileWrapper
{
public:
    FILE* f = NULL;
    std::string filename;

    virtual ~CFileWrapper(void);
    void openfile(const std::string& filename, bool write);
    void closefile(void);
    bool good(void);
    bool eof(void);
};

class FileReader : virtual public FileWrapper
{
public:
    virtual std::string getline(void) =0;
    virtual ~FileReader(void);
};

class UncompressedFileReader : public CFileWrapper, virtual public FileReader {
    public:
    UncompressedFileReader(void);
    UncompressedFileReader(const std::string& fn);
    std::string getline(void);
};

#ifdef HAVE_ZLIB

class GZFileWrapper : virtual public FileWrapper {
public:
    gzFile gzf;
    std::string filename;
    GZFileWrapper(void);
    void openfile(const std::string& filename, bool write);
    void closefile(void);
    bool good(void);
    bool eof(void);    
};

class GZFileReader : public GZFileWrapper, public FileReader
{

public:
    GZFileReader(void);
    GZFileReader(const std::string& filename);

    std::string getline(void);
    

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