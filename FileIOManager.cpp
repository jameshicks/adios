#include "FileIOManager.hpp"

UncompressedFile::UncompressedFile(void) { return; }

bool UncompressedFile::good(void) { return !(feof(f) || ferror(f)); }

bool UncompressedFile::eof(void) { return feof(f); }

void UncompressedFile::openfile(const std::string& fn, bool write)
{
    filename = fn;
    f = (filename.compare("-") == 0) ? stdout : fopen(filename.c_str(), write ? "w" : "r");

    if ((f == NULL) || feof(f) || ferror(f)) {
        throw std::invalid_argument("Couldn't open file: " + filename);
    }
}

void UncompressedFile::closefile(void)
{

    if (f == stdout || !f) { return; } // Dont close stdout or noop if already closed
    int r = fclose(f);
    if (r != 0) {
        puts(strerror(errno));
        throw std::runtime_error("Couldn't close file: " + filename);
    }
    f = NULL;
}

UncompressedFile::UncompressedFile(const std::string& filename)
{
    openfile(filename, false);
}

UncompressedFile::~UncompressedFile(void)
{
    closefile();
}

std::string UncompressedFile::getline(void)
{
    std::string line;
    char readbuf[READBUF_SIZE];

    size_t readlen = 0;
    while (fgets(readbuf, sizeof(readbuf), f) != NULL) {
        int e = ferror(f);

        readlen = strlen(readbuf);
        line.append(readbuf);
        if ((readbuf[readlen - 1] == '\n') || !good()) { break; }
    }
    if (line.back() == '\n') { line.pop_back(); }
    return line;

}

#ifdef HAVE_ZLIB
GZFile::GZFile(void) { return; }
GZFile::GZFile(const std::string filename)
{
    openfile(filename, false);
}

void GZFile::openfile(const std::string& fn, bool write)
{
    filename = fn;
    gzf = gzopen(filename.c_str(), "r");
    if (gzf == NULL) { throw std::invalid_argument("Couldn't open file: " + filename); }
}

void GZFile::closefile(void)
{
    int success = gzclose(gzf);
    if (success != Z_OK) { throw std::runtime_error("Couldn't close file: " + filename); }
}

bool GZFile::good(void)
{
    int errcode = 0;
    gzerror(gzf, &errcode);
    return !(eof() || (errcode != Z_OK));
}

bool GZFile::eof(void) { return gzeof(gzf); }

std::string GZFile::getline(void)
{
    std::string line;
    char readbuf[READBUF_SIZE];
    size_t readlen = 0;

    while (gzgets(gzf, readbuf, sizeof(readbuf))) {
        readlen = strlen(readbuf);
        line.append(readbuf);
        if ((readbuf[readlen - 1] == '\n') || !good()) { break; }
    }

    if (line.back() == '\n') { line.pop_back(); }
    return line;

}



#endif

DelimitedFileWriter::DelimitedFileWriter(const std::string& fn, char delimiter)
{
    openfile(fn, true);
    delim = delimiter;
}

void DelimitedFileWriter::write(const std::string& s)
{
    int res = 0;
    res = fputs(s.c_str(), f);
    if (res == EOF) { throw std::runtime_error("Couldn't write to file"); }
}

void DelimitedFileWriter::write(char d)
{
    int res = 0;
    res = fputc(d, f);
    if (res == EOF) { throw std::runtime_error("Couldn't write to file"); }
}

void DelimitedFileWriter::writetoks(const std::vector<std::string>& toks)
{
    for (size_t i = 0; i < toks.size(); ++i) {
        write(toks[i]);
        write((i != toks.size() - 1) ? delim : '\n');
    }
}

Logstream::Logstream(const std::string& fn) {
    logfile = std::ofstream(fn);
}