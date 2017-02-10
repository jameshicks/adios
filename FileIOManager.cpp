#include "FileIOManager.hpp"

bool CFileWrapper::good(void) { return !(feof(f) || ferror(f)); }
bool CFileWrapper::eof(void) { return feof(f); }

void CFileWrapper::openfile(const std::string& fn, bool write) {
    filename = fn;
    f = (filename.compare("-") == 0) ? stdout : fopen(filename.c_str(), write ? "w" : "r");
    if (f == NULL) throw std::invalid_argument("Couldn't open file: " + filename);
}

void CFileWrapper::closefile(void) {

    if (f == stdout || !f) return; // Dont close stdout or noop if already closed
    int r = fclose(f);
    if (r == EOF) {
        puts(strerror(errno));
        throw std::runtime_error("Couldn't close file: " + filename);
    }
    f = NULL;
}

CFileWrapper::~CFileWrapper(void) { closefile(); }

FileReader::FileReader(const std::string& fn) {
    openfile(fn, false);
}

std::string FileReader::getline(void) {
    std::string line;
    char readbuf[READBUF_SIZE];
    size_t readlen = 0;
    while (fgets(readbuf, sizeof(readbuf), f)) {
        readlen = strlen(readbuf);
        line.append(readbuf);
        if ((readbuf[readlen-1] == '\n') || !good()) break;
    } 
    if (line.back() == '\n') line.pop_back();
    return line;

}

#ifdef HAVE_ZLIB
GZFileReader::GZFileReader(const std::string& filename) {
    openfile(filename, false);
}

void GZFileReader::openfile(const std::string& fn, bool write) {
    filename = fn;
    gzf = gzopen(filename.c_str(), "r");
    if (gzf == NULL) throw std::invalid_argument("Couldn't open file: " + filename);
}

void GZFileReader::closefile(void) {
    int success = gzclose(gzf);
    if (success != Z_OK) throw std::runtime_error("Couldn't close file: " + filename);
}

std::string GZFileReader::getline(void) {
    std::string line;
    char readbuf[READBUF_SIZE];
    size_t readlen = 0;
    
    while (gzgets(gzf, readbuf, sizeof(readbuf))) {
        readlen = strlen(readbuf);
        line.append(readbuf);
        if ((readbuf[readlen-1] == '\n') || !good()) break;
    } 
    
    if (line.back() == '\n') line.pop_back();
    return line;

}

bool GZFileReader::good(void){
    int errcode = 0;
    gzerror(gzf, &errcode);
    return !(eof() || (errcode != Z_OK));
}

bool GZFileReader::eof(void) { return gzeof(gzf); }

#endif

DelimitedFileWriter::DelimitedFileWriter(const std::string& fn, char delimiter) { 
    openfile(fn, true);
    delim = delimiter;
}

void DelimitedFileWriter::write(const std::string& s) {
    int res = 0;
    res = fputs(s.c_str(), f);
    if (res == EOF) throw std::runtime_error("Couldn't write to file");
}

void DelimitedFileWriter::write(char d) {
    int res = 0;
    res = fputc(d, f);
    if (res == EOF) throw std::runtime_error("Couldn't write to file");
}

void DelimitedFileWriter::writetoks(const std::vector<std::string>& toks) {
    for (size_t i = 0; i < toks.size(); ++i) {
        write(toks[i]);
        write((i != toks.size()-1) ? delim : '\n');
    }
}