#include "FileIOManager.hpp"

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