#include "FileIOManager.hpp"

void CFileWrapper::openfile(const std::string& fn, const std::string& mode) {
    filename = fn;
    f = fopen(filename.c_str(), mode.c_str());
    if (f == NULL) throw std::invalid_argument("Couldn't open file: " + filename);
}

void CFileWrapper::closefile(void) {
    int r = fclose(f);
    if (r == EOF) throw std::runtime_error("Couldn't close file" + filename);
}

CFileWrapper::~CFileWrapper(void) { closefile(); }

FileReader::FileReader(const std::string& fn) {
    openfile(fn, "r");
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