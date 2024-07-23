#ifndef PARSER_H
#define PARSER_H

#include <utils/__types__.h>

#include <filesystem>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

namespace fs = std::filesystem;

class Parser {
   public:
    Parser();
    ~Parser() = default;

    void saveData(const std::string &filename, const array &data);

    array readData(const std::string &filename);

   private:
    fs::path __parent_path__;
    fs::path __data_path__;
};

#endif  // PARSER_H