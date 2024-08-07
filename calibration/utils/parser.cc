#include <utils/parser.h>

Parser::Parser() {
    __parent_path__ = fs::current_path();
    __data_path__ = __parent_path__ / "data";

    if (!fs::exists(__data_path__)) {
        fs::create_directories(__data_path__);
    }
}

void Parser::saveData(const std::string& newFile, const array& arr) {
    fs::path __new_file_path__ = __data_path__ / newFile;
    std::ofstream file(__new_file_path__);

    Eigen::IOFormat csvFmt(Eigen::FullPrecision, Eigen::DontAlignCols, ", ", "\n");

    if (file.is_open()) {
        file << arr.format(csvFmt);
        file.close();
    } else {
        std::cerr << "Failed to open file at " << __new_file_path__ << std::endl;
    }
}

array Parser::readData(const std::string& file) {
    fs::path __file_path__ = __data_path__ / file;
    std::ifstream dataFile(__file_path__);

    if (!dataFile.is_open()) {
        std::cerr << "Failed to open file" << std::endl;
        return array();
    }

    std::vector<double> matrixEntries;
    std::string matrixRowStr;
    std::string matrixEntriesStr;
    int matrixRowCount = 0;

    while (getline(dataFile, matrixRowStr)) {
        std::stringstream matrixRowStrStream(matrixRowStr);
        while (getline(matrixRowStrStream, matrixEntriesStr, ',')) {
            matrixEntries.push_back(std::stod(matrixEntriesStr)); 
        }
        matrixRowCount++;
    }

    if (matrixRowCount == 0) {
        std::cerr << "No data found in file: " << __file_path__ << std::endl;
        return array();
    }

    return Eigen::Map<Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic, Eigen::RowMajor>>(
        matrixEntries.data(), matrixRowCount, matrixEntries.size() / matrixRowCount);
}