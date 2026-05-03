// src/io/topology/TopParserStringUtils.cpp

#include "TopParserStringUtils.hpp"
#include "TopParserMain.hpp"
#include "TopParserUtilities.hpp"
#include <chrono>
#include <filesystem>
#include <fstream>
#include <random>
#include <stdexcept>
#include <system_error>

namespace pygcmc {
namespace io {

namespace {

class TemporaryTopologyFile {
public:
    TemporaryTopologyFile(const std::string& stem, const std::string& suffix, const std::string& content)
        : directory_(createTemporaryDirectory(stem)), file_(directory_ / (stem + suffix)) {
        std::ofstream out(file_);
        if (!out) {
            throw std::runtime_error("Failed to create temporary file for topology parsing");
        }
        out << content;
    }

    ~TemporaryTopologyFile() {
        std::error_code ec;
        std::filesystem::remove_all(directory_, ec);
    }

    const std::filesystem::path& path() const {
        return file_;
    }

private:
    static std::filesystem::path createTemporaryDirectory(const std::string& stem) {
        const auto base = std::filesystem::temp_directory_path();
        std::random_device random_device;

        for (int attempt = 0; attempt < 100; ++attempt) {
            const auto now = std::chrono::steady_clock::now().time_since_epoch().count();
            const auto token = std::to_string(now) + "-" + std::to_string(random_device()) + "-" + std::to_string(attempt);
            auto candidate = base / ("pygcmc-" + stem + "-" + token);

            std::error_code ec;
            if (std::filesystem::create_directory(candidate, ec)) {
                return candidate;
            }
        }

        throw std::runtime_error("Failed to create unique temporary directory for topology parsing");
    }

    std::filesystem::path directory_;
    std::filesystem::path file_;
};

} // namespace

model::Topology TopParserStringUtils::parse_string(const std::string& top_str) {
    // Handle empty string case - first trim whitespace
    const std::string trimmed_str = TopParserUtilities::trim(top_str);
    if (trimmed_str.empty()) {
        return model::Topology();
    }

    TemporaryTopologyFile temp_file("topology", ".top", top_str);
    return TOPParser::parse_file(temp_file.path().string());
}

} // namespace io
} // namespace pygcmc
