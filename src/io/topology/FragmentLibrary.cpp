#include "FragmentLibrary.hpp"
#include <fstream>
#include <sstream>

namespace pygcmc {
namespace io {
namespace topology {

void FragmentLibrary::addTemplate(const TemplateData& t) {
    byName_[t.name] = t;
    if (t.typeId >= 0) byType_[t.typeId] = t.name;
}

bool FragmentLibrary::loadFromITP(const std::string& path, const std::string& name, int typeId) {
    // Minimal stub: create placeholder if file exists
    std::ifstream f(path);
    if (!f.good()) return false;

    TemplateData t;
    t.name = name;
    t.typeId = typeId;
    // NOTE: Real parser to extract atoms/bonds. For now: empty atoms (caller may fill).
    addTemplate(t);
    return true;
}

bool FragmentLibrary::loadFromDirectory(const std::string& /*dir*/) {
    // Optional: not implemented in stub
    return false;
}

const FragmentLibrary::TemplateData* FragmentLibrary::get(const std::string& name) const {
    auto it = byName_.find(name);
    return it == byName_.end() ? nullptr : &it->second;
}

const FragmentLibrary::TemplateData* FragmentLibrary::getByType(int typeId) const {
    auto it = byType_.find(typeId);
    if (it == byType_.end()) return nullptr;
    return get(it->second);
}

} // namespace topology
} // namespace io
} // namespace pygcmc