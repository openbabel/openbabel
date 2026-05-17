/**********************************************************************
ketformat.cpp - Read & write Ketcher KET (JSON) chemical structure files.

The KET format is Ketcher's native serialization format, also produced and
consumed by Indigo. It is a JSON document with a top-level "root" object that
references molecules, R-groups, reactions, monomers, templates, and a variety
of meta objects (arrows, plus signs, text, simple shapes, images).

This format supports:
  - molecules (atoms, bonds, S-groups, R-group attachment sites)
  - R-groups
  - reactions assembled from arrow + plus + molecule positions
  - reaction agents (components placed near the arrow but not on either side)
  - verbatim preservation of unknown top-level sections (monomers,
    macromolecule templates, monomer shapes, V2 rich text, images, etc.)
    so round-trips don't drop Ketcher-only data.

For information about the spec, see https://github.com/epam/Indigo and the
KET format reference in this repository (ket-format.md).

This file is part of the Open Babel project.
For more information, see <http://openbabel.org/>

This program is free software; you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation version 2 of the License.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
***********************************************************************/

#include <openbabel/atom.h>
#include <openbabel/alias.h>
#include <openbabel/babelconfig.h>
#include <openbabel/bond.h>
#include <openbabel/elements.h>
#include <openbabel/generic.h>
#include <openbabel/json.h>
#include <openbabel/mol.h>
#include <openbabel/obconversion.h>
#include <openbabel/obiter.h>
#include <openbabel/obmolecformat.h>
#include <openbabel/reactionfacade.h>
#include <openbabel/stereo/stereo.h>

#include <algorithm>
#include <array>
#include <cmath>
#include <cstdint>
#include <map>
#include <memory>
#include <set>
#include <sstream>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

namespace OpenBabel {

namespace {

// ============================================================================
//  Constants
// ============================================================================

// OBPairData attribute names used to persist KET-specific atom/bond/mol data
// that has no native OpenBabel representation.
constexpr const char *kAttrKetAlias              = "_ket_alias";
constexpr const char *kAttrKetCip                = "_ket_cip";
constexpr const char *kAttrKetStereoLabel        = "_ket_stereoLabel";
constexpr const char *kAttrKetAttachmentPoints   = "_ket_attachmentPoints";
constexpr const char *kAttrKetMapping            = "_ket_mapping";
constexpr const char *kAttrKetInvRet             = "_ket_invRet";
constexpr const char *kAttrKetExactChangeFlag    = "_ket_exactChangeFlag";
constexpr const char *kAttrKetQueryProperties    = "_ket_queryProperties";
constexpr const char *kAttrKetHCount             = "_ket_hCount";
constexpr const char *kAttrKetRingBondCount      = "_ket_ringBondCount";
constexpr const char *kAttrKetSubstitutionCount  = "_ket_substitutionCount";
constexpr const char *kAttrKetUnsaturatedAtom    = "_ket_unsaturatedAtom";
constexpr const char *kAttrKetSelected           = "_ket_selected";
constexpr const char *kAttrKetExplicitValence    = "_ket_explicitValence";
constexpr const char *kAttrKetImplicitHCountSet  = "_ket_implicitHCountSet";

constexpr const char *kAttrKetBondTopology       = "_ket_bond_topology";
constexpr const char *kAttrKetBondCenter         = "_ket_bond_center";
constexpr const char *kAttrKetBondCustomQuery    = "_ket_bond_customQuery";
constexpr const char *kAttrKetBondType           = "_ket_bond_type";  // raw KET integer (4-10)
constexpr const char *kAttrKetBondStereobox      = "_ket_bond_stereobox";

constexpr const char *kAttrKetVersion            = "_ket_version";
constexpr const char *kAttrKetPassthrough        = "_ket_passthrough";
constexpr const char *kAttrKetStereoFlagPos      = "_ket_stereoFlagPosition";
constexpr const char *kAttrKetRGroupNumber       = "_ket_rgroupNumber";
// Set on the first atom of each parsed molecule fragment so the writer can
// reuse the molecule's original $ref name (mol0, mol1, ...) on round-trip.
constexpr const char *kAttrKetMolRef             = "_ket_molRef";

// Bond type integers from the KET spec.
enum KetBondType : int {
    KET_BOND_SINGLE             = 1,
    KET_BOND_DOUBLE             = 2,
    KET_BOND_TRIPLE             = 3,
    KET_BOND_AROMATIC           = 4,
    KET_BOND_SINGLE_OR_DOUBLE   = 5,
    KET_BOND_SINGLE_OR_AROMATIC = 6,
    KET_BOND_DOUBLE_OR_AROMATIC = 7,
    KET_BOND_ANY                = 8,
    KET_BOND_COORDINATION       = 9,
    KET_BOND_HYDROGEN           = 10
};

// Bond stereo integers (Biovia convention).
enum KetBondStereo : int {
    KET_STEREO_NONE     = 0,
    KET_STEREO_UP       = 1,
    KET_STEREO_CISTRANS = 3,
    KET_STEREO_EITHER   = 4,
    KET_STEREO_DOWN     = 6
};

// ============================================================================
//  Helper utilities
// ============================================================================

inline bool getMemberDouble(const rapidjson::Value &v, const char *name, double &out)
{
    if (!v.HasMember(name)) return false;
    const auto &m = v[name];
    if (m.IsNumber()) { out = m.GetDouble(); return true; }
    return false;
}

inline bool getMemberInt(const rapidjson::Value &v, const char *name, int &out)
{
    if (!v.HasMember(name)) return false;
    const auto &m = v[name];
    if (m.IsInt())    { out = m.GetInt();    return true; }
    if (m.IsUint())   { out = static_cast<int>(m.GetUint());   return true; }
    if (m.IsNumber()) { out = static_cast<int>(m.GetDouble()); return true; }
    return false;
}

inline bool getMemberBool(const rapidjson::Value &v, const char *name, bool &out)
{
    if (!v.HasMember(name)) return false;
    const auto &m = v[name];
    if (m.IsBool()) { out = m.GetBool(); return true; }
    return false;
}

inline bool getMemberString(const rapidjson::Value &v, const char *name, std::string &out)
{
    if (!v.HasMember(name)) return false;
    const auto &m = v[name];
    if (m.IsString()) { out = m.GetString(); return true; }
    return false;
}

// Stringify a sub-document as a compact JSON string (used for raw passthrough).
std::string stringifyValue(const rapidjson::Value &v)
{
    rapidjson::StringBuffer buf;
    rapidjson::Writer<rapidjson::StringBuffer> writer(buf);
    v.Accept(writer);
    return std::string(buf.GetString(), buf.GetSize());
}

// Set or replace a string OBPairData attribute on a generic OBBase target.
void setPairData(OBBase *target, const char *attr, const std::string &value)
{
    if (auto *existing = target->GetData(attr)) {
        if (auto *p = dynamic_cast<OBPairData *>(existing)) {
            p->SetValue(value);
            return;
        }
        target->DeleteData(existing);
    }
    auto *pd = new OBPairData();
    pd->SetAttribute(attr);
    pd->SetValue(value);
    pd->SetOrigin(fileformatInput);
    target->SetData(pd);
}

void setPairDataInt(OBBase *target, const char *attr, int value)
{
    setPairData(target, attr, std::to_string(value));
}

bool getPairDataString(OBBase *target, const char *attr, std::string &out)
{
    auto *gd = target->GetData(attr);
    if (!gd) return false;
    auto *p = dynamic_cast<OBPairData *>(gd);
    if (!p) return false;
    out = p->GetValue();
    return true;
}

bool getPairDataInt(OBBase *target, const char *attr, int &out)
{
    std::string s;
    if (!getPairDataString(target, attr, s)) return false;
    try { out = std::stoi(s); return true; }
    catch (...) { return false; }
}

bool ketVersionMajor(const std::string &version, int &major)
{
    if (version.empty()) return false;
    size_t dot = version.find('.');
    const std::string majorPart =
        (dot == std::string::npos) ? version : version.substr(0, dot);
    try {
        major = std::stoi(majorPart);
        return true;
    } catch (...) {
        return false;
    }
}

bool ketVersionIsSupported(const std::string &version)
{
    int major = 0;
    return !ketVersionMajor(version, major) || major <= 2;
}

// Convert an isotope keyword to (element, mass) for D and T.
bool decodeIsotopeShorthand(const std::string &label, int &elem, int &iso)
{
    if (label == "D") { elem = OBElements::Hydrogen; iso = 2; return true; }
    if (label == "T") { elem = OBElements::Hydrogen; iso = 3; return true; }
    return false;
}

// ============================================================================
//  Reader state
// ============================================================================

// Geometric placement of a molecule fragment relative to an arrow.
struct ComponentBounds {
    int firstAtomIdx = 0;   // 1-based first atom in the merged OBMol
    int lastAtomIdx  = 0;   // 1-based last atom (inclusive)
    double centroidX = 0.0;
    double centroidY = 0.0;
};

struct Vec3d { double x = 0.0, y = 0.0, z = 0.0; };

struct ArrowGeometry {
    Vec3d tail;
    Vec3d head;
    std::string mode;
    bool valid = false;
};

struct PlusGeometry {
    Vec3d pos;
};

// One entry in the original `root.nodes` array — captured so the writer can
// re-emit nodes in the same order, with the same shape (ref vs inline).
struct NodeOrderEntry {
    enum class Kind { Ref, Inline };
    Kind kind = Kind::Ref;
    // For Ref: the referenced name (e.g. "mol0"). For Inline: raw JSON of the
    // inline meta object.
    std::string value;
};

// Tracks unknown top-level objects and source-document layout so the writer
// can reconstruct the file with maximum fidelity.
struct KetPassthrough {
    // Ordered list of root.nodes entries — both known refs (mol/rgroup) and
    // every inline meta object. The writer walks this to preserve original
    // ordering of nodes (z-order, grouping).
    std::vector<NodeOrderEntry> nodeOrder;
    // Ordered list of root.templates ref names.
    std::vector<std::string> templateOrder;
    // Raw JSON of root.connections (the spec calls these monomer/polymer
    // connectivity; we don't translate them).
    std::string connections;
    // Raw JSON of root.annotation, if any.
    std::string documentAnnotation;
    // unknownRefs: ordered list of (refName, rawJSON) for $ref nodes whose
    // top-level body we don't natively translate (monomer, monomerTemplate,
    // ambiguousMonomer, monomerShape, etc.). Their refName entries appear in
    // `nodeOrder` (or `templateOrder` for templates); the body lives here.
    std::vector<std::pair<std::string, std::string>> unknownRefs;
    // unknownTopMembers: top-level members other than root / ket_version /
    // any known ref. Echoed verbatim on output.
    std::vector<std::pair<std::string, std::string>> unknownTopMembers;

    bool empty() const
    {
        return nodeOrder.empty() && templateOrder.empty() &&
               connections.empty() && documentAnnotation.empty() &&
               unknownRefs.empty() && unknownTopMembers.empty();
    }

    std::string serialize() const
    {
        rapidjson::Document d(rapidjson::kObjectType);
        auto &al = d.GetAllocator();

        // nodeOrder: array of [kind, value] pairs (kind = "ref" or "inline").
        rapidjson::Value nodes(rapidjson::kArrayType);
        for (const auto &e : nodeOrder) {
            rapidjson::Value pair(rapidjson::kArrayType);
            pair.PushBack(rapidjson::Value(
                e.kind == NodeOrderEntry::Kind::Ref ? "ref" : "inline",
                al).Move(), al);
            pair.PushBack(rapidjson::Value(e.value.c_str(), al).Move(), al);
            nodes.PushBack(pair, al);
        }
        d.AddMember("nodeOrder", nodes, al);

        auto pushStringArray = [&](const char *key,
                                   const std::vector<std::string> &items) {
            rapidjson::Value arr(rapidjson::kArrayType);
            for (const auto &s : items)
                arr.PushBack(rapidjson::Value(s.c_str(), al).Move(), al);
            d.AddMember(rapidjson::Value(key, al), arr, al);
        };
        auto pushStringPairs =
            [&](const char *key,
                const std::vector<std::pair<std::string, std::string>> &items) {
                rapidjson::Value arr(rapidjson::kArrayType);
                for (const auto &kv : items) {
                    rapidjson::Value pair(rapidjson::kArrayType);
                    pair.PushBack(
                        rapidjson::Value(kv.first.c_str(), al).Move(), al);
                    pair.PushBack(
                        rapidjson::Value(kv.second.c_str(), al).Move(), al);
                    arr.PushBack(pair, al);
                }
                d.AddMember(rapidjson::Value(key, al), arr, al);
            };
        pushStringArray("templateOrder", templateOrder);
        pushStringPairs("refs", unknownRefs);
        pushStringPairs("topMembers", unknownTopMembers);
        if (!connections.empty())
            d.AddMember("connections",
                        rapidjson::Value(connections.c_str(), al), al);
        if (!documentAnnotation.empty())
            d.AddMember("annotation",
                        rapidjson::Value(documentAnnotation.c_str(), al), al);

        rapidjson::StringBuffer buf;
        rapidjson::Writer<rapidjson::StringBuffer> writer(buf);
        d.Accept(writer);
        return std::string(buf.GetString(), buf.GetSize());
    }

    static KetPassthrough deserialize(const std::string &raw)
    {
        KetPassthrough p;
        if (raw.empty()) return p;
        rapidjson::Document d;
        if (d.Parse(raw.c_str()).HasParseError() || !d.IsObject())
            return p;

        if (d.HasMember("nodeOrder") && d["nodeOrder"].IsArray()) {
            for (auto &v : d["nodeOrder"].GetArray()) {
                if (v.IsArray() && v.Size() == 2 &&
                    v[0].IsString() && v[1].IsString()) {
                    NodeOrderEntry e;
                    e.kind = (std::string(v[0].GetString()) == "ref")
                                 ? NodeOrderEntry::Kind::Ref
                                 : NodeOrderEntry::Kind::Inline;
                    e.value = v[1].GetString();
                    p.nodeOrder.push_back(std::move(e));
                }
            }
        }
        auto readStringArray = [&](const char *key,
                                   std::vector<std::string> &dst) {
            if (!d.HasMember(key) || !d[key].IsArray()) return;
            for (auto &v : d[key].GetArray())
                if (v.IsString()) dst.emplace_back(v.GetString());
        };
        auto readStringPairs =
            [&](const char *key,
                std::vector<std::pair<std::string, std::string>> &dst) {
                if (!d.HasMember(key) || !d[key].IsArray()) return;
                for (auto &v : d[key].GetArray()) {
                    if (v.IsArray() && v.Size() == 2 &&
                        v[0].IsString() && v[1].IsString())
                        dst.emplace_back(v[0].GetString(), v[1].GetString());
                }
            };
        readStringArray("templateOrder", p.templateOrder);
        readStringPairs("refs", p.unknownRefs);
        readStringPairs("topMembers", p.unknownTopMembers);
        if (d.HasMember("connections") && d["connections"].IsString())
            p.connections = d["connections"].GetString();
        if (d.HasMember("annotation") && d["annotation"].IsString())
            p.documentAnnotation = d["annotation"].GetString();
        return p;
    }
};

// ============================================================================
//  Reader
// ============================================================================

class KetReader {
public:
    explicit KetReader(OBMol *mol, OBConversion *conv) : mol_(mol)
    {
        expandAliases_ = conv && conv->IsOption("a", OBConversion::INOPTIONS);
    }

    bool parse(const rapidjson::Document &doc, std::string &errorMsg)
    {
        if (!doc.IsObject() || !doc.HasMember("root") || !doc["root"].IsObject()) {
            errorMsg = "KET document must have a 'root' object";
            return false;
        }

        if (doc.HasMember("ket_version") && doc["ket_version"].IsString()) {
            const std::string version = doc["ket_version"].GetString();
            if (!ketVersionIsSupported(version)) {
                errorMsg = "unsupported KET major version: " + version;
                return false;
            }
            setPairData(mol_, kAttrKetVersion, version);
        }

        const auto &root = doc["root"];

        // Pre-scan root.nodes to build the preserved nodeOrder and split refs
        // from inline meta. The full original order is restored verbatim on
        // write.
        std::vector<std::string> moleculeRefsInOrder;
        std::unordered_map<std::string, ComponentBounds> molBounds;
        std::unordered_set<std::string> consumedTopMembers;
        consumedTopMembers.insert("ket_version");
        consumedTopMembers.insert("root");

        // Templates section: capture refs in order so the writer can preserve
        // them. The referenced top-level objects are recorded as
        // `unknownRefs` (we don't translate monomer / ambiguous /
        // monomerGroup templates into atoms).
        if (root.HasMember("templates") && root["templates"].IsArray()) {
            for (const auto &t : root["templates"].GetArray()) {
                if (!t.IsObject() || !t.HasMember("$ref") || !t["$ref"].IsString())
                    continue;
                const std::string refName = t["$ref"].GetString();
                passthrough_.templateOrder.emplace_back(refName);
                if (doc.HasMember(refName.c_str()) &&
                    doc[refName.c_str()].IsObject()) {
                    passthrough_.unknownRefs.emplace_back(
                        refName, stringifyValue(doc[refName.c_str()]));
                    consumedTopMembers.insert(refName);
                }
            }
        }

        if (root.HasMember("connections"))
            passthrough_.connections = stringifyValue(root["connections"]);
        if (root.HasMember("annotation"))
            passthrough_.documentAnnotation = stringifyValue(root["annotation"]);

        mol_->BeginModify();

        if (root.HasMember("nodes") && root["nodes"].IsArray()) {
            for (const auto &node : root["nodes"].GetArray()) {
                if (!node.IsObject()) continue;

                if (node.HasMember("$ref") && node["$ref"].IsString()) {
                    const std::string refName = node["$ref"].GetString();
                    NodeOrderEntry e;
                    e.kind = NodeOrderEntry::Kind::Ref;
                    e.value = refName;
                    passthrough_.nodeOrder.push_back(std::move(e));

                    if (!doc.HasMember(refName.c_str())) continue;
                    consumedTopMembers.insert(refName);
                    const auto &top = doc[refName.c_str()];
                    if (!top.IsObject() || !top.HasMember("type") ||
                        !top["type"].IsString()) {
                        passthrough_.unknownRefs.emplace_back(
                            refName, stringifyValue(top));
                        continue;
                    }
                    const std::string nodeType = top["type"].GetString();
                    if (nodeType == "molecule") {
                        ComponentBounds cb;
                        const int firstIdx =
                            static_cast<int>(mol_->NumAtoms()) + 1;
                        parseMolecule(top, cb);
                        molBounds.emplace(refName, cb);
                        moleculeRefsInOrder.emplace_back(refName);
                        // Mark this molecule's first atom with its original
                        // $ref so the writer can preserve the name.
                        if (mol_->NumAtoms() >= static_cast<unsigned int>(firstIdx)) {
                            OBAtom *firstAtom = mol_->GetAtom(firstIdx);
                            if (firstAtom)
                                setPairData(firstAtom, kAttrKetMolRef, refName);
                        }
                    } else if (nodeType == "rgroup") {
                        parseRGroup(refName, top);
                    } else {
                        passthrough_.unknownRefs.emplace_back(
                            refName, stringifyValue(top));
                    }
                } else if (node.HasMember("type") && node["type"].IsString()) {
                    NodeOrderEntry e;
                    e.kind = NodeOrderEntry::Kind::Inline;
                    e.value = stringifyValue(node);
                    passthrough_.nodeOrder.push_back(std::move(e));

                    const std::string nodeType = node["type"].GetString();
                    if (nodeType == "arrow")
                        parseArrow(node);
                    else if (nodeType == "plus")
                        parsePlus(node);
                    // For other inline meta (text, image, simpleObject,
                    // multi-tailed-arrow) the raw JSON in nodeOrder is the
                    // only representation we keep.
                }
            }
        }

        // Anything left at the top level that isn't `root`/`ket_version`/an
        // already-consumed ref becomes a passthrough top-member so the
        // writer can echo it back.
        for (auto it = doc.MemberBegin(); it != doc.MemberEnd(); ++it) {
            std::string name = it->name.GetString();
            if (consumedTopMembers.count(name)) continue;
            passthrough_.unknownTopMembers.emplace_back(
                name, stringifyValue(it->value));
        }

        mol_->EndModify();

        // After EndModify, OB has set up the topology but we need to redo
        // a few things that EndModify wipes (selected flags, etc.) — handled
        // via OBPairData attributes so they survive.

        // Aromatic bonds: preserve KET aromaticity rather than forcing a
        // Kekule form chosen by Open Babel.
        if (needsKekulization_) {
            mol_->SetAromaticPerceived();
            FOR_BONDS_OF_MOL(b, mol_) {
                if (b->IsAromatic()) {
                    b->GetBeginAtom()->SetAromatic();
                    b->GetEndAtom()->SetAromatic();
                }
            }
        }

        for (const auto &alias : aliasesToExpand_)
            alias.first->Expand(*mol_, alias.second->GetIdx());

        // Flush accumulated S-group JSON to a single OBPairData (avoiding
        // O(N^2) re-parse+re-serialize that the naive append would cause).
        flushSGroupAccumulator();

        // Reaction assembly: if exactly one arrow was found, assign roles by
        // projecting each component's centroid onto the arrow axis.
        if (!arrows_.empty())
            assembleReaction(molBounds, moleculeRefsInOrder);

        // Dimension detection — KET coords are 2D conventionally (z = 0).
        bool any3D = false;
        FOR_ATOMS_OF_MOL(a, mol_) {
            if (std::fabs(a->z()) > 1e-9) { any3D = true; break; }
        }
        mol_->SetDimension(any3D ? 3 : 2);

        // Perceive stereo from 2D wedges/hashes if present.
        if (!any3D)
            StereoFrom2D(mol_);
        else
            StereoFrom3D(mol_);

        if (!passthrough_.empty())
            setPairData(mol_, kAttrKetPassthrough, passthrough_.serialize());

        return true;
    }

private:
    OBMol *mol_;
    bool expandAliases_ = false;
    bool needsKekulization_ = false;
    std::vector<std::pair<AliasData *, OBAtom *>> aliasesToExpand_;
    std::vector<ArrowGeometry> arrows_;
    std::vector<PlusGeometry> pluses_;
    KetPassthrough passthrough_;
    // Accumulates S-group entries (with remapped atom/bond indices) and is
    // serialized to OBPairData once via flushSGroupAccumulator().
    rapidjson::Document sgroupAccumulator_{rapidjson::kArrayType};

    // ----- Atom / bond parsing -------------------------------------------

    void parseMolecule(const rapidjson::Value &node, ComponentBounds &cb)
    {
        const int atomOffset = static_cast<int>(mol_->NumAtoms());
        const int bondOffset = static_cast<int>(mol_->NumBonds());

        cb.firstAtomIdx = atomOffset + 1;

        const bool hasAtoms =
            node.HasMember("atoms") && node["atoms"].IsArray();

        // Track positions so atomMap[i] yields the 1-based OBAtom index.
        std::vector<int> atomMap;
        if (hasAtoms) atomMap.reserve(node["atoms"].Size());

        double sumX = 0.0, sumY = 0.0;
        int posCount = 0;
        const rapidjson::SizeType atomCount =
            hasAtoms ? node["atoms"].Size() : 0;
        for (rapidjson::SizeType i = 0; i < atomCount; ++i) {
            const auto &a = node["atoms"][i];
            OBAtom *atom = parseAtom(a);
            atomMap.push_back(atom->GetIdx());
            if (a.HasMember("location") && a["location"].IsArray() &&
                a["location"].Size() >= 2) {
                sumX += a["location"][0].GetDouble();
                sumY += a["location"][1].GetDouble();
                ++posCount;
            }
        }
        if (posCount > 0) {
            cb.centroidX = sumX / posCount;
            cb.centroidY = sumY / posCount;
        }
        cb.lastAtomIdx = static_cast<int>(mol_->NumAtoms());

        if (node.HasMember("bonds") && node["bonds"].IsArray()) {
            for (const auto &b : node["bonds"].GetArray())
                parseBond(b, atomMap);
        }

        if (node.HasMember("sgroups") && node["sgroups"].IsArray()) {
            for (const auto &s : node["sgroups"].GetArray())
                parseSGroup(s, atomMap, bondOffset);
        }

        // Per-molecule passthrough fields (properties, highlight, stereo
        // flag position) are stored on the FIRST atom of this molecule —
        // not the OBMol — so they stay tied to their component when the
        // OBMol contains multiple molecules.
        OBAtom *firstAtom = (cb.firstAtomIdx > 0 &&
                             mol_->NumAtoms() >=
                                 static_cast<unsigned int>(cb.firstAtomIdx))
                                ? mol_->GetAtom(cb.firstAtomIdx)
                                : nullptr;

        if (firstAtom && node.HasMember("properties") &&
            node["properties"].IsArray()) {
            setPairData(firstAtom, "_ket_properties",
                        stringifyValue(node["properties"]));
            // Also surface each entry as OBPairData on the OBMol so other
            // OpenBabel APIs can find them.
            for (const auto &p : node["properties"].GetArray()) {
                std::string key, val;
                if (getMemberString(p, "key", key) &&
                    getMemberString(p, "value", val)) {
                    setPairData(mol_, key.c_str(), val);
                }
            }
        }

        if (firstAtom && node.HasMember("stereoFlagPosition") &&
            node["stereoFlagPosition"].IsObject()) {
            setPairData(firstAtom, kAttrKetStereoFlagPos,
                        stringifyValue(node["stereoFlagPosition"]));
        }

        if (firstAtom && node.HasMember("highlight"))
            setPairData(firstAtom, "_ket_highlight",
                        stringifyValue(node["highlight"]));
    }

    OBAtom *parseAtom(const rapidjson::Value &a)
    {
        OBAtom *atom = mol_->NewAtom();

        // Element/label resolution.
        int elem = 0;
        std::string label;
        if (a.HasMember("type") && a["type"].IsString()) {
            const std::string atomType = a["type"].GetString();
            if (atomType == "rg-label") {
                // R-group site: store as a dummy atom and stash the R-group number(s).
                atom->SetAtomicNum(0);
                if (a.HasMember("$refs") && a["$refs"].IsArray() &&
                    a["$refs"].Size() > 0) {
                    std::string ref = a["$refs"][0].GetString();
                    if (ref.rfind("rg-", 0) == 0) {
                        try {
                            int n = std::stoi(ref.substr(3));
                            setPairDataInt(atom, kAttrKetRGroupNumber, n);
                        } catch (...) {}
                    }
                    // Preserve all refs JSON for multi-refs cases.
                    setPairData(atom, "_ket_rgroup_refs",
                                stringifyValue(a["$refs"]));
                }
                if (a.HasMember("attachmentOrder"))
                    setPairData(atom, "_ket_attachmentOrder",
                                stringifyValue(a["attachmentOrder"]));
            } else if (atomType == "atom-list") {
                // Query: pick the first listed element so we have a valid atom,
                // but record the whole list as a passthrough.
                atom->SetAtomicNum(0);
                if (a.HasMember("elements") && a["elements"].IsArray() &&
                    a["elements"].Size() > 0 && a["elements"][0].IsString()) {
                    int el = OBElements::GetAtomicNum(a["elements"][0].GetString());
                    if (el > 0) atom->SetAtomicNum(el);
                }
                setPairData(atom, "_ket_atomList", stringifyValue(a));
            }
        } else if (a.HasMember("label") && a["label"].IsString()) {
            label = a["label"].GetString();
            int iso = 0;
            if (decodeIsotopeShorthand(label, elem, iso)) {
                atom->SetAtomicNum(elem);
                atom->SetIsotope(iso);
            } else if (label == "*" || label == "A" || label == "Q" ||
                       label == "X" || label == "M" || label == "AH" ||
                       label == "QH" || label == "XH" || label == "MH" ||
                       label == "R") {
                atom->SetAtomicNum(0);
                setPairData(atom, "_ket_label", label);
            } else {
                elem = OBElements::GetAtomicNum(label.c_str());
                atom->SetAtomicNum(elem);
                if (elem == 0)
                    setPairData(atom, "_ket_label", label);
            }
        } else {
            atom->SetAtomicNum(0);
        }

        // Position.
        if (a.HasMember("location") && a["location"].IsArray() &&
            a["location"].Size() >= 3) {
            double x = a["location"][0].GetDouble();
            double y = a["location"][1].GetDouble();
            double z = a["location"][2].GetDouble();
            atom->SetVector(x, y, z);
        }

        // Charge, isotope, radical.
        int charge = 0;
        if (getMemberInt(a, "charge", charge))      atom->SetFormalCharge(charge);
        int isotope = 0;
        if (getMemberInt(a, "isotope", isotope) && isotope > 0)
            atom->SetIsotope(isotope);
        int radical = 0;
        if (getMemberInt(a, "radical", radical) && radical > 0)
            atom->SetSpinMultiplicity(static_cast<short>(radical));

        // Implicit hydrogens.
        int implH = -1;
        if (getMemberInt(a, "implicitHCount", implH) && implH >= 0) {
            atom->SetImplicitHCount(static_cast<unsigned int>(implH));
            setPairData(atom, kAttrKetImplicitHCountSet, "true");
        }

        // Explicit valence — KET uses 15 to mean zero.
        int evalence = -1;
        if (getMemberInt(a, "explicitValence", evalence) && evalence >= 0)
            setPairDataInt(atom, kAttrKetExplicitValence, evalence);

        // Reaction-mapping fields.
        int mapping = 0;
        if (getMemberInt(a, "mapping", mapping) && mapping != 0) {
            setPairDataInt(atom, kAttrKetMapping, mapping);
            // OpenBabel atom map class lives in OBPairInteger "Atom Class".
            auto *pi = new OBPairInteger();
            pi->SetAttribute("Atom Class");
            pi->SetValue(mapping);
            pi->SetOrigin(fileformatInput);
            atom->SetData(pi);
        }
        int invRet = 0;
        if (getMemberInt(a, "invRet", invRet) && invRet != 0)
            setPairDataInt(atom, kAttrKetInvRet, invRet);
        bool exactChangeFlag = false;
        if (getMemberBool(a, "exactChangeFlag", exactChangeFlag) && exactChangeFlag)
            setPairData(atom, kAttrKetExactChangeFlag, "true");

        // Stereo / CIP / alias.
        std::string s;
        if (getMemberString(a, "stereoLabel", s))
            setPairData(atom, kAttrKetStereoLabel, s);
        if (getMemberString(a, "cip", s))
            setPairData(atom, kAttrKetCip, s);
        if (getMemberString(a, "alias", s)) {
            setPairData(atom, kAttrKetAlias, s);
            auto *ad = new AliasData();
            ad->SetAlias(s);
            ad->SetOrigin(fileformatInput);
            atom->SetData(ad);
            if (expandAliases_ && atom->GetAtomicNum() == 0)
                aliasesToExpand_.push_back(std::make_pair(ad, atom));
        }

        // Attachment points (bitmask).
        int aps = 0;
        if (getMemberInt(a, "attachmentPoints", aps) && aps != 0)
            setPairDataInt(atom, kAttrKetAttachmentPoints, aps);

        // Query bits — preserve verbatim.
        int rbCount = 0;
        if (getMemberInt(a, "ringBondCount", rbCount))
            setPairDataInt(atom, kAttrKetRingBondCount, rbCount);
        int subCount = 0;
        if (getMemberInt(a, "substitutionCount", subCount))
            setPairDataInt(atom, kAttrKetSubstitutionCount, subCount);
        bool unsat = false;
        if (getMemberBool(a, "unsaturatedAtom", unsat) && unsat)
            setPairData(atom, kAttrKetUnsaturatedAtom, "true");
        int hCount = 0;
        if (getMemberInt(a, "hCount", hCount))
            setPairDataInt(atom, kAttrKetHCount, hCount);
        bool selected = false;
        if (getMemberBool(a, "selected", selected) && selected)
            setPairData(atom, kAttrKetSelected, "true");

        if (a.HasMember("queryProperties"))
            setPairData(atom, kAttrKetQueryProperties,
                        stringifyValue(a["queryProperties"]));

        return atom;
    }

    void parseBond(const rapidjson::Value &b, const std::vector<int> &atomMap)
    {
        if (!b.HasMember("atoms") || !b["atoms"].IsArray() ||
            b["atoms"].Size() < 2)
            return;

        const int aIdx = b["atoms"][0].GetInt();
        const int bIdx = b["atoms"][1].GetInt();
        if (aIdx < 0 || bIdx < 0 ||
            aIdx >= static_cast<int>(atomMap.size()) ||
            bIdx >= static_cast<int>(atomMap.size()))
            return;

        const int beginIdx = atomMap[aIdx];
        const int endIdx   = atomMap[bIdx];

        int ketType = KET_BOND_SINGLE;
        getMemberInt(b, "type", ketType);

        int order = 1;
        int flags = 0;
        switch (ketType) {
        case KET_BOND_SINGLE:               order = 1; break;
        case KET_BOND_DOUBLE:               order = 2; break;
        case KET_BOND_TRIPLE:               order = 3; break;
        case KET_BOND_AROMATIC:
            order = 1;
            flags |= OBBond::Aromatic;
            needsKekulization_ = true;
            break;
        default:
            // Query / coordination / hydrogen bonds: represent as single-order
            // and stash the original KET type for round-trip.
            order = 1;
            break;
        }

        int stereo = 0;
        getMemberInt(b, "stereo", stereo);
        switch (stereo) {
        case KET_STEREO_UP:       flags |= OBBond::Wedge;       break;
        case KET_STEREO_DOWN:     flags |= OBBond::Hash;        break;
        case KET_STEREO_EITHER:   flags |= OBBond::WedgeOrHash; break;
        case KET_STEREO_CISTRANS: flags |= OBBond::CisOrTrans;  break;
        default: break;
        }

        if (!mol_->AddBond(beginIdx, endIdx, order, flags))
            return;

        OBBond *bond = mol_->GetBond(beginIdx, endIdx);
        if (!bond) return;

        if (ketType == KET_BOND_AROMATIC ||
            (ketType != KET_BOND_SINGLE && ketType != KET_BOND_DOUBLE &&
             ketType != KET_BOND_TRIPLE))
            setPairDataInt(bond, kAttrKetBondType, ketType);

        int topology = 0;
        if (getMemberInt(b, "topology", topology) && topology != 0)
            setPairDataInt(bond, kAttrKetBondTopology, topology);
        int center = 0;
        if (getMemberInt(b, "center", center) && center != 0)
            setPairDataInt(bond, kAttrKetBondCenter, center);
        int stereobox = 0;
        if (getMemberInt(b, "stereobox", stereobox) && stereobox != 0)
            setPairDataInt(bond, kAttrKetBondStereobox, stereobox);

        std::string s;
        if (getMemberString(b, "customQuery", s))
            setPairData(bond, kAttrKetBondCustomQuery, s);
        if (getMemberString(b, "cip", s))
            setPairData(bond, kAttrKetCip, s);
        bool sel = false;
        if (getMemberBool(b, "selected", sel) && sel)
            setPairData(bond, kAttrKetSelected, "true");
    }

    void parseSGroup(const rapidjson::Value &s,
                     const std::vector<int> &atomMap,
                     int bondOffset)
    {
        // OpenBabel doesn't have a native S-group container that survives
        // round-trip cleanly, so we serialize each S-group as a JSON blob
        // attached to the OBMol. The writer reads these back verbatim.
        //
        // We accumulate into a single in-memory rapidjson Array and flush it
        // to OBPairData once at end-of-parse — repeatedly parsing and
        // re-serializing the stored JSON for every S-group would be O(N^2)
        // in the number of S-groups.
        if (!s.IsObject()) return;

        auto &al = sgroupAccumulator_.GetAllocator();
        rapidjson::Value copy(s, al);  // deep copy with our allocator

        auto remap = [](rapidjson::Value &arr, const std::vector<int> &mapping,
                        int offset) {
            if (!arr.IsArray()) return;
            for (auto &v : arr.GetArray()) {
                if (!v.IsInt()) continue;
                const int idx = v.GetInt();
                if (idx >= 0 && idx < static_cast<int>(mapping.size()))
                    v.SetInt(mapping[idx] - 1 + offset);
            }
        };

        if (copy.HasMember("atoms")) remap(copy["atoms"], atomMap, 0);
        if (copy.HasMember("bonds")) {
            auto &arr = copy["bonds"];
            if (arr.IsArray()) {
                for (auto &v : arr.GetArray()) {
                    if (v.IsInt()) v.SetInt(v.GetInt() + bondOffset);
                }
            }
        }

        sgroupAccumulator_.PushBack(copy, al);
    }

    void flushSGroupAccumulator()
    {
        if (sgroupAccumulator_.Empty()) return;
        rapidjson::StringBuffer buf;
        rapidjson::Writer<rapidjson::StringBuffer> w(buf);
        sgroupAccumulator_.Accept(w);
        setPairData(mol_, "_ket_sgroup",
                    std::string(buf.GetString(), buf.GetSize()));
    }

    void parseRGroup(const std::string &refName, const rapidjson::Value &node)
    {
        // R-group payloads are full molecule fragments. We store them as
        // passthrough so the writer can echo them — translating them into
        // OpenBabel's mol-of-fragments model would require structural
        // changes outside this format's scope.
        passthrough_.unknownRefs.emplace_back(refName, stringifyValue(node));
    }

    // ----- Meta parsing --------------------------------------------------

    void parseArrow(const rapidjson::Value &node)
    {
        ArrowGeometry ag;
        if (!node.HasMember("data") || !node["data"].IsObject()) return;
        const auto &data = node["data"];
        if (!data.HasMember("pos") || !data["pos"].IsArray() ||
            data["pos"].Size() < 2)
            return;
        auto readVec = [](const rapidjson::Value &v, Vec3d &out) {
            if (!v.IsObject()) return;
            getMemberDouble(v, "x", out.x);
            getMemberDouble(v, "y", out.y);
            getMemberDouble(v, "z", out.z);
        };
        readVec(data["pos"][0], ag.tail);
        readVec(data["pos"][1], ag.head);
        getMemberString(data, "mode", ag.mode);
        ag.valid = true;
        arrows_.push_back(ag);
    }

    void parsePlus(const rapidjson::Value &node)
    {
        PlusGeometry pg;
        if (node.HasMember("location") && node["location"].IsArray() &&
            node["location"].Size() >= 2) {
            pg.pos.x = node["location"][0].GetDouble();
            pg.pos.y = node["location"][1].GetDouble();
            if (node["location"].Size() >= 3)
                pg.pos.z = node["location"][2].GetDouble();
        }
        pluses_.push_back(pg);
    }

    // ----- Reaction assembly ---------------------------------------------

    // Compute per-component centroids for all connected fragments of the
    // assembled molecule, then assign reaction roles based on each centroid's
    // projection onto the arrow axis. Anything close to the arrow's
    // perpendicular axis is tagged as an AGENT.
    void assembleReaction(const std::unordered_map<std::string, ComponentBounds> &molBounds,
                          const std::vector<std::string> &orderedRefs)
    {
        if (arrows_.empty()) return;
        // Use the first arrow for left/right partitioning. (Multi-tailed and
        // pathway reactions are out of scope; they round-trip as inline meta.)
        const ArrowGeometry &arrow = arrows_.front();
        const double midX = 0.5 * (arrow.tail.x + arrow.head.x);
        const double midY = 0.5 * (arrow.tail.y + arrow.head.y);
        const double dx = arrow.head.x - arrow.tail.x;
        const double dy = arrow.head.y - arrow.tail.y;
        const double axisLen = std::sqrt(dx * dx + dy * dy);
        if (axisLen < 1e-9) return;
        const double ux = dx / axisLen;
        const double uy = dy / axisLen;
        // Perpendicular distance threshold — molecules near the arrow axis
        // (within half the arrow length perpendicular distance) are agents.
        const double agentBand = std::max(0.5 * axisLen, 1.0);

        mol_->SetIsReaction();
        OBReactionFacade facade(mol_);
        unsigned int compId = 1;

        // Process components in the order they appear in root.nodes for stable
        // numbering, since OBReactionFacade preserves per-atom data.
        for (const auto &refName : orderedRefs) {
            auto it = molBounds.find(refName);
            if (it == molBounds.end()) continue;
            const auto &cb = it->second;
            const double cx = cb.centroidX;
            const double cy = cb.centroidY;

            // Vector from arrow midpoint to centroid.
            const double rx = cx - midX;
            const double ry = cy - midY;
            // Project onto arrow direction; sign indicates side.
            const double projected = rx * ux + ry * uy;
            // Perpendicular distance from arrow axis.
            const double perp = std::fabs(rx * (-uy) + ry * ux);

            OBReactionRole role;
            if (perp > agentBand)
                role = AGENT;
            else if (projected < 0.0)
                role = REACTANT;
            else
                role = PRODUCT;

            for (int ai = cb.firstAtomIdx; ai <= cb.lastAtomIdx; ++ai) {
                OBAtom *atom = mol_->GetAtom(ai);
                if (!atom) continue;
                facade.SetRole(atom, role);
                facade.SetComponentId(atom, compId);
            }
            ++compId;
        }
    }
};

// ============================================================================
//  Writer helpers
// ============================================================================

struct OutComponent {
    std::string refName;
    std::vector<int> atomIdx;
    double centroidX = 0.0;
    double centroidY = 0.0;
    OBReactionRole role = NO_REACTIONROLE;
};

void writeAtomEntry(rapidjson::Value &arr,
                    OBAtom *atom,
                    rapidjson::Document::AllocatorType &al)
{
    rapidjson::Value obj(rapidjson::kObjectType);
    const int anum = atom->GetAtomicNum();
    const unsigned int iso = atom->GetIsotope();
    std::string label;

    // Check for explicit KET label passthrough (covers query atoms / pseudo).
    if (auto *labelData = atom->GetData("_ket_label")) {
        if (auto *pd = dynamic_cast<OBPairData *>(labelData))
            label = pd->GetValue();
    }
    if (label.empty()) {
        if (anum == OBElements::Hydrogen && iso == 2) label = "D";
        else if (anum == OBElements::Hydrogen && iso == 3) label = "T";
        else if (anum > 0) label = OBElements::GetSymbol(anum);
        else label = "*";  // dummy / R-site placeholder
    }

    // Query atom-list: replay the originally-captured type/elements/notList
    // from `_ket_atomList`; we still emit location/charge/etc through the
    // regular path below.
    std::string atomListRaw;
    const bool isAtomList =
        getPairDataString(atom, "_ket_atomList", atomListRaw);
    if (isAtomList) {
        rapidjson::Document tmp;
        if (!tmp.Parse(atomListRaw.c_str()).HasParseError() && tmp.IsObject()) {
            obj.AddMember("type", "atom-list", al);
            if (tmp.HasMember("elements") && tmp["elements"].IsArray()) {
                rapidjson::Value v;
                v.CopyFrom(tmp["elements"], al);
                obj.AddMember("elements", v, al);
            }
            if (tmp.HasMember("notList") && tmp["notList"].IsBool() &&
                tmp["notList"].GetBool())
                obj.AddMember("notList", true, al);
        }
    }
    // R-site: emit as rg-label. Prefer the preserved multi-ref JSON when
    // available; otherwise fall back to a single rg-N ref.
    std::string rgRefsRaw;
    const bool hasRgRefs =
        !isAtomList && getPairDataString(atom, "_ket_rgroup_refs", rgRefsRaw);
    int rgNum = 0;
    const bool hasSingleRgNum =
        !isAtomList && getPairDataInt(atom, kAttrKetRGroupNumber, rgNum) &&
        rgNum > 0;

    if (hasRgRefs) {
        rapidjson::Document tmp;
        if (!tmp.Parse(rgRefsRaw.c_str()).HasParseError() && tmp.IsArray()) {
            obj.AddMember("type", "rg-label", al);
            rapidjson::Value v;
            v.CopyFrom(tmp, al);
            obj.AddMember("$refs", v, al);
        }
    } else if (hasSingleRgNum) {
        obj.AddMember("type", "rg-label", al);
        rapidjson::Value refs(rapidjson::kArrayType);
        refs.PushBack(rapidjson::Value(
                          ("rg-" + std::to_string(rgNum)).c_str(), al), al);
        obj.AddMember("$refs", refs, al);
    } else if (!isAtomList) {
        obj.AddMember("label", rapidjson::Value(label.c_str(), al), al);
    }

    // attachmentOrder is only meaningful for rg-label atoms.
    if (hasRgRefs || hasSingleRgNum) {
        std::string attOrder;
        if (getPairDataString(atom, "_ket_attachmentOrder", attOrder)) {
            rapidjson::Document tmp;
            if (!tmp.Parse(attOrder.c_str()).HasParseError() && tmp.IsArray()) {
                rapidjson::Value v;
                v.CopyFrom(tmp, al);
                obj.AddMember("attachmentOrder", v, al);
            }
        }
    }

    // Coordinates (always 3 components).
    rapidjson::Value loc(rapidjson::kArrayType);
    loc.PushBack(atom->x(), al);
    loc.PushBack(atom->y(), al);
    loc.PushBack(atom->z(), al);
    obj.AddMember("location", loc, al);

    const int charge = atom->GetFormalCharge();
    if (charge != 0) obj.AddMember("charge", charge, al);

    if (iso > 0 && !(anum == OBElements::Hydrogen && (iso == 2 || iso == 3)))
        obj.AddMember("isotope", static_cast<int>(iso), al);

    const int radical = atom->GetSpinMultiplicity();
    if (radical > 0) obj.AddMember("radical", radical, al);

    // Pass-through KET attributes.
    int v = 0;
    if (getPairDataInt(atom, kAttrKetMapping, v) && v != 0)
        obj.AddMember("mapping", v, al);
    else if (auto *pi = dynamic_cast<OBPairInteger *>(atom->GetData("Atom Class"))) {
        obj.AddMember("mapping", pi->GetGenericValue(), al);
    }
    if (getPairDataInt(atom, kAttrKetInvRet, v) && v != 0)
        obj.AddMember("invRet", v, al);
    std::string ecf;
    if (getPairDataString(atom, kAttrKetExactChangeFlag, ecf) && ecf == "true")
        obj.AddMember("exactChangeFlag", true, al);

    std::string s;
    if (getPairDataString(atom, kAttrKetStereoLabel, s))
        obj.AddMember("stereoLabel", rapidjson::Value(s.c_str(), al), al);
    if (getPairDataString(atom, kAttrKetCip, s))
        obj.AddMember("cip", rapidjson::Value(s.c_str(), al), al);
    if (getPairDataString(atom, kAttrKetAlias, s))
        obj.AddMember("alias", rapidjson::Value(s.c_str(), al), al);

    if (getPairDataInt(atom, kAttrKetAttachmentPoints, v) && v != 0)
        obj.AddMember("attachmentPoints", v, al);
    if (getPairDataInt(atom, kAttrKetExplicitValence, v) && v >= 0)
        obj.AddMember("explicitValence", v, al);

    if (getPairDataInt(atom, kAttrKetRingBondCount, v))
        obj.AddMember("ringBondCount", v, al);
    if (getPairDataInt(atom, kAttrKetSubstitutionCount, v))
        obj.AddMember("substitutionCount", v, al);
    if (getPairDataString(atom, kAttrKetUnsaturatedAtom, s) && s == "true")
        obj.AddMember("unsaturatedAtom", true, al);
    if (getPairDataInt(atom, kAttrKetHCount, v))
        obj.AddMember("hCount", v, al);

    std::string explicitImplicitH;
    if (getPairDataString(atom, kAttrKetImplicitHCountSet, explicitImplicitH) ||
        atom->GetImplicitHCount() > 0) {
        obj.AddMember("implicitHCount",
                      static_cast<int>(atom->GetImplicitHCount()), al);
    }

    if (getPairDataString(atom, kAttrKetSelected, s) && s == "true")
        obj.AddMember("selected", true, al);

    // queryProperties is preserved as a raw JSON sub-document.
    std::string qp;
    if (getPairDataString(atom, kAttrKetQueryProperties, qp)) {
        rapidjson::Document tmp;
        if (!tmp.Parse(qp.c_str()).HasParseError()) {
            rapidjson::Value v2;
            v2.CopyFrom(tmp, al);
            obj.AddMember("queryProperties", v2, al);
        }
    }

    arr.PushBack(obj, al);
}

void writeBondEntry(rapidjson::Value &arr,
                    OBBond *bond,
                    int beginLocal, int endLocal,
                    rapidjson::Document::AllocatorType &al)
{
    rapidjson::Value obj(rapidjson::kObjectType);

    int ketType = 0;
    if (getPairDataInt(bond, kAttrKetBondType, ketType) && ketType > 0) {
        obj.AddMember("type", ketType, al);
    } else {
        // Map OB order + aromatic flag back to KET type.
        if (bond->IsAromatic())
            obj.AddMember("type", KET_BOND_AROMATIC, al);
        else
            obj.AddMember("type", static_cast<int>(bond->GetBondOrder()), al);
    }

    rapidjson::Value atoms(rapidjson::kArrayType);
    atoms.PushBack(beginLocal, al);
    atoms.PushBack(endLocal, al);
    obj.AddMember("atoms", atoms, al);

    // Stereo flag.
    int stereo = 0;
    if (bond->IsWedge())       stereo = KET_STEREO_UP;
    else if (bond->IsHash())   stereo = KET_STEREO_DOWN;
    else if (bond->IsWedgeOrHash()) stereo = KET_STEREO_EITHER;
    else if (bond->IsCisOrTrans()) stereo = KET_STEREO_CISTRANS;
    if (stereo != 0)
        obj.AddMember("stereo", stereo, al);

    int v = 0;
    if (getPairDataInt(bond, kAttrKetBondTopology, v) && v != 0)
        obj.AddMember("topology", v, al);
    if (getPairDataInt(bond, kAttrKetBondCenter, v) && v != 0)
        obj.AddMember("center", v, al);
    if (getPairDataInt(bond, kAttrKetBondStereobox, v) && v != 0)
        obj.AddMember("stereobox", v, al);

    std::string s;
    if (getPairDataString(bond, kAttrKetBondCustomQuery, s))
        obj.AddMember("customQuery", rapidjson::Value(s.c_str(), al), al);
    if (getPairDataString(bond, kAttrKetCip, s))
        obj.AddMember("cip", rapidjson::Value(s.c_str(), al), al);
    if (getPairDataString(bond, kAttrKetSelected, s) && s == "true")
        obj.AddMember("selected", true, al);

    arr.PushBack(obj, al);
}

void writeMoleculeBody(rapidjson::Document &doc,
                       OBMol *mol,
                       const OutComponent &c,
                       rapidjson::Document::AllocatorType &al)
{
    rapidjson::Value molObj(rapidjson::kObjectType);
    molObj.AddMember("type", "molecule", al);

    // Build a local-index map: source OBMol idx (1-based) -> local index 0..N.
    std::unordered_map<int, int> globalToLocal;
    globalToLocal.reserve(c.atomIdx.size());
    for (size_t i = 0; i < c.atomIdx.size(); ++i)
        globalToLocal.emplace(c.atomIdx[i], static_cast<int>(i));

    rapidjson::Value atomsArr(rapidjson::kArrayType);
    for (int gidx : c.atomIdx) {
        OBAtom *atom = mol->GetAtom(gidx);
        if (!atom) continue;
        writeAtomEntry(atomsArr, atom, al);
    }
    molObj.AddMember("atoms", atomsArr, al);

    // Emit bonds and build a global→local bond-index map so we can rewrite
    // S-group bond references back to per-molecule-local KET indices.
    rapidjson::Value bondsArr(rapidjson::kArrayType);
    std::unordered_map<int, int> bondGlobalToLocal;
    FOR_BONDS_OF_MOL(b, mol) {
        const int beginG = b->GetBeginAtomIdx();
        const int endG   = b->GetEndAtomIdx();
        const auto itB = globalToLocal.find(beginG);
        const auto itE = globalToLocal.find(endG);
        if (itB == globalToLocal.end() || itE == globalToLocal.end())
            continue;
        const int localIdx = static_cast<int>(bondsArr.Size());
        bondGlobalToLocal.emplace(static_cast<int>(b->GetIdx()), localIdx);
        writeBondEntry(bondsArr, &*b, itB->second, itE->second, al);
    }
    if (bondsArr.Size() > 0)
        molObj.AddMember("bonds", bondsArr, al);

    // Recover S-groups from passthrough storage on the source molecule.
    // The reader stored atom indices as 0-based global, and bond indices as
    // 0-based global; we remap both back to local indices for this component.
    std::string sgRaw;
    if (getPairDataString(mol, "_ket_sgroup", sgRaw)) {
        rapidjson::Document sgDoc;
        if (!sgDoc.Parse(sgRaw.c_str()).HasParseError() && sgDoc.IsArray()) {
            rapidjson::Value sgArr(rapidjson::kArrayType);
            for (auto &item : sgDoc.GetArray()) {
                if (!item.IsObject() || !item.HasMember("atoms") ||
                    !item["atoms"].IsArray())
                    continue;
                // Only include S-groups whose atoms are all in this component.
                bool allHere = true;
                for (auto &v : item["atoms"].GetArray()) {
                    if (!v.IsInt()) { allHere = false; break; }
                    if (!globalToLocal.count(v.GetInt() + 1)) {
                        allHere = false; break;
                    }
                }
                if (!allHere) continue;

                rapidjson::Value entry;
                entry.CopyFrom(item, al);

                if (entry.HasMember("atoms") && entry["atoms"].IsArray()) {
                    for (auto &v : entry["atoms"].GetArray()) {
                        if (!v.IsInt()) continue;
                        const int gidx = v.GetInt() + 1;  // 1-based global
                        v.SetInt(globalToLocal[gidx]);    // 0-based local
                    }
                }
                if (entry.HasMember("bonds") && entry["bonds"].IsArray()) {
                    rapidjson::Value remapped(rapidjson::kArrayType);
                    for (auto &v : entry["bonds"].GetArray()) {
                        if (!v.IsInt()) continue;
                        auto it = bondGlobalToLocal.find(v.GetInt());
                        if (it != bondGlobalToLocal.end())
                            remapped.PushBack(it->second, al);
                    }
                    entry["bonds"] = remapped;
                }
                sgArr.PushBack(entry, al);
            }
            if (sgArr.Size() > 0)
                molObj.AddMember("sgroups", sgArr, al);
        }
    }

    // Per-molecule passthrough (stereoFlagPosition / highlight / properties)
    // lives on the FIRST atom of this component, not the OBMol — otherwise a
    // multi-molecule OBMol would echo the same blocks into every component.
    OBAtom *firstAtom =
        c.atomIdx.empty() ? nullptr : mol->GetAtom(c.atomIdx.front());

    if (firstAtom) {
        std::string sfp;
        if (getPairDataString(firstAtom, kAttrKetStereoFlagPos, sfp)) {
            rapidjson::Document tmp;
            if (!tmp.Parse(sfp.c_str()).HasParseError()) {
                rapidjson::Value v;
                v.CopyFrom(tmp, al);
                molObj.AddMember("stereoFlagPosition", v, al);
            }
        }
        std::string highlightRaw;
        if (getPairDataString(firstAtom, "_ket_highlight", highlightRaw)) {
            rapidjson::Document tmp;
            if (!tmp.Parse(highlightRaw.c_str()).HasParseError()) {
                rapidjson::Value v;
                v.CopyFrom(tmp, al);
                molObj.AddMember("highlight", v, al);
            }
        }
        std::string propsRaw;
        if (getPairDataString(firstAtom, "_ket_properties", propsRaw)) {
            rapidjson::Document tmp;
            if (!tmp.Parse(propsRaw.c_str()).HasParseError()) {
                rapidjson::Value v;
                v.CopyFrom(tmp, al);
                molObj.AddMember("properties", v, al);
            }
        }
    }

    doc.AddMember(rapidjson::Value(c.refName.c_str(), al), molObj, al);
}

}  // anonymous namespace

// ============================================================================
//  Format class
// ============================================================================

class KETFormat : public OBMoleculeFormat {
public:
    KETFormat()
    {
        OBConversion::RegisterFormat("ket", this, "chemical/x-indigo-ket");
        OBConversion::RegisterOptionParam("a", this, 0, OBConversion::INOPTIONS);
    }

    const char *Description() override
    {
        return "Ketcher KET JSON format\n"
               "Native format used by Ketcher and Indigo for chemical structures\n"
               "and reactions.\n\n"
               "This implementation reads and writes molecules (atoms, bonds,\n"
               "S-groups, R-group attachment sites) and full reactions assembled\n"
               "from arrow + plus + molecule geometry. Unknown sections such as\n"
               "macromolecule monomers, monomer templates, V2 rich text, and\n"
               "images are preserved verbatim across round-trips.\n\n"
               "Write Options, e.g. -xm\n"
               "  m  minified output (no indentation)\n\n"
               "Read Options, e.g. -ia\n"
               "  a  expand KET atom aliases where chemically meaningful\n\n";
    }

    const char *SpecificationURL() override
    {
        return "https://github.com/epam/Indigo";
    }

    const char *GetMIMEType() override
    {
        return "chemical/x-indigo-ket";
    }

    unsigned int Flags() override { return READONEONLY | WRITEONEONLY; }

    bool ReadMolecule(OBBase *pOb, OBConversion *pConv) override;
    bool WriteMolecule(OBBase *pOb, OBConversion *pConv) override;
};

KETFormat theKETFormat;

bool KETFormat::ReadMolecule(OBBase *pOb, OBConversion *pConv)
{
    OBMol *pmol = pOb->CastAndClear<OBMol>();
    if (pmol == nullptr) return false;

    std::istream &ifs = *pConv->GetInStream();
    if (!ifs.good()) return false;

    // KET documents are read once per stream. After a successful parse the
    // stream is fully consumed so the next call sees peek() == EOF.
    if (ifs.peek() == EOF) return false;

    rapidjson::Document doc;
    rapidjson::IStreamWrapper isw(ifs);
    doc.ParseStream(isw);
    if (doc.HasParseError()) {
        std::stringstream msg;
        msg << "KET JSON parse error at offset " << doc.GetErrorOffset()
            << ": " << rapidjson::GetParseError_En(doc.GetParseError());
        obErrorLog.ThrowError("KETFormat", msg.str(), obError);
        return false;
    }

    KetReader reader(pmol, pConv);
    std::string err;
    if (!reader.parse(doc, err)) {
        obErrorLog.ThrowError("KETFormat", err, obError);
        return false;
    }
    return true;
}

namespace {

// Decide each connected component's $ref name, preferring the name preserved
// from the original document (stored as `_ket_molRef` on the first atom).
// Falls back to "molN" for fresh components or those whose preserved name
// would collide with one already in use.
void assignComponentRefNames(OBMol *pmol, std::vector<OutComponent> &components)
{
    std::unordered_set<std::string> used;

    // First pass: claim preserved names where they don't collide.
    for (auto &c : components) {
        if (c.atomIdx.empty()) continue;
        std::string preserved;
        OBAtom *first = pmol->GetAtom(c.atomIdx.front());
        if (first && getPairDataString(first, kAttrKetMolRef, preserved) &&
            !preserved.empty() && !used.count(preserved)) {
            c.refName = preserved;
            used.insert(preserved);
        }
    }

    // Second pass: anything still un-named gets the next available molN.
    int counter = 0;
    auto nextFree = [&]() {
        for (;;) {
            std::string candidate = "mol" + std::to_string(counter++);
            if (!used.count(candidate)) {
                used.insert(candidate);
                return candidate;
            }
        }
    };
    for (auto &c : components) {
        if (c.refName.empty() && !c.atomIdx.empty())
            c.refName = nextFree();
    }
}

// Detect whether the preserved nodeOrder already contains any object with the
// given inline type (e.g. "arrow", "plus"), so we don't synthesize duplicates.
bool nodeOrderHasInlineType(const KetPassthrough &pt, const char *type)
{
    for (const auto &e : pt.nodeOrder) {
        if (e.kind != NodeOrderEntry::Kind::Inline) continue;
        rapidjson::Document tmp;
        if (tmp.Parse(e.value.c_str()).HasParseError()) continue;
        if (tmp.IsObject() && tmp.HasMember("type") &&
            tmp["type"].IsString() &&
            std::string(tmp["type"].GetString()) == type)
            return true;
    }
    return false;
}

// True if the version string's major component is > 1.
bool majorVersionExceedsOne(const std::string &version)
{
    int major = 0;
    return ketVersionMajor(version, major) && major > 1;
}

}  // anonymous namespace

bool KETFormat::WriteMolecule(OBBase *pOb, OBConversion *pConv)
{
    OBMol *pmol = dynamic_cast<OBMol *>(pOb);
    if (pmol == nullptr) return false;
    std::ostream &ofs = *pConv->GetOutStream();

    rapidjson::Document doc(rapidjson::kObjectType);
    rapidjson::Document::AllocatorType &al = doc.GetAllocator();

    // ----- Step 1: collect components ------------------------------------
    std::vector<OutComponent> components;
    const bool isReaction = pmol->IsReaction();

    if (isReaction) {
        OBReactionFacade facade(pmol);
        auto collectRole = [&](OBReactionRole role) {
            // Group atoms by component id; OBReactionFacade::GetComponent
            // enumerates compIds for a given role, but we need the original
            // (1-based) atom indices in pmol so we can later pull bonds.
            std::map<unsigned int, std::vector<int>> byCompId;
            FOR_ATOMS_OF_MOL(a, pmol) {
                if (facade.GetRole(&*a) != role) continue;
                const unsigned int cid = facade.GetComponentId(&*a);
                byCompId[cid].push_back(a->GetIdx());
            }
            for (auto &kv : byCompId) {
                OutComponent c;
                c.role = role;
                c.atomIdx = std::move(kv.second);
                double sumX = 0.0, sumY = 0.0;
                for (int idx : c.atomIdx) {
                    OBAtom *a = pmol->GetAtom(idx);
                    sumX += a->x();
                    sumY += a->y();
                }
                if (!c.atomIdx.empty()) {
                    c.centroidX = sumX / c.atomIdx.size();
                    c.centroidY = sumY / c.atomIdx.size();
                }
                components.push_back(std::move(c));
            }
        };
        collectRole(REACTANT);
        collectRole(AGENT);
        collectRole(PRODUCT);
    } else {
        const int natoms = static_cast<int>(pmol->NumAtoms());
        std::vector<int> comp(natoms + 1, 0);
        int next = 0;
        for (int seed = 1; seed <= natoms; ++seed) {
            if (comp[seed] != 0) continue;
            ++next;
            std::vector<int> stack{seed};
            comp[seed] = next;
            std::vector<int> members;
            while (!stack.empty()) {
                const int cur = stack.back();
                stack.pop_back();
                members.push_back(cur);
                OBAtom *a = pmol->GetAtom(cur);
                FOR_NBORS_OF_ATOM(nbr, a) {
                    const int nb = nbr->GetIdx();
                    if (comp[nb] == 0) {
                        comp[nb] = next;
                        stack.push_back(nb);
                    }
                }
            }
            std::sort(members.begin(), members.end());
            OutComponent c;
            c.atomIdx = std::move(members);
            double sumX = 0.0, sumY = 0.0;
            for (int idx : c.atomIdx) {
                OBAtom *a = pmol->GetAtom(idx);
                sumX += a->x();
                sumY += a->y();
            }
            if (!c.atomIdx.empty()) {
                c.centroidX = sumX / c.atomIdx.size();
                c.centroidY = sumY / c.atomIdx.size();
            }
            components.push_back(std::move(c));
        }
        // If pmol has no atoms we leave components empty so root.nodes is
        // empty as well — the document is still structurally valid and we
        // never emit a $ref that can't be resolved to a top-level molecule.
    }

    assignComponentRefNames(pmol, components);

    // ----- Step 2: load passthrough -------------------------------------
    KetPassthrough pt;
    {
        std::string ptRaw;
        if (getPairDataString(pmol, kAttrKetPassthrough, ptRaw))
            pt = KetPassthrough::deserialize(ptRaw);
    }

    // refName → index into components, used to splice fresh molecule refs at
    // their original positions. Empty components (no atoms) are excluded so
    // we never emit a $ref that resolves to no top-level body.
    std::unordered_map<std::string, size_t> componentByRef;
    for (size_t i = 0; i < components.size(); ++i) {
        if (components[i].refName.empty() || components[i].atomIdx.empty())
            continue;
        componentByRef.emplace(components[i].refName, i);
    }

    // ----- Step 3: build root + nodes ------------------------------------
    std::string version;
    getPairDataString(pmol, kAttrKetVersion, version);
    if (majorVersionExceedsOne(version))
        doc.AddMember("ket_version",
                      rapidjson::Value(version.c_str(), al), al);

    rapidjson::Value root(rapidjson::kObjectType);
    rapidjson::Value nodes(rapidjson::kArrayType);

    std::unordered_set<std::string> emittedRefs;

    // Walk the original nodes order, preserving layout and ordering.
    for (const auto &entry : pt.nodeOrder) {
        if (entry.kind == NodeOrderEntry::Kind::Inline) {
            rapidjson::Document tmp;
            if (tmp.Parse(entry.value.c_str()).HasParseError()) continue;
            rapidjson::Value v;
            v.CopyFrom(tmp, al);
            nodes.PushBack(v, al);
        } else {
            // Ref node — emit `{"$ref": refName}` either way. The actual
            // body comes from either our current components (preferred) or
            // from passthrough.unknownRefs.
            rapidjson::Value refObj(rapidjson::kObjectType);
            refObj.AddMember("$ref",
                             rapidjson::Value(entry.value.c_str(), al), al);
            nodes.PushBack(refObj, al);
            emittedRefs.insert(entry.value);
        }
    }

    // Append components added programmatically (or originally first-time
    // written) that weren't part of the preserved nodeOrder.
    for (const auto &c : components) {
        if (c.refName.empty() || c.atomIdx.empty()) continue;
        if (emittedRefs.count(c.refName)) continue;
        rapidjson::Value refObj(rapidjson::kObjectType);
        refObj.AddMember("$ref",
                         rapidjson::Value(c.refName.c_str(), al), al);
        nodes.PushBack(refObj, al);
        emittedRefs.insert(c.refName);
    }

    // Synthesize reaction meta only if the molecule is a reaction AND we
    // don't already have arrows/pluses from passthrough.
    if (isReaction) {
        std::vector<const OutComponent *> reactants, products;
        for (const auto &c : components) {
            if (c.role == REACTANT) reactants.push_back(&c);
            else if (c.role == PRODUCT) products.push_back(&c);
        }
        if (!reactants.empty() && !products.empty()) {
            auto byX = [](const OutComponent *a, const OutComponent *b) {
                return a->centroidX < b->centroidX;
            };
            std::sort(reactants.begin(), reactants.end(), byX);
            std::sort(products.begin(), products.end(), byX);

            const bool haveArrow = nodeOrderHasInlineType(pt, "arrow") ||
                                   nodeOrderHasInlineType(pt, "multi-tailed-arrow");
            const bool havePlus  = nodeOrderHasInlineType(pt, "plus");

            const double tailX = reactants.back()->centroidX;
            const double headX = products.front()->centroidX;
            double sumYr = 0.0, sumYp = 0.0;
            for (auto *p : reactants) sumYr += p->centroidY;
            for (auto *p : products)  sumYp += p->centroidY;
            const double midY = 0.5 * (sumYr / reactants.size() +
                                       sumYp / products.size());
            const double gap = std::max(0.5, 0.1 * std::fabs(headX - tailX));

            if (!haveArrow) {
                rapidjson::Value arrow(rapidjson::kObjectType);
                arrow.AddMember("type", "arrow", al);
                rapidjson::Value data(rapidjson::kObjectType);
                data.AddMember("mode", "open-angle", al);
                rapidjson::Value pos(rapidjson::kArrayType);
                for (double x : {tailX + gap, headX - gap}) {
                    rapidjson::Value p(rapidjson::kObjectType);
                    p.AddMember("x", x, al);
                    p.AddMember("y", midY, al);
                    p.AddMember("z", 0, al);
                    pos.PushBack(p, al);
                }
                data.AddMember("pos", pos, al);
                arrow.AddMember("data", data, al);
                nodes.PushBack(arrow, al);
            }

            if (!havePlus) {
                auto emitPluses =
                    [&](const std::vector<const OutComponent *> &list) {
                        for (size_t i = 0; i + 1 < list.size(); ++i) {
                            const double px = 0.5 * (list[i]->centroidX +
                                                     list[i + 1]->centroidX);
                            const double py = 0.5 * (list[i]->centroidY +
                                                     list[i + 1]->centroidY);
                            rapidjson::Value plus(rapidjson::kObjectType);
                            plus.AddMember("type", "plus", al);
                            rapidjson::Value loc(rapidjson::kArrayType);
                            loc.PushBack(px, al);
                            loc.PushBack(py, al);
                            loc.PushBack(0, al);
                            plus.AddMember("location", loc, al);
                            plus.AddMember(
                                "prop",
                                rapidjson::Value(rapidjson::kObjectType), al);
                            nodes.PushBack(plus, al);
                        }
                    };
                emitPluses(reactants);
                emitPluses(products);
            }
        }
    }

    root.AddMember("nodes", nodes, al);

    // ----- Step 4: connections, templates, annotation --------------------
    if (!pt.connections.empty()) {
        rapidjson::Document tmp;
        if (!tmp.Parse(pt.connections.c_str()).HasParseError()) {
            rapidjson::Value v;
            v.CopyFrom(tmp, al);
            root.AddMember("connections", v, al);
        }
    }

    rapidjson::Value templates(rapidjson::kArrayType);
    for (const auto &name : pt.templateOrder) {
        rapidjson::Value refObj(rapidjson::kObjectType);
        refObj.AddMember("$ref", rapidjson::Value(name.c_str(), al), al);
        templates.PushBack(refObj, al);
    }
    if (templates.Size() > 0)
        root.AddMember("templates", templates, al);

    if (!pt.documentAnnotation.empty()) {
        rapidjson::Document tmp;
        if (!tmp.Parse(pt.documentAnnotation.c_str()).HasParseError()) {
            rapidjson::Value v;
            v.CopyFrom(tmp, al);
            root.AddMember("annotation", v, al);
        }
    }

    doc.AddMember("root", root, al);

    // ----- Step 5: emit molecule bodies ----------------------------------
    for (const auto &c : components) {
        if (c.atomIdx.empty()) continue;
        writeMoleculeBody(doc, pmol, c, al);
    }

    // ----- Step 6: emit preserved top-level objects ---------------------
    // unknownRefs cover monomers, monomer templates, rgroups stored verbatim,
    // etc. We dedupe against doc members already emitted (molecule bodies).
    std::unordered_set<std::string> docMemberNames;
    for (auto it = doc.MemberBegin(); it != doc.MemberEnd(); ++it)
        docMemberNames.insert(it->name.GetString());

    for (const auto &kv : pt.unknownRefs) {
        if (docMemberNames.count(kv.first)) continue;
        rapidjson::Document tmp;
        if (tmp.Parse(kv.second.c_str()).HasParseError()) continue;
        rapidjson::Value v;
        v.CopyFrom(tmp, al);
        doc.AddMember(rapidjson::Value(kv.first.c_str(), al), v, al);
        docMemberNames.insert(kv.first);
    }
    for (const auto &kv : pt.unknownTopMembers) {
        if (docMemberNames.count(kv.first)) continue;
        rapidjson::Document tmp;
        if (tmp.Parse(kv.second.c_str()).HasParseError()) continue;
        rapidjson::Value v;
        v.CopyFrom(tmp, al);
        doc.AddMember(rapidjson::Value(kv.first.c_str(), al), v, al);
        docMemberNames.insert(kv.first);
    }

    // ----- Step 7: serialize --------------------------------------------
    rapidjson::OStreamWrapper osw(ofs);
    if (pConv->IsOption("m", pConv->OUTOPTIONS)) {
        rapidjson::Writer<rapidjson::OStreamWrapper> writer(osw);
        doc.Accept(writer);
    } else {
        rapidjson::PrettyWriter<rapidjson::OStreamWrapper> writer(osw);
        writer.SetIndent(' ', 2);
        doc.Accept(writer);
    }
    return true;
}

}  // namespace OpenBabel
