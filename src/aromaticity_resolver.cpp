/*
 * AromaticityResolver - Implementation
 *
 * Port of OpenChemLib's RingCollection aromaticity detection to C++.
 */

#include <openbabel/aromaticity_resolver.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <algorithm>
#include <iostream>

static constexpr bool debug = false;

namespace OpenChemLib {

// ============================================================================
// Constructor & Destructor
// ============================================================================

AromaticityResolver::AromaticityResolver(OpenBabel::OBMol& mol)
    : mMol(mol) {

    // Initialize feature arrays
    mAtomRingFeatures.resize(mol.NumAtoms(), 0);
    mBondRingFeatures.resize(mol.NumBonds(), 0);

    // Find rings using OpenBabel's SSSR
    findSmallRings();

    // Build annelation map (detect shared bonds between rings)
    buildAnnelationMap();

    // Determine aromaticity for all rings
    determineAromaticity();

    // Propagate features to atoms and bonds
    propagateRingFeatures();
}

AromaticityResolver::~AromaticityResolver() {
    // Nothing to clean up (using std::vector which manages its own memory)
}

// ============================================================================
// Public Query Methods - Bonds
// ============================================================================

bool AromaticityResolver::isAromaticBond(const OpenBabel::OBBond* bond) const {
    int bondIdx = toBondIndex(bond);
    if (bondIdx < 0 || bondIdx >= (int)mBondRingFeatures.size()) {
        return false;
    }
    return (mBondRingFeatures[bondIdx] & FEATURES_AROMATIC) != 0;
}

bool AromaticityResolver::isDelocalizedBond(const OpenBabel::OBBond* bond) const {
    int bondIdx = toBondIndex(bond);
    if (bondIdx < 0 || bondIdx >= (int)mBondRingFeatures.size()) {
        return false;
    }
    return (mBondRingFeatures[bondIdx] & FEATURES_DELOCALIZED) != 0;
}

bool AromaticityResolver::isHeteroAromaticBond(const OpenBabel::OBBond* bond) const {
    int bondIdx = toBondIndex(bond);
    if (bondIdx < 0 || bondIdx >= (int)mBondRingFeatures.size()) {
        return false;
    }
    return (mBondRingFeatures[bondIdx] & FEATURES_HETERO_AROMATIC) != 0;
}

// ============================================================================
// Public Query Methods - Atoms
// ============================================================================

bool AromaticityResolver::isAromaticAtom(const OpenBabel::OBAtom* atom) const {
    int atomIdx = toAtomIndex(atom);
    if (atomIdx < 0 || atomIdx >= (int)mAtomRingFeatures.size()) {
        return false;
    }
    return (mAtomRingFeatures[atomIdx] & FEATURES_AROMATIC) != 0;
}

bool AromaticityResolver::isDelocalizedAtom(const OpenBabel::OBAtom* atom) const {
    int atomIdx = toAtomIndex(atom);
    if (atomIdx < 0 || atomIdx >= (int)mAtomRingFeatures.size()) {
        return false;
    }
    return (mAtomRingFeatures[atomIdx] & FEATURES_DELOCALIZED) != 0;
}

bool AromaticityResolver::isHeteroAromaticAtom(const OpenBabel::OBAtom* atom) const {
    int atomIdx = toAtomIndex(atom);
    if (atomIdx < 0 || atomIdx >= (int)mAtomRingFeatures.size()) {
        return false;
    }
    return (mAtomRingFeatures[atomIdx] & FEATURES_HETERO_AROMATIC) != 0;
}

// ============================================================================
// Public Query Methods - Rings
// ============================================================================

int AromaticityResolver::getRingCount() const {
    return mRingAtoms.size();
}

bool AromaticityResolver::isAromaticRing(int ringNo) const {
    if (ringNo < 0 || ringNo >= (int)mIsAromaticRing.size()) {
        return false;
    }
    return mIsAromaticRing[ringNo];
}

bool AromaticityResolver::isDelocalizedRing(int ringNo) const {
    if (ringNo < 0 || ringNo >= (int)mIsDelocalizedRing.size()) {
        return false;
    }
    return mIsDelocalizedRing[ringNo];
}

int AromaticityResolver::getRingSize(int ringNo) const {
    if (ringNo < 0 || ringNo >= (int)mRingAtoms.size()) {
        return 0;
    }
    return mRingAtoms[ringNo].size();
}

const std::vector<int>& AromaticityResolver::getRingAtoms(int ringNo) const {
    static const std::vector<int> empty;
    if (ringNo < 0 || ringNo >= (int)mRingAtoms.size()) {
        return empty;
    }
    return mRingAtoms[ringNo];
}

const std::vector<int>& AromaticityResolver::getRingBonds(int ringNo) const {
    static const std::vector<int> empty;
    if (ringNo < 0 || ringNo >= (int)mRingBonds.size()) {
        return empty;
    }
    return mRingBonds[ringNo];
}

// ============================================================================
// Ring Detection
// ============================================================================

void AromaticityResolver::findSmallRings() {
    // Use OpenBabel's SSSR (Smallest Set of Smallest Rings)
    mMol.FindSSSR();

    std::vector<OpenBabel::OBRing*>& sssr = mMol.GetSSSR();

    for (OpenBabel::OBRing* obRing : sssr) {
        int ringSize = obRing->Size();

        // Only process small rings (up to MAX_SMALL_RING_SIZE)
        if (ringSize > MAX_SMALL_RING_SIZE) {
            continue;
        }

        // Limit number of rings to prevent explosion with complex structures
        if ((int)mRingAtoms.size() >= MAX_SMALL_RING_COUNT) {
            break;
        }

        // Extract atom indices (convert from 1-based to 0-based)
        std::vector<int> atomIndices;
        atomIndices.reserve(ringSize);
        for (int i = 0; i < ringSize; i++) {
            int atomIdx = obRing->_path[i] - 1;  // Convert to 0-based
            atomIndices.push_back(atomIdx);
        }

        // Extract bond indices
        std::vector<int> bondIndices;
        bondIndices.reserve(ringSize);
        for (int i = 0; i < ringSize; i++) {
            int atom1Idx = atomIndices[i];
            int atom2Idx = atomIndices[(i + 1) % ringSize];

            OpenBabel::OBAtom* atom1 = mMol.GetAtom(atom1Idx + 1);  // OB uses 1-based
            OpenBabel::OBAtom* atom2 = mMol.GetAtom(atom2Idx + 1);

            OpenBabel::OBBond* bond = mMol.GetBond(atom1, atom2);
            if (bond) {
                bondIndices.push_back(toBondIndex(bond));
            }
        }

        // Skip ring if any bond was not found (malformed/stale molecule data)
        if (bondIndices.size() != (size_t)ringSize) {
            continue;
        }

        // Store ring data
        mRingAtoms.push_back(atomIndices);
        mRingBonds.push_back(bondIndices);
    }

    // Initialize aromaticity arrays
    mIsAromaticRing.resize(mRingAtoms.size(), false);
    mIsDelocalizedRing.resize(mRingAtoms.size(), false);
    mAromaticityHandled.resize(mRingAtoms.size(), false);
    mHeteroPosition.resize(mRingAtoms.size(), -1);
}

// ============================================================================
// Annelation Map Building
// ============================================================================

void AromaticityResolver::buildAnnelationMap() {
    int numRings = getRingCount();

    // Initialize annelation map - for each ring, for each bond position, which other ring shares it (-1 if none)
    mAnnelatedRing.resize(numRings);
    for (int ringNo = 0; ringNo < numRings; ringNo++) {
        int ringSize = mRingBonds[ringNo].size();
        mAnnelatedRing[ringNo].resize(ringSize, -1);
    }

    // Find shared bonds
    for (int ring1 = 0; ring1 < numRings; ring1++) {
        for (size_t pos1 = 0; pos1 < mRingBonds[ring1].size(); pos1++) {
            int bond1 = mRingBonds[ring1][pos1];

            // Check if this bond is in any other ring
            for (int ring2 = ring1 + 1; ring2 < numRings; ring2++) {
                for (size_t pos2 = 0; pos2 < mRingBonds[ring2].size(); pos2++) {
                    int bond2 = mRingBonds[ring2][pos2];

                    if (bond1 == bond2) {
                        // Found shared bond!
                        mAnnelatedRing[ring1][pos1] = ring2;
                        mAnnelatedRing[ring2][pos2] = ring1;
                    }
                }
            }
        }
    }
}

// ============================================================================
// Aromaticity Determination
// ============================================================================

void AromaticityResolver::determineAromaticity() {
    int numRings = getRingCount();

    // Use iterative approach to handle complex annelation dependencies
    // Keep trying until no more rings become aromatic
    int maxPasses = 10;  // Prevent infinite loops

    for (int pass = 0; pass < maxPasses; pass++) {
        bool anyNewAromatic = false;

        // Try to determine aromaticity for all unhandled rings
        for (int ringNo = 0; ringNo < numRings; ringNo++) {
            if (!mAromaticityHandled[ringNo]) {
                bool wasAromatic = mIsAromaticRing[ringNo];
                determineAromaticityForRing(ringNo);

                // If ring became aromatic, mark as handled
                if (mIsAromaticRing[ringNo]) {
                    mAromaticityHandled[ringNo] = true;
                    if (!wasAromatic) {
                        anyNewAromatic = true;
                    }
                }
                // If ring is definitely not aromatic (no unhandled annelation), mark as handled
                else {
                    // Check if this ring has unhandled annelated rings
                    bool hasUnhandledAnnelation = false;
                    for (size_t i = 0; i < mAnnelatedRing[ringNo].size(); i++) {
                        int annelated = mAnnelatedRing[ringNo][i];
                        if (annelated != -1 && !mAromaticityHandled[annelated]) {
                            hasUnhandledAnnelation = true;
                            break;
                        }
                    }

                    if (!hasUnhandledAnnelation) {
                        mAromaticityHandled[ringNo] = true;
                    }
                }
            }
        }

        // If no new aromatic rings were found, we're done
        if (!anyNewAromatic) {
            // Mark all remaining rings as handled (not aromatic)
            for (int ringNo = 0; ringNo < numRings; ringNo++) {
                mAromaticityHandled[ringNo] = true;
            }
            break;
        }
    }
}

bool AromaticityResolver::determineAromaticityForRing(int ringNo) {
    const std::vector<int>& ringAtoms = mRingAtoms[ringNo];
    const std::vector<int>& ringBonds = mRingBonds[ringNo];
    int ringSize = ringAtoms.size();

    // Check if all atoms qualify as aromatic
    for (int atomIdx : ringAtoms) {
        if (!qualifiesAsAromaticAtom(atomIdx)) {
            if (debug && ringNo == 0) {
                std::cerr << "Ring 0: Atom " << atomIdx << " doesn't qualify as aromatic (atomicNo="
                          << mMol.GetAtom(atomIdx + 1)->GetAtomicNum() << ")" << std::endl;
            }
            return true;  // Not aromatic, but successfully determined
        }
    }

    // Build bond sequence (binary pattern of pi bonds)
    // This now handles annelated rings by treating shared aromatic bonds as pi bonds
    int aromaticButNotDelocalizedSequence = 0;
    int bondSequence = buildBondSequence(ringNo, aromaticButNotDelocalizedSequence);


    bool hasDelocalizationLeak = false;
    bool isAromatic = checkAromaticPattern(bondSequence, ringSize, ringNo,
                                          aromaticButNotDelocalizedSequence,
                                          hasDelocalizationLeak);


    if (isAromatic) {
        mIsAromaticRing[ringNo] = true;
        if (!hasDelocalizationLeak) {
            mIsDelocalizedRing[ringNo] = true;
        }
    }

    return true;
}

void AromaticityResolver::propagateRingFeatures() {
    int numRings = getRingCount();

    for (int ringNo = 0; ringNo < numRings; ringNo++) {
        if (!mIsAromaticRing[ringNo]) {
            continue;
        }

        bool isDelocalized = mIsDelocalizedRing[ringNo];
        bool hasHetero = mHeteroPosition[ringNo] >= 0;

        // Propagate to atoms
        for (int atomIdx : mRingAtoms[ringNo]) {
            mAtomRingFeatures[atomIdx] |= FEATURES_AROMATIC;
            if (isDelocalized) {
                mAtomRingFeatures[atomIdx] |= FEATURES_DELOCALIZED;
            }
            if (hasHetero) {
                mAtomRingFeatures[atomIdx] |= FEATURES_HETERO_AROMATIC;
            }
        }

        // Propagate to bonds
        for (int bondIdx : mRingBonds[ringNo]) {
            // Skip azulene bridge bonds
            if (mAzuleneBridgeBonds.find(bondIdx) != mAzuleneBridgeBonds.end()) {
                continue;
            }

            mBondRingFeatures[bondIdx] |= FEATURES_AROMATIC;
            if (isDelocalized) {
                mBondRingFeatures[bondIdx] |= FEATURES_DELOCALIZED;
            }
            if (hasHetero) {
                mBondRingFeatures[bondIdx] |= FEATURES_HETERO_AROMATIC;
            }
        }
    }
}

// ============================================================================
// Helper Methods
// ============================================================================

bool AromaticityResolver::qualifiesAsAromaticAtom(int atomIdx) const {
    OpenBabel::OBAtom* atom = mMol.GetAtom(atomIdx + 1);  // Convert to 1-based
    if (!atom) return false;

    int atomicNum = atom->GetAtomicNum();

    // Carbon, nitrogen, oxygen, sulfur, selenium, tellurium can be aromatic
    // Also boron (5) and phosphorus (15) in special cases
    return (atomicNum == 5 ||  // B
            atomicNum == 6 ||  // C
            atomicNum == 7 ||  // N
            atomicNum == 8 ||  // O
            atomicNum == 15 || // P
            atomicNum == 16 || // S
            atomicNum == 34 || // Se
            atomicNum == 52);  // Te
}

bool AromaticityResolver::qualifiesAsPiBond(int bondIdx) const {
    OpenBabel::OBBond* bond = mMol.GetBond(bondIdx);
    if (!bond) return false;

    // Pi bonds are double or triple bonds
    int bondOrder = bond->GetBondOrder();
    return bondOrder >= 2;
}

int AromaticityResolver::buildBondSequence(int ringNo, int& aromaticButNotDelocalizedSequence) const {
    const std::vector<int>& ringBonds = mRingBonds[ringNo];
    int bondSequence = 0;
    aromaticButNotDelocalizedSequence = 0;

    for (size_t i = 0; i < ringBonds.size(); i++) {
        int bondIdx = ringBonds[i];
        bondSequence <<= 1;
        aromaticButNotDelocalizedSequence <<= 1;

        if (qualifiesAsPiBond(bondIdx)) {
            // Regular pi bond (double/triple)
            bondSequence |= 1;
            if (debug) std::cerr << "  Bond " << i << ": pi bond (order >= 2)" << std::endl;
        }
        // TODO: else if (includeTautomericBonds && qualifiesAsAmideTypeBond(bondIdx))
        else {
            // Check if bond is shared with an annelated aromatic ring
            int annelated = mAnnelatedRing[ringNo][i];
            if (annelated != -1) {
                if (mAromaticityHandled[annelated]) {
                    if (mIsAromaticRing[annelated]) {
                        // Shared bond with aromatic ring - treat as pi bond!
                        bondSequence |= 1;
                        if (!mIsDelocalizedRing[annelated]) {
                            aromaticButNotDelocalizedSequence |= 1;
                        }
                        if (debug) std::cerr << "  Bond " << i << ": shared with aromatic ring " << annelated << std::endl;
                    }
                } else {
                    if (debug) std::cerr << "  Bond " << i << ": shared with unhandled ring " << annelated << std::endl;
                }
            } else {
                if (debug) std::cerr << "  Bond " << i << ": sigma bond (no annelation)" << std::endl;
            }
        }
    }

    if (debug) {
        std::cerr << "Final bondSequence: " << bondSequence << " (binary: ";
        for (int b = ringBonds.size()-1; b >= 0; b--) {
            std::cerr << ((bondSequence >> b) & 1);
        }
        std::cerr << ")" << std::endl;
    }

    return bondSequence;
}

bool AromaticityResolver::checkAromaticPattern(
    int bondSequence,
    int ringSize,
    int ringNo,
    int aromaticButNotDelocalizedSequence,
    bool& hasDelocalizationLeak) const {

    const std::vector<int>& ringAtoms = mRingAtoms[ringNo];
    hasDelocalizationLeak = false;

    switch (ringSize) {
    case 6:
        // 6-membered rings: require alternating pattern 010101 or 101010
        hasDelocalizationLeak = true;
        if ((bondSequence & 21) == 21) {   // 010101 = 21 in decimal
            if ((aromaticButNotDelocalizedSequence & 21) == 0) {
                hasDelocalizationLeak = false;
            }
            return true;
        }
        if ((bondSequence & 42) == 42) {   // 101010 = 42 in decimal
            if ((aromaticButNotDelocalizedSequence & 42) == 0) {
                hasDelocalizationLeak = false;
            }
            return true;
        }
        break;

    case 5:
        // 5-membered rings: require specific patterns with heteroatom
        // Pattern: alternating bonds with one heteroatom contributing lone pair
        {
            const int patterns[] = {10, 5, 18, 9, 20};  // 01010, 00101, 10010, 01001, 10100
            hasDelocalizationLeak = true;

            for (int position = 0; position < 5; position++) {
                if ((bondSequence & patterns[position]) == patterns[position]) {
                    OpenBabel::OBAtom* atom = mMol.GetAtom(ringAtoms[position] + 1);
                    int atomicNum = atom->GetAtomicNum();
                    int charge = atom->GetFormalCharge();

                    switch (atomicNum) {
                    case 6:  // Carbon
                        if (charge == -1) {  // Cyclopentadienyl anion
                            mHeteroPosition[ringNo] = position;
                            if ((aromaticButNotDelocalizedSequence & patterns[position]) == 0) {
                                hasDelocalizationLeak = false;
                            }
                            return true;
                        }
                        break;
                    case 7:  // Nitrogen
                        if (charge <= 0) {  // Pyrrole-type nitrogen
                            mHeteroPosition[ringNo] = position;
                            return true;
                        }
                        break;
                    case 8:  // Oxygen (furan)
                        mHeteroPosition[ringNo] = position;
                        return true;
                    case 16:  // Sulfur (thiophene)
                    case 34:  // Selenium
                    case 52:  // Tellurium
                        if (atom->GetTotalDegree() == 2) {
                            mHeteroPosition[ringNo] = position;
                            return true;
                        }
                        break;
                    }
                }
            }
        }
        break;

    case 3:
        // 3-membered rings: cyclopropenium cation (C+) or boron anion (B)
        {
            const int patterns[] = {2, 1, 4};  // 010, 001, 100
            hasDelocalizationLeak = true;

            for (int position = 0; position < 3; position++) {
                if ((bondSequence & patterns[position]) == patterns[position]) {
                    OpenBabel::OBAtom* atom = mMol.GetAtom(ringAtoms[position] + 1);
                    int atomicNum = atom->GetAtomicNum();
                    int charge = atom->GetFormalCharge();

                    if ((atomicNum == 6 && charge == 1) ||      // C+ (cyclopropenium)
                        (atomicNum == 5 && charge == 0)) {      // B (borane)
                        mHeteroPosition[ringNo] = position;
                        if ((aromaticButNotDelocalizedSequence & patterns[position]) == 0) {
                            hasDelocalizationLeak = false;
                        }
                        return true;
                    }
                }
            }
        }
        break;

    case 7:
        // 7-membered rings: tropylium cation
        {
            const int patterns[] = {42, 21, 74, 37, 82, 41, 84};
            // 0101010, 0010101, 1001010, 0100101, 1010010, 0101001, 1010100
            hasDelocalizationLeak = true;

            for (int position = 0; position < 7; position++) {
                if ((bondSequence & patterns[position]) == patterns[position]) {
                    OpenBabel::OBAtom* atom = mMol.GetAtom(ringAtoms[position] + 1);
                    int atomicNum = atom->GetAtomicNum();
                    int charge = atom->GetFormalCharge();

                    if ((atomicNum == 6 && charge == 1) ||      // C+ (tropylium)
                        (atomicNum == 5 && charge == 0)) {      // B
                        mHeteroPosition[ringNo] = position;
                        if ((aromaticButNotDelocalizedSequence & patterns[position]) == 0) {
                            hasDelocalizationLeak = false;
                        }
                        return true;
                    }

                    // TODO: Check for azulene-like annelated 7+5 ring systems
                    // This is complex and rarely needed, skipping for now
                }
            }
        }
        break;
    }

    return false;
}

} // namespace OpenChemLib

