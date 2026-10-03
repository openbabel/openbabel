/**
 * \file clogp.cpp
 * \brief Plugin to calculate cLogP (atom-increment model) in OpenBabel
 *
 * This implementation is based on OpenChemLib (https://github.com/actelion/openchemlib),
 * an open source Java-based chemistry library providing cheminformatics core functionality.
 * The cLogP calculation algorithm and constants have been adapted from OpenChemLib's
 * AtomTypeCalculator.java implementation.
 */
#include <openbabel/atom.h>
#include <openbabel/descriptor.h>
#include <openbabel/mol.h>
#include <openbabel/obiter.h>
#include <openbabel/ring.h>
#include <unordered_map>

// Include Java-compatible aromaticity detection
#include <openbabel/aromaticity_resolver.h>

// ============================================================================
// AMPHOLYTIC DETECTION - Full Java-compatible implementation
// Ported from AtomFunctionAnalyzer.java
// ============================================================================

namespace OpenChemLib {

using namespace OpenBabel;

// Helper: Count double bonds to N, O, or S (getFakeOxoCount)
static int getFakeOxoCount(OBMol& mol, const OBAtom* atom) {
    int count = 0;
    FOR_BONDS_OF_ATOM(bond, atom) {
        if (bond->GetBondOrder() == 2) {
            OBAtom* nbr = bond->GetNbrAtom(const_cast<OBAtom*>(atom));
            const unsigned int atomicNo = nbr->GetAtomicNum();
            if (atomicNo == 7 || atomicNo == 8 || atomicNo == 16) {
                ++count;
            }
        }
    }
    return count;
}

// Helper: Check if atom is stabilized by adjacent carbonyl/sulfonyl (isStabilized)
static bool isStabilized(OBMol& mol, const OBAtom* atom, const bool twice) {
    bool alreadyFound = false;
    FOR_BONDS_OF_ATOM(bond, atom) {
        if (!bond->IsAromatic() && bond->GetBondOrder() == 1) {
            OBAtom* conn = bond->GetNbrAtom(const_cast<OBAtom*>(atom));
            if (!conn->IsAromatic()) {
                const unsigned int connAtomicNo = conn->GetAtomicNum();
                const int fakeOxoCount = getFakeOxoCount(mol, conn);

                if ((connAtomicNo == 6 && fakeOxoCount == 1) ||
                    (connAtomicNo == 16 && fakeOxoCount == 2)) {
                    if (alreadyFound || !twice)
                        return true;
                    alreadyFound = true;
                }
            }
        }
    }
    return false;
}

// Helper: Check for vinylog fake oxo pattern (isVinylogFakeOxo)
static bool isVinylogFakeOxo(OBMol& mol, const OBAtom* atom) {
    FOR_BONDS_OF_ATOM(bond, atom) {
        if (bond->GetBondOrder() != 1) {
            OBAtom* conn = bond->GetNbrAtom(const_cast<OBAtom*>(atom));
            FOR_BONDS_OF_ATOM(bond2, conn) {
                if (bond2->GetBondOrder() == 1) {
                    OBAtom* conn2 = bond2->GetNbrAtom(conn);
                    if (getFakeOxoCount(mol, conn2) != 0)
                        return true;
                }
            }
        }
    }
    return false;
}

// Helper: Get atom pi electrons
static int getAtomPi(const OBAtom* atom) {
    int pi = 0;
    FOR_BONDS_OF_ATOM(bond, atom) {
        const unsigned int order = bond->GetBondOrder();
        if (order > 1) {
            pi += (order - 1);
        }
    }
    return pi;
}

// Helper: Count non-H neighbors
static int getConnAtoms(const OBAtom* atom) {
    int count = 0;
    FOR_NBORS_OF_ATOM(nbr, atom) {
        if (nbr->GetAtomicNum() != 1) {
            ++count;
        }
    }
    return count;
}

// MAIN: isAcidicOxygen - Java-compatible implementation
bool isAcidicOxygen(OBMol& mol, const OBAtom* atom, const bool considerCharge = true) {
    if (atom->GetAtomicNum() != 8)
        return false;

    if (considerCharge && atom->GetFormalCharge() != 0)
        return false;

    const int connCount = getConnAtoms(atom);
    if (connCount != 1)
        return false;

    // Check bond order (should be single bond)
    OBBond* bond = nullptr;
    FOR_BONDS_OF_ATOM(b, atom) {
        OBAtom* nbr = b->GetNbrAtom(const_cast<OBAtom*>(atom));
        if (nbr->GetAtomicNum() != 1) {
            bond = b;
            break;
        }
    }

    if (!bond || bond->GetBondOrder() != 1)
        return false;

    OBAtom* connAtom = bond->GetNbrAtom(const_cast<OBAtom*>(atom));

    // COOH pattern
    if (connAtom->GetAtomicNum() == 6) {
        FOR_NBORS_OF_ATOM(nbr, connAtom) {
            if (nbr->GetIdx() == atom->GetIdx())
                continue;

            if (nbr->GetAtomicNum() != 8)
                continue;

            OBBond* bond2 = mol.GetBond(connAtom, nbr);
            if (bond2 && bond2->GetBondOrder() == 2)
                return true;
        }
    }
    // (N+)-OH pattern
    else if (connAtom->GetAtomicNum() == 7) {
        if (connAtom->GetFormalCharge() == 1)
            return true;
    }
    // Sulfonic/sulfuric acid: S(=O)(=O)OH
    else if (connAtom->GetAtomicNum() == 16) {
        int nDoubleBondedO = 0;
        FOR_NBORS_OF_ATOM(nbr, connAtom) {
            if (nbr->GetIdx() == atom->GetIdx())
                continue;

            if (nbr->GetAtomicNum() != 8)
                continue;

            OBBond* bond2 = mol.GetBond(connAtom, nbr);
            if (bond2 && bond2->GetBondOrder() == 2)
                ++nDoubleBondedO;
        }

        if (nDoubleBondedO == 2)
            return true;
    }

    return false;
}

// MAIN: isBasicNitrogen - Full Java-compatible implementation (200+ lines)
bool isBasicNitrogen(OBMol& mol, const OBAtom* atom, const AromaticityResolver* aromaticity, const bool considerCharge = true) {
    if (atom->GetAtomicNum() != 7)
        return false;

    if (considerCharge && atom->GetFormalCharge() != 0)
        return false;

    const int connAtoms = getConnAtoms(atom);
    const int atomPi = getAtomPi(atom);

    if (connAtoms + atomPi > 3)
        return false;

    // ===== AROMATIC NITROGEN HANDLING =====
    if (aromaticity && aromaticity->isAromaticAtom(atom)) {
        if (atomPi != 1)
            return false;  // pyrrol type

        // Count rings containing this atom (up to size 7)
        int ringCount = 0;
        for (int size = 3; size <= 7; size++) {
            if (atom->IsInRingSize(size))
                ringCount++;
        }

        if (ringCount != 1)
            return false;

        // Find the aromatic ring containing this nitrogen
        const std::vector<OBRing*>& rings = mol.GetSSSR();
        for (OBRing* ring : rings) {
            if (!ring->IsAromatic())
                continue;

            if (!ring->IsInRing(atom->GetIdx()))
                continue;

            const int ringSize = ring->Size();
            if (ringSize != 5 && ringSize != 6)
                continue;

            const std::vector<int>& ringAtoms = ring->_path;

            // Find nitrogen index in ring
            int nIndex = -1;
            for (size_t i = 0; i < ringAtoms.size(); i++) {
                if (ringAtoms[i] == (int)atom->GetIdx()) {
                    nIndex = i;
                    break;
                }
            }

            if (nIndex == -1)
                continue;

            int enablerCount = 0;

            // Setup ortho/para and meta positions
            std::vector<int> opi, mi;  // ortho-para influences, meta influences

            if (ringSize == 5) {
                opi.push_back(ringAtoms[(nIndex - 1 < 0) ? nIndex + 4 : nIndex - 1]);
                opi.push_back(ringAtoms[(nIndex - 4 < 0) ? nIndex + 1 : nIndex - 4]);
                mi.push_back(ringAtoms[(nIndex - 2 < 0) ? nIndex + 3 : nIndex - 2]);
                mi.push_back(ringAtoms[(nIndex - 3 < 0) ? nIndex + 2 : nIndex - 3]);
            }

            if (ringSize == 6) {
                opi.push_back(ringAtoms[(nIndex - 1 < 0) ? nIndex + 5 : nIndex - 1]);
                opi.push_back(ringAtoms[(nIndex - 3 < 0) ? nIndex + 3 : nIndex - 3]);
                opi.push_back(ringAtoms[(nIndex - 5 < 0) ? nIndex + 1 : nIndex - 5]);
                mi.push_back(ringAtoms[(nIndex - 2 < 0) ? nIndex + 4 : nIndex - 2]);
                mi.push_back(ringAtoms[(nIndex - 4 < 0) ? nIndex + 2 : nIndex - 4]);
            }

            // Count other aromatic nitrogens in ring (decrease basicity)
            for (int atomIdx : ringAtoms) {
                OBAtom* ringAtom = mol.GetAtom(atomIdx);
                if (atomIdx != (int)atom->GetIdx() &&
                    ringAtom->GetAtomicNum() == 7 &&
                    getAtomPi(ringAtom) == 1) {
                    enablerCount--;
                }
            }

            // Check ortho/para positions for electron-donating groups
            for (int opiIdx : opi) {
                OBAtom* opiAtom = mol.GetAtom(opiIdx);

                // Find exocyclic atom (first non-aromatic bond from opiAtom)
                OBAtom* exoCyclicAtom = nullptr;

                FOR_BONDS_OF_ATOM(bond, opiAtom) {
                    if (!aromaticity->isAromaticBond(bond)) {
                        exoCyclicAtom = bond->GetNbrAtom(opiAtom);
                        break;
                    }
                }

                if (exoCyclicAtom) {
                    // Check for NH2 group (electron-donating)
                    if (exoCyclicAtom->GetAtomicNum() == 7 &&
                        getAtomPi(exoCyclicAtom) == 0 &&
                        (getConnAtoms(exoCyclicAtom) + getAtomPi(exoCyclicAtom) <= 3) &&
                        !isStabilized(mol, exoCyclicAtom, false)) {
                        enablerCount++;
                        continue;
                    }

                    // Check for OH group (electron-donating)
                    if (exoCyclicAtom->GetAtomicNum() == 8 &&
                        getConnAtoms(exoCyclicAtom) == 1) {
                        enablerCount += 2;
                        continue;
                    }

                    // Check if the substituent is itself in an aromatic ring
                    // containing an electron-withdrawing nitrogen (decreases basicity)
                    if (aromaticity->isAromaticAtom(exoCyclicAtom)) {
                        for (OBRing* ring2 : rings) {
                            if (ring2->IsAromatic() && ring2->IsInRing(exoCyclicAtom->GetIdx())) {
                                const std::vector<int>& ratom = ring2->_path;
                                for (const int ratomIdx : ratom) {
                                    const OBAtom* ra = mol.GetAtom(ratomIdx);
                                    if (ra->GetAtomicNum() == 7 && getAtomPi(ra) == 1) {
                                        --enablerCount;
                                        break;
                                    }
                                }
                                break;
                            }
                        }
                    }
                }
            }

            // Check meta positions for electron-withdrawing groups
            for (int miIdx : mi) {
                OBAtom* miAtom = mol.GetAtom(miIdx);

                // Find exocyclic atom
                OBAtom* exoCyclicAtom = nullptr;
                FOR_BONDS_OF_ATOM(bond, miAtom) {
                    if (!aromaticity->isAromaticBond(bond)) {
                        exoCyclicAtom = bond->GetNbrAtom(miAtom);
                        break;
                    }
                }

                if (miAtom->GetAtomicNum() == 6) {
                    if (exoCyclicAtom && getFakeOxoCount(mol, exoCyclicAtom) != 0)
                        enablerCount--;
                }
                else if (miAtom->GetAtomicNum() == 7) {
                    if (getAtomPi(miAtom) == 0 &&
                        (!exoCyclicAtom ||
                         (!aromaticity->isAromaticAtom(exoCyclicAtom) &&
                          getFakeOxoCount(mol, exoCyclicAtom) == 0))) {
                        enablerCount++;  // imidazole type
                    }
                }
            }

            return enablerCount > 0;
        }

        return false;
    }

    // ===== NITRILE CHECK =====
    if (atomPi > 1)
        return false;  // nitrile

    // ===== IMINE HANDLING (atomPi == 1) =====
    if (atomPi == 1) {
        OBAtom* imineC = nullptr;
        int supporterCount = 0;

        FOR_BONDS_OF_ATOM(bond, atom) {
            OBAtom* conn = bond->GetNbrAtom(const_cast<OBAtom*>(atom));
            const int bondOrder = bond->GetBondOrder();

            if (bondOrder == 2) {
                if (conn->GetAtomicNum() != 6)
                    return false;  // N=O, N=N, N=S, etc.
                imineC = conn;
                continue;
            }

            if (conn->GetAtomicNum() == 8)  // C=N-O
                return false;

            if (conn->GetAtomicNum() == 7) {
                --supporterCount;
                if (isStabilized(mol, conn, false))  // C=N-N-C=O
                    --supporterCount;
                continue;
            }

            if (aromaticity && aromaticity->isAromaticAtom(conn))
                --supporterCount;
        }

        if (!imineC)
            return false;

        int aromaticNeighborCount = 0;
        FOR_BONDS_OF_ATOM(bond, imineC) {
            if (bond->GetBondOrder() == 1) {
                OBAtom* conn = bond->GetNbrAtom(imineC);

                if (getFakeOxoCount(mol, conn) != 0)
                    return false;

                if (aromaticity && aromaticity->isAromaticAtom(conn))
                    ++aromaticNeighborCount;

                if (conn->GetAtomicNum() == 7 && !isStabilized(mol, conn, true))
                    ++supporterCount;

                if (conn->GetAtomicNum() == 8 || conn->GetAtomicNum() == 16)
                    --supporterCount;  // S-C=N or O-C=N mostly below pKa=7
            }
        }

        if (aromaticNeighborCount == 2)  // two phenyl substituents reduce pKa to 6
            --supporterCount;

        return (supporterCount >= 0);
    }

    // ===== SATURATED AMINE HANDLING (atomPi == 0) =====
    FOR_BONDS_OF_ATOM(bond, atom) {
        OBAtom* conn = bond->GetNbrAtom(const_cast<OBAtom*>(atom));

        if (aromaticity && aromaticity->isAromaticAtom(conn))
            return false;

        if (conn->GetAtomicNum() != 6)
            return false;

        if (getFakeOxoCount(mol, conn) != 0)
            return false;

        if (getAtomPi(conn) != 0 && isVinylogFakeOxo(mol, conn))
            return false;
    }

    return true;
}

}  // namespace OpenChemLib

// ============================================================================
// END AMPHOLYTIC DETECTION
// ============================================================================

namespace {

constexpr uint64_t ATOM_FLAG_COUNT = 15;
constexpr uint64_t CONN_FLAG_COUNT = 11;
constexpr int MAX_NEIGHBORS = 4;  // Maximum neighbors to consider for atom type

constexpr uint64_t cPropertiesAtomSmallRing = 0x00000001;
constexpr uint64_t cPropertiesAtomCharged = 0x00001000;
constexpr uint64_t cPropertiesConnBondOrder = 0x00000020;
constexpr uint64_t cPropertiesConnAtomTypeSimple = 0x00000040;
constexpr uint64_t cPropertiesConnAtomAromatic = 0x00000800;

// Mode for cLogP: Use the EXACT same mode as OpenChemLib
// constexpr uint64_t MODE_CLOGP = 0x00001861;
constexpr uint64_t MODE_CLOGP =
    cPropertiesAtomSmallRing | cPropertiesConnBondOrder |
    cPropertiesConnAtomTypeSimple | cPropertiesConnAtomAromatic |
    cPropertiesAtomCharged;

// Atom flags and shifts
constexpr uint64_t ATOM_FLAG_SMALLRING = 0x00000040;   // bit 6
constexpr uint64_t ATOM_FLAG_CHARGED = 0x00002000;     // bit 13
constexpr uint64_t ATOM_FLAG_AMPHOLYTIC = 0x00004000;  // bit 14

// Connection flags and shifts
constexpr uint64_t CONN_SHIFT_BONDORDER = 4;
constexpr uint64_t CONN_FLAG_AROMATIC = 0x00000200;

static const int8_t cAtomicNoCode[] = {
    -1, -1, -1, 0,  0,  1,  2,  //  H  ,He ,Li ,Be ,B  ,C  ,
    3,  4,  5,  -1, 0,  0,      //  N , O  ,F  ,Ne ,Na ,Mg ,
    0,  6,  7,  8,  9,  -1,     //  Al ,Si ,P  ,S  ,Cl ,Ar ,
    0,  0,  10, 10, 10, 10,     //  K  ,Ca ,Sc ,Ti ,V  ,Cr ,
    10, 10, 10, 10, 10, 10,     //  Mn ,Fe ,Co ,Ni ,Cu ,Zn ,
    1,  11, 11, 12, 13, -1,     //  Ga ,Ge ,As ,Se ,Br ,Kr ,
    0,  0,  10, 10, 10, 10,     //  Rb ,Sr ,Y  ,Zr ,Nb ,Mo ,
    10, 10, 10, 10, 10, 10,     //  Tc ,Ru ,Rh ,Pd ,Ag ,Cd ,
    0,  0,  0,  11, 14, -1,     //  In ,Sn ,Sb ,Te ,I  ,Xe ,
    0,  0,  15, 15, 15, 15,     //  Cs ,Ba ,La ,Ce ,Pr ,Nd ,
    15, 15, 15, 15, 15, 15,     //  Pm ,Sm ,Eu ,Gd ,Tb ,Dy ,
    15, 15, 15, 15, 15, 15,     //  Ho ,Er ,Tm ,Yb ,Lu ,Hf ,
    10, 10, 10, 10, 10, 10,     //  Ta ,W , Re ,Os ,Ir ,Pt ,
    10, 10, 1,  1,  1,  1,      //  Au ,Hg ,Tl ,Pb ,Bi ,Po ,
    -1, -1, -1, -1, 15, 15,     //  At ,Rn ,Fr ,Ra ,Ac ,Th ,
    15, 15, 15, 15, 15, 15,     //  Pa ,U , Np ,Pu ,Am ,Cm ,
    15, 15, 15, 15, 15, 15,     //  Bk ,Cf ,Es ,Fm ,Md ,No ,
    15, -1, -1, -1, -1, -1,     //  Lr ,?? ,?? ,?? ,?? ,?? ,
    -1, -1, -1, -1, -1, -1,     //  ?? ,?? ,?? ,?? ,?? ,?? ,
    -1, -1, -1, -1, -1, -1,     //  ?? ,?? ,?? ,?? ,?? ,?? ,
    -1, -1, -1, -1, -1, -1,     //  ?? ,?? ,?? ,?? ,?? ,?? ,
    -1, -1, -1, -1, -1, -1,     //  ?? ,?? ,?? ,?? ,?? ,?? ,
    -1, -1, -1, -1, -1, -1,     //  ?? ,?? ,?? ,?? ,?? ,?? ,
    -1, -1, -1, -1, -1, -1,     //  ?? ,?? ,?? ,R1 ,R2 ,R3 ,
    -1, -1, -1, -1, -1, -1,     //  A  ,A1 ,A2 ,A3 ,?? ,?? ,
    -1, -1, -1, -1, -1, -1,     //  D  ,T  ,X  ,R  ,H2 ,H+
    -1, -1, -1, -1, -1, -1,     //  Nnn,HYD,Pol,?? ,?? ,?? ,
    -1, -1, -1, -1, -1, -1,     //  ?? ,?? ,?? ,?? ,?? ,?? ,
    -1, -1, -1, -1, -1, -1,     //  ?? ,?? ,Ala,Arg,Asn,Asp,
    -1, -1, -1, -1, -1, -1,     //  Cys,Gln,Glu,Gly,His,Ile,
    -1, -1, -1, -1, -1, -1,     //  Leu,Lys,Met,Phe,Pro,Ser,
    -1, -1, -1, -1};            //  Thr,Trp,Tyr,Val,

static const int8_t cSimpleAtomicNoCode[] = {
    -1, -1, -1, 0,  0,  0,  2,  //  H  ,He ,Li ,Be ,B  ,C  ,
    5,  5,  5,  -1, 0,  0,      //  N , O  ,F  ,Ne ,Na ,Mg ,
    0,  0,  9,  9,  9,  -1,     //  Al ,Si ,P  ,S  ,Cl ,Ar ,
    0,  0,  0,  0,  0,  0,      //  K  ,Ca ,Sc ,Ti ,V  ,Cr ,
    0,  0,  0,  0,  0,  0,      //  Mn ,Fe ,Co ,Ni ,Cu ,Zn ,
    0,  0,  0,  9,  9,  -1,     //  Ga ,Ge ,As ,Se ,Br ,Kr ,
    0,  0,  0,  0,  0,  0,      //  Rb ,Sr ,Y  ,Zr ,Nb ,Mo ,
    0,  0,  0,  0,  0,  0,      //  Tc ,Ru ,Rh ,Pd ,Ag ,Cd ,
    0,  0,  0,  0,  9,  -1,     //  In ,Sn ,Sb ,Te ,I  ,Xe ,
    0,  0,  0,  0,  0,  0,      //  Cs ,Ba ,La ,Ce ,Pr ,Nd ,
    0,  0,  0,  0,  0,  0,      //  Pm ,Sm ,Eu ,Gd ,Tb ,Dy ,
    0,  0,  0,  0,  0,  0,      //  Ho ,Er ,Tm ,Yb ,Lu ,Hf ,
    0,  0,  0,  0,  0,  0,      //  Ta ,W , Re ,Os ,Ir ,Pt ,
    0,  0,  0,  0,  0,  0,      //  Au ,Hg ,Tl ,Pb ,Bi ,Po ,
    -1, -1, -1, -1, 0,  0,      //  At ,Rn ,Fr ,Ra ,Ac ,Th ,
    0,  0,  0,  0,  0,  0,      //  Pa ,U , Np ,Pu ,Am ,Cm ,
    0,  0,  0,  0,  0,  0,      //  Bk ,Cf ,Es ,Fm ,Md ,No ,
    0,  -1, -1, -1, -1, -1,     //  Lr ,?? ,?? ,?? ,?? ,?? ,
    -1, -1, -1, -1, -1, -1,     //  ?? ,?? ,?? ,?? ,?? ,?? ,
    -1, -1, -1, -1, -1, -1,     //  ?? ,?? ,?? ,?? ,?? ,?? ,
    -1, -1, -1, -1, -1, -1,     //  ?? ,?? ,?? ,?? ,?? ,?? ,
    -1, -1, -1, -1, -1, -1,     //  ?? ,?? ,?? ,?? ,?? ,?? ,
    -1, -1, -1, -1, -1, -1,     //  ?? ,?? ,?? ,?? ,?? ,?? ,
    -1, -1, -1, -1, -1, -1,     //  ?? ,?? ,?? ,R1 ,R2 ,R3 ,
    -1, -1, -1, -1, -1, -1,     //  A  ,A1 ,A2 ,A3 ,?? ,?? ,
    -1, -1, -1, -1, -1, -1,     //  D  ,T  ,X  ,R  ,H2 ,H+
    -1, -1, -1, -1, -1, -1,     //  Nnn,HYD,Pol,?? ,?? ,?? ,
    -1, -1, -1, -1, -1, -1,     //  ?? ,?? ,?? ,?? ,?? ,?? ,
    -1, -1, -1, -1, -1, -1,     //  ?? ,?? ,Ala,Arg,Asn,Asp,
    -1, -1, -1, -1, -1, -1,     //  Cys,Gln,Glu,Gly,His,Ile,
    -1, -1, -1, -1, -1, -1,     //  Leu,Lys,Met,Phe,Pro,Ser,
    -1, -1, -1, -1};            //  Thr,Trp,Tyr,Val,
}

static const long ATOM_TYPES[] = {
    0x80002L,           0x80004L,           0x90002L,
    0x90003L,           0x90004L,           0x90005L,
    0x90008L,           0x90009L,           0x9000dL,
    0x9000eL,           0x92003L,           0x92004L,
    0x92008L,           0x94003L,           0xa8002L,
    0xa8003L,           0xa8004L,           0xa8009L,
    0xa800eL,           0xaa004L,           0xc8002L,
    0xc8003L,           0xc8004L,           0xc8005L,
    0xc800eL,           0xca004L,           0x100004L,
    0x110002L,          0x110003L,          0x110004L,
    0x110008L,          0x128004L,          0x12a003L,
    0x148004L,          0x148008L,          0x190002L,
    0x190003L,          0x1090002L,         0x1090003L,
    0x1090004L,         0x1090005L,         0x1090008L,
    0x1090009L,         0x109000dL,         0x109000eL,
    0x1092004L,         0x10a8002L,         0x10a8004L,
    0x10a8009L,         0x10aa004L,         0x48080002L,
    0x48080008L,        0x48080048L,        0x48090002L,
    0x48090003L,        0x48090004L,        0x48090008L,
    0x4809000cL,        0x48090042L,        0x48090043L,
    0x48090044L,        0x48090048L,        0x48092003L,
    0x48094003L,        0x48094043L,        0x54090002L,
    0x54090003L,        0x54090004L,        0x54090008L,
    0x54090042L,        0x54090043L,        0x54090044L,
    0x54090048L,        0x54092003L,        0x540a8002L,
    0x540a8008L,        0x540a8042L,        0x64080004L,
    0x64090002L,        0x64090003L,        0x64090004L,
    0x64090008L,        0x64090042L,        0x64090043L,
    0x64090044L,        0x64090048L,        0x640a8002L,
    0x640a8003L,        0x640a8004L,        0x640a8042L,
    0x640c8002L,        0x640c8008L,        0x640c800aL,
    0x640c8042L,        0x88090002L,        0x88090003L,
    0x88090042L,        0x88090043L,        0x880a8002L,
    0x880a8003L,        0x880a8042L,        0x880a8043L,
    0x880a8044L,        0x880c8002L,        0x880c8003L,
    0x880c8042L,        0x880c8043L,        0x880cc043L,
    0x94090002L,        0x94090003L,        0x94090042L,
    0x94090043L,        0x94090044L,        0x94090048L,
    0x940a8002L,        0x940a8003L,        0x940a8042L,
    0x940a8043L,        0x94112003L,        0x9412a003L,
    0xa4128002L,        0xc8090002L,        0xc80aa003L,
    0xd4090002L,        0xd40a8002L,        0xd40c8002L,
    0x809010042L,       0x809010043L,       0x809010044L,
    0x809010048L,       0x809012043L,       0x809014043L,
    0x815010042L,       0x815010043L,       0x815010044L,
    0x815010048L,       0x815028042L,       0x815028043L,
    0x815028044L,       0x815028048L,       0x825010042L,
    0x825010043L,       0x825028042L,       0x825028043L,
    0x825028044L,       0x848090002L,       0x848090003L,
    0x848090004L,       0x848090008L,       0x84809000cL,
    0x848090042L,       0x848090043L,       0x848090044L,
    0x848090048L,       0x84809004cL,       0x8480a8002L,
    0x8480a8003L,       0x8480a8004L,       0x8480a8042L,
    0x8480a8048L,       0x8480c8002L,       0x8480c8003L,
    0x8480c8004L,       0x8480c8008L,       0x8480c8042L,
    0x8480c8043L,       0x8480c8044L,       0x8480c8048L,
    0x848110002L,       0x848110003L,       0x848110042L,
    0x848110043L,       0x848128002L,       0x848128003L,
    0x848128043L,       0x848190002L,       0x8481a8002L,
    0x849090002L,       0x849090003L,       0x849090004L,
    0x849090008L,       0x849090042L,       0x849090043L,
    0x849090044L,       0x849090048L,       0x854090002L,
    0x854090004L,       0x854090042L,       0x8540a8002L,
    0x854110002L,       0x854128002L,       0x855090002L,
    0x855090003L,       0x24048080043L,     0x24048090002L,
    0x24048090003L,     0x24048090042L,     0x24048090043L,
    0x24048092003L,     0x24048092043L,     0x24048094003L,
    0x24048094043L,     0x2a048090002L,     0x2a048090003L,
    0x2a048090042L,     0x2a048090043L,     0x2a054090002L,
    0x2a054090003L,     0x2a054090042L,     0x2a0540a8002L,
    0x2a0540a8042L,     0x32048090002L,     0x32048090003L,
    0x32048090042L,     0x32048090043L,     0x32054090002L,
    0x32054090003L,     0x32054090042L,     0x320540a8002L,
    0x32064090002L,     0x320640c8002L,     0x44048090002L,
    0x44048090042L,     0x44048090043L,     0x44048092003L,
    0x44054090002L,     0x44054090042L,     0x44054090043L,
    0x44054092003L,     0x44054092043L,     0x440540a8002L,
    0x440540a8042L,     0x44064090002L,     0x44064090042L,
    0x440640a8042L,     0x440640c8002L,     0x440640c8042L,
    0x4a048090002L,     0x4a048090008L,     0x4a048090042L,
    0x4a048090043L,     0x4a048090048L,     0x4a054090002L,
    0x4a054090042L,     0x4a054092003L,     0x4a0540a8002L,
    0x4a0540a8007L,     0x4a0540a8042L,     0x4a0540a8048L,
    0x4a0540aa003L,     0x4a064090002L,     0x4a064090042L,
    0x4a0640a8002L,     0x4a0640a8042L,     0x4a0640c8042L,
    0x52048090042L,     0x52054090002L,     0x52054090042L,
    0x520540a8002L,     0x520540a8042L,     0x52064090042L,
    0x520640a8002L,     0x520640a8042L,     0x404808090042L,
    0x404808090043L,    0x404808092043L,    0x4048080a8042L,
    0x4048080a8043L,    0x4048080aa043L,    0x4048080c8042L,
    0x404809010042L,    0x404809010043L,    0x40a808090042L,
    0x40a808090043L,    0x40a808092043L,    0x40a8080a8042L,
    0x40a8080a8043L,    0x40a8080aa043L,    0x40a8080c8042L,
    0x40a809010042L,    0x40a809010043L,    0x40a814090042L,
    0x40a814090043L,    0x40a8140a8042L,    0x40a8140c8042L,
    0x40a815010042L,    0x40a815010048L,    0x40a815028042L,
    0x412808090042L,    0x4128080a8042L,    0x4128080c8042L,
    0x412809010042L,    0x412814090042L,    0x4128140a8042L,
    0x4128140c8042L,    0x412815010042L,    0x412825010042L,
    0x424048090002L,    0x424048090003L,    0x424048090042L,
    0x424048090043L,    0x424054090002L,    0x424054090003L,
    0x424054090042L,    0x424054090043L,    0x4240540a8002L,
    0x4240540a8042L,    0x424064090002L,    0x424064090003L,
    0x424064090042L,    0x424064090043L,    0x4240640a8003L,
    0x4240640a8042L,    0x4240640c8003L,    0x424088090002L,
    0x424088090042L,    0x4240880a8002L,    0x4240880a8042L,
    0x4240880aa003L,    0x4240880c8002L,    0x4240880c8008L,
    0x4240880c8042L,    0x424094090002L,    0x424094090008L,
    0x424094090042L,    0x424094090043L,    0x424094090048L,
    0x4240940a8002L,    0x4240940a8008L,    0x4240940a8042L,
    0x4240940aa003L,    0x4240940c8002L,    0x4240940c8042L,
    0x4240a40a8002L,    0x4240a40a8042L,    0x424809010042L,
    0x424809010043L,    0x424815010042L,    0x424815010043L,
    0x424815012043L,    0x424815028042L,    0x424815028043L,
    0x424825010042L,    0x424825028042L,    0x424848090002L,
    0x424848090003L,    0x424848090042L,    0x424848090043L,
    0x4248480a8002L,    0x4248480a8042L,    0x424848110002L,
    0x424848110042L,    0x424848128002L,    0x424848128008L,
    0x424848128042L,    0x424849090002L,    0x42a048090002L,
    0x42a048090042L,    0x42a054090002L,    0x42a054090042L,
    0x42a064090042L,    0x42a094090002L,    0x42a0940a8002L,
    0x42a0a40a8002L,    0x42a809010042L,    0x42a815010042L,
    0x42a815028042L,    0x42a848090002L,    0x42a848110042L,
    0x42a848128002L,    0x42a849090002L,    0x12024048080043L,
    0x12024048090001L,  0x12024048090002L,  0x12024048090042L,
    0x12024048092003L,  0x12024048092043L,  0x15024048090002L,
    0x15024048090042L,  0x15024048092003L,  0x15024048092043L,
    0x1502a048090002L,  0x1502a048090042L,  0x1502a054090002L,
    0x1502a054090042L,  0x1502a0540a8002L,  0x19024048090002L,
    0x19024048090003L,  0x19024048090042L,  0x19024048090043L,
    0x1902a048090002L,  0x1902a054090002L,  0x1902a0540a8002L,
    0x19032048090002L,  0x19032048090042L,  0x19032054090002L,
    0x19032064090002L,  0x190320640a8002L,  0x190320640c8002L,
    0x22032048090003L,  0x2502a048090007L,  0x2502a054090007L,
    0x2502a0540a8007L,  0x2502a0540a8047L,  0x25032054090007L,
    0x250320540a8007L,  0x250320640a8007L,  0x250320640a804aL,
    0x2504a048090008L,  0x2504a048090048L,  0x2504a054090008L,
    0x2504a0540a8008L,  0x2504a0540a8048L,  0x2902a054090007L,
    0x2902a0540a8007L,  0x290320540a8007L,  0x202404064090043L,
    0x212024048090002L, 0x212024048090042L, 0x212024048092003L,
    0x21202a048090002L, 0x21202a048090042L, 0x21202a054090002L,
    0x21202a054090042L, 0x21202a0540a8002L, 0x2120320640c8002L,
    0x21204a0540a8007L, 0x21204a094090008L, 0x21204a094090048L,
    0x21204a0940a8008L, 0x21204a0940a8048L, 0x2120520540a8007L,
    0x212424048090002L, 0x212424048090042L, 0x212424054090002L,
    0x212424054090042L, 0x212424094128008L, 0x212424094128048L,
    0x2124248480a8002L, 0x2124248480a8042L, 0x215024048090002L,
    0x215024048090042L};

static const float INCREMENTS[] = {
    0.6967f,  0.0f,     0.4886f,  -0.4727f, -0.0749f, 0.6262f,  0.2735f,
    0.57f,    0.701f,   0.9534f,  -3.6435f, -2.1509f, 0.4975f,  -2.1995f,
    -0.2809f, -0.826f,  -0.1785f, -1.6203f, -1.096f,  -0.362f,  0.1395f,
    -0.2975f, -1.2908f, 1.0162f,  -1.3825f, -2.5384f, 0.3317f,  0.4292f,
    -0.5824f, -0.1834f, 0.1306f,  -0.5015f, 0.6262f,  -0.5258f, 0.4244f,
    -0.161f,  -0.2778f, 0.5111f,  -0.4357f, -0.1041f, 0.3424f,  -0.0615f,
    0.6035f,  0.7227f,  0.4346f,  -1.6821f, -0.331f,  -0.498f,  -1.4915f,
    -0.3651f, 0.4597f,  0.3384f,  0.6633f,  0.4544f,  0.1597f,  0.6339f,
    0.3504f,  0.0449f,  0.342f,   0.2611f,  0.4046f,  0.5219f,  -1.065f,
    -1.2314f, -1.8023f, -0.3632f, -0.4108f, 0.3057f,  -0.1456f, -0.2713f,
    -0.5193f, 0.4526f,  0.5539f,  0.8374f,  -0.707f,  -0.4881f, -0.41f,
    0.0f,     0.1479f,  0.3448f,  0.4298f,  0.5579f,  -0.1265f, -0.0425f,
    0.0767f,  0.6635f,  -0.3812f, -0.8368f, 1.0287f,  -0.1021f, 0.3587f,
    -0.5945f, 0.1692f,  -0.1218f, 0.3283f,  0.2239f,  0.2043f,  0.059f,
    -0.4835f, 0.6165f,  -0.4011f, 0.5578f,  -0.2164f, -0.0175f, 0.2981f,
    0.11f,    0.2715f,  -0.2995f, -0.467f,  0.1566f,  0.0468f,  -0.1321f,
    1.3686f,  0.0f,     -0.4116f, 1.0186f,  -0.3935f, 0.5223f,  0.3685f,
    0.2577f,  1.5193f,  0.2705f,  0.3791f,  0.012f,   -0.3397f, 0.1483f,
    0.2766f,  0.3593f,  0.7715f,  0.315f,   -1.6185f, 0.1889f,  -0.2652f,
    -0.0965f, 0.4202f,  0.1871f,  -0.3684f, -0.0778f, 0.8943f,  0.3694f,
    0.2879f,  0.4489f,  -0.2601f, 0.4771f,  0.1923f,  0.4381f,  0.1695f,
    0.4525f,  0.3352f,  0.1583f,  0.4036f,  -0.048f,  0.5023f,  -0.2649f,
    0.7691f,  -0.3552f, 1.0301f,  -0.1141f, -0.5932f, 0.1749f,  0.1313f,
    -0.1804f, 0.3994f,  0.2291f,  0.3169f,  0.3599f,  -0.0039f, -0.2956f,
    0.4409f,  -0.1609f, 0.3775f,  -0.1346f, 0.2839f,  0.5129f,  0.1266f,
    0.4294f,  0.2806f,  0.4907f,  0.354f,   0.2192f,  0.1565f,  0.6935f,
    0.3618f,  0.6735f,  0.5778f,  -0.5636f, 0.5569f,  0.3038f,  -0.3276f,
    -0.6992f, 0.0103f,  -0.4659f, 0.3818f,  0.3317f,  0.1841f,  0.7071f,
    0.1227f,  0.7949f,  -0.6592f, -1.3148f, -0.4067f, -0.1316f, -0.4925f,
    -0.0929f, -0.4352f, -0.2207f, -0.996f,  -0.7238f, -0.5469f, -1.2939f,
    -0.0136f, 0.0657f,  0.7189f,  0.057f,   0.6619f,  -0.6381f, -0.8073f,
    0.2355f,  0.3048f,  -0.0199f, -0.0752f, 0.4461f,  0.1559f,  1.1168f,
    -1.8039f, 0.2365f,  -0.2206f, 0.448f,   -1.134f,  -0.5653f, -0.4053f,
    -0.1361f, 0.2199f,  0.0536f,  -0.021f,  0.6985f,  0.9643f,  -0.4152f,
    -1.0369f, -0.183f,  0.5883f,  -0.2918f, -0.5294f, -0.6541f, 0.9473f,
    -0.1906f, -0.8484f, -0.3457f, 0.9541f,  1.4231f,  -0.7924f, -0.602f,
    0.08f,    -0.2596f, 0.8382f,  -0.4416f, -0.3704f, -0.7487f, -0.1079f,
    -0.2992f, -0.3277f, 0.0251f,  -0.9188f, 0.1094f,  0.8231f,  -3.2333f,
    0.035f,   0.3818f,  -0.738f,  0.2791f,  0.3206f,  0.5662f,  -0.3784f,
    0.4032f,  -1.7948f, -0.1554f, 0.3785f,  0.0534f,  -0.1653f, -0.0987f,
    -0.1005f, -0.6461f, 0.8035f,  -0.2405f, -0.1238f, -0.3576f, 0.0961f,
    -0.6401f, 0.2029f,  0.2359f,  0.4951f,  0.1921f,  -0.3745f, 0.3463f,
    0.2899f,  -0.1533f, -0.417f,  0.377f,   0.6998f,  0.594f,   0.5912f,
    -0.5571f, 0.0238f,  -0.2475f, 0.0307f,  -0.3875f, -0.7437f, 0.5144f,
    0.0057f,  0.7655f,  0.172f,   -2.5624f, -0.3066f, 0.3647f,  0.1727f,
    -0.0329f, -0.1893f, 0.0702f,  -1.2454f, 0.1496f,  -1.3825f, 0.4146f,
    -0.2668f, -0.1106f, 0.0362f,  -0.3189f, -0.7278f, -0.0894f, -0.2277f,
    -0.2394f, 1.0112f,  -0.2962f, 0.7776f,  0.2945f,  -0.2234f, 0.2764f,
    0.8011f,  -0.1744f, 0.1581f,  -1.7552f, -0.3848f, 0.5993f,  0.5268f,
    -0.0417f, 0.4733f,  -0.3401f, -0.145f,  0.7088f,  -0.1318f, 0.0426f,
    -0.5028f, 0.3832f,  -0.0118f, -0.4358f, 0.3749f,  -0.1203f, -0.5648f,
    -0.1972f, -0.8769f, -0.3675f, -0.2004f, -0.607f,  -0.1857f, 0.3468f,
    -0.3624f, 0.5358f,  -0.3701f, 0.1336f,  0.9545f,  0.1139f,  -0.1699f,
    0.3317f,  0.2891f,  0.2613f,  -0.0344f, -1.9498f, -2.0249f, -0.6005f,
    -0.6258f, -1.2349f, 0.328f,   -0.5434f, -0.7712f, -0.9057f, -0.1668f,
    -0.9905f, -0.0372f, -1.1638f, 0.1262f,  -0.5248f, -0.1538f, -0.3682f,
    0.3249f,  0.065f,   0.0511f,  -0.4607f, 0.2231f,  0.2822f,  0.1397f,
    -0.4938f, 0.3948f,  -0.4075f, -0.6411f, -0.0091f, -0.1333f, -0.5192f,
    -0.1661f, 0.3317f,  -0.0707f, 0.4806f,  0.3828f,  0.2229f,  0.616f,
    0.3371f,  0.1884f,  0.1381f,  -1.4915f, 0.2833f,  -0.1226f, -3.9189f,
    -0.4592f, -0.3435f, -0.6654f, -0.5056f, -0.8631f, 0.1536f,  -0.6427f,
    -0.0884f, -0.0471f, 0.1106f,  0.3821f,  -0.2392f, -0.4051f, 0.0891f,
    -0.6972f, -0.4699f, 0.0922f,  0.0806f,  -0.6774f, -0.0622f, -0.93f,
    0.1337f};

namespace OpenBabel {

// Forward declarations for helper functions
static uint64_t ComputeAtomType(OBMol& mol, OBAtom* atom, const OpenChemLib::AromaticityResolver* aromaticity);
static const std::unordered_map<uint64_t, double>& GetIncrementMap();

class OBAPI OBDescriptorCLogP : public OBDescriptor {
 public:
  OBDescriptorCLogP(const char* ID) : OBDescriptor(ID) {}
  const char* Description() override {
    return "Calculated logP via atom-type increments";
  }
  double Predict(OBBase* pOb, std::string* param = nullptr) override {
    OBMol* mol = dynamic_cast<OBMol*>(pOb);
    if (!mol) {
      return std::numeric_limits<double>::quiet_NaN();
    }

    // Create Java-compatible aromaticity resolver
    OpenChemLib::AromaticityResolver aromaticity(*mol);

    const auto& incrMap = GetIncrementMap();
    double clogp = 0.0;

    // Sum per-atom increments
    for (unsigned int i = 1; i <= mol->NumAtoms(); ++i) {
      OBAtom* atom = mol->GetAtom(i);
      if (atom->GetAtomicNum() == 1)  // skip hydrogens
        continue;
      const uint64_t type = ComputeAtomType(*mol, atom, &aromaticity);

      const auto it = incrMap.find(type);
      if (it != incrMap.end()) {
        clogp += it->second;
      }
      // Note: Missing atom types are silently ignored (same behavior as Java)
    }

    return clogp;
  }
};

// Initialize increment map from parallel arrays.
// Lambda-initialized static is constructed exactly once (C++11 thread-safe).
static const std::unordered_map<uint64_t, double>& GetIncrementMap() {
  static const std::unordered_map<uint64_t, double> map = []() {
    constexpr size_t n = sizeof(INCREMENTS) / sizeof(INCREMENTS[0]);
    std::unordered_map<uint64_t, double> m;
    m.reserve(n);
    for (size_t i = 0; i < n; ++i) {
      m[ATOM_TYPES[i]] = INCREMENTS[i];
    }
    return m;
  }();
  return map;
}

static uint64_t ComputeAtomType(OBMol& mol, OBAtom* atom, const OpenChemLib::AromaticityResolver* aromaticity) {
  const uint64_t mode = MODE_CLOGP;

  // Get neighbor count (equivalent to Java's mol.getConnAtoms(atom))
  int neighborCount = 0;
  FOR_NBORS_OF_ATOM(nbr, atom) {
    if (nbr->GetAtomicNum() != 1) {  // skip explicit hydrogen
      neighborCount++;
    }
  }

  // Initialize neighbor type array
  std::vector<uint64_t> neighborType(neighborCount, 0);
  int actualNeighborCount = 0;

  // Process each neighbor
  FOR_NBORS_OF_ATOM(nbr, atom) {
    if (nbr->GetAtomicNum() == 1)  // skip explicit hydrogen
      continue;

    uint64_t connType = 0;

    // Bond order (equivalent to Java's mol.getConnBondOrder(atom, i))
    if ((mode & cPropertiesConnBondOrder) != 0) {
      // Find the bond between current atom and neighbor
      OBBond* bond = mol.GetBond(atom, nbr);
      if (bond) {
        uint64_t connBondOrder = bond->GetBondOrder();
        // Handle aromatic bonds like Java does (use Java-compatible aromaticity)
        if (connBondOrder < 3 && aromaticity && aromaticity->isAromaticBond(bond)) {
          connBondOrder = 0;
        }

        connType |= (connBondOrder << CONN_SHIFT_BONDORDER);
      }
    }

    // Atomic number code (equivalent to Java's cSimpleAtomicNoCode[mol.getAtomicNo(connAtom)])
    if ((mode & cPropertiesConnAtomTypeSimple) != 0) {
      constexpr size_t kCodeSize = sizeof(cSimpleAtomicNoCode) / sizeof(cSimpleAtomicNoCode[0]);
      if (nbr->GetAtomicNum() >= kCodeSize || cSimpleAtomicNoCode[nbr->GetAtomicNum()] == -1) {
        continue;
      }
      connType += cSimpleAtomicNoCode[nbr->GetAtomicNum()];
    }

    // Aromatic flag (equivalent to Java's mol.isAromaticAtom(connAtom))
    if ((mode & cPropertiesConnAtomAromatic) != 0) {
      if (aromaticity && aromaticity->isAromaticAtom(nbr)) {
        connType |= CONN_FLAG_AROMATIC;
      }
    }

    // Insertion sort to maintain sorted order (equivalent to Java's insertion
    // sort)
    int index = 0;
    while (index < actualNeighborCount && connType < neighborType[index]) {
      index++;
    }

    // Shift elements to make room
    for (int j = actualNeighborCount; j > index; --j) {
      neighborType[j] = neighborType[j - 1];
    }

    neighborType[index] = connType;
    ++actualNeighborCount;
  }

  // Limit to maximum neighbors
  if (actualNeighborCount > MAX_NEIGHBORS) {
    actualNeighborCount = MAX_NEIGHBORS;
  }

  // Build atom type by packing neighbors
  uint64_t atomType = 0;
  for (int i = 0; i < actualNeighborCount; ++i) {
    atomType <<= CONN_FLAG_COUNT;
    atomType |= neighborType[i];
  }

  // Final shift for atom flags
  atomType <<= ATOM_FLAG_COUNT;

  // Add central atom's atomic number code
  constexpr size_t kCodeSize = sizeof(cAtomicNoCode) / sizeof(cAtomicNoCode[0]);
  if (atom->GetAtomicNum() >= kCodeSize || cAtomicNoCode[atom->GetAtomicNum()] == -1) {
    return 0;
  }
  atomType |= cAtomicNoCode[atom->GetAtomicNum()];

  // Small ring flag (equivalent to Java's mol.isSmallRingAtom(atom))
  if ((mode & cPropertiesAtomSmallRing) != 0) {
    // Check if atom is in small ring (7 or less)
    if (atom->IsInRingSize(3) || atom->IsInRingSize(4) || atom->IsInRingSize(5) || atom->IsInRingSize(6) || atom->IsInRingSize(7)) {
      atomType |= ATOM_FLAG_SMALLRING;
    }
  }

  // Charged flag (equivalent to Java's AtomFunctionAnalyzer.hasUnbalancedAtomCharge(mol, atom))
  if ((mode & cPropertiesAtomCharged) != 0) {
    // Java logic from AtomFunctionAnalyzer.hasUnbalancedAtomCharge:
    // Returns true if atom charge is unbalanced (not neutralized by neighboring charges)
    int atomCharge = atom->GetFormalCharge();
    bool hasUnbalancedCharge = false;

    if (atomCharge != 0) {
      hasUnbalancedCharge = true;  // Assume unbalanced by default

      // Calculate sum of charges from connected atoms
      int sumChargeConnectedAtoms = 0;
      FOR_NBORS_OF_ATOM(nbr, atom) {
        sumChargeConnectedAtoms += nbr->GetFormalCharge();
      }

      // Check if charge is balanced by neighbors
      if (abs(atomCharge) <= abs(sumChargeConnectedAtoms)) {
        // If signs are opposite, charge is balanced
        if ((atomCharge > 0) != (sumChargeConnectedAtoms > 0)) {
          hasUnbalancedCharge = false;
        }
      }

      if (hasUnbalancedCharge) {
        atomType |= ATOM_FLAG_CHARGED;
      }
    }

    // Ampholytic nitrogen check (Java: isBasicNitrogen + isAcidicOxygen)
    // ENABLED: Full Java-compatible implementation ported from AtomFunctionAnalyzer.java
    if (OpenChemLib::isBasicNitrogen(mol, atom, aromaticity)) {
      // Check if molecule has any acidic oxygen
      bool hasAcidic = false;
      for (unsigned int i = 1; i <= mol.NumAtoms(); i++) {
        OBAtom* testAtom = mol.GetAtom(i);
        if (OpenChemLib::isAcidicOxygen(mol, testAtom)) {
          hasAcidic = true;
          break;
        }
      }

      if (hasAcidic) {
        atomType |= ATOM_FLAG_AMPHOLYTIC;  // Set ampholytic flag
      }
    }
  }

  return atomType;
}

OBDescriptorCLogP clogP("clogP");
}  // namespace OpenBabel
