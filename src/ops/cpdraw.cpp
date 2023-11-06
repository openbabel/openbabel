/**********************************************************************
cpdraw.cpp - An OBOp for Cp (cyclopentadienyl) detection and svg depiction

Copyright (C) 2023 by Jesus N. M.

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
#include <openbabel/babelconfig.h>
#include <iostream>
#include<openbabel/op.h>
#include<openbabel/mol.h>
#include<openbabel/bond.h>
#include <openbabel/obiter.h>
#include <openbabel/elements.h>
#include <openbabel/oberror.h>
#include <openbabel/ring.h>
#include <openbabel/cpcomplex.h>
#include <openbabel/generic.h>
#include <openbabel/obmolecformat.h>
#include "../formats/smilesformat.cpp"

using namespace std;

namespace OpenBabel
{
    class OpCpDraw : public OBOp
    {
    public:
        OpCpDraw(const char* ID) : OBOp(ID, false) {};
        const char* Description() {
            return
                "Cp estructure-like detection and svg depiction.\n"
                "Works after 'gen2D' plugin generates 2D coordinates for the atoms"
                "Detects and modify carbons (in Cp ligands) coordinates and bonds"
                "Jesus N. M.";
        }
        virtual bool WorksWith(OBBase* pOb) const { return dynamic_cast<OBMol*>(pOb) != nullptr; }
        virtual bool Do(OBBase* pOb, const char* OptionText = nullptr, OpMap* pOptions = nullptr, OBConversion* pConv = nullptr);
        //! \return If @p bond is likely to be a cp-bond like
        bool isCpBond(OBBond* bond); //Metodo que comprueba si a priori podría ser un enlace tipo Cp
        //! Finds the ring of which the carbon with idx @p carbonIdx is a part of, among the rings of @p rlist (obtained from a SSSR perspective), and stores it in @p result.
        //! \returns whether it was found or not
        bool FindRingWithCarbon(vector<OBRing*>& rlist, int carbonIdx, OBRing*& result);
        //! Canonize the input SMILES and identify blocks
        void CanonizeOgm(OBMol* mol, OBConversion* pConv); 


    };

    /////////////////////////////////////////////////////////////////
    OpCpDraw theOpCpDraw("CpDraw"); //Global instance

    /////////////////////////////////////////////////////////////////
    bool OpCpDraw::Do(OBBase* pOb, const char* OptionText, OpMap* pOptions, OBConversion* pConv)
    {
        OBMol* pmol = dynamic_cast<OBMol*>(pOb);
        if (!pmol)
            return false;

        //Sacamos el SMILES Canonico e identificamos los bloques (no se modifica nada internamente, es solamente para identificar bloques)
        CanonizeOgm(pmol, pConv);
        //cout << "CanSmiles: " << pmol->GetCanSmiles() << "\n";
        


        /* ------------------- Algoritmo de deteccion de Cp ----------------------- */
        vector<int> test(pmol->NumBonds(),0);
        vector<pair<int,int>> cpBonds(pmol->NumBonds(), std::make_pair(-1,-1)); //Por defecto lo ponemos a -1. con esto comprobamos tb que si el valor es -1, no es un cpbond
        OBAtom* begin;
        OBAtom* end;
        int contador = 0;
        int atomNMetal = -1; 
        bool foundOne = false;
        bool cpTest = false;
        FOR_BONDS_OF_MOL(b,pmol) {
            begin = b->GetBeginAtom();
            end = b->GetEndAtom();
            if (isCpBond(&*b)) {
                if (!foundOne) 
                    foundOne = true;
                cpBonds[b->GetIdx()] = make_pair(begin->GetIdx(), end->GetIdx());
            }

            //Debug
            contador++;
            
        }
        

        //Si ha detectado aunque sea 1 posible cpBond seguimos con el algoritmo. Si no, acabamos.
        if (foundOne) {

            std::vector<std::vector<pair<int, int>>> individualCpBonds;
            std::vector<pair<int, int>> temp_cp;
            int usedCarbons = 0;
            BranchBlock* newBranch = nullptr, * currentBranch = nullptr;
            unsigned int current_carbon = 0; //current carbon
            unsigned int current_metal = 0, new_metal = 0;
            for (int i = 0; i < cpBonds.size(); i++) {
                if (cpBonds[i].first != -1) {
                    //Sacamos el idx del carbono 
                    if (pmol->GetBond(i)->GetBeginAtom()->IsOgmMetal()) {
                        current_carbon = pmol->GetBond(i)->GetEndAtom()->GetIdx();
                        new_metal = pmol->GetBond(i)->GetBeginAtom()->GetIdx();
                    }
                    else {
                        current_carbon = pmol->GetBond(i)->GetBeginAtom()->GetIdx();
                        new_metal = pmol->GetBond(i)->GetEndAtom()->GetIdx();
                    }   
                    newBranch = pmol->FindBranch(current_carbon);
                    if (newBranch) {
                        if (!current_metal) 
                            current_metal = new_metal;

                        if (!currentBranch) //Si es el inicio del bloque que tratamos, lo asignamos a current 
                            currentBranch = newBranch;
                        if (newBranch == currentBranch) { //Si es la misma branch, es que seguimos identificando el mismo cp que la iteracion anterior
                            //temp_cp.push_back(cpBonds[i]);
                            usedCarbons++;
                        }
                        else {
                            //Ha encontrado un bloque pero es distinto que el anterior, por lo que ya hemos acabado con el cp anterior 
                            //(puede que lo haya encontrado entero, en cuyo caso lo metemos definitaivamente a individualCps; 
                            //o puede que no todos los carbonos del bloque se hayan usado, por lo que no lo metemos)

                            if (usedCarbons == currentBranch->Size() && usedCarbons >= 3) { //todo el bloque usado, y minimo 3 carbonos (para que haga un poligono cerrado, de momento los compuestos con deticity no los contemplo)
                                for (int j = 0; j < currentBranch->Size(); j++) {
                                    temp_cp.push_back(make_pair(currentBranch->GetAtomIdx(j),current_metal));
                                }
                                individualCpBonds.push_back(temp_cp);
                                if (current_metal != new_metal)
                                    current_metal = new_metal;
                            }

                            //Reseteamos el contador para el siguiente bloque
                            usedCarbons = 0;
                            currentBranch = nullptr;
                            temp_cp.clear();
                            i--; //Volvemos a la interacion anterior para que no se salte el carbono que ha hecho cambiar de bloque
                        }
                    }
                    else {//No me ha encontrado una branch en la que esté el carbono

                    }
                }
                else {
                }
                //Si es el ultimo elemento del vector, comprobamos si nos dejamos algun bloque sin insertar
                if (i == cpBonds.size() - 1)
                    if (currentBranch && (usedCarbons == currentBranch->Size())) { //todo el bloque usado
                        for (int j = 0; j < currentBranch->Size(); j++) {
                            temp_cp.push_back(make_pair(currentBranch->GetAtomIdx(j), current_metal));
                        }
                        individualCpBonds.push_back(temp_cp);
                    }
            }



            //Algunas comprobaciones de seguridad
            // -    Comprobar que el metal tiene al menos el tantos bonds como carbonos detectamos en los ciclos
            // -    Que no me detecte carbonos aislados (que quizas sean parte de un anillo, pero el anillo entero no sea Cp) 
            //      Cojo un carbono de los que tengo marcado un bondCp, y compruebo si todos los demas carbonos del mismo anillo son tambien parte de un bondCp
            //Esto es para que si despues del bucle, tenemos bonds que tienen pinta de Cp, pero luego vemos que el metal en cuestion no tiene enlaces suficientes, algo ha ido mal, o realmente no eran Cp bonds
            OBAtom* atom;
            int carbonIdx = 0;
            vector<OBRing*> rlist;
            vector<OBRing*>::iterator itr;
            OBAtomIterator it;
            OBRing* ringCarbon = nullptr;
            vector<int> rpath;
            bool goodInsert = true;


            /*------------------ Para cada uno de los cp individuales detectados -----------------*/
            int ncps = individualCpBonds.size();
            for (int icp = 0; icp < individualCpBonds.size(); icp++) {
                goodInsert = true;
                CpComplex* cp = new CpComplex;
                unsigned int cpMetalIdx;
                std::vector<pair<int, int>>cpIndvBonds = individualCpBonds[icp];

                (pmol->GetAtom(cpIndvBonds[0].first)->IsOgmMetal()) ? //no sabemos si el metal es el first o el second
                    cpMetalIdx = pmol->GetAtom(cpIndvBonds[0].first)->GetIdx() :
                    cpMetalIdx = pmol->GetAtom(cpIndvBonds[0].second)->GetIdx();
                cp->SetMetalIdx(cpMetalIdx);

                int idxInsert;
                for (int i = 0; i < cpIndvBonds.size(); i++) {
                    (pmol->GetAtom(cpIndvBonds[i].first)->IsOgmMetal()) ?
                        idxInsert = pmol->GetAtom(cpIndvBonds[i].second)->GetIdx() :
                        idxInsert = pmol->GetAtom(cpIndvBonds[i].first)->GetIdx();
                    cp->AddIdxCarbon(idxInsert);
                    atom = pmol->GetAtom(idxInsert);
                    cp->AddCpAtom(atom);
                }

                //Descartar carbonos aislados para no modificar anillos erroneos
                //Sacamos la lista de anillos totales de la molecula

                // Este check no me hace falta ya realmente despues de usar los bloques
                if (!pmol->HasSSSRPerceived())
                    pmol->FindSSSR();

                rlist = pmol->GetSSSR();

                //De momento, solo puedo comprobar esto con 1 anillo. Tenemos el mismo problema de separar los Cp. En cpBonds tengo todos los posibles enlaces: no se cuando salgo de comprobar 1 y me meto en otro 
                //Esto cuando tenga la separacion de Cps, tendre que cambiarlo
                map<int, bool> atomVisited; //map<key:idx_carbon, bool:checked_good>
                for (auto id : cp->GetIdxCarbons()) //Relleno el map con los mismos idx de los carbonos del cp
                    atomVisited.emplace(id, 0);
                vector<unsigned int> cppath = cp->GetIdxCarbons();
                for (int i = 0; i < cp->GetCarbonsSize() && goodInsert; i++) {
                    carbonIdx = cp->GetCarbonIdx(i);

                    if (atomVisited.at(carbonIdx) == 0) {
                        /*Esto deberia hacer que devolviera una lista de obbrings, por si el carbono está en varios anillos distintos, y luego comprobar todos los anillos.
                        Simplemente viendo si todos los carbonos del _path están en el Cp, seria valido. Si el Cp tiene menos, el Cp está mal detectado y no se insertará en la molecula. Lo mismo si tiene más.
                        Esto no funciona del todo para todas las moleculas. Mol3 da problemas*/
                        if (FindRingWithCarbon(rlist, carbonIdx, ringCarbon)) {
                            rpath = ringCarbon->_path;

                            for (int j = 0; j < rpath.size(); j++) {
                                if (pmol->GetAtom(rpath[j])->GetAtomicNum() == 6) { //Si es un carbono (lo compruebo, porque en los ciclos de cp, normalmente se dividen en ciclos mas pequeños de 3 (C-C-M), y para saltarme el metal)
                                    auto it = find(cppath.begin(), cppath.end(), rpath[j]);
                                    if (it != cppath.end()) {
                                        //Me ha encontrado el idx del rpath[i] en el cp. Aniado ese atomo como visto
                                        if (atomVisited[*it] == 0) atomVisited.at(*it) = 1;
                                    }
                                    else { //No me ha encontrado el idx de un rpath donde el carbono que estabamos comprobando si estaba, pero rpth[i] no está, por lo que no son el mismo anillo. Abortamos
                                        goodInsert = false;
                                        break;
                                    }
                                }
                            }
                            //Si los encuentra todos, ese anillo contiene todos los carbonos 
                        }

                    }

                }

                //Si es valido lo podemos insertar en la molecula.
                if (goodInsert) {
                    cp->SetParent(pmol);
                    pmol->AddCpComplex(*cp); //Posible memory leak
                    //Sabiendo que es valido, marcamos todos los atomos del cp con un flag especial (servirá mas tarde en el dibujado para detectar estos atomos en concreto)
                    for (atom = cp->BeginAtomCp(it); atom; atom = cp->NextAtomCp(it)) {
                        atom->SetInCp();
                    }
                }

                //cout << "\n\n";


                //Usamos el Cp (si existe) que hemos aniadido a la molecula (es distinto del que hemos creado antes, ya que internamnete al añadirlo, se copia en uno nuevo)
                cp = pmol->GetCpComplex(icp + 1);
                if (cp) {
                    double offsetX = 0.0, offsetY = 0.0;
                    OBAtom* atomMetal;
                    atomMetal = pmol->GetAtom(cp->GetMetalIdx());


                    /* ------------------ Proceso de creacion del pentagono en perspectiva ---------------------*/
                    /* 1. Creo el dummy y lo coloco en la posicion que toque segun el numero de Cps que haya detectado (el primero tiene que ir arriba)
                    *  2. El radio, cojo uno estandar y que me lo haga siempre del mismo tamaño mejor. Para el cp_centroid en la formula de mover los carbonos, usaria la posicion del dummy
                    *  3. Pongo los carbonos en pentagono regular, y luego en perspectiva
                    *  4. Creo el circulo del cp
                    */

                    //Calculamos la posicion en la que deberia ir el dummy en base a los cp encontrados
                    vector3 dummy_coords;
                    double bondLengthSum = 0.0;
                    OBBondIterator j;
                    for (OBBond* bond = pmol->BeginBond(j); bond; bond = pmol->NextBond(j))
                        bondLengthSum += bond->GetLength();
                    const double averageBondLength = bondLengthSum / pmol->NumBonds();
                    double dummy_radius = averageBondLength * 1.2; //medida arbitraria que me parece buena
                    double alphaD, _xD, _yD;
                    alphaD = 2 * M_PI * icp / ncps;
                    _xD = atomMetal->GetX() + dummy_radius * cos(alphaD + M_PI/2); //Para usar un origen distinto de (0,0), sumamos (x,y). Para empezar la division de angulos desde arriba, sumamos PI/2
                    _yD = atomMetal->GetY() + dummy_radius * sin(alphaD + M_PI/2);
                    dummy_coords.Set(_xD, _yD, 0.0);


                    // Create bond dummy atom to center metal
                    OBAtom* dummy;
                    dummy = pmol->NewAtom();
                    dummy->SetAtomicNum(0);
                    dummy->SetVector(dummy_coords);

                    cp->SetDummyIdx(dummy->GetIdx());
                    pmol->AddBond(atomMetal->GetIdx(), dummy->GetIdx(), 1);
                    //cout << "Dummy bond created \n"; //Debug

                    cp->SetCentroid(dummy_coords);
                    

                    vector3 centroidCp;
                    centroidCp = cp->GetCentroid();


                    //Movemos los carbonos en forma de poligono regular en base al centro calculado de centroidCp con un radio de centroide-carbono
                    int n_lados = cp->GetCarbonsSize();
                    double radioPentagono = 1.0;
                    double alpha, _x, _y, _z;
                    for (int i = 0; i < n_lados; i++) {
                        alpha = 2 * M_PI * i / n_lados;
                        _x = centroidCp.x() + radioPentagono * cos(alpha); //Para usar un origen distinto de (0,0), sumamos (x,y)
                        _y = centroidCp.y() + radioPentagono * sin(alpha);
                        pmol->GetAtom(cp->GetCarbonIdx(i))->SetVector(_x, _y, 0.0);
                    }


                    //Mover de nuevo los carbonos para darle perspectiva al pentagono
                    //Aplicamos una rotacion en un eje paralelo al Eje de coordenadas X. Dicho eje quiero que pase justamente por el centro del poligono para que la rotacion del poligono sea sobre su centro (centroide) (por lo que utilizamos las coordenadas del centroide en la formula)
                    OBAtom* _atom;
                    double alphaR = deg2rads(PERSPECTIVE_DEG);
                    for (int i = 0; i < n_lados; i++) {
                        _atom = pmol->GetAtom(cp->GetCarbonIdx(i));
                        _y = _atom->GetY() * cos(alphaR) - _atom->GetZ() * sin(alphaR) + (centroidCp.GetY() * (1.0 - cos(alphaR)) + _atom->GetZ() * sin(alphaR));
                        _z = _atom->GetY() * sin(alphaR) + _atom->GetZ() * cos(alphaR) + (centroidCp.GetZ() * (1.0 - cos(alphaR)) - _atom->GetY() * sin(alphaR));
                        _atom->SetVector(_atom->GetX(), _y, _z);                            // Al ser rotacion en eje X, la x queda igual
                    }


                    //Creacion del circulo central y su rotacion
                    vector3 _coord;
                    int circle_sides = 30;
                    radioPentagono = radioPentagono * 0.6;
                    double _cx = 0.0, _cy = 0.0, _cz = 0.0;
                    double __cy, __cz;
                    for (int i = 0; i < circle_sides; i++) {
                        alpha = 2 * M_PI * i / circle_sides;
                        _cx = centroidCp.x() + radioPentagono * cos(alpha); //Para usar un origen distinto de (0,0), sumamos (x,y)
                        _cy = centroidCp.y() + radioPentagono * sin(alpha);

                        //lo rotamos directamente antes de meterlo
                        __cy = _cy * cos(alphaR) - _cz * sin(alphaR) + (centroidCp.GetY() * (1.0 - cos(alphaR)) + _cz * sin(alphaR));
                        __cz = _cy * sin(alphaR) + _cz * cos(alphaR) + (centroidCp.GetZ() * (1.0 - cos(alphaR)) - _cz * sin(alphaR));
                        _coord.Set(_cx, __cy, __cz);
                        cp->AddCircleCoord(_coord);
                    } 
                }
            }

             /*------ Quitar los bonds M - C / C - M -----*/
            OBBond* bondToDelete;
            OBAtom* metal, * carbono;

            for (int i = 1; i <= pmol->GetCpSize(); i++) {
                CpComplex* cp = pmol->GetCpComplex(i);

                for (int i = 0; i < cp->GetCarbonsSize(); i++) {
                    metal = pmol->GetAtom(cp->GetMetalIdx());
                    carbono = pmol->GetAtom(cp->GetCarbonIdx(i));
                    bondToDelete = pmol->GetBond(metal, carbono); //Me da igual quien sea el begin o el end, el metodo saca el bond entre los atomos especificados, pero no es estricto en ese sentido
                    pmol->DeleteBond(bondToDelete);
                }


               /* cout << "Cp Bonds deleted. New molecule bond list: \n"; //Debug
                FOR_BONDS_OF_MOL(b, pmol) {
                    begin = b->GetBeginAtom();
                    end = b->GetEndAtom();
                    cout << "[" << begin->GetIdx() << "]" << OBElements::GetSymbol(begin->GetAtomicNum()) << "-"
                        << "[" << end->GetIdx() << "]" << OBElements::GetSymbol(end->GetAtomicNum()) << "\n";
                }*/
            }
        }

        return true;
    }    


    bool OpCpDraw::isCpBond(OBBond* bond)
    {
        if (bond->GetBondOrder() != 1)
            return false;

        OBAtom* M = nullptr, * C = nullptr; //M for metal and C for carbon

        OBAtom* begin = bond->GetBeginAtom();
        if (begin->IsOgmMetal()) //Si es un Ogm metal
            M = begin;
        if (begin->GetAtomicNum() == 6)
            C = begin;

        OBAtom* end = bond->GetEndAtom();
        if (end->IsOgmMetal())
            M = end;
        if (end->GetAtomicNum() == 6)
            C = end;

        if (!M || !C)
            return false;

        return C->IsInRing();
    }

    bool OpCpDraw::FindRingWithCarbon(vector<OBRing*> &rlist, int carbonIdx, OBRing* &result) //Esto puede resultar en error, porque devuelve el 1º anillo que encuentre que cumpla la condicion. Pero el mismo carbono puede estar en varios anillos, y no ser el 1º que encuentre el que queremos comprobar
    {
        vector<OBRing*>::iterator itr;
        for (itr = rlist.begin(); itr != rlist.end(); ++itr)
            if ((*itr)->IsInRing(carbonIdx)){
                result = (*itr);
                return (true);
            }
        return (false);
    }

    void OpCpDraw::CanonizeOgm(OBMol* pmol, OBConversion* pConv)
    {
        OBBitVec fragatoms(pmol->NumAtoms());

        OBPairData* dp = (OBPairData*)pmol->GetData("SMILES_Fragment");
        const char* ppF = pConv->IsOption("F");
        if (dp) {
            fragatoms.FromString(dp->GetValue(), pmol->NumAtoms());
        }
        else if (ppF) { // Use info from option "F"
            fragatoms.FromString(ppF, pmol->NumAtoms());
        }
        // If no "SMILES_Fragment" data, fill the entire OBBitVec
        // with 1's so that the SMILES will be for the whole molecule.
        else {
            FOR_ATOMS_OF_MOL(a, *pmol)
            {
                fragatoms.SetBitOn(a->GetIdx());
            }
        }

        std::string buffer;
        buffer.reserve(1000);

        if (pmol->NumAtoms() > 0 || pmol->IsReaction()) {
            CreateCansmiString(*pmol, buffer, fragatoms, pConv);
        }
    }

}//namespace
