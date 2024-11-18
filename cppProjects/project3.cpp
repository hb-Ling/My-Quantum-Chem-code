/* UTF-8 */
/* This cpp file is used for practice.*/

/* This code use direct SCF, so it will be difficult to reach convergence limit
For better code implementation, you can check project3_DIIS.cpp, which supports DIIS algorithm and accelerates convergence

*/

/* This code will:
1. Read and process all electronic integrals calculated by other programs.
2. Orthogonalization of the Basis Set
3. Create inital guess density matrix
4. Perform SCF calculations and obtain single point energy

*/

#include <vector>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <cmath>
#include <algorithm>
#include <format>
#include <stdlib.h>
#include <F:/cpp/eigen-3.4.0/Eigen/Dense>
#include <F:/cpp/eigen-3.4.0/Eigen/Eigenvalues>
#define M_PI 3.141592654
#define BOHR_TO_A 0.52918
#define RAD_TO_ANG 180 / M_PI
#define AMU_TO_G 1.6605e-24
#define H_PLANCK 6.626e-34
#define HZ_TO_CM
#define _LIMIT 10e-3
#define _C 2.99792458e8 // light speed
#define HARTREE_TO_J 4.36e-18
#define _judgeBond 1.4 // use to judge whether two atoms can form a bond.
using namespace std;

const vector<string> ATOMS = {"H", "He", "Li", "Be", "B", "C", "N", "O", "F", "Ne", "Na", "Mg", "Al", "Si", "P", "S", "Cl"};
const double Atom_mass[] = {1.008, 4.003, 6.94, 9.01, 10.81, 12.01, 14.01, 15.99, 19.00};          // in amu
const double Radius[] = {0.586, 0.529, 2.419, 1.814, 1.587, 1.436, 1.341, 1.247, 1.209, 2.0, 2.0}; // covalent radius in Bohr
vector<double> cross(vector<double> &v1, vector<double> &v2)                                       // return v1 X v2
{
    assert(v1.size() == 3 && v2.size() == 3);
    vector<double> ans(3, 0);
    ans[0] = v1[1] * v2[2] - v1[2] * v2[1];
    ans[1] = v1[2] * v2[0] - v1[0] * v2[2];
    ans[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return ans;
}
Eigen::MatrixXf sort_columns(Eigen::MatrixXf init_matrix, vector<int> &eigen_order)
{
    Eigen::MatrixXf sorted = init_matrix;
    assert(eigen_order.size() == init_matrix.cols());
    for (int col = 0; col < init_matrix.cols(); col++)
    {
        for (int row = 0; row < init_matrix.rows(); row++)
            sorted(row, eigen_order[col]) = init_matrix(row, col);
    }
    return sorted;
}
int index_map(int miu, int v)
{
    if (miu >= v)
        return miu * (miu + 1) / 2 + v;
    else
        return v * (v + 1) / 2 + miu;
}

class Molecular
{
private:
    int atoms_num = 0;     // molecular contains how many atoms --> at least one
    int orbital_num = 0;   // number of basis function
    int *atom_list = NULL; // array of Z
    int electron = 0;
    double E_elec = 0;          // obtain from SCF
    double V_NN = 0;            // Nuclear Repulsion energy, fixed under B-O appromixation
    double (*Coords)[3] = NULL; // Z_atom (x,y,z)     ; a 2-D array with size [atoms_num][3]
    double (*Forces)[3] = NULL; // Z_atom (fx,fy,fz)  ; a 2-D array with size [atoms_num][3]
    const char *basis = NULL;

public:
    Molecular(void)
    {
    }
    ~Molecular() // Destructor
    {
        atoms_num = 0;
        V_NN = 0;
        orbital_num = 0;
        electron = 0;
        E_elec = 0;
        V_NN = 0;
        basis = NULL;
        delete Coords;
        delete atom_list;
        delete Forces;
    }

    void build_Coords(const char *filename)
    {
        ifstream fin;
        fin.open(filename, ios::in);
        if (!fin.is_open())
        {
            cout << "Error: Could not open the xyz file! Please ensure the file path is vaild." << endl;
            return;
        }
        char buf[512] = {0};
        bool start = true;
        int count = 0;
        while (fin.getline(buf, sizeof(buf)))
        {
            if (start == true)
            {
                int length;
                sscanf(buf, "%d", &length);
                atoms_num = length;
                start = false;
                Coords = new double[atoms_num][3];
                atom_list = new int[atoms_num];
                continue;
            }
            sscanf(buf, "%d%lf%lf%lf", &atom_list[count], &Coords[count][0], &Coords[count][1], &Coords[count][2]);
            count += 1;
            assert(count <= atoms_num);
        }
        fin.close();
        int e_num = 0;
        for (int i = 0; i < atoms_num; i++)
            e_num += atom_list[i];
        electron = e_num;
        return;
    }
    void print_Coords()
    {
        for (int i = 0; i < atoms_num; i++)
            printf("%d%12.6lf%12.6lf%12.6lf\n", atom_list[i], Coords[i][0], Coords[i][1], Coords[i][2]);
    }
    void write_Coords(ofstream &fout) // reference
    {
        fout << "Coordinates of Molecule (in Angstrom)\n"
             << endl;
        fout << format("{:^8s}{:^8s}{:>12s}{:>12s}{:>12s}", "Index", "Atom", "Coord x", "Coord y", "Coord z") << endl;
        for (int i = 0; i < atoms_num; i++)
            fout << format("{:^8d}{:^8s}{:>12.6f}{:>12.6f}{:>12.6f}",
                           i + 1, ATOMS[atom_list[i] - 1], Coords[i][0] * BOHR_TO_A, Coords[i][1] * BOHR_TO_A, Coords[i][2] * BOHR_TO_A)
                 << endl;
        fout << endl;
    }
    int get_num()
    {
        return atoms_num;
    }
    int get_atlist(int index)
    {
        return atom_list[index];
    }
    double get_Ele()
    {
        return E_elec;
    }
    double get_VNN()
    {
        return V_NN;
    }
    void read_VNN(const char *filename)
    {
        ifstream enu;
        enu.open(filename, ios::in);
        if (!enu.is_open())
        {
            cout << "Error: Could not open the enuc file! Please ensure the file path is vaild." << endl;
            return;
        }
        char buf[32] = {0};
        enu.getline(buf, sizeof(buf));
        sscanf(buf, "%lf", &V_NN);
        enu.close();
    }
    void get_orbital_num(string &Basis) // return number of basis function (contract Gaussian function)
    {                                   // Cartisean GF --> 6 functions for d (xx,yy,zz,yx,yz,zx)
        int sum = 0;
        if (Basis == "STO-3G")
        {
            for (int i = 0; i < atoms_num; i++)
            {
                if (atom_list[i] == 1)
                    sum += 1; // 1s
                if (atom_list[i] >= 3 && atom_list[i] <= 9)
                    sum += 5; // 1s + 2s + 3*2p
            }
        }
        if (Basis == "DZ") // Double Zeta
        {
            for (int i = 0; i < atoms_num; i++)
            {
                if (atom_list[i] == 1)
                    sum += 2; // 2 function for 1s
                if (atom_list[i] >= 3 && atom_list[i] <= 9)
                    sum += 10; // 2 for 1s, 2 for 2s, 6 for 2p
            }
        }
        if (Basis == "DZP") // Double Zeta with Polarization
        {
            for (int i = 0; i < atoms_num; i++)
            {
                if (atom_list[i] == 1)
                {
                    sum += 2; // 2 functions for 1s
                    sum += 3; // p ploarization for H
                }
                if (atom_list[i] >= 3 && atom_list[i] <= 9)
                {
                    sum += 10; // 2x(1s + 2s + 3*2p)
                    sum += 6;  // 6d ploarization for p orbital
                }
            }
        }
        orbital_num = sum;
        basis = Basis.c_str();
        return;
    }
    int get_orbital_num(void)
    {
        return orbital_num;
    }
    Eigen::MatrixXf read_1eint(const char *Workdir, const char *Basis) // return H_core = T + V
    {
        string workdir = Workdir;
        string basis = Basis;
        assert(orbital_num >= 1);
        Eigen::MatrixXf H_core(orbital_num, orbital_num);

        ifstream T, V;
        string T_file = workdir + basis + "/t.dat";
        string V_file = workdir + basis + "/v.dat";
        // read T
        T.open(T_file.c_str(), ios::in);
        if (!T.is_open())
        {
            cout << "Can't open T integral file." << endl;
            return {};
        }
        char buf[512] = {0};
        while (T.getline(buf, sizeof(buf)))
        {
            int i, j;
            double inte;
            sscanf(buf, "%d%d%lf", &i, &j, &inte);
            H_core(i - 1, j - 1) = inte;
        }
        T.close();
        // read V
        V.open(V_file.c_str(), ios::in);
        if (!V.is_open())
        {
            cout << "Can't open V integral file." << endl;
            return {};
        }
        while (V.getline(buf, sizeof(buf)))
        {
            int i, j;
            double inte;
            sscanf(buf, "%d%d%lf", &i, &j, &inte);
            H_core(i - 1, j - 1) += inte;
        }
        V.close();
        for (int i = 0; i < orbital_num; i++)
        {
            for (int j = i + 1; j < orbital_num; j++) // j > i
                H_core(i, j) = H_core(j, i);
        }
        return H_core;
    }
    Eigen::MatrixXf read_overlap_int(const char *Workdir, const char *Basis)
    {
        string workdir = Workdir;
        string basis = Basis;
        assert(orbital_num >= 1);
        string S_file = workdir + basis + "/s.dat";
        Eigen::MatrixXf S(orbital_num, orbital_num);
        ifstream fin;
        fin.open(S_file.c_str(), ios::in);
        if (!fin.is_open())
        {
            cout << "Can't open overlap inte file!";
            return {};
        }
        char buf[512] = {0};
        while (fin.getline(buf, sizeof(buf)))
        {
            int i, j;
            double inte;
            sscanf(buf, "%d%d%lf", &i, &j, &inte);
            S(i - 1, j - 1) = inte;
            S(j - 1, i - 1) = inte;
        }
        fin.close();
        return S;
    }
    vector<double> read_rep_int(const char *Workdir, const char *Basis) // µνλσ ≡ µν(µν+1)/2+λσ.
    {
        string workdir = Workdir;
        string basis = Basis;
        vector<double> rep;
        rep.resize(index_map(index_map(orbital_num, orbital_num), index_map(orbital_num, orbital_num)) + 10);
        ifstream fin;
        string r_file = workdir + basis + "/eri.dat";
        fin.open(r_file.c_str(), ios::in);
        if (!fin.is_open())
        {
            cout << "Can't open repulsion inte file." << endl;
            return {};
        }
        char buf[512] = {0};
        while (fin.getline(buf, sizeof(buf)))
        {
            int miu, v, lambda, sigma;
            double inte;
            sscanf(buf, "%d%d%d%d%lf", &miu, &v, &lambda, &sigma, &inte);
            rep[index_map(index_map(miu, v), index_map(lambda, sigma))] = inte;
            /*
            rep[index_map(index_map(v, miu), index_map(lambda, sigma))] = inte;
            rep[index_map(index_map(miu, v), index_map(sigma, lambda))] = inte;
            rep[index_map(index_map(v, miu), index_map(sigma, lambda))] = inte;
            rep[index_map(index_map(lambda, sigma), index_map(miu, v))] = inte;
            rep[index_map(index_map(lambda, sigma), index_map(v, miu))] = inte;
            rep[index_map(index_map(sigma, lambda), index_map(miu, v))] = inte;
            rep[index_map(index_map(sigma, lambda), index_map(v, miu))] = inte;
            */
        }
        fin.close();
        return rep;
    }
    Eigen::MatrixXf ortho_overlap(Eigen::MatrixXf &init_S) // return S^-0.5; Lowdin Diag
    {
        Eigen::MatrixXf ortho_S(orbital_num, orbital_num);
        Eigen::EigenSolver<Eigen::MatrixXf> es(init_S);
        Eigen::MatrixXf eigenvectors(orbital_num, orbital_num), eigenvalues(orbital_num, orbital_num);
        eigenvectors = es.eigenvectors().real();
        for (int i = 0; i < orbital_num; i++)
        {
            for (int j = 0; j < orbital_num; j++)
            {
                if (i == j)
                {
                    assert(real(es.eigenvalues()[i]) > 0);
                    eigenvalues(i, j) = 1 / sqrt(real(es.eigenvalues()[i]));
                }
                else
                    eigenvalues(i, j) = 0;
            }
        }
        ortho_S = eigenvectors * eigenvalues * eigenvectors.transpose();
        return ortho_S;
    }
    Eigen::MatrixXf initial_guess_D(Eigen::MatrixXf &S_sym, Eigen::MatrixXf &H_core, bool is_RHF = true) // return inital Density Matrix
    {
        Eigen::MatrixXf F_init(orbital_num, orbital_num), C_init(orbital_num, orbital_num), D_init(orbital_num, orbital_num);
        F_init = S_sym.transpose() * H_core * S_sym;
        Eigen::EigenSolver<Eigen::MatrixXf> es(F_init);
        vector<int> order(orbital_num, 0);
        for (int i = 0; i < orbital_num; i++)
        {
            int count = 0;
            for (int j = 0; j < orbital_num; j++)
            {
                if (i == j)
                    continue;
                if (es.eigenvalues()[i].real() - es.eigenvalues()[j].real() > 0)
                    count++;
            }
            order[i] = count;
        }
        Eigen::MatrixXf sorted_eigenvector = sort_columns(es.eigenvectors().real(), order);
        // sort eigenvectors based on the order of eigenvalue

        C_init = S_sym * sorted_eigenvector; // C_init = C_0, the coffeciant matrix

        // find index corresponding to an occupied MO

        int occupied_orbital = 0;
        if (is_RHF == true)
            occupied_orbital = electron / 2;
        else
            cout << "Sorry, program can't support UHF calculations right now." << endl;

        for (int miu = 0; miu < orbital_num; miu++)
        {
            for (int v = 0; v < orbital_num; v++)
            {
                D_init(miu, v) = 0;
                for (int i = 0; i < occupied_orbital; i++)
                    D_init(miu, v) += C_init(miu, i) * C_init(v, i);
            }
        }

        double E_init = 0;
        for (int i = 0; i < orbital_num; i++)
        {
            for (int j = 0; j < orbital_num; j++)
            {
                E_init += D_init(i, j) * (H_core(i, j) + F_init(i, j));
            }
        }
        E_elec = E_init;
        /*
                for (int i = 0; i < orbital_num; i++)
                {
                    for (int j = 0; j < orbital_num; j++)
                        cout << format("{:>14.8f}", D_init(i, j));
                    cout << endl;
                }
        */
        return D_init;
    }
    Eigen::MatrixXf build_Fock_matrix(Eigen::MatrixXf &H_core, Eigen::MatrixXf &D, vector<double> &rep)
    {
        Eigen::MatrixXf Fock = H_core;
        for (int miu = 1; miu <= orbital_num; miu++)
        {
            for (int v = 1; v <= orbital_num; v++)
            {
                for (int lambda = 1; lambda <= orbital_num; lambda++)
                {
                    for (int sigma = 1; sigma <= orbital_num; sigma++)
                    {
                        int index1 = index_map(index_map(miu, v), index_map(lambda, sigma));
                        int index2 = index_map(index_map(miu, lambda), index_map(v, sigma));
                        Fock(miu - 1, v - 1) += D(lambda - 1, sigma - 1) * (2 * rep[index1] - rep[index2]);
                    }
                }
            }
        }
        return Fock;
    }
    Eigen::MatrixXf build_Density_matrix(Eigen::MatrixXf &Fock, Eigen::MatrixXf &S_sym, bool is_RHF = true)
    {
        Eigen::MatrixXf Density(orbital_num, orbital_num);
        Eigen::MatrixXf F_sym = S_sym.transpose() * Fock * S_sym;
        Eigen::EigenSolver<Eigen::MatrixXf> es(F_sym);
        Eigen::MatrixXf eigenvectors(orbital_num, orbital_num);
        vector<int> order(orbital_num, 0);
        for (int i = 0; i < orbital_num; i++)
        {
            int count = 0;
            for (int j = 0; j < orbital_num; j++)
            {
                if (i == j)
                    continue;
                if (es.eigenvalues()[i].real() - es.eigenvalues()[j].real() > 0)
                    count++;
            }
            order[i] = count;
        }
        Eigen::MatrixXf sorted_eigenvector = sort_columns(es.eigenvectors().real(), order);
        Eigen::MatrixXf Coef_matrix = S_sym * sorted_eigenvector;
        int occupied_orbital = 0;
        if (is_RHF == true)
            occupied_orbital = electron / 2;

        for (int miu = 0; miu < orbital_num; miu++)
        {
            for (int v = 0; v < orbital_num; v++)
            {
                Density(miu, v) = 0;
                for (int m = 0; m < occupied_orbital; m++)
                    Density(miu, v) += Coef_matrix(miu, m) * Coef_matrix(v, m);
            }
        }
        return Density;
    }
    double cal_SCF_Energy(Eigen::MatrixXf &Density, Eigen::MatrixXf &H_core, Eigen::MatrixXf &Fock) // return value will include V_NN
    {
        double E_SCF = 0;
        for (int miu = 0; miu < orbital_num; miu++)
        {
            for (int v = 0; v < orbital_num; v++)
                E_SCF += Density(miu, v) * (H_core(miu, v) + Fock(miu, v));
        }
        E_elec = E_SCF;
        return E_SCF + V_NN;
    }
    double cal_RMSD(Eigen::MatrixXf &D_1, Eigen::MatrixXf &D_2)
    {
        double rmsd = 0;
        for (int i = 0; i < orbital_num; i++)
            for (int j = 0; j < orbital_num; j++)
                rmsd += pow((D_1(i, j) - D_2(i, j)), 2);
        return sqrt(rmsd);
    }
    Eigen::MatrixXf get_Coef_matrix(Eigen::MatrixXf &Fock, Eigen::MatrixXf &S_sym)
    {
        Eigen::MatrixXf F_sym = S_sym.transpose() * Fock * S_sym;
        Eigen::EigenSolver<Eigen::MatrixXf> es(F_sym);
        Eigen::MatrixXf eigenvectors(orbital_num, orbital_num);
        vector<int> order(orbital_num, 0);
        for (int i = 0; i < orbital_num; i++)
        {
            int count = 0;
            for (int j = 0; j < orbital_num; j++)
            {
                if (i == j)
                    continue;
                if (es.eigenvalues()[i].real() - es.eigenvalues()[j].real() > 0)
                    count++;
            }
            order[i] = count;
        }
        Eigen::MatrixXf sorted_eigenvector = sort_columns(es.eigenvectors().real(), order);
        return S_sym * sorted_eigenvector;
    }
};

int main(int argc, char *argv[]) // command: ./project3 {-molecule} {-basis set} -{SCF}; default parameter is: h2o STO-3G tight
{
    // file processing;
    string Workdir = "F:/cpp/ProgrammingProjects-master/Project#03/input/h2o/";
    string output = "F:/cpp/ProgrammingProjects-master/Project#03/input/h2o.log";
    string basis = "STO-3G";
    string inter = "h2o";
    double E_limit = 1e-6, RMSD_limit = 1e-6;
    int SCF_LIMIT = 100; // max iter num for normal limit
    for (int i = 0; i < argc; i++)
    {
        if (argc <= 1) // default arg.
            break;
        if (i == 1)
        {
            inter = argv[1];
            Workdir = ("F:/cpp/ProgrammingProjects-master/Project#03/input/" + inter + "/").c_str();
            cout << format("Set Workdir to {:s}.", Workdir.c_str()) << endl;
        }
        if (i == 2)
        {
            basis = argv[2];
            cout << "Basis set is " << basis << endl;
        }
        if (i == 3)
        {
            cout << argv[3] << endl;
            if (argv[3] == "tight")
            {
                E_limit = 1e-12;
                RMSD_limit = 1e-12;
                SCF_LIMIT = 200;
                cout << "Converge limit of SCF is set to tight." << endl;
            }
            if (argv[3] == "normal")
            {
                cout << "Converge limit of SCF is set to normal." << endl;
                continue;
            }
            if (argv[3] == "loose")
            {
                cout << "Converge limit of SCF is set to loose." << endl;
                E_limit = 1e-3;
                RMSD_limit = 1e-3;
                SCF_LIMIT = 60;
            }
        }
    }
    output = ("F:/cpp/ProgrammingProjects-master/Project#03/input/" + inter + "/" + basis + "/" + inter + ".log").c_str();
    ofstream fout;
    fout.open(output, ios::out);
    if (!fout.is_open())
    {
        cout << "Error: Could not write the log file! Please ensure the file path is vaild." << endl;
        return 0;
    }
    fout << format("Workdir: {:s}.", Workdir.c_str()) << endl;
    fout << format("Output file path: {:s}.", output.c_str()) << endl;

    // initalize molecule
    Molecular mol;
    mol.build_Coords((Workdir + basis + "/geom.dat").c_str());
    mol.get_orbital_num(basis);

    // read integral
    mol.read_VNN((Workdir + basis + "/enuc.dat").c_str());
    Eigen::MatrixXf H_core = mol.read_1eint(Workdir.c_str(), basis.c_str());
    Eigen::MatrixXf S = mol.read_overlap_int(Workdir.c_str(), basis.c_str());
    vector<double> rep = mol.read_rep_int(Workdir.c_str(), basis.c_str());

    // preparation
    Eigen::MatrixXf S_sym = mol.ortho_overlap(S);
    Eigen::MatrixXf D_init = mol.initial_guess_D(S_sym, H_core);
    double pre_Ele = mol.get_Ele();

    fout << "The first Density Matrix" << endl;
    for (int i = 0; i < mol.get_orbital_num(); i++)
    {
        for (int j = 0; j < mol.get_orbital_num(); j++)
        {
            fout << format("{:>14.8f}", D_init(i, j)) << "  ";
        }
        fout << endl;
    }
    fout << endl;
    // SCF Procedure
    double rmsd = 1, delta_E = 1, E_tot = 0; // E_tot = E_ele + V_NN
    int iter = 1;
    Eigen::MatrixXf Fock = mol.build_Fock_matrix(H_core, D_init, rep);
    fout << "The first Fock Matrix" << endl;
    for (int i = 0; i < mol.get_orbital_num(); i++)
    {
        for (int j = 0; j < mol.get_orbital_num(); j++)
        {
            fout << format("{:>14.8f}", Fock(i, j)) << "  ";
        }
        fout << endl;
    }

    fout << endl;
    fout << format("{:^4s}{:^18s}{:^18s}{:^18s}", "Iter", "Total Energy", "Delta_E", "RMS(D)") << endl;
    fout << format("{:>4d}{:>18.12f}", 0, mol.get_Ele() + mol.get_VNN()) << endl;

    Eigen::MatrixXf D = mol.build_Density_matrix(Fock, S_sym);
    Eigen::MatrixXf pre_D = D_init;
    E_tot = mol.cal_SCF_Energy(D, H_core, Fock);
    delta_E = mol.get_Ele() - pre_Ele;
    rmsd = mol.cal_RMSD(D, pre_D);
    fout << format("{:>4d}{:>18.12f}{:>18.12f}{:>18.12f}", iter, mol.get_Ele() + mol.get_VNN(), delta_E, rmsd) << endl;
    bool SCF_sign = true;
    while (rmsd > RMSD_limit || fabs(delta_E) > E_limit) // update Fock and Density matrix.
    {
        iter++;
        if (iter > SCF_LIMIT)
        {
            fout << "\nAlready reached step limit, now exiting SCF procedure...\n";
            SCF_sign = false;
            break;
        }
        pre_Ele = mol.get_Ele();
        pre_D = D;
        Eigen::MatrixXf pre_Fock = Fock;
        if (iter > 10 && iter <= 20) // open dumping method
            Fock = 0.8 * pre_Fock + 0.2 * mol.build_Fock_matrix(H_core, D, rep);
        else if (iter > 20)
            Fock = 0.99 * pre_Fock + 0.01 * mol.build_Fock_matrix(H_core, D, rep);
        else
            Fock = mol.build_Fock_matrix(H_core, D, rep);

        D = mol.build_Density_matrix(Fock, S_sym);
        E_tot = mol.cal_SCF_Energy(D, H_core, Fock);
        delta_E = mol.get_Ele() - pre_Ele;
        rmsd = mol.cal_RMSD(D, pre_D);

        fout << format("{:>4d}{:>18.12f}{:>18.12f}{:>18.12f}", iter, E_tot, delta_E, rmsd) << endl;
    }
    if (SCF_sign == true)
    {
        cout << "\nSuccessfully complete." << endl;
        fout << "\nSuccessfully complete." << endl;
    }
    else
    {
        cout << "\nSCF completed without convergence." << endl;
        fout << "\nSCF completed without convergence." << endl;
    }
    fout.close();
    cout << "\nProgram Normally Terminated.\n";
    return 0;
}