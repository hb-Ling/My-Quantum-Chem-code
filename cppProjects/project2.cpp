/* UTF-8 */
/* This cpp file is used for practice.*/
/*This code will:
1. Read the geometry (like project1).
2. Read the Hessian matrix, then mass weighted it.
3. Diagonalize the mass-weighted Hessian, then obtain harmonical vibration freq.
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

vector<double> cross(vector<double> v1, vector<double> v2) // return v1 X v2
{
    assert(v1.size() == 3 && v2.size() == 3);
    vector<double> ans(3, 0);
    ans[0] = v1[1] * v2[2] - v1[2] * v2[1];
    ans[1] = v1[2] * v2[0] - v1[0] * v2[2];
    ans[2] = v1[0] * v2[1] - v1[1] * v2[0];
    return ans;
}

class Molecular
{
private:
    int atoms_num = 0; // molecular contains how many atoms --> at least one
    int *atom_list = NULL;
    double (*Coords)[3] = NULL; // Z_atom (x,y,z)     ; a 2-D array with size [atoms_num][3]
    double (*Forces)[3] = NULL; // Z_atom (fx,fy,fz)  ; a 2-D array with size [atoms_num][3]
public:
    Molecular(void)
    {
    }
    ~Molecular() // Destructor
    {
        atoms_num = 0;
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
    vector<vector<double>> get_distance(void) // build distance martix. (use this when frequently requiring distance value)
    {
        vector<vector<double>> dis; // vector<vector<double>>(atoms_num, vector(atoms_num, 0));
        dis.resize(atoms_num);
        for (int i = 0; i < atoms_num; i++)
            dis[i].resize(atoms_num);
        for (int i = 0; i < atoms_num; i++)
        {
            for (int j = 0; j < atoms_num; j++)
            {
                if (i == j)
                    dis[i][i] = 0;
                else if (i > j)
                    dis[i][j] = dis[j][i];
                else
                    dis[i][j] = sqrt(pow((Coords[i][0] - Coords[j][0]), 2) + pow((Coords[i][1] - Coords[j][1]), 2) + pow((Coords[i][2] - Coords[j][2]), 2));
            }
        }
        return dis;
    }
    double get_distance(int i, int j) // calculate distance between atom i and atom j.
    {
        if (i >= atoms_num || j >= atoms_num)
        {
            cout << "From function get_distance():\nError! Invaild atom label." << endl;
            return -1;
        }
        return sqrt(pow((Coords[i][0] - Coords[j][0]), 2) + pow((Coords[i][1] - Coords[j][1]), 2) + pow((Coords[i][2] - Coords[j][2]), 2));
    }

    vector<double> get_center_of_mass() // return coordinate [x,y,z] in Bohr
    {
        if (atoms_num < 1)
        {
            cout << "Error! Invalid Molecule." << endl;
            return {};
        }
        vector<double> ans(3, 0);
        for (int a = 0; a < 3; a++)
        {
            double sum_mass = 0;
            for (int i = 0; i < atoms_num; i++)
            {
                ans[a] += Atom_mass[atom_list[i] - 1] * Coords[i][a];
                sum_mass += Atom_mass[atom_list[i] - 1];
            }
            ans[a] /= sum_mass;
        }
        return ans;
    }
    void normalize_coords(vector<double> center) // translate the molecular and set center of mass to (0,0,0).
    {
        if (atoms_num < 1)
            return;
        if (center.size() != 3)
            return;
        for (int i = 0; i < atoms_num; i++)
        {
            Coords[i][0] -= center[0];
            Coords[i][1] -= center[1];
            Coords[i][2] -= center[2];
        }
        return;
    }

    vector<vector<double>> read_Hessian(const char *input) // return a 3N X 3N matrix
    {
        if (atoms_num < 1)
        {
            cout << "Error from read_Hessian():\nInvalid Molecular! Too few atoms." << endl;
            return {};
        }
        int N = atoms_num; // 3N X 3N matrix
        vector<vector<double>> hessian;
        hessian.resize(3 * N);
        for (int i = 0; i < 3 * N; i++)
            hessian[i].resize(3 * N);

        ifstream fin;
        fin.open(input, ios::in); // file format: (3N)^2 X 3
        if (!fin.is_open())
        {
            cout << "Error: Could not open the Hessian file! Please ensure the file path is vaild." << endl;
            return {};
        }
        char buf[512] = {0};
        bool start = true;
        int count = 0, i = 0, j = 0;
        while (fin.getline(buf, sizeof(buf)))
        {
            if (start == true)
            {
                int length;
                sscanf(buf, "%d", &length);
                assert(atoms_num == length);
                start = false;
                continue;
            }
            if (count == N)
            {
                count = 0;
                j = 0;
                i += 1;
            }
            double a, b, c;
            sscanf(buf, "%lf%lf%lf", &a, &b, &c);
            hessian[i][j] = a;
            hessian[i][j + 1] = b;
            hessian[i][j + 2] = c;
            count += 1;
            j += 3;
        }
        fin.close();
        return hessian;
    }
    vector<vector<double>> mass_weighted(vector<vector<double>> &hessian) // return mass_weighted hessian
    {
        for (int i = 0; i < 3 * atoms_num; i++)
        {
            for (int j = 0; j < 3 * atoms_num; j++)
            {
                hessian[i][j] = hessian[i][j] / sqrt(Atom_mass[atom_list[i / 3] - 1] * Atom_mass[atom_list[j / 3] - 1]);
            }
        }
        return hessian;
    }
    vector<double> get_freq(vector<vector<double>> &hessian)
    {
        vector<double> eigens(3 * atoms_num, 0);
        Eigen::MatrixXf Hess(3 * atoms_num, 3 * atoms_num);
        for (int i = 0; i < 3 * atoms_num; i++)
            for (int j = 0; j < 3 * atoms_num; j++)
                Hess(i, j) = hessian[i][j];
        Eigen::EigenSolver<Eigen::MatrixXf> es(Hess);
        for (int i = 0; i < 3 * atoms_num; i++)
        {
            eigens[i] = real(es.eigenvalues()[i]);
        }
        sort(eigens.begin(), eigens.end(), greater<double>());
        return eigens;
    }
};

int main(int argc, char *argv[])
{
    const char *input = "F:/cpp/ProgrammingProjects-master/Project#02/input/h2o";
    const char *output = NULL;
    for (int i = 0; i < argc; i++)
    {
        if (argc <= 1) // default arg.
            continue;
        if (i == 1) // read input file name, or use default setting.
        {
            input = argv[1];
        }
    }
    string inter = input;
    string geom_file = inter + "_geom.txt";
    string hess_file = inter + "_hessian.txt";
    output = (inter + ".log").c_str();
    ofstream fout;
    fout.open(output, ios::out);
    if (!fout.is_open())
    {
        cout << "Error: Could not write the xyz file! Please ensure the file path is vaild." << endl;
        return 0;
    }

    Molecular mol;
    mol.build_Coords(geom_file.c_str());
    vector<vector<double>> hess;
    hess.resize(mol.get_num());
    for (int i = 0; i < mol.get_num(); i++)
        hess[i].resize(mol.get_num());
    hess = mol.read_Hessian(hess_file.c_str());
    hess = mol.mass_weighted(hess);

    fout << "Mass_weighted Hessian (Force Constant) Matrix: (in atom units)" << endl;
    for (int i = 0; i < 3 * mol.get_num(); i++)
    {
        for (int j = 0; j < 3 * mol.get_num(); j++)
        {
            fout << format("{:>12.6f}", hess[i][j]);
        }
        fout << endl;
    }
    fout << endl;
    vector<double> eigens(3 * mol.get_num(), 0);
    eigens = mol.get_freq(hess);
    fout << "Eigenvalues of Mass-weighted Hessian Matrix (in hartree/amu-bohr^2):" << endl;
    for (int i = 0; i < 3 * mol.get_num(); i++)
        fout << format(" {:<3d}: {:>12.6f}", i, eigens[i]) << endl;
    fout << "Frequency (cm^-1):" << endl;
    for (int i = 0; i < 3 * mol.get_num(); i++)
    {
        if (eigens[i] < 0)
            fout << format(" {:<3d}:{:>12.3f} (Imaginary Freqency)", i, -sqrt(fabs(eigens[i]) * HARTREE_TO_J * 1e3 / (AMU_TO_G * pow(BOHR_TO_A * 1e-10, 2))) / (_C * 100 * 2 * M_PI)) << endl;
        else
            fout << format(" {:<3d}:{:>12.3f}", i, sqrt(eigens[i] * HARTREE_TO_J * 1e3 / (AMU_TO_G * pow(BOHR_TO_A * 1e-10, 2))) / (_C * 100 * 2 * M_PI)) << endl;
    }
    fout.close();
    cout << "\nSuccessful Terminated.\n";
    return 0;
}