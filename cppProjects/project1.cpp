/* UTF-8 */
/* This cpp file is used for practice.*/

#include <vector>
#include <fstream>
#include <iostream>
#include <assert.h>
#include <cmath>
#include <algorithm>
#include <format>
#include <F:/cpp/eigen-3.4.0/Eigen/Dense>
#include <F:/cpp/eigen-3.4.0/Eigen/Eigenvalues>
#define M_PI 3.141592654
#define BOHR_TO_A 0.52918
#define RAD_TO_ANG 180 / M_PI
#define AMU_TO_G 1.6605e-24
#define H_PLANCK 6.626e-34
#define _LIMIT 10e-3
#define _C 2.99792458e8 // light speed
#define _judgeBond 1.4  // use to judge whether two atoms can form a bond.
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
    // double (*Forces)[3] = NULL; // Z_atom (fx,fy,fz)  ; a 2-D array with size [atoms_num][3]
public:
    Molecular(void)
    {
    }
    ~Molecular() // Destructor
    {
        atoms_num = 0;
        delete Coords;
        delete atom_list;
        // delete Forces;
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

    double get_angle(int i, int j, int k, bool only_bond = false) // calculate angle between atom i-j-k, assert j is center atom.
    {
        if (i >= atoms_num || j >= atoms_num || k >= atoms_num)
        {
            cout << "From function get_angle()\nError! Invaild atom label" << endl;
            return -1;
        }
        double dist_a = get_distance(i, j);
        double dist_b = get_distance(j, k);
        if (only_bond == true)
        {
            if (dist_a > _judgeBond * (Radius[atom_list[i] - 1] + Radius[atom_list[j] - 1]))
                return -1;
            if (dist_b > _judgeBond * (Radius[atom_list[j] - 1] + Radius[atom_list[k] - 1]))
                return -1;
        }

        double vec_ji[3] = {0}, vec_jk[3] = {0};
        for (int a = 0; a < 3; a++)
        {
            vec_ji[a] = -(Coords[j][a] - Coords[i][a]) / dist_a;
            vec_jk[a] = -(Coords[j][a] - Coords[k][a]) / dist_b;
        }
        double dot = vec_ji[0] * vec_jk[0] + vec_ji[1] * vec_jk[1] + vec_ji[2] * vec_jk[2];
        if (fabs(dot - 1) < _LIMIT)
            return 0;
        if (fabs(dot + 1) < _LIMIT)
            return M_PI;
        return acos(dot);
    }
    double get_plane_angle(int i, int j, int k, int l) // return out_of_plane_angle i-j-k-l
    {
        if (i == j || i == k || i == l || j == k || j == l || k == l)
        {
            cout << "From function get_plane_angle():\nError! Repeat atom label" << endl;
            return -1;
        }
        if (i >= atoms_num || j >= atoms_num || k >= atoms_num || l >= atoms_num)
        {
            cout << "From function get_plane_angle():\nError! Invalid atom label" << endl;
            return -1;
        }

        double dist_ij = get_distance(i, j), dist_jk = get_distance(k, j), dist_kl = get_distance(k, l);
        if (dist_ij > _judgeBond * (Radius[atom_list[i] - 1] + Radius[atom_list[j] - 1]))
            return -1;
        if (dist_jk > _judgeBond * (Radius[atom_list[k] - 1] + Radius[atom_list[j] - 1]))
            return -1;
        if (dist_kl > _judgeBond * (Radius[atom_list[k] - 1] + Radius[atom_list[l] - 1]))
            return -1;

        double dist_ki = get_distance(k, i);
        vector<double> vec_kl(3, 0), vec_kj(3, 0), vec_ki(3, 0), vec_med(3, 0);
        for (int a = 0; a < 3; a++)
        {
            vec_kl[a] = -(Coords[k][a] - Coords[l][a]) / dist_kl;
            vec_kj[a] = -(Coords[k][a] - Coords[j][a]) / dist_jk;
            vec_ki[a] = -(Coords[k][a] - Coords[i][a]) / dist_ki;
        }
        vec_med = cross(vec_kj, vec_kl);
        return asin((vec_med[0] * vec_ki[0] + vec_med[1] * vec_ki[1] + vec_med[2] * vec_ki[2]) / sin(get_angle(j, k, l)));
    }
    double get_torsion_angle(int i, int j, int k, int l) // return torsion angle i-j-k-l
    {
        if (i == j || i == k || i == l || j == k || j == l || k == l)
        {
            cout << "From function get_torsion_angle():\nError! Repeat atom label" << endl;
            return -1;
        }
        if (i >= atoms_num || j >= atoms_num || k >= atoms_num || l >= atoms_num)
        {
            cout << "From function get_torsion_angle():\nError! Invalid atom label" << endl;
            return -1;
        }

        double dist_ij = get_distance(i, j), dist_jk = get_distance(k, j), dist_kl = get_distance(k, l);
        if (dist_ij > _judgeBond * (Radius[atom_list[i] - 1] + Radius[atom_list[j] - 1]))
            return -1;
        if (dist_jk > _judgeBond * (Radius[atom_list[k] - 1] + Radius[atom_list[j] - 1]))
            return -1;
        if (dist_kl > _judgeBond * (Radius[atom_list[k] - 1] + Radius[atom_list[l] - 1]))
            return -1;

        vector<double> v1(3, 0), v2(3, 0), vec_ij(3, 0), vec_jk(3, 0), vec_kl(3, 0);
        for (int a = 0; a < 3; a++)
        {
            vec_ij[a] = -(Coords[i][a] - Coords[j][a]) / dist_ij;
            vec_jk[a] = -(Coords[j][a] - Coords[k][a]) / dist_jk;
            vec_kl[a] = -(Coords[k][a] - Coords[l][a]) / dist_kl;
        }
        v1 = cross(vec_ij, vec_jk);
        v2 = cross(vec_jk, vec_kl);
        double dot = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
        double divd = (sin(get_angle(i, j, k)) * sin(get_angle(j, k, l)));
        if (divd < 10e-8)
        {
            cout << "From function get_torsion_angle():\nError! Div with 0! Please check the geometry." << endl;
            return -1;
        }
        if (fabs(dot / divd - 1) < _LIMIT)
            return 0;
        if (fabs(dot / divd + 1) < _LIMIT)
            return M_PI;
        return acos(dot / divd);
    }
    void write_dist(ofstream &fout)
    {
        fout << "Distance List (in Angstrom)" << endl;
        for (int i = 0; i < atoms_num - 1; i++)
            for (int j = i + 1; j < atoms_num; j++)
                fout << format("{:>2s}{:<2d}-{:>2s}{:<2d}:{:8.2f}",
                               ATOMS[atom_list[i] - 1], i + 1, ATOMS[atom_list[j] - 1], j + 1, get_distance(i, j) * BOHR_TO_A)
                     << endl;
        fout << endl;
    }
    void write_angle(ofstream &fout)
    {
        fout << "Bond Angle List (in degree)" << endl;
        for (int i = 0; i < atoms_num - 1; i++)
        {
            for (int k = i + 1; k < atoms_num; k++)
            {
                for (int j = 0; j < atoms_num; j++)
                {
                    if (j == i || j == k)
                        continue;
                    double angles = get_angle(i, j, k, true);
                    if (angles < 0)
                        continue;
                    fout << format("{:>2s}{:<2d}-{:>2s}{:<2d}-{:>2s}{:<2d}:{:8.2f}",
                                   ATOMS[get_atlist(i) - 1], i + 1, ATOMS[get_atlist(j) - 1], j + 1, ATOMS[get_atlist(k) - 1], k + 1,
                                   angles * RAD_TO_ANG)
                         << endl;
                }
            }
        }
        fout << endl;
    }
    void write_plane_angle(ofstream &fout)
    {
        fout << "Plane Angle List (in degree)" << endl;
        for (int i = 0; i < atoms_num - 1; i++)
        {
            for (int l = i + 1; l < atoms_num; l++)
            {
                for (int j = 0; j < atoms_num; j++)
                {
                    if (j == i || j == l)
                        continue;
                    for (int k = 0; k < atoms_num; k++)
                    {
                        if (j == k || k == i || k == l)
                            continue;
                        double angles = get_plane_angle(i, j, k, l);
                        if (angles < 0)
                            continue;
                        fout << format("{:>2s}{:<2d}-{:>2s}{:<2d}-{:>2s}{:<2d}-{:>2s}{:<2d}:{:8.2f}",
                                       ATOMS[atom_list[i] - 1], i + 1, ATOMS[atom_list[j] - 1], j + 1, ATOMS[atom_list[k] - 1], k + 1,
                                       ATOMS[atom_list[l] - 1], l + 1, angles * RAD_TO_ANG)
                             << endl;
                    }
                }
            }
        }
        fout << endl;
    }
    void write_torsion_angle(ofstream &fout)
    {
        fout << "Torsional Angle List (in degree)" << endl;
        for (int i = 0; i < atoms_num - 1; i++)
        {
            for (int l = i + 1; l < atoms_num; l++)
            {
                for (int j = 0; j < atoms_num; j++)
                {
                    if (j == i || j == l)
                        continue;
                    for (int k = 0; k < atoms_num; k++)
                    {
                        if (j == k || k == i || k == l)
                            continue;
                        double angles = get_torsion_angle(i, j, k, l);
                        if (angles < 0)
                            continue;
                        fout << format("{:>2s}{:<2d}-{:>2s}{:<2d}-{:>2s}{:<2d}-{:>2s}{:<2d}:{:8.2f}",
                                       ATOMS[atom_list[i] - 1], i + 1, ATOMS[atom_list[j] - 1], j + 1, ATOMS[atom_list[k] - 1], k + 1,
                                       ATOMS[atom_list[l] - 1], l + 1, angles * RAD_TO_ANG)
                             << endl;
                    }
                }
            }
        }
        fout << endl;
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
    vector<vector<double>> get_inertia_tensor()
    { // return a 3x3 vector
        vector<vector<double>> ans;
        ans.resize(3);
        for (int i = 0; i < 3; i++)
            ans[i].resize(3);

        for (int i = 0; i < 3; i++) // 0 for x; 1 for y; 2 for z;
        {
            for (int j = 0; j < 3; j++)
            {
                if (i == j)
                {
                    for (int x = 0; x < atoms_num; x++)
                    {
                        for (int y = 0; y < 3; y++)
                        {
                            if (y == i)
                                continue;
                            ans[i][i] += Atom_mass[atom_list[x] - 1] * Coords[x][y] * Coords[x][y];
                            // printf("Mass: %12.6lfCoords: %12.6lf\n", Atom_mass[atom_list[x] - 1], Coords[x][y]);
                        }
                    }
                }
                else if (i < j)
                    for (int x = 0; x < atoms_num; x++)
                    {
                        ans[i][j] += Atom_mass[atom_list[x] - 1] * Coords[x][i] * Coords[x][j];
                    }
                else
                    ans[i][j] = ans[j][i];
            }
        }
        return ans;
    }
    vector<double> get_principle_moment(vector<vector<double>> tensor)
    {
        vector<double> ans;
        ans.resize(3);
        if (tensor.size() != 3)
        {
            cout << "Error! Invalid tensor matrix." << endl;
            return {};
        }
        Eigen::Matrix3d Tensor;
        for (int i = 0; i < 3; i++)
            for (int j = 0; j < 3; j++)
                Tensor(i, j) = tensor[i][j];
        Eigen::EigenSolver<Eigen::Matrix3d> es(Tensor);
        // Eigen::Vector3cd eigenvalues = es.eigenvalues();
        ans[0] = real(es.eigenvalues()[0]);
        ans[1] = real(es.eigenvalues()[1]);
        ans[2] = real(es.eigenvalues()[2]);
        sort(ans.begin(), ans.end());
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
    string judge_symm(vector<double> principle, double converge_limit = 0.01)
    {

        if (principle.size() != 3)
            return "Error Type.";
        if (atoms_num == 1)
            return "Single Atom";
        if (atoms_num == 2)
            return "Diatomic";
        if (fabs(principle[1] - principle[0]) < converge_limit || fabs(principle[1] - principle[2]) < converge_limit)
        {
            if (fabs(principle[2] - principle[0]) < converge_limit)
                return "Spherical top";
            else
                return "Symmetric top";
        }
        return "Asymmetric top";
    }
};

int main(int argc, char *argv[])
{
    const char *input = "F:/cpp/ProgrammingProjects-master/Project#01/input/benzene.dat";
    const char *output = "F:/cpp/ProgrammingProjects-master/Project#01/input/benzene.log";
    for (int i = 0; i < argc; i++)
    {
        if (i == 0)
            continue; // default arg.
        if (i == 1)   // read input file name, or use default setting.
        {
            input = argv[1];
            string intermediate = argv[1];
            int length = intermediate.size();
            int index = intermediate.find(".dat");
            if (index == string::npos)
            {
                cout << "Invalid File name!" << endl;
                return 0;
            }
            intermediate.replace(index, index + 4, ".log");
            output = intermediate.c_str();
        }
    }

    ofstream fout;
    fout.open(output, ios::out);
    if (!fout.is_open())
    {
        cout << "Error: Could not write the xyz file! Please ensure the file path is vaild." << endl;
        return 0;
    }
    Molecular mol;
    mol.build_Coords(input);
    // print coords, distances and angles
    mol.write_Coords(fout);
    mol.write_dist(fout);
    mol.write_angle(fout);
    mol.write_plane_angle(fout);
    mol.write_torsion_angle(fout);

    vector<double> center = mol.get_center_of_mass();
    fout << endl;
    fout << format("Center of mass (in Angstrom) = {:<12.6f}{:<12.6f}{:<12.6f}", center[0] * BOHR_TO_A, center[1] * BOHR_TO_A, center[2] * BOHR_TO_A) << endl;
    fout << endl;
    mol.normalize_coords(center);

    fout << "Tensor Matrix:" << endl;
    vector<vector<double>> tensor = mol.get_inertia_tensor();
    for (int i = 0; i < 3; i++)
        fout << format("{:>18.6f}{:>18.6f}{:>18.6f}", tensor[i][0], tensor[i][1], tensor[i][2]) << endl;
    fout << endl;

    vector<double> principle = mol.get_principle_moment(tensor);
    fout << format("Principle Moment I (Amu.Bohr^2):\n{:<18.6f}{:<18.6f}{:<18.6f}\n", principle[0], principle[1], principle[2]) << endl;

    fout << format("I in g*cm^2:\n{:<18.6e}{:<18.6e}{:<18.6e}\n", principle[0] * AMU_TO_G * pow(BOHR_TO_A * 1.0e-8, 2),
                   principle[1] * AMU_TO_G * pow(BOHR_TO_A * 1.0e-8, 2), principle[2] * AMU_TO_G * pow(BOHR_TO_A * 1.0e-8, 2))
         << endl;
    fout << format("Molecule Type: {:s}\n", (mol.judge_symm(principle)).c_str()) << endl;

    vector<double> Rotation(3, 0);
    for (int a = 0; a < 3; a++)
    {
        Rotation[a] = H_PLANCK * 1.0e7 * 1.0e-6 / (8 * M_PI * M_PI * principle[a] * AMU_TO_G * pow(BOHR_TO_A * 1.0e-8, 2));
        fout << format("Rotational Constant {:c}: {:<12.6f} MHz", 'A' + a, Rotation[a])
             << format("  or  {:<12.6f} cm^-1", Rotation[a] * 10e4 / _C) << endl;
    }
    fout << endl;
    fout << "Normal Termination." << endl;
    // fout << format("Time usage: {:6.2f} ms\n", 0.00) << endl;
    fout.close();
    printf("Successful Terminated.\n");
    return 0;
}