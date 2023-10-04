// This C++ code can be used to run Monte Carlo simulations of a lattice-based mesoscale model introduced in:
// - Long Li, Jinglei Hu, Bartosz Rozycki, Fan Song, Nano Letters 20(1):722-728 (2020)
// - Long Li, Jinglei Hu, Xinghua Shi, Bartosz Rozycki, Fan Song, Soft Matter 17(7):1912-1920 (2021)
//
// Authors: Lukasz Milewski & Bartosz Rozycki
// Institute of Physics
// Polish Academy of Sciences
// October 2023
//
// The development of this code was supported by the National Science Centre via grant number 2021/40/Q/NZ1/00017.
//

#include <iostream>
#include <fstream>
#include <random>
#include <cmath>
#include <string>
#include <cstdlib>
using namespace std;

#define K 50
// K = linear extension of the square lattice

random_device rd;
mt19937 gen(rd());
int seed = 1;
uniform_real_distribution<> dis(-0.5, 0.5);
uniform_real_distribution<> dist(0.0, 1.0);
// random number generator

int bc[K][2];
// auxiliary matrix for periodic boundary conditions

// arrays L, M, N, PNUM, RNUM below: [K][K][0] refers to the upper membrane, [K][K][1] refers to the lower membrane

double L[K][K][2];
// location of the two membranes relative to the reference plane z=0

int M[K][K][2];
// array with values 0 or 1 to indicate absence or presence of protein particles

int N[K][K][2];
// array with values 0 or 1 to indicate absence or presence of lipid rafts

int PNUM[K][K][2];
// indexes of the receptors and ligands

int RNUM[K][K][2];
// indexes of rafts on the two membranes

vector<int> particle_x, particle_y;
// Cartesian coordinates of all the protein particles

vector<int> raft_x, raft_y;
// Cartesian coordinates of all th lipid rafts

double Cp = 0.05;
// area concentration of the protein particles (in units of 1/a^2)

double X = 0.1;
// area fraction of the lipid rafts

int N0 = int(round(Cp * K * K));
// number of receptors

int N1 = int(round(Cp * K * K));
// number of ligands

int N2 = int(round(X * K * K));
// number of lipid rafts on the upper membrane

int N3 = int(round(X * K * K));
// number of lipid rafts on the lower membrane

double a = 10.0;
// lattice constant (in units of nm)

double kappa = 10.0;
// bending rigidity modulus (in kT units)

double U = 1.6;
// energy of contact between rafts (in kT units)

double U_a = 3.0;
// energy of association of a protein particle with a raft (in kT units)

double U_b = 6.0;
// energy of receptor-ligand binding (in kT units)

double lmax = 4.0;
// maximal variation of membrane vertical displacements (in units of nm)

double l_c = 15.0;
// length of the extracellular domains of the receptor-ligand complex (in units of nm)

double l_b = 1.0;
// width of the receptor-ligand binding potential (in units of nm)

double l_minus = l_c - 0.5 * l_b;
// lower end of the receptor-ligand binding range (in units of nm)

double l_plus = l_c + 0.5 * l_b;
// upper end of the receptor-ligand binding range (in units of nm)

double pos0 = 0.5 * l_c;
double pos1 = -0.5 * l_c;
// initial positions of the membranes (planar at the start of the simulation)

int nxi = 100;
// number of MC cycles between recording observables

int nlmax = 1000;
// number of MC cycles between computing the acceptance rate during the equilibration

int neq = 10000;
// number of MC cycles for equilibration

int ncyc = 20000;
// number of MC cycles used for data acquisition

int res = 0;
// set it to 1 if you read in an initial configuration from a restart file

int nsnapshot = 5000;
// number of MC cycles between recording system configurations

string header, filename, data_out, snapname;
// output file names and header

void arguments_values(char *argv[]);
// input parameters

void initialise();
// initialization of arrays and vectors, and placing particles and rafts randomly on the membranes

void save_config(ofstream &, int);
// saving configurations to file

void acceptance_rate_try(int, int);
// calculate acceptance rate

double averageL();
// calculate the lattice-averaged distance between the two membranes

double calc_xi();
// calculate the instantaneous relative roughness of the membranes

int how_many_pairs();
// calculate the number or receptor-ligand complexes

double energy_interactions_move(int);
// calculate the receptor-ligand interaction energy for protein particle trial moves and receptor-ligand pair trial moves

double energy_interactions_bending(double, double);
// calculate change in the receptor-ligand interaction energy for membranes bending trial moves

double energy_bending(int, int, int, double, double);
// calculate change in the membrane bending energy

double energy_r_r(int, int);
// calculate change in the energy of contacts between lipid rafts

double energy_r_p_particle(int, int);
// calculate the raft-particle association energy for protein particle trial moves

double energy_r_p_raft(int, int);
// calculate the raft-particle association energy for raft trial moves

void dynamics_particle();
// protein particle trial moves

void dynamics_pairs();
// receptor-ligand pair trial moves

void dynamics_rafts();
// lipid raft trial moves

int dynamics_bending();
// membrane bending trial moves

double total_me_energy(int);
// calculate the total energy of membrane bending

double total_r_p_energy();
// calculate the total energy of raft-particle association

double total_r_r_energy();
// calculate the total energy of raft-raft contacts

double total_RL_energy(int);
// calculate the total energy of receptor-ligand interactions

int main(int argc, char *argv[])
{
	srand(seed);
	// feed the random number generator with a seed
	
	if (argc == 12)
	{
		arguments_values(argv);
		// input parameters
	}
	else
	{
		cout << "default values\n";
	}
	initialise();
	if (N0 > K * K or N1 > K * K or N0 < 0 or N1 < 0 or ncyc < neq or l_c < 0 or l_b < 0 or lmax < 0 or kappa < 0)
	{
		cerr << "Wrong parameters! Program has been terminated" << endl;
		exit(1);
	}
	ofstream xi_out;
	// stream for data output file

	ofstream snapshot;
	// stream for configuration output file

	xi_out.open(data_out, ios_base::trunc);
	// open data output file

	snapshot.open(snapname, ios_base::trunc);
	// open configuration output file

	xi_out << header << "\n";
	xi_out << "# column 1: time" << endl;
	xi_out << "# column 2: average membrane separation" << endl;
	xi_out << "# column 3: membrane roughness" << endl;
	xi_out << "# column 4: number of R-L bonds" << endl;
	xi_out << "# column 5: number of free particles" << endl;
	xi_out << "# column 6: membrane bending energy" << endl;
	xi_out << "# column 7: raft contact energy" << endl;
	xi_out << "# column 8: particle-raft association energy" << endl;
	xi_out << "# column 9: energy of R-L bonds" << endl;
	xi_out << "# column 10: total energy" << endl;

	double xi, avl, acc_rate;
	// xi = membrane roughness
	// avl = lattice-averaged distance between the membranes
	// acc_rate = acceptance rate

	double num = 0;
	double E_me, E_rr, E_rp, E_RL, E_total;
	// E_me = total energy of membrane bending
	// E_rr = total energy of raft-raft interactions
	// E_rp = total energy of raft-particle association
	// E_RL = total energy of receptor-ligand binding
	// E_total = total energy of the simulation system

	int pairs = 0;
	// number of receptor-ligands complexes
	
	int free = 0;
	// number of free receptors

	int ntry = 0;
	// total number of trial moves for acceptance rate calculation

	int nacc = 0;
	// number of accepted trial moves for acceptance rate calculation

	for (int nmc = 1; nmc <= ncyc; nmc++)
	{
	// loop over MC cycles
		for (int kdyn = 0; kdyn < N2 + N3; kdyn++)
		{
			dynamics_rafts();
			// raft trial moves repeated N2+N3 times
		}
		for (int kdyn = 0; kdyn < N0 + N1; kdyn++)
		{
			dynamics_particle();
			// protein particle trial moves repeated N2+N3 times

			dynamics_pairs();
			// receptor-ligand pairs trial moves repeated N2+N3 times
		}
		for (int ix = 0; ix < K; ix++)
		{
			for (int iy = 0; iy < K; iy++)
			{
				for (int k = 0; k < 2; k++)
				{
					nacc = nacc + dynamics_bending();
					// membrane bending trial moves repeated 2*K*K times
					// dynamics_bending() returns 1 if trial move is accepted
					// dynamics_bending() returns 0 if trial move is rejected
					ntry++;
				}
			}
		}
		if (nmc % nxi == 0 and nmc > neq)
		{
			avl = averageL();
			xi = calc_xi();
			pairs = how_many_pairs();
			free = N0 - pairs;
			E_me = total_me_energy(0) + total_me_energy(1);
			E_rr = total_r_r_energy();
			E_rp = total_r_p_energy();
			E_RL = total_RL_energy(pairs);
			E_total = E_me + E_rr + E_rp + E_RL;
			xi_out << nmc - neq << " " << avl << " " << xi << " " << pairs << " " << N0 - pairs << " " << E_me << " "
				<< E_rr << " " << E_rp << " " << E_RL << " " << E_total << endl;
			// write data to output file
		}
		if (nmc % nsnapshot == 0)
		{
			save_config(snapshot, nmc);
			// write configuration to file
		}
		if (nmc % nlmax == 0 and nmc < neq)
		{
			acceptance_rate_try(nacc, ntry);
			// compute acceptance rate and adjust lmax

			nacc = 0;
			ntry = 0;
		}
	}
	snapshot.close();
	// close configuration output file

	xi_out.close();
	// close data output file

	return 0;
}

void arguments_values(char *argv[])
// parameters given as input when program starts
{
	Cp = stod(argv[1]);
	N0 = int(round(Cp * K * K));
	N1 = int(round(Cp * K * K));
	X = stod(argv[2]);
	N2 = int(round(X * K * K));
	N3 = int(round(X * K * K));
	kappa = stod(argv[3]);
	U = stod(argv[4]);
	U_a = stod(argv[5]);
	U_b = stod(argv[6]);
	neq = stoi(argv[7]);
	ncyc = stoi(argv[8]);
	nxi = stoi(argv[9]);
	nsnapshot = stoi(argv[10]);
	seed = stoi(argv[11]);
}

void initialise()
// initialization of arrays and vectors
// initial random distribution of protein particles and rafts randomly on the membranes
{
	header = "# 3-dimensional membrane model \n# square lattice with " 
		 + to_string(K * K) + " sites \n"
		 + "# Cp = " + to_string(Cp) + "\n"
		 + "# x = " + to_string(X) + "\n"
		 + "# N0 = " + to_string(N0) + "\n"
		 + "# N1 = " + to_string(N1) + "\n"
		 + "# N2 = " + to_string(N2) + "\n"
		 + "# N3 = " + to_string(N3) + "\n"
		 + "# a = " + to_string(a) + "\n"
		 + "# kappa = " + to_string(kappa) + "\n"
		 + "# l_max0 = " + to_string(lmax) + "\n"
		 + "# U = " + to_string(U) + "\n"
		 + "# U_a = " + to_string(U_a) + "\n"
		 + "# U_b = " + to_string(U_b) + "\n"
		 + "# l_c = " + to_string(l_c) + "\n"
		 + "# l_b = " + to_string(l_b) + "\n"
		 + "# nequil = " + to_string(neq) + "\n"
		 + "# ncycle = " + to_string(ncyc) + "\n"
		 + "# nxi = " + to_string(nxi) + "\n"
		 + "# nlmax = " + to_string(nlmax) + "\n"
		 + "# nsnapshot = " + to_string(nsnapshot) + "\n"
		 + "# RNG seed = " + to_string(seed);
	
	filename = "K_" + to_string(K)
		+ "_cp_" + to_string(Cp)
		+ "_x_" + to_string(X)
		+ "_U_" + to_string(U)
		+ "_U_a_" + to_string(U_a)
		+ "_U_b_" + to_string(U_b)
		+ "_" + to_string(int(seed));

	data_out = filename + "_data.out";
	snapname = filename + "_snapshot.txt";
	
	particle_x.resize(N0 + N1);
	// set the size of the vector with the x-coordinates of protein particles
	
	particle_y.resize(N0 + N1);
	// set the size of the vector with the y-coordinates of protein particles
	
	raft_x.resize(N2 + N3);
	// set the size of the vector with the x-coordinates of raft particles
	
	raft_y.resize(N2 + N3);
	// set the size of the vector with the y-coordinates of raft particles
	
	for (int i = 0; i < K; i++)
	{
		bc[i][0] = i + 1;
		bc[i][1] = i - 1;
	}
	bc[K - 1][0] = 0;
	bc[0][1] = K - 1;
	// bc matrix is used to implement periodic boundary conditions
	
	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < K; j++)
		{
			for (int k = 0; k < 2; k++)
			{
				M[i][j][k] = 0;
				PNUM[i][j][k] = -1;
				N[i][j][k] = 0;
				RNUM[i][j][k] = -1;
				// at start, membranes with no receptors, ligands or rafts
			}
		}
	}
	if (res != 1)
	// randomly distribute protein particles and rafts on the membranes

	{
		for (int i = 0; i < K; i++)
		{
			for (int j = 0; j < K; j++)
			{
				L[i][j][0] = pos0;
				L[i][j][1] = pos1;
				// planar membranes located at z = pos0 (upper membrane) and z = pos1 (lower membrane)
			}
		}
		int aa = 0;
		// counter
		
		int rx, ry;
		while (aa < N0)
		// distribute N0 receptors numbered from 0 to N0-1
		{
			rx = rand() % K;
			ry = rand() % K;
			// randomly choose x- and y-coordinates to place a receptor on the upper membrane
			if (M[rx][ry][0] == 0)
			{
				M[rx][ry][0] = 1;
				// if the chosen membrane patch is not occupied by another receptor, placed the receptor with index aa here
			
				particle_x[aa] = rx;
				particle_y[aa] = ry;
				// store the x- and y-coordinates of the receptor with index aa
			
				PNUM[rx][ry][0] = aa;
				// receptor index aa stored in particle index array
				aa++;
			}
		}
		aa = 0;
		// reset counter
		
		while (aa < N1)
		// distribute N1 ligands numbered from N0 to N0+N1-1
		{
			rx = rand() % K;
			ry = rand() % K;
			// randomly choose x- and y-coordinates to place a ligand on the lower membrane
			if (M[rx][ry][1] == 0)
			{
				M[rx][ry][1] = 1;
				// if the chosen membrane patch is not occupied by another ligand, placed the ligand with index aa here
			
				particle_x[aa + N0] = rx;
				particle_y[aa + N0] = ry;
				// store the x- and y-coordinates of the ligand with index aa
				
				PNUM[rx][ry][1] = aa;
				// ligand index aa stored in particle index array
				aa++;
			}
		}
		aa = 0;
		//reset counter
		
		while (aa < N2)
		// distribute N2 rafts numbered from 0 to N2-1 on the upper membrane
		{
			rx = rand() % K;
			ry = rand() % K;
			// randomly choose x- and y-coordinates to place a raft on the upper membrane
			if (N[rx][ry][0] == 0)
			{
				N[rx][ry][0] = 1;
				// if the chosen membrane patch is not occupied by another raft, placed the raft with index aa here

				raft_x[aa] = rx;
				raft_y[aa] = ry;
				// store the x- and y-coordinates of the raft with index aa

				RNUM[rx][ry][0] = aa;
				// raft index aa stored in raft index array
				aa++;
			}
		}
		aa = 0;
		//reset counter

		while (aa < N3)
		// distribute N3 rafts numbered from N2 to N2+N3-1 on the lower membrane
		{
			rx = rand() % K;
			ry = rand() % K;
			// randomly choose x- and y-coordinates to place a raft on the lower membrane
			if (N[rx][ry][1] == 0)
			{
				N[rx][ry][1] = 1;
				// if the chosen membrane patch is not occupied by another raft, placed the raft with index aa here

				raft_x[aa + N2] = rx;
				raft_y[aa + N2] = ry;
				// store the x- and y-coordinates of the raft with index aa

				RNUM[rx][ry][1] = aa;
				// raft index aa stored in raft index array
				aa++;
			}
		}
	}
}

void save_config(ofstream &fout, int t)
// save configurations to file
{
	int rn, ln, lun, lln;
	// nr = index of a receptor
	// nl = index of a ligand
	// lun = index of an upper membrane raft
	// lln = index of a lower membrane raft

	fout << "# MC step: " << t << endl;
	for (int x = 0; x < K; x++)
	{
		for (int y = 0; y < K; y++)
		{
			rn = PNUM[x][y][0] + 1;
			ln = PNUM[x][y][1] + 1;
			lun = RNUM[x][y][0] + 1;
			lln = RNUM[x][y][1] + 1;
			// 1 is added to have indexes as natural numbers larger than 0

			fout << x << " " << y << " " << L[x][y][0] << " " << L[x][y][1] << " " << rn << " " << ln << " " << lun
				 << " " << lln << endl;
			// eight columns: 
			// 1. x coordinate
			// 2. y coordinate
			// 3. local position of the upper membrane at this lattice site
			// 4. local position of the lower membrane at this lattice site
			// 5. receptor index (0 if none) at this lattice site
			// 6. ligand index (0 if none) at this lattice site
			// 7. raft index (0 if none) at this lattice site
			// 8. raft index (0 if none) at this lattice site
		}
	}
}

void acceptance_rate_try(int nacc, int ntry)
// calculate the acceptance rate for membrane bending trial moves and adjust lmax
{
	double a, b, r;
	a = nacc;
	b = ntry;
	r = a / b;
	// acceptance rate = number of accepted trial moves divided by total number of trial moves 
	
	if (r > 0.5)
	//if acceptance rate higher than 50%
	{
		lmax = lmax * 1.05;
	}
	if (r < 0.3)
	//if acceptance rate lower than 30%
	{
		lmax = lmax * 0.95;
	}
}

double averageL()
// lattice-averaged distance between the membranes
{
	double ll;
	// lattice-averaged distance between the membranes

	double li;
	// local distance between the membranes at a given lattice site

	double sum = 0.0;
	// sum of li

	for (int ix = 0; ix < K; ix++)
	{
		for (int iy = 0; iy < K; iy++)
		{
			li = L[ix][iy][0] - L[ix][iy][1];
			// local distance between the membranes = difference between the local positions

			sum = sum + li;
		}
	}
	ll = sum / (K * K);
	return ll;
}

double calc_xi()
// instantaneous relative roughness of the membranes
{
	double xi;
	// membrane roughness
	
	double li;
	// local distance between the membranes
	
	double sum, sumsqr, mean, meansqr;
	sum = 0.0;
	sumsqr = 0.0;

	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < K; j++)
		{
			li = L[i][j][0] - L[i][j][1];
			// local distance between the membranes = difference between the local positions

			sum = sum + li;
			sumsqr = sumsqr + li * li;
		}
	}
	mean = sum / (K * K);
	meansqr = sumsqr / (K * K);
	xi = sqrt(meansqr - (mean * mean));
	return xi;
}

int how_many_pairs()
// number of receptor-ligand complexes
{
	int pairs = 0;
	double li;
	// local distance between the membranes

	for (int i = 0; i < K; i++)
	{
		for (int j = 0; j < K; j++)
		{
			if (M[i][j][0] == 1 and M[i][j][1] == 1)
			// there is a receptor on the upper membrane and a ligand on the lower membrane
			{
				li = L[i][j][0] - L[i][j][1];
				// local distance between the membranes = difference between the local positions
				
				if (li > l_minus and li < l_plus)
				// specified distance range required for the receptor-ligand binding 
				{
					pairs++;
				}
			}
		}
	}
	return pairs;
}

double energy_interactions_move(int nr)
// receptor-ligand interaction energy for trial moves of protein particles and their pairs 
// nr = particle index
{
	double E = 0;
	// receptor-ligand interaction energy

	double li;
	// local distance between the membranes

	int x, y;
	// x- and y-coordinates

	int n, m;
	// receptor absence (n=0) or presence (n=1)
	// ligand absence (m=0) or presence (m=1)

	x = particle_x[nr];
	y = particle_y[nr];
	// x- and y-coordinates of particle with index nr

	n = M[x][y][0];
	// receptor absence (n=0) or presence (n=1) at this patch of the upper membrane

	m = M[x][y][1];
	// ligand absence (m=0) or presence (m=1) at this patch of the lower membrane

	li = L[x][y][0] - L[x][y][1];
	// local distance between the membranes

	if (n == 1 and m == 1)
	// there is a receptor on the upper membrane and a ligand on the lower membrane
	{
		if (li > l_minus and li < l_plus)
		// specified distance range required for the receptor-ligand binding
		{
			E = -U_b;
			// receptor-ligand interaction energy < 0
		}
	}
	return E;
}

double energy_interactions_bending(double li_old, double li_new)
// change in the receptor-ligand interaction energy for trial moves of membrane bending
// li_old and li_new, respectively, denote the local distance between the membranes before and after the trial move
{
	double E_old = 0;
	// receptor-ligand energy before the membrane bending trial move

	double E_new = 0;
	// receptor-ligand energy after the membrane bending trial move
	
	if (li_new > l_minus and li_new < l_plus)
	// specified distance range required for the receptor-ligand binding
	{
		E_new = -U_b;
	}
	if (li_old > l_minus and li_old < l_plus)
	// specified distance range required for the receptor-ligand binding
	{
		E_old = -U_b;
	}
	return E_new - E_old;
	// receptor-ligand interaction energy change
}

double energy_bending(int x, int y, int k, double l_old, double l_new)
// change in the membrane bending energy for membrane trial moves 
// x and y are the coordinates of the selected site
// k is the membrane index: k=0 for the upper membrane, k=1 for the lower membrane
// l_old is local membrane position at the selected site before the trial move
// l_new is local membrane position at the selected site after the trial move
{
	vector<double> l(13);
	// local positions of the given site and of its 12 neighbors

	int x_right, x_left, x_right_right, x_left_left, y_up, y_down, y_up_up, y_down_down;
	// x- and y-coordinates of neighbors

	double E;
	// membrane bending energy change

	double ka = kappa / (a * a);
	// constant in an energy change equation

	// below: periodic boundary conditions are implemented using the bc matrix

	x_right = bc[x][0];
	x_left = bc[x][1];
	y_up = bc[y][1];
	y_down = bc[y][0];
	// nearest neighbors of the given site

	x_right_right = bc[x_right][0];
	x_left_left = bc[x_left][1];
	y_up_up = bc[y_up][1];
	y_down_down = bc[y_down][0];
	// next-nearest neighbors of the given site

	l[0] = L[x][y][k];
	// local position at the given site

	l[1] = L[x_right][y][k];
	l[2] = L[x][y_down][k];
	l[3] = L[x_left][y][k];
	l[4] = L[x][y_up][k];
	// local positions at its nearest neighbors

	l[5] = L[x_right][y_down][k];
	l[6] = L[x_left][y_down][k];
	l[7] = L[x_left][y_up][k];
	l[8] = L[x_right][y_up][k];
	// local positions at its next-nearest neighbors

	l[9] = L[x_right_right][y][k];
	l[10] = L[x][y_down_down][k];
	l[11] = L[x_left_left][y][k];
	l[12] = L[x][y_up_up][k];
	// local positions at its next-next nearest neighbors

	/* top view on the membrane:
	..................l[12].................
	..........l[7].....l[4]....l[8].........
	.l[11]....l[3].....l[0]....l[1]....l[9].
	..........l[6].....l[2]....l[5].........
	..................l[10].................
	*/
	
	E = ka * (10 * ((l_new * l_new) - (l_old * l_old)) + (l_new - l_old)
			* (-8 * (l[1] + l[2] + l[3] + l[4]) + 2 * (l[5] + l[6] + l[7] + l[8]) + (l[9] + l[10] + l[11] + l[12])));
	// membrane bending energy change 
	// calculated according to the Helfrich theory of membrane elasticity
	// T. R. Weikl & R. Lipowsky. Membrane adhesion and domain formation. In: Advances in planar lipid bilayers and liposomes 5, 63-127 (2006).

	return E;
}

void dynamics_particle()
// protein particle trial moves
{
	int r;
	// protein particle index

	int k;
	// membrane index: k=0 for the upper membrane, k=1 for the lower membrane

	int x_old, y_old;
	// x- and y- coordinates before the particle trial move

	int x_new, y_new;
	// x- and y-coordinates after the particle trial move

	int move;
	// direction of movement (right, left, up or down)

	double E_old, E_int_old;
	// energies before the trial move

	double E_new, E_int_new;
	// energies after the trial move

	double deltaE;
	// energy change

	r = rand() % (N0 + N1);
	// particle index selected randomly

	if (r < N0)
	// the selected particle is a receptor
	{
		k = 0;
		// the selected particle will be moved on the upper membrane
	}
	else
	// the selected particle is a ligand
	{
		k = 1;
		// the selected particle will be moved on the lower membrane
	}

	E_int_old = energy_interactions_move(r);
	// energy of interactions between the selected particle and its potential binding partner on the apposing membrane before the trial move

	E_old = energy_r_p_particle(r, k) + E_int_old;
	// energy of interactions between the selected particle and lipid rafts on the selected membrane before the trial move

	x_old = particle_x[r];
	y_old = particle_y[r];
	// x- and y-coordinates of the selected particle before the trial move

	move = rand() % 4;
	// randomly selected direction of movement

	// below: periodic boundary conditions are implemented using the bc matrix

	if (move == 0)
	{
		x_new = bc[x_old][0];
		y_new = y_old;
		// move to the right
	}
	if (move == 1)
	{
		x_new = bc[x_old][1];
		y_new = y_old;
		// move to the left
	}
	if (move == 2)
	{
		y_new = bc[y_old][0];
		x_new = x_old;
		// move up
	}
	if (move == 3)
	{
		y_new = bc[y_old][1];
		x_new = x_old;
		// move downw
	}
	if (M[x_new][y_new][k] == 0)
	// proceed only if no other protein particle is located at the membrane patch with coordinates (x_new, y_new)
	{
		
		M[x_old][y_old][k] = 0;
		PNUM[x_old][y_old][k] = -1;
		// the selected particle is removed from the membrane patch with coordinates (x_old, y_old)

		M[x_new][y_new][k] = 1;
		PNUM[x_new][y_new][k] = r;
		// the selected particle is placed on the membrane patch with coordinates (x_new, y_new)

		particle_x[r] = x_new;
		particle_y[r] = y_new;
		// the x- and y-coordinates of the moved particle are updated

		E_int_new = energy_interactions_move(r);
		// energy of interactions between the selected particle and its potential binding partner on the apposing membrane after the trial move

		E_new = energy_r_p_particle(r, k) + E_int_new;
		// energy of interactions between the selected particle and lipid rafts on the selected membrane after the trial move

		deltaE = E_new - E_old;
		// total energy change for the trial move = the raft-particle association energy change plus the receptor-ligand interaction energy change

		if (deltaE > 0)
		{
			if (dist(gen) > exp(-deltaE))
			// Metropolis criterion
			{
				// reject the trial move
				M[x_new][y_new][k] = 0;
				M[x_old][y_old][k] = 1;
				PNUM[x_new][y_new][k] = -1;
				PNUM[x_old][y_old][k] = r;
				particle_x[r] = x_old;
				particle_y[r] = y_old;
				// arrays and vectors restored to the values prior to the trial move
			}
		}
		if (E_int_new == -U_b and E_int_old == -U_b)
		// non-physical moves in which a bound protein particle directly hops into a new complex at a neighboring site are rejected
		{
			// reject the trial move
			M[x_new][y_new][k] = 0;
			M[x_old][y_old][k] = 1;
			PNUM[x_new][y_new][k] = -1;
			PNUM[x_old][y_old][k] = r;
			particle_x[r] = x_old;
			particle_y[r] = y_old;
			// arrays and vectors restored to the values prior to the trial move
		}
	}
}

void dynamics_pairs()
// receptor-ligand pair trial moves

{
	int r;
	// protein particle index

	int x_old, y_old;
	// x- and y- coordinates before the particle pair trial move

	int x_new, y_new;
	// x- and y- coordinates after the particle pair trial move

	int i0;
	// receptor index
	// receptors are numbered from 0 to N0-1

	int i1;
	// ligand index
	// ligands are numbered from N0 to N0+N1-1

	int move;
	// direction of movement (right, left, up or downw)

	double E_old;
	//energy before the trial move

	double E_new;
	//energy after the trial move

	double deltaE;
	// total energy change for the particle pair trial move

	r = rand() % (N0 + N1);
	// particle index selected randomly

	x_old = particle_x[r];
	y_old = particle_y[r];
	// x- and y-coordinates of the selected particle before the trial move

	i0 = PNUM[x_old][y_old][0];
	// receptor index
	// i0 = -1 if no receptor is located at the site with coordinates (x_old, y_old)

	i1 = PNUM[x_old][y_old][1];
	// ligand index
	// i1 = -1 if no ligand is located at the site with coordinates (x_old, y_old)

	if (i0 >= 0 and i0 < N0 and i1 >= N0 and i1 < N0 + N1)
	// proceed only if there is both a receptor and a ligand at the site with coordinates (x_old, y_old)
	{
		E_old = energy_r_p_particle(i0, 0) + energy_r_p_particle(i1, 1) + energy_interactions_move(r);
		// energy before the trial move

		move = rand() % 4;
		// randomly selected direction of movement

		// below: periodic boundary conditions are implemented using the bc matrix

		if (move == 0)
		{
			x_new = bc[x_old][0];
			y_new = y_old;
			// move to the right
		}
		if (move == 1)
		{
			x_new = bc[x_old][1];
			y_new = y_old;
			// move to the left
		}
		if (move == 2)
		{
			x_new = x_old;
			y_new = bc[y_old][0];
			// move up
		}
		if (move == 3)
		{
			x_new = x_old;
			y_new = bc[y_old][1];
			// move down
		}
		if (M[x_new][y_new][0] == 0 and M[x_new][y_new][1] == 0)
		// proceed only if no receptor and no ligand is located at the site with coordinates (x_new, y_new)
		{
			M[x_old][y_old][0] = 0;
			M[x_old][y_old][1] = 0;
			PNUM[x_old][y_old][0] = -1;
			PNUM[x_old][y_old][1] = -1;
			// the selected receptor-ligand pair is removed from the site with coordinates (x_old, y_old)

			M[x_new][y_new][0] = 1;
			M[x_new][y_new][1] = 1;
			PNUM[x_new][y_new][0] = i0;
			PNUM[x_new][y_new][1] = i1;
			// the selected receptor-ligand pair is placed on the site with coordinates (x_old, y_old)

			particle_x[i0] = x_new;
			particle_y[i0] = y_new;
			particle_x[i1] = x_new;
			particle_y[i1] = y_new;
			// the x- and y-coordinates of the moved particle pair are updated

			E_new = energy_r_p_particle(i0, 0) + energy_r_p_particle(i1, 1) + energy_interactions_move(r);
			// energy after the trial move

			deltaE = E_new - E_old;
			// total energy change consists of three terms: 
			// 1) raft-receptor association energy change
			// 2) raft-ligand association energy change
			// 3) receptor-ligand interactions energy change

			if (deltaE > 0)
			{
				if (dist(gen) > exp(-deltaE))
				//Metropolis criterion
				{
					// reject the trial move
					M[x_new][y_new][0] = 0;
					M[x_new][y_new][1] = 0;
					M[x_old][y_old][0] = 1;
					M[x_old][y_old][1] = 1;
					PNUM[x_new][y_new][0] = -1;
					PNUM[x_new][y_new][1] = -1;
					PNUM[x_old][y_old][0] = i0;
					PNUM[x_old][y_old][1] = i1;
					particle_x[i0] = x_old;
					particle_y[i0] = y_old;
					particle_x[i1] = x_old;
					particle_y[i1] = y_old;
					// arrays and vectors restored to the values prior to the trial move
				}
			}
		}
	}
}

int dynamics_bending()
// membrane bending trial move
{
	int x;
	int y;
	// lattice site coordinates

	int k;
	// membrane index: k=0 for the upper membrane, k=1 for the lower membrane
	
	x = rand() % K;
	y = rand() % K;
	k = rand() % 2;
	// randomly choose a lattice site and one of the two membranes

	double l_old;
	// membrane position at the site with coordinates (x,y) before the trial move

	double l_new;
	// membrane position at the site with coordinates (x,y) after the trial move

	double li_old;
	// local distance between the membranes before the trial move

	double li_new;
	// local distance between the membranes after the trial move

	double delta_lk;
	// change of membrane position at the site with coordinates (x,y)

	double deltaE;
	// total energy change for the membrane bending trial move

	double p;
	// random real number from 0 to 1 for implementation of the Metropolis criterion

	int iacc = 0;
	// counter 
	// iacc = 1 if the trial move is accepted
	// iacc = 0 if the trial move is rejected

	delta_lk = lmax * dis(gen);
	// change of the membrane local position
	// dis(gen) is a random number with a uniform distribution in the interval from -0.5 to 0.5

	l_old = L[x][y][k];
	// membrane local position before the trial move

	l_new = l_old + delta_lk;
	// membrane local position after the trial move

	if (k == 0)
	// the upper membrane is chosen
	{
		li_old = l_old - L[x][y][1];
		li_new = l_new - L[x][y][1];
		// local distance between the membranes before and after the trial move
	}
	else
	// the lower membrane is chosen
	{
		li_old = L[x][y][0] - l_old;
		li_new = L[x][y][0] - l_new;
		// local distance between the membranes before and after the trial move
	}
	if (li_new > 0)
	// proceed only if the upper membrane is above the lower membrane after the trial move
	{
		L[x][y][k] = l_new;
		// membrane position at the site with coordinates (x,y) is changed from l_old to l_new

		deltaE = energy_bending(x, y, k, l_old, l_new);
		// membrane bending energy change

		if (M[x][y][0] == 1 and M[x][y][1] == 1)
		// check if there is both a receptor and a ligand at the chosen site
		// change in the local distance between the membranes can influence the receptor-ligand interaction energy
		{

			deltaE = deltaE + energy_interactions_bending(li_old, li_new);
			// receptor-ligand interaction energy change
		}
		iacc = 1;
		// iacc = 1 because the trial move is accepted now

		if (deltaE > 0)
		{
			p = dist(gen) * exp(deltaE);
			if (p > 1)
			// Metropolis criterion
			{
				// reject the trial move

				L[x][y][k] = l_old;
				// membrane position at the site with coordinates (x,y) is changed back from l_new to l_old

				iacc = 0;
				// iacc = 0 because the trial move is rejected now
			}
		}
	}
	return iacc;
	// iacc is returned and used of acceptance rate calculation
}

double energy_r_r(int i, int k)
// change in the energy of interactions between lipid rafts
// i = raft index
// k = membrane index: k=0 for the upper membrane, k=1 for the lower membrane
{
	int xp, yp;
	// x- and y-coordinates of the selected raft

	int x_right, x_left, y_up, y_down;
	// x- and y-coordinates of nearest neighbor sites

	int nn;
	// number of nearest neighbor rafts

	double E;
	xp = raft_x[i];
	yp = raft_y[i];
	// x- and y-coordinates of the selected raft

	// below: periodic boundary conditions are implemented using the bc matrix to determine the coordinates of the nearest neighbor sites

	x_right = bc[xp][0];
	x_left = bc[xp][1];
	y_up = bc[yp][0];
	y_down = bc[yp][1];

	nn = N[x_right][yp][k] + N[x_left][yp][k] + N[xp][y_up][k] + N[xp][y_down][k];
	// number of nearest neighbor rafts

	E = -U * nn;
	// every adjacent raft adds -U to the raft-raft contact energy

	return E;
}

double energy_r_p_particle(int i, int k)
// raft-particle association energy for protein particle trial moves
// i = protein particle index
// k = membrane index: k=0 for the upper membrane, k=1 for the lower membrane
{
	int x, y;
	// x- and y-coordinates of the selected protein particle

	double E = 0;
	// raft-particle association energy

	x = particle_x[i];
	y = particle_y[i];

	if (M[x][y][k] == 1 and N[x][y][k] == 1)
	// raft-particle association energy = -U_a if there is both a protein particle and a raft on the same membrane patch
	// otherwise, raft-particle association energy = 0
	{
		E = -U_a;
	}
	return E;
}

double energy_r_p_raft(int i, int k)
// raft-particle association energy for raft trial moves
// i = raft index
// k = membrane index: k=0 for the upper membrane, k=1 for the lower membrane
{
	int x, y;
	// x- and y-coordinates of the selected raft

	double E = 0.0;
	// raft-particle association energy

	x = raft_x[i];
	y = raft_y[i];
	// x- and y-coordinates of the selected raft

	if (M[x][y][k] == 1 and N[x][y][k] == 1)
	// raft-particle association energy = -U_a if there is both a protein particle and a raft on the same membrane patch
	// otherwise, raft-particle association energy = 0
	{
		E = -U_a;
	}
	return E;
}

void dynamics_rafts()
// lipid raft trial moves
{
	int r;
	// raft index

	int k;
	// membrane index: k=0 for the upper membrane, k=1 for the lower membrane

	int x_old, y_old;
	// x- and y-coordinates before the trial move

	int x_new, y_new;
	// x- and y-coordinates after the trial move

	int move;
	// direction of movement

	double E_old, E_int_old;
	// energies before the trial move

	double E_new, E_int_new;
	// energies after the trial move

	double deltaE;
	// energy change

	r = rand() % (N2 + N3);
	// raft index is selected on random

	if (r < N2)
	// upper membrane raft is selected
	{
		k = 0;
		// the selected raft will be attempted to move on the upper membrane
	}
	else
	// lower membrane raft is selected
	{
		k = 1;
		// the selected raft will be attempted to move on the lower membrane
	}

	E_old = energy_r_r(r, k) + energy_r_p_raft(r, k);
	// energy before the trial move

	x_old = raft_x[r];
	y_old = raft_y[r];
	// coordinates of the raft before the trial move

	move = rand() % 4;
	// randomly selected direction of movement

	// below: periodic boundary conditions are implemented using the bc matrix

	if (move == 0)
	{
		x_new = bc[x_old][0];
		y_new = y_old;
		// move to the left
	}
	if (move == 1)
	{
		x_new = bc[x_old][1];
		y_new = y_old;
		// move to the right
	}
	if (move == 2)
	{
		x_new = x_old;
		y_new = bc[y_old][0];
		// move up
	}
	if (move == 3)
	{
		x_new = x_old;
		y_new = bc[y_old][1];
		// move down
	}
	if (N[x_new][y_new][k] == 0)
	// proceed only if the membrane patch with coordinates (x_new, y_new) is not occupied by another raft
	{
		N[x_old][y_old][k] = 0;
		RNUM[x_old][y_old][k] = -1;
		// the selected raft is removed from the membrane patch with coordinates (x_old, y_old)
		
		N[x_new][y_new][k] = 1;
		RNUM[x_new][y_new][k] = r;
		// the selected raft is placed on the membrane patch with coordinates (x_new, y_new)
		
		raft_x[r] = x_new;
		raft_y[r] = y_new;
		// the x- and y-coordinates of the moved raft are updated
		
		E_new = energy_r_r(r, k) + energy_r_p_raft(r, k);
		// energy after the trial move
		
		deltaE = E_new - E_old;
		// total energy change consists of two terms: 
		// 1) raft-raft contact energy change 
		// 2) raft-particle association energy change
		
		if (deltaE > 0)
		{
			if (dist(gen) > exp(-deltaE))
			// Metropolis criterion
			{
				// reject the trial move
				N[x_new][y_new][k] = 0;
				N[x_old][y_old][k] = 1;
				RNUM[x_new][y_new][k] = -1;
				RNUM[x_old][y_old][k] = r;
				raft_x[r] = x_old;
				raft_y[r] = y_old;
				// the selected raft is moved back from the site with coordinates (x_new, y_new) to the site with coordinates (x_old, y_old)
				// arrays and vectors are restored to the values prior to the trial move
			}
		}
	}
}

double total_me_energy(int m)
// total energy of membrane bending
// calculated according to the Helfrich theory of membrane elasticity
// m = membrane index: m=0 for the upper membrane, m=1 for the lower membrane
{
	int x_right, x_left, y_up, y_down;
	// x- and y-coordinates of adjacent sites

	double E;
	// total energy of membrane bending

	double DL;
	// discretized Laplacian

	double sumDL2 = 0.0;
	// auxiliary variable
	
	for (int x = 0; x < K; x++)
	{
		for (int y = 0; y < K; y++)
		// loop over all sites with coordinates (x,y)
		{
			x_right = bc[x][0];
			// x-coordinate of the nearest neighbor to the right

			x_left = bc[x][1];
			// x-coordinate of the nearest neighbor to the left

			y_up = bc[y][1];
			// y-coordinate of the nearest neighbor above

			y_down = bc[y][0];
			// y-coordinate of the nearest neighbor below
			
			DL = L[x_right][y][m] + L[x][y_down][m] + L[x_left][y][m] + L[x][y_up][m] - 4 * L[x][y][m];
			// discretized Laplacian at a given site

			sumDL2 = sumDL2 + DL * DL;
			// sum of (DL)^2 over all sites
		}
	}
	E = sumDL2 * kappa / (2 * a * a);
	// total energy of membrane bending 

	return E;
}

double total_r_p_energy()
// total energy of raft-particle association
{
	double E = 0;
	// total energy of raft-particle association

	int nrp = 0;
	// number of raft-particle pairs
	
	for (int x = 0; x < K; x++)
	{
		for (int y = 0; y < K; y++)
		//loop over all sites
		{
			if (M[x][y][0] == 1 and N[x][y][0] == 1)
			// upper membrane
			{
				nrp++;
			}
			if (M[x][y][1] == 1 and N[x][y][1] == 1)
			// lower membrane
			{
				nrp++;
			}
			// nrp is increased by 1 if there is both a particle and a raft at the same membrane patch
		}
	}
	E = -U_a * nrp;
	// any raft-particle pair adds -U_a to the total energy of raft-particle association

	return E;
}

double total_r_r_energy()
// total energy of contacts between rafts
{
	int x_right, x_left, y_up, y_down;
	// x- and y-coordinates of adjacent sites
	
	double nn0;
	// number of rafts neighboring to a raft with coordinates (x,y) on the upper membrane

	double nn1;
	// number of rafts neighboring to a raft with coordinates (x,y) on the lower membrane
	
	double sumnn0 = 0;
	// sum of nn0 over all sites

	double sumnn1 = 0;
	// sum of nn1 over all sites
	
	double E;
	// total energy of contacts between rafts

	for (int x = 0; x < K; x++)
	{
		for (int y = 0; y < K; y++)
		// loop over all sites with coordinates (x,y)
		{
			x_right = bc[x][0];
			// x-coordinate of the nearest neighbor to the right

			x_left = bc[x][1];
			// x-coordinate of the nearest neighbor to the left

			y_up = bc[y][0];
			// y-coordinate of the nearest neighbor above

			y_down = bc[y][1];
			// y-coordinate of the nearest neighbor below

			if (N[x][y][0] == 1)
			// there is a raft at site with coordinates (x,y) on the upper membrane
			{
				nn0 = N[x_right][y][0] + N[x_left][y][0] + N[x][y_up][0] + N[x][y_down][0];
				// number of rafts neighboring to the raft with coordinates (x,y) on the upper membrane
			}
			else
			// there is no raft at site with coordinates (x,y) on the upper membrane
			{
				nn0 = 0;
			}
			if (N[x][y][1] == 1)
			// there is a raft at site with coordinates (x,y) on the lower membrane
			{
				nn1 = N[x_right][y][1] + N[x_left][y][1] + N[x][y_up][1] + N[x][y_down][1];
				// number of rafts neighboring to the raft with coordinates (x,y) on the lower membrane
			}
			else
			// there is no raft at site with coordinates (x,y) on the lower membrane
			{
				nn1 = 0;
			}
			sumnn0 = sumnn0 + nn0;
			// sum of nn0 over all sites

			sumnn1 = sumnn1 + nn1;
			// sum of nn1 over all sites
		}
	}
	E = -U * 0.5 * (sumnn0 + sumnn1);
	// sum over all nearest-neighbor raft pairs, multiplied by 0.5 to avoid double counting
	// any raft-raft contact adds -U to the total energy of raft interactions

	return E;
}

double total_RL_energy(int pairs)
// total energy of receptor-ligand interactions
// pairs = number of receptor-ligand complexes
{
	return -U_b * pairs;
	// any receptor-ligand complex adds -U_b to the total energy of receptor-ligand interactions
}


