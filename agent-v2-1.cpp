//
// Created by FeiChenyi on 5/28/18.
//

#include "agent-v2-1.h"

int main() {
    cout << "\n\n\n\n\n\n\n";
    cout << "Program running\n";

    // The program input.
    int trial;
    srand ((trial+1)*time(NULL));
    trial = 1;
    string my_index = "v2";
    string trials = "4";
    string repeats = "8";


    NU_1 = NU_1*.1;
    NU_0 = NU_0*1;
    A_1 = A_1*1.;
    // A_0 = A_0*1000.;
    Fr0 = 0.01;
    wFr = 2.0;

    // trial = atoi(argv[1]);  //Optional input for seed
/*
	string my_index = argv[1];
	Lf = atof(argv[2])/10.0;
	E_1 = atof(argv[3]);
	NU_1 = atof(argv[4]) * (E_0) / 10.0; // surface drag coefficient


	// The ambient drag is chosen to be a fixed ratio of the surface drag.
	NU_0 = NU_1 / 100; // ambient (positional) drag
*/
    string my_cd = "/Users/FeiChenyi/Google Drive/Ongoing Works/biofilm-agent/";
    string my_output = "output";


    my_output = my_output + "/" + my_index + "/" + trials;
    // mkdir(my_output.c_str(), 0700);

/*	my_output = my_output + "/" + argv[2];
	mkdir(my_output.c_str(), 0700);

	my_output = my_output + "/" + argv[3];
	mkdir(my_output.c_str(), 0700);

	my_output = my_output + "/" + argv[4];
	mkdir(my_output.c_str(), 0700);*/


    my_name = my_output + "/pack " + repeats;

    cout << my_name << endl;


    // Initializations:
    cells.clear();
    cells.reserve(10000);


    // Create first cell:
    double ti = 0.0; // initial angle
    double nxi = cos(ti);
    double nyi = sin(ti);
    double nzi = 0.;
    if (xy == 1) { nzi = 0; }
    if (xz == 1) { nyi = 0; }
    double nm = sqrt(nxi*nxi + nyi*nyi + nzi*nzi);

    double init_xj = 8.*R;
    double init_yj = 8.*R;
    double tj = 0.; // initial angle
    double nxj = cos(tj);
    double nyj = sin(tj);
    double nzj = 0.;
    if (xy == 1) { nzj = 0; }
    if (xz == 1) { nyj = 0; }
    double nmj = sqrt(nxj*nxj + nyj*nyj + nzj*nzj);

    double init_l = L0;
    double init_z = -pow(A_1,0.6666666666666666)/
                    (2.*pow(2,0.3333333333333333)*pow(E_1,0.6666666666666666)) + R + nzi*init_l/2.0;

    if (E_1 == 0) {
        init_z = R;
    }


    ofstream myfile;

    string my_mono_name = my_cd + my_name+"-line.txt";
    myfile.open(my_mono_name);

    cout << "Declaring new cell" << endl;

    Cell * mb = new Cell(0, 0, init_z, nxi/nm, nyi/nm, nzi/nm, init_l);
    Cell * mr = new Cell(init_xj, init_yj, init_z, nxj/nmj, nyj/nmj, nzj/nmj, init_l);
    mb->set_lin("bW");
    mr->set_lin("rA");

    cout << "Pushing new cell" << endl;

    cells.push_back(mb);
    cells.push_back(mr);

    myfile << 0 << " " << cells.size() << endl;
    for(int j = 0; j <cells.size(); j++)
    {
        myfile << cells[j]->get_x() << " " << cells[j]->get_y() << " " << cells[j]->get_z() << " "
               << cells[j]->get_nx() << " " << cells[j]->get_ny() << " " << cells[j]->get_nz() << " "
               << cells[j]->get_l() << " " << cells[j]->get_ai() << " " << 1 << " " << 0 << " " << 0 << " " << 0 << " " << cells[j]->get_lin() << endl;
    }



    cout << "Calling grow" << endl;

    simple_grow(T, myfile);


    myfile.close();

    cout << "Done" << endl;

    return 0;
}